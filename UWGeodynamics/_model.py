from __future__ import print_function
import os
import sys
import numpy as np
import underworld as uw
import underworld.function as fn
import UWGeodynamics.shapes as shapes
import UWGeodynamics.surfaceProcesses as surfaceProcesses
import json
from json import JSONEncoder
from UWGeodynamics import rcParams
from .scaling import Dimensionalize
from .scaling import nonDimensionalize as nd
from .scaling import UnitRegistry as u
from .lithopress import LithostaticPressure
from ._utils import PressureSmoother, PassiveTracers
from ._rheology import ViscosityLimiter
from ._material import Material
from ._plots import Plots
from ._visugrid import Visugrid
from ._velocity_boundaries import VelocityBCs
from ._thermal_boundaries import TemperatureBCs
from ._mesh_advector import _mesh_advector
from ._frictional_boundary import FrictionBoundaries
from collections import OrderedDict


_attributes_to_save = {
    "elementRes": lambda x: tuple(int(val) for val in x.split(",")) ,
    "minCoord": lambda x: tuple([val if u.Quantity(val).dimensionless else u.Quantity(val) for val in x.split(",")]),
    "maxCoord": lambda x: tuple([val if u.Quantity(val).dimensionless else u.Quantity(val) for val in x.split(",")]),
    "name": lambda x: str(x),
    "gravity": lambda x: tuple([val if u.Quantity(val).dimensionless else u.Quantity(val) for val in x.split(",")]),
    "periodic": lambda x: tuple(bool(val) for val in x.split(",")),
    "elementType": lambda x: str(x),
    "Tref": lambda x: x if u.Quantity(x).dimensionless else u.Quantity(x)
    }

class _ModelEncoder(JSONEncoder):

    def default(self, obj):
        model = {}
        for attribute in _attributes_to_save:
            val = obj[attribute]
            if isinstance(val, (list, tuple)):
                model[attribute] = ", ".join([str(v) for v in val]) 
            else:
                model[attribute] = str(val)
        
        materials = []
        for material in obj["materials"]:
            string = material._json_()
            materials.append(eval(string))

        model["materials"] = materials

        return model


class Model(Material):
    """ This class provides the main UWGeodynamics Model

    Attributes
    ----------

    materialField:
    pressureField:
    velocityField:
    temperature:
    tractionField:

    """

    def __init__(self, elementRes, minCoord, maxCoord,
                 name=None, gravity=None,
                 periodic=None, elementType=None,
                 Tref=273.15 * u.degK):

        """
        Parameters
        ----------
        elementRes:
            Resolution of the mesh in number of elements for each axis (degree
            of freedom)
        minCoord:
            Minimum coordinates for each axis.
        maxCoord:
            Maximum coordinates for each axis.
        name:
            The Model name.
        gravity:
            Acceleration due to gravity vector.
        periodic:
            Mesh periodicity.
        elementType:
            Type of finite element used (Only Q1/dQ0 are currently supported

        Returns
        --------
        Model Class Object

        Examples
        --------
        >>> import UWGeodynamics as GEO
        >>> u = GEO.UnitRegistry
        >>> Model = Model = GEO.Model(elementRes=(64, 64),
                            minCoord=(-1. * u.meter, -50. * u.centimeter),
                            maxCoord=(1. * u.meter, 50. * u.centimeter))

        """

        super(Model, self).__init__()
        
        # Process __init__ arguments
        if not name:
            self.name = rcParams["model.name"]
        else:
            self.name = name

        self.minCoord = minCoord
        self.maxCoord = maxCoord
        self.top = maxCoord[-1]
        self.bottom = minCoord[-1]
        self.dimension = dim = len(maxCoord)
        
        if not gravity:
            self.gravity = [0.0 for val in range(self.dimension)]
            self.gravity[-1] = -1.0 * rcParams["gravity"]
            self.gravity = tuple(self.gravity)
        else:
            self.gravity = gravity
        
        self.Tref = Tref

        if not elementType:
            self.elementType = rcParams["element.type"]
        else:
            self.elementType = elementType

        self.elementRes = elementRes
        self.outputDir = rcParams["output.directory"]
        
        # Compute model dimensions
        self.length = maxCoord[0] - minCoord[0]
        self.height = maxCoord[-1] - minCoord[-1]

        if dim == 3:
            self.width = maxCoord[1] - minCoord[1]

        if periodic:
            self.periodic = periodic
        else:
            periodic = tuple([False for val in range(dim)])
            self.periodic = periodic

        # Get non-dimensional extents along each axis
        minCoord = tuple([nd(val) for val in self.minCoord])
        maxCoord = tuple([nd(val) for val in self.maxCoord])

        # Initialize model mesh
        self.mesh = uw.mesh.FeMesh_Cartesian(elementType=self.elementType,
                                             elementRes=self.elementRes,
                                             minCoord=minCoord,
                                             maxCoord=maxCoord,
                                             periodic=self.periodic)

        # Add common mesh variables
        # Temperature Field is initialised to None
        self.temperature = None
        self.pressureField = uw.mesh.MeshVariable(mesh=self.mesh.subMesh,
                                                  nodeDofCount=1)
        self.velocityField = uw.mesh.MeshVariable(mesh=self.mesh,
                                                  nodeDofCount=self.mesh.dim)
        self.tractionField = uw.mesh.MeshVariable(mesh=self.mesh,
                                                  nodeDofCount=self.mesh.dim)
        self._strainRateField = uw.mesh.MeshVariable(mesh=self.mesh,
                                                     nodeDofCount=1)

        # Initialise fields to 0.
        self.pressureField.data[...] = 0.
        self.velocityField.data[...] = 0.
        self.tractionField.data[...] = 0.

        # symmetric component of the gradient of the flow velocityField.
        self.strainRate = fn.tensor.symmetric(self.velocityField.fn_gradient)
        self._strainRate_2ndInvariant = fn.tensor.second_invariant(
            self.strainRate
        )

        # Create the material swarm
        self.swarm = uw.swarm.Swarm(mesh=self.mesh, particleEscape=True)
        self._swarmLayout = uw.swarm.layouts.GlobalSpaceFillerLayout(
            swarm=self.swarm,
            particlesPerCell=rcParams["swarm.particles.per.cell"])

        self.swarm.populate_using_layout(layout=self._swarmLayout)

        self.population_control = uw.swarm.PopulationControl(
            self.swarm,
            aggressive=rcParams["popcontrol.aggressive"],
            splitThreshold=rcParams["popcontrol.split.threshold"],
            maxSplits=rcParams["popcontrol.max.splits"],
            particlesPerCell=rcParams["popcontrol.particles.per.cell"])

        # timing and checkpointing
        self.checkpointID = 0
        self.time = 0.0 * u.megayears
        self.step = 0
        self._dt = None

        # viscosity limiter
        self.minViscosity = rcParams["minimum.viscosity"]
        self.maxViscosity = rcParams["maximum.viscosity"]
        self._viscosity_limiter = ViscosityLimiter(self.minViscosity,
                                                   self.maxViscosity)

        self._defaultMaterial = self.index
        self.materials = [self]
        self._material_drawing_order = None

        # Create a series of aliases for the boundary sets
        self._left_wall = self.mesh.specialSets["MinI_VertexSet"]
        self._right_wall = self.mesh.specialSets["MaxI_VertexSet"]
        
        if self.mesh.dim == 2:
            self._top_wall = self.mesh.specialSets["MaxJ_VertexSet"]
            self._bottom_wall = self.mesh.specialSets["MinJ_VertexSet"]
        else:
            self._front_wall = self.mesh.specialSets["MinJ_VertexSet"]
            self._back_wall = self.mesh.specialSets["MaxJ_VertexSet"]
            self._top_wall = self.mesh.specialSets["MaxK_VertexSet"]
            self._bottom_wall = self.mesh.specialSets["MinK_VertexSet"]

        # Boundary Conditions
        self.velocityBCs = None
        self.temperatureBCs = None
        self.frictionalBCs = None
        self._isostasy = None
        self.surfaceProcesses = None
        self._surfaceProcesses = None

        self.pressSmoother = PressureSmoother(self.mesh, self.pressureField)

        # Passive Tracers
        self.passive_tracers = []

        # Visualisation
        self.plot = Plots(self)
        self._visugrid = None

        # Mesh advector
        self._advector = None
        
        # Initialise remaining attributes
        self._default_strain_rate = 1e-30 / u.second
        self._solution_exist = False
        self._isYielding = None
        self._temperatureDot = None
        self._temperature = None
        self._outputDir = None
        self.DiffusivityFn = None
        self.HeatProdFn = None
        self._buoyancyFn = None

        self._initialize()

    def _initialize(self):

        self.swarm_advector = uw.systems.SwarmAdvector(
            swarm=self.swarm,
            velocityField=self.velocityField,
            order=2
        )

        # Add Common Swarm Variables
        self.materialField = self.swarm.add_variable(dataType="int", count=1)
        self.plasticStrain = self.swarm.add_variable(dataType="double",
                                                     count=1)
        self._viscosityField = self.swarm.add_variable(dataType="double",
                                                       count=1)
        self._densityField = self.swarm.add_variable(dataType="double",
                                                     count=1)
        self.meltField = self.swarm.add_variable(dataType="double", count=1)

        # Initialise materialField to Model material
        self.materialField.data[:] = self.index
        
        # Initialise remaining fields to 0.
        self.plasticStrain.data[:] = 0.0
        self._viscosityField.data[:] = 0.
        self._densityField.data[:] = 0.
        self.meltField.data[:] = 0.

        # Create a bunch of tools to project swarmVariable onto the mesh
        self._projMaterialField = uw.mesh.MeshVariable(mesh=self.mesh,
                                                       nodeDofCount=1)

        self._materialFieldProjector = uw.utils.MeshVariable_Projection(
            self._projMaterialField, self.materialField, type=0)

        self._projViscosityField = uw.mesh.MeshVariable(
            mesh=self.mesh, nodeDofCount=1)

        self._viscosityFieldProjector = uw.utils.MeshVariable_Projection(
            self._projViscosityField, self._viscosityField, type=0)

        self._projPlasticStrain = uw.mesh.MeshVariable(mesh=self.mesh,
                                                       nodeDofCount=1)
        self._plasticStrainProjector = uw.utils.MeshVariable_Projection(
            self._projPlasticStrain, self.plasticStrain, type=0)

        self._projDensityField = uw.mesh.MeshVariable(mesh=self.mesh,
                                                      nodeDofCount=1)
        self._densityFieldProjector = uw.utils.MeshVariable_Projection(
            self._projDensityField, self._densityField, type=0)

    def __getitem__(self, name):
        return self.__dict__[name]

    def _repr_html_(self):
        return _model_html_repr(self) 
    
    @property
    def x(self):
        return fn.input()[0]

    @property
    def y(self):
        return fn.input()[1]

    @property
    def z(self):
        return fn.input()[2]

    @property
    def outputDir(self):
        """ Output Directory """
        return self._outputDir

    @outputDir.setter
    def outputDir(self, value):
        if uw.rank() == 0:
            if not os.path.exists(value):
                os.makedirs(value)
        self._outputDir = value

    def restart(self, restartDir=None, step=None):
        """ Restart a Model

        parameters:
        -----------
            restartDir: (string)
                Directory that contains the outputs of the model
                you want to restart
            step: (int)
                Step from which you want to restart the model.

        """
        if not restartDir:
            restartDir = self.outputDir
        if not step:
            step = max([int(os.path.splitext(filename)[0].split("-")[-1])
                        for filename in os.listdir(restartDir) if "-" in
                        filename])

        self.checkpointID = step
        self.mesh.load(os.path.join(restartDir, "mesh.h5"))
        self.swarm = uw.swarm.Swarm(mesh=self.mesh, particleEscape=True)
        self.swarm.load(os.path.join(restartDir, 'swarm-%s.h5' % step))
        self._initialize()
        self.materialField.load(os.path.join(restartDir,
                                             "materialField-%s.h5" % step))
        self.temperature.load(os.path.join(restartDir,
                                           'temperature-%s.h5' % step))
        self.pressureField.load(os.path.join(restartDir,
                                             'pressureField-%s.h5' % step))
        self.plasticStrain.load(os.path.join(restartDir,
                                             'plasticStrain-%s.h5' % step))
        self.velocityField.load(os.path.join(restartDir,
                                             'velocityField-%s.h5' % step))

    @property
    def projMaterialField(self):
        """ Material field projected on the mesh """
        self._materialFieldProjector.solve()
        return self._projMaterialField

    @property
    def projPlasticStrain(self):
        """ Plastic Strain Field projected on the mesh """
        self._plasticStrainProjector.solve()
        return self._projPlasticStrain

    @property
    def strainRateField(self):
        """ Strain Rate Field """
        self._strainRateField.data[:] = self._strainRate_2ndInvariant.evaluate(
            self.mesh)
        return self._strainRateField

    @property
    def projViscosityField(self):
        """ Viscosity Field projected on the mesh """
        self.viscosityField.data[...] = self._viscosityFn.evaluate(self.swarm)
        self._viscosityFieldProjector.solve()
        return self._projViscosityField

    @property
    def viscosityField(self):
        """ Viscosity Field on particles """
        self._viscosityField.data[:] = self._viscosityFn.evaluate(self.swarm)
        return self._viscosityField

    @property
    def projDensityField(self):
        """ Density Field projected on the mesh """
        self.densityField.data[...] = self._densityFn.evaluate(self.swarm)
        self._densityFieldProjector.solve()
        return self._projDensityField

    @property
    def densityField(self):
        """ Density Field on particles """
        self._densityField.data[:] = self._densityFn.evaluate(self.swarm)
        return self._densityField

    @property
    def surfaceProcesses(self):
        """ Surface processes handler """
        return self._surfaceProcesses

    @surfaceProcesses.setter
    def surfaceProcesses(self, value):
        self._surfaceProcesses = value
        if isinstance(value, surfaceProcesses.Badlands):
            self._surfaceProcesses.Model = self

    def set_temperatureBCs(self, left=None, right=None, top=None, bottom=None,
                           front=None, back=None,
                           indexSets=None, materials=None):
        """ Set Model thermal conditions

        Parameters:
            left:
                Temperature or flux along the left wall.
                Default is 'None'
            right:
                Temperature or flux along the right wall.
                Default is 'None'
            top:
                Temperature or flux along the top wall.
                Default is 'None'
            bottom:
                Temperature or flux along the bottom wall.
                Default is 'None'
            front:
                Temperature or flux along the front wall.
                Default is 'None'
            back:
                Temperature or flux along the back wall.
                Default is 'None'
            indexSets: (set, temperature)
                underworld mesh index set and associate temperature
            materials:
                list of materials for which temperature need to be
                fixed. [(material, temperature)]

        """

        if not self.temperature:
            self.temperature = uw.mesh.MeshVariable(mesh=self.mesh,
                                                    nodeDofCount=1)
            self._temperatureDot = uw.mesh.MeshVariable(mesh=self.mesh,
                                                        nodeDofCount=1)
            self.temperature.data[...] = nd(self.Tref)
            self._temperatureDot.data[...] = 0.

        self.temperatureBCs = TemperatureBCs(self, left=left, right=right,
                                             top=top, bottom=bottom,
                                             back=back, front=front,
                                             indexSets=indexSets,
                                             materials=materials)
        return self.temperatureBCs.get_conditions()

    @property
    def temperature(self):
        """ Temperature Field """
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        self._temperature = value

    @property
    def _temperatureBCs(self):
        if not self.temperatureBCs:
            raise ValueError("Set Boundary Conditions")
        return self.temperatureBCs.get_conditions()

    @property
    def _advdiffSystem(self):
        """ Advection Diffusion System """

        DiffusivityMap = {}
        for material in self.materials:
            DiffusivityMap[material.index] = nd(material.diffusivity)

        self.DiffusivityFn = fn.branching.map(fn_key=self.materialField,
                                              mapping=DiffusivityMap)

        HeatProdMap = {}
        for material in self.materials:
            if all([material.density,
                    material.capacity,
                    material.radiogenicHeatProd]):

                HeatProdMap[material.index] = (nd(material.radiogenicHeatProd) /
                                               (nd(material.density) *
                                                nd(material.capacity)))
            else:
                HeatProdMap[material.index] = 0.

            # Melt heating
            if material.latentHeatFusion and self._dt:
                dynamicHeating = self._get_dynamic_heating(material)
                HeatProdMap[material.index] += dynamicHeating

        self.HeatProdFn = fn.branching.map(fn_key=self.materialField,
                                           mapping=HeatProdMap)

        obj = uw.systems.AdvectionDiffusion(
            self.temperature,
            self._temperatureDot,
            velocityField=self.velocityField,
            fn_diffusivity=self.DiffusivityFn,
            fn_sourceTerm=self.HeatProdFn,
            conditions=[self._temperatureBCs]
        )
        return obj

    @property
    def _stokes(self):
        """ Stokes solver """

        gravity = tuple([nd(val) for val in self.gravity])
        self._buoyancyFn = self._densityFn * gravity

        if any([material.viscosity for material in self.materials]):

            stokes_object = uw.systems.Stokes(velocityField=self.velocityField,
                                              pressureField=self.pressureField,
                                              conditions=self._velocityBCs,
                                              fn_viscosity=self._viscosityFn,
                                              fn_bodyforce=self._buoyancyFn,
                                              fn_one_on_lambda=None)

            solver = uw.systems.Solver(stokes_object)
            solver.set_inner_method(rcParams["solver"])
            solver.set_penalty(rcParams["penalty"])

        return solver

    def _init_melt_fraction(self):
        """ Initialize the Melt Fraction Field """

        # Initialize Melt Fraction to material property
        meltFractionMap = {}
        for material in self.materials:
            if material.meltFraction:
                meltFractionMap[material.index] = material.meltFraction

        if meltFractionMap:
            InitFn = fn.branching.map(fn_key=self.materialField,
                                      mapping=meltFractionMap, fn_default=0.0)
            self.meltField.data[:] = InitFn.evaluate(self.swarm)

    def set_velocityBCs(self, left=None, right=None, top=None, bottom=None,
                        front=None, back=None, indexSets=None):
        """ Set Model kinematic conditions

        Parameters:
            left:
                Velocity along the left wall.
                Default is 'None'
            right:
                Velocity along the right wall.
                Default is 'None'
            top:
                Velocity along the top wall.
                Default is 'None'
            bottom:
                Velocity along the bottom wall.
                Default is 'None'
            front:
                Velocity along the front wall.
                Default is 'None'
            back:
                Velocity along the back wall.
                Default is 'None'
            indexSets: (set, velocity)
                underworld mesh index set and associate velocity
        """
        self._velocityBCs_saved_args = locals()
        del(self._velocityBCs_saved_args["self"])

        self.velocityBCs = VelocityBCs(self, left=left,
                                       right=right, top=top,
                                       bottom=bottom, front=front,
                                       back=back, indexSets=indexSets)
        return self.velocityBCs.get_conditions()

    @property
    def _velocityBCs(self):
        """ Retrieve kinematic boundary conditions """
        if not self.velocityBCs:
            raise ValueError("Set Boundary Conditions")
        return self.velocityBCs.get_conditions()

    def add_material(self, material=None, shape=None, name="unknown", reset=False):
        """ Add Material to the Model

        Parameters:
        -----------
            material:
                An UWGeodynamics material. If None the material is
                initialized to the global properties.
            shape:
                Shape of the material. See UWGeodynamics.shape
            name:
                Material name
            reset: (bool)
                Reset the material Field before adding the new
                material. Default is False.
        """

        if reset:
            self.materialField.data[...] = 0
            self.materials = [self]

        if material is not None:
            mat = material
            mat.name = name
        else:
            mat = Material()

            mat.name = name

            # Initialize some properties to Model property
            mat.diffusivity = self.diffusivity
            mat.capacity = self.capacity
            mat.thermalExpansivity = self.thermalExpansivity
            mat.radiogenicHeatProd = self.radiogenicHeatProd

        if isinstance(shape, shapes.Layer):
            mat.top = shape.top
            mat.bottom = shape.bottom

            if self.mesh.dim == 3:
                shape.minY = self.minCoord[1]
                shape.maxY = self.maxCoord[1]

        mat.shape = shape
        mat.indices = self._get_material_indices(mat)
        self.materials.reverse()
        self.materials.append(mat)
        self.materials.reverse()
        self._fill_model()

        return mat

    def _fill_model(self):
        """ Initialize the Material Field from the list of materials"""

        conditions = [(obj.shape.fn, obj.index)
                      for obj in self.materials if obj.shape is not None]

        conditions.append((True, self._defaultMaterial))
        vals = fn.branching.conditional(conditions).evaluate(self.swarm)
        self.materialField.data[:] = vals

    @property
    def material_drawing_order(self):
        """ Order of material drawing. Overlapping shapes reset the material"""
        return self._material_drawing_order

    @material_drawing_order.setter
    def material_drawing_order(self, value):
        if value:
            self.materials = value
            self._material_drawing_order = value
            self._fill_model()

    @property
    def _densityFn(self):
        densityMap = {}
        for material in self.materials:

            if self.temperature:
                densityMap[material.index] = nd(material.density) * (1.0 -
                                             nd(material.thermalExpansivity) *
                                             (self.temperature - nd(self.Tref)))
            else:
                densityMap[material.index] = nd(material.density)

            if material.meltExpansion:
                fact = material.meltExpansion * self.meltField
                densityMap[material.index] = (densityMap[material.index] *
                                              (1.0 - fact))

        return fn.branching.map(fn_key=self.materialField, mapping=densityMap)

    def set_frictional_boundary(self, right=False, left=False,
                                top=False, bottom=False,
                                thickness=2, friction=0.0):
        """ Set Frictional Boundary conditions

        Frictional boundaries are implemented as a thin layer of frictional
        material along the chosen boundaries.

        Parameters:
        -----------

            right: (bool)
                if True create a frictional boundary.
            left: (bool)
                if True create a frictional boundary.
            top: (bool)
                if True create a frictional boundary.
            bottom: (bool)
                if True create a frictional boundary.
            thickness: (int)
                thickness of the boundaries (in number of
                elements)
            friction: (float)
                Friction coefficient at the boundaries.

        Returns:
        --------
            Underworld mesh variable that maps the boundaries.
            (Values are set to 1 if the mesh element applies
            a frictional condition, 0 otherwise).
        """

        self.frictionalBCs = FrictionBoundaries(self, right=right,
                                                left=left,
                                                top=top,
                                                bottom=bottom,
                                                thickness=thickness,
                                                friction=friction)

        return self.frictionalBCs

    @property
    def _viscosityFn(self):
        """ Create the Viscosity Function """

        ViscosityMap = {}
        BGViscosityMap = {}

        # Viscous behavior
        for material in self.materials:
            if material.viscosity:
                ViscosityHandler = material.viscosity
                ViscosityHandler.pressureField = self.pressureField
                ViscosityHandler.strainRateInvariantField = (
                    self._strainRate_2ndInvariant)
                ViscosityHandler.temperatureField = self.temperature
                ViscosityHandler.viscosityLimiter = self._viscosity_limiter
                ViscosityMap[material.index] = ViscosityHandler.muEff

        # Melt Modifier
        for material in self.materials:
            if material.viscosity and material.viscosityChange > 1.0:
                X1 = material.viscosityChangeX1
                X2 = material.viscosityChangeX2
                change = (1.0 + (material.viscosityChange - 1.0) /
                          (X2 - X1) * (self.meltField - X1))
                conditions = [(self.meltField < X1, 1.0),
                              (self.meltField > X2, material.viscosityChange),
                              (True, change)]
                ViscosityMap[material.index] *= fn.branching.conditional(
                    conditions)

        # Plasticity
        PlasticityMap = {}
        for material in self.materials:

            if material.plasticity:

                YieldHandler = material.plasticity
                YieldHandler.pressureField = self.pressureField
                YieldHandler.plasticStrain = self.plasticStrain
                if self.mesh.dim == 2:
                    yieldStress = YieldHandler._get_yieldStress2D()
                if self.mesh.dim == 3:
                    yieldStress = YieldHandler._get_yieldStress3D()
                eijdef = nd(self._default_strain_rate)
                eij = fn.branching.conditional(
                    [(self._strainRate_2ndInvariant < sys.float_info.epsilon,
                      eijdef),
                     (True, self._strainRate_2ndInvariant)])
                muEff = 0.5 * yieldStress / eij
                muEff = self._viscosity_limiter.apply(muEff)
                PlasticityMap[material.index] = muEff

            if self.frictionalBCs is not None:
                from copy import copy

                # Only affect plastic materials
                if material.plasticity:

                    YieldHandler = copy(material.plasticity)
                    YieldHandler.frictionCoefficient = (
                        self.frictionalBCs.friction)
                    YieldHandler.pressureField = self.pressureField
                    YieldHandler.plasticStrain = self.plasticStrain
                    if self.mesh.dim == 2:
                        yieldStress = YieldHandler._get_yieldStress2D()
                    if self.mesh.dim == 3:
                        yieldStress = YieldHandler._get_yieldStress3D()
                    eijdef = nd(self._default_strain_rate)
                    eij = fn.branching.conditional(
                        [(self._strainRate_2ndInvariant <
                          sys.float_info.epsilon,
                          eijdef),
                         (True, self._strainRate_2ndInvariant)])
                    muEff = 0.5 * yieldStress / eij
                    muEff = self._viscosity_limiter.apply(muEff)

                    conditions = [(self.frictionalBCs._mask == 1, muEff),
                                  (True, PlasticityMap[material.index])]

                    PlasticityMap[material.index] = fn.branching.conditional(
                        conditions
                    )

        # Combine rheologies
        EffViscosityMap = {}
        PlasticMap = {}
        for material in self.materials:
            idx = material.index
            if material.viscosity and material.plasticity:
                EffViscosityMap[idx] = fn.misc.min(PlasticityMap[idx],
                                                   ViscosityMap[idx])
                BGViscosityMap[idx] = ViscosityMap[idx]
                PlasticMap[idx] = 0.
            elif material.viscosity:
                EffViscosityMap[idx] = ViscosityMap[idx]
                BGViscosityMap[idx] = ViscosityMap[idx]
                PlasticMap[idx] = 0.
            elif material.plasticity:
                EffViscosityMap[idx] = PlasticityMap[idx]
                BGViscosityMap[idx] = PlasticityMap[idx]
                PlasticMap[idx] = 1.0

        viscosityFn = fn.branching.map(fn_key=self.materialField,
                                       mapping=EffViscosityMap)

        backgroundViscosityFn = fn.branching.map(fn_key=self.materialField,
                                                 mapping=BGViscosityMap)

        isPlastic = fn.branching.map(fn_key=self.materialField,
                                     mapping=PlasticMap)
        yieldConditions = [(viscosityFn < backgroundViscosityFn, 1.0),
                           (isPlastic > 0.5, 1.0),
                           (True, 0.0)]

        # Do not yield at the very first solve
        if self._solution_exist:
            self._isYielding = (fn.branching.conditional(yieldConditions) *
                               self._strainRate_2ndInvariant)
        else:
            self._isYielding = fn.misc.constant(0.0)

        return viscosityFn

    @property
    def _yieldStressFn(self):
        """ Calculate Yield stress function from viscosity and strain rate"""
        eij = self._strainRate_2ndInvariant
        eijdef = nd(self._default_strain_rate)
        return 2.0 * self._viscosityFn * fn.misc.max(eij, eijdef)

    def solve_temperature_steady_state(self):
        """ Solve for steady state temperature

        Returns:
        --------
            Updated temperature Field
        """

        if self.materials:

            DiffusivityMap = {}
            for material in self.materials:
                DiffusivityMap[material.index] = nd(material.diffusivity)

            self.DiffusivityFn = fn.branching.map(fn_key=self.materialField,
                                                  mapping=DiffusivityMap)

            HeatProdMap = {}
            for material in self.materials:

                if all([material.density,
                        material.capacity,
                        material.radiogenicHeatProd]):

                    HeatProdMap[material.index] = (
                        nd(material.radiogenicHeatProd) /
                        (nd(material.density) *
                         nd(material.capacity))
                    )

                else:
                    HeatProdMap[material.index] = 0.

            self.HeatProdFn = fn.branching.map(fn_key=self.materialField,
                                               mapping=HeatProdMap)
        else:
            self.DiffusivityFn = fn.misc.constant(nd(self.diffusivity))
            self.HeatProdFn = fn.misc.constant(nd(self.radiogenicHeatProd))

        conditions = self._temperatureBCs

        heatequation = uw.systems.SteadyStateHeat(
            temperatureField=self.temperature,
            fn_diffusivity=self.DiffusivityFn,
            fn_heating=self.HeatProdFn,
            conditions=conditions
        )

        heatsolver = uw.systems.Solver(heatequation)
        heatsolver.solve(nonLinearIterate=True)

        return self.temperature

    def get_lithostatic_pressureField(self):
        """ Calculate the lithostatic Pressure Field

        Returns:
        --------
            Pressure Field, array containing the pressures
            at the bottom of the Model
        """

        gravity = np.abs(nd(self.gravity[-1]))  # Ugly!!!!!
        lithoPress = LithostaticPressure(self.mesh, self._densityFn, gravity)
        self.pressureField.data[:], LPressBot = lithoPress.solve()
        self.pressSmoother.smooth()

        return self.pressureField, LPressBot

    def _calibrate_pressureField(self):
        """ Pressure Calibration callback function """

        surfaceArea = uw.utils.Integral(fn=1.0, mesh=self.mesh,
                                        integrationType='surface',
                                        surfaceIndexSet=self._top_wall)

        surfacepressureFieldIntegral = uw.utils.Integral(
            fn=self.pressureField,
            mesh=self.mesh,
            integrationType='surface',
            surfaceIndexSet=self._top_wall
        )

        area, _ = surfaceArea.evaluate()
        p0, _ = surfacepressureFieldIntegral.evaluate()
        offset = p0 / area
        self.pressureField.data[:] -= offset
        self.pressSmoother.smooth()

        for material in self.materials:
            if material.viscosity:
                material.viscosity.firstIter = False

    def _get_material_indices(self, material):
        """ Get mesh indices of a Material

        Parameter:
        ----------
            material:
                    The Material of interest
        Returns;
        --------
            underworld IndexSet

        """
        nodes = np.arange(0, self.mesh.nodesLocal)[self.projMaterialField.evaluate(self.mesh.data[0:self.mesh.nodesLocal])[:, 0] == material.index]
        return uw.mesh.FeMesh_IndexSet(self.mesh, topologicalIndex=0, size=self.mesh.nodesGlobal, fromObject=nodes)

    def solve(self):
        """ Solve Stokes """
        self._stokes.solve(nonLinearIterate=True,
                          callback_post_solve=self._calibrate_pressureField,
                          nonLinearTolerance=rcParams["nonlinear.tolerance"])
        self._solution_exist = True

    def init_model(self, temperature=True, pressureField=True):
        """ Initialize the Temperature Field as steady state,
            Initialize the Pressure Field as Lithostatic

        Parameters:
        -----------
            temperature: (bool) default to True
            pressure: (bool) default to True
        """

        # Init Temperature Field
        if self.temperature and temperature:
            self.solve_temperature_steady_state()

        # Init pressureField Field
        if self.pressureField and pressureField:
            self.get_lithostatic_pressureField()
        return

    def run_for(self, duration=None, checkpoint=None, timeCheckpoints=None):
        """ Run the Model

        Parameters:
        -----------

            duration: duration of the Model in Time units or nb of steps
            checkpoint: checkpoint interval
            timeCheckpoints: Additional checkpoints
        """

        step = self.step
        time = nd(self.time)
        units = duration.units
        duration = time + nd(duration)

        next_checkpoint = None

        if timeCheckpoints:
            timeCheckpoints = [nd(val) for val in timeCheckpoints]

        if checkpoint:
            next_checkpoint = time + nd(checkpoint)
        
        #pbar = tqdm(total=duration-time)
        while time < duration:
            self.solve()

            # Whats the longest we can run before reaching the end of the model
            # or a checkpoint?
            # Need to generalize that
            dt = self.swarm_advector.get_max_dt()

            if self.temperature:
                dt = min(dt, self._advdiffSystem.get_max_dt())

            if checkpoint:
                dt = min(dt, next_checkpoint - time)

            self._dt = min(dt, duration - time)

            uw.barrier()

            self._update()

            step += 1
            self.time += Dimensionalize(self._dt, units)
            time += self._dt

            if time == next_checkpoint:
                self.checkpointID += 1
                self.checkpoint()
                self.output_glucifer_figures(self.checkpointID)
                next_checkpoint += nd(checkpoint)

            #pbar.update(self._dt)
            if checkpoint or step % 1 == 0:
                if uw.rank() == 0:
                    print("Time: ", str(self.time.to(units)),
                          'dt:', str(Dimensionalize(self._dt, units)))
        return 1

    def _update(self):
        """ Update Function

        The following function processes the mesh and swarm variables
        between two solves. It takes care of mesh, swarm advection and
        update the fields according to the Model state.

        """

        dt = self._dt
        # Increment plastic strain
        plasticStrainIncrement = dt * self._isYielding.evaluate(self.swarm)
        self.plasticStrain.data[:] += plasticStrainIncrement

        if any([material.melt for material in self.materials]):
            # Calculate New meltField
            meltFraction = self._get_melt_fraction()
            self.meltField.data[:] = meltFraction.evaluate(self.swarm)

        # Solve for temperature
        if self.temperature:
            self._advdiffSystem.integrate(dt)

        # Integrate Swarms in time
        self.swarm_advector.integrate(dt, update_owners=True)

        if self.passive_tracers:
            for tracers in self.passive_tracers:
                tracers.integrate(dt)

        if self._advector:
            self._advector.advect_mesh(dt)

        # Do pop control
        self.population_control.repopulate()

        if self.surfaceProcesses:
            self.surfaceProcesses.solve(dt)

        if self._isostasy:
            self._isostasy.solve()

        if self._visugrid:
            self._visugrid.advect(dt)

    def mesh_advector(self, axis):
        """ Initialize the mesh advector

        Parameters:
        -----------
            axis:
                list of axis (or degree of freedom) along which the
                mesh is allowed to deform
        """
        self._advector = _mesh_advector(self, axis)

    def add_passive_tracers(self, name=None, vertices=None,
                            particleEscape=True):
        """ Add a swarm of passive tracers to the Model

        Parameters:
        -----------
            name:
                Name of the swarm of tracers.
            vertices:
                Numpy array that contains the coordinates of the tracers.
            particleEscape: (bool)
                Allow or prevent tracers from escaping the boundaries of the
                Model (default to True)
        """
        tracers = PassiveTracers(self.mesh,
                                 self.velocityField,
                                 name=name,
                                 vertices=vertices,
                                 particleEscape=particleEscape)

        self.passive_tracers.append(tracers)

        return tracers

    def _get_melt_fraction(self):
        """ Melt Fraction function

        Returns:
        -------
            Underworld function that calculates the Melt fraction on the
            particles.

        """
        meltMap = {}
        for material in self.materials:
            if material.melt:
                T_s = material.solidus.temperature(self.pressureField)
                T_l = material.liquidus.temperature(self.pressureField)
                T_ss = (self.temperature - 0.5 * (T_s + T_l)) / (T_l - T_s)
                value = (0.5 + T_ss + (T_ss * T_ss - 0.25) *
                         (0.4256 + 2.988 * T_ss))
                conditions = [((-0.5 < T_ss) & (T_ss < 0.5),
                               fn.misc.min(value, material.meltFractionLimit)),
                              (True, 0.0)]
                meltMap[material.index] = fn.branching.conditional(conditions)

        return fn.branching.map(fn_key=self.materialField,
                                mapping=meltMap, fn_default=0.0)

    def _get_dynamic_heating(self, material):
        """ Calculate additional heating source due to melt

        Returns:
        --------
            Underworld function

        """

        ratio = material.latentHeatFusion / material.capacity

        if not ratio.dimensionless:
            raise ValueError("Unit Error in either Latent Heat Fusion or Capacity (Material: "+material.name)
        ratio = ratio.magnitude

        dF = (self._get_melt_fraction() - self.meltField) / self._dt
        return (ratio * dF) * self.temperature

    @property
    def _lambdaFn(self):
        """ Initialize compressibility """

        materialMap = {}
        if any([material.compressibility for material in self.materials]):
            for material in self.materials:
                if material.compressibility:
                    materialMap[material.index] = material.compressibility

            return uw.function.branching.map(fn_key=self.materialField,
                                             mapping=materialMap,
                                             fn_default=0.0)
        return

    def checkpoint(self, variables=None, checkpointID=None):
        """ Do a checkpoint (Save fields)

        Parameters:
        -----------
            variables:
                list of fields/variables to save
            checkpointID:
                checkpoint ID.

        """
        if not variables:
            variables = rcParams["default.outputs"]

        if not checkpointID:
            checkpointID = self.checkpointID

        for variable in variables:
            if variable == "temperature" and not self.temperature:
                continue
            self._save_field(variable, checkpointID)

        if checkpointID > 2:
            for field in rcParams["swarm.variables"]:
                try:
                    os.remove(os.path.join(self.outputDir,
                                           field + "-%s" % checkpointID - 2))
                except:
                    pass

    def output_glucifer_figures(self, step):
        """ Output glucifer Figures to the gldb store """

        import glucifer
        GluciferStore = glucifer.Store(os.path.join(self.outputDir,
                                       "glucifer"))
        GluciferStore.step = step

        pressure = self.plot.pressureField(store=GluciferStore, show=False)
        pressure.save()

        temperature = self.plot.temperature(store=GluciferStore, show=False)
        temperature.save()

        velocity = self.plot.velocityField(store=GluciferStore, show=False)
        velocity.save()

        strainrate = self.plot.strainRate(store=GluciferStore, show=False)
        strainrate.save()

        material = self.plot.material(projected=True, store=GluciferStore,
                                      show=False)
        material.save()

        strain = self.plot.plastic_strain(projected=True, store=GluciferStore,
                                          show=False)
        strain.save()

        density = self.plot.density(projected=True, store=GluciferStore,
                                    show=False)
        density.save()

        viscosity = self.plot.viscosity(projected=True, store=GluciferStore,
                                        show=False)
        viscosity.save()

    def add_visugrid(self, elementRes, minCoord=None, maxCoord=None):
        """ Add a tracking grid to the Model

        This is essentially a lagrangian grid that deforms with the materials.

        Parameters:
        -----------
            elementRes:
                Grid resolution in number of elements along each axis (x, y, z).
            minCoord:
                Minimum coordinate for each axis.
                Used to define the extent of the grid,
            maxCoord:
                Maximum coordinate for each axis.
                Used to define the extent of the grid,
        """

        if not maxCoord:
            maxCoord = self.maxCoord

        if not minCoord:
            minCoord = self.minCoord

        self._visugrid = Visugrid(self, elementRes, minCoord, maxCoord,
                                  self.velocityField)

    def _save_field(self, field, checkpointID, units=None):
        """ Save Field """

        if field in rcParams["mesh.variables"]:
            if not units:
                units = rcParams[field + ".SIunits"]
            mH = self.mesh.save(os.path.join(self.outputDir, "mesh.h5"),
                                units=u.kilometers)
            file_prefix = os.path.join(self.outputDir,
                                       field + '-%s' % checkpointID)
            obj = getattr(self, field)
            handle = obj.save('%s.h5' % file_prefix, units=units)
            obj.xdmf('%s.xdmf' % file_prefix, handle, field, mH, 'mesh',
                     modeltime=self.time.magnitude)

        elif field in rcParams["swarm.variables"]:
            if not units:
                units = rcParams[field + ".SIunits"]
            sH = self.swarm.save(os.path.join(self.outputDir,
                                 'swarm-%s.h5' % checkpointID),
                                 units=u.kilometers)
            file_prefix = os.path.join(self.outputDir,
                                       field + '-%s' % checkpointID)
            obj = getattr(self, field)
            handle = obj.save('%s.h5' % file_prefix)
            obj.xdmf('%s.xdmf' % file_prefix, handle, field, sH, 'swarm',
                     modeltime=self.time.magnitude)
        else:
            raise ValueError(field, ' is not a valid variable name \n')

    def save(self):
        with open("Model.json", "w") as f:
            json.dump(self, f, sort_keys=True, indent=4, cls=_ModelEncoder)


def load_model(filename):
    with open(filename, "r") as f:
        model = json.load(f)

    for key, val in _attributes_to_save.iteritems():
        model[key] = _attributes_to_save[key](model[key])

    return Model(**model)


_html_global = OrderedDict()
_html_global["Number of Elements"] = "elementRes"
_html_global["length"] = "length"
_html_global["width"] = "width"
_html_global["height"] = "height"


def _model_html_repr(Model):
    header = "<table>"
    footer = "</table>"
    html = ""
    for key, val in _html_global.iteritems():
        value = Model.__dict__.get(val)
        html += "<tr><td>{0}</td><td>{1}</td></tr>".format(key, value)

    return header + html + footer      

