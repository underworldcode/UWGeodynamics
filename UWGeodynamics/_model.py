from __future__ import print_function
import os
import json
import json_encoder
from collections import Iterable
from collections import OrderedDict
import numpy as np
import underworld as uw
import underworld.function as fn
import UWGeodynamics.shapes as shapes
import UWGeodynamics.surfaceProcesses as surfaceProcesses
from . import scaling_coefficients
from . import rcParams, uwgeodynamics_fname
from .scaling import Dimensionalize
from .scaling import nonDimensionalize as nd
from .scaling import UnitRegistry as u
from .lithopress import LithostaticPressure
from ._utils import PressureSmoother, PassiveTracers
from ._rheology import ViscosityLimiter, StressLimiter
from ._material import Material
from ._plots import Plots
from ._visugrid import Visugrid
from ._velocity_boundaries import VelocityBCs
from ._thermal_boundaries import TemperatureBCs
from ._mesh_advector import _mesh_advector
from ._frictional_boundary import FrictionBoundaries
from ._mesh import FeMesh_Cartesian
from ._swarm import Swarm
from ._meshvariable import MeshVariable
from ._swarmvariable import SwarmVariable
from scipy import interpolate
from six import string_types, iteritems


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
                 name=None, gravity=None, periodic=None, elementType=None,
                 temperatureBCs=None, velocityBCs=None, materials=list(),
                 outputDir=None, frictionalBCs=None, surfaceProcesses=None,
                 isostasy=None, visugrid=None, advector=None):

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

        if not elementType:
            self.elementType = rcParams["element.type"]
        else:
            self.elementType = elementType

        self.elementRes = elementRes

        if outputDir:
            self.outputDir = outputDir
        else:
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
        self.mesh = FeMesh_Cartesian(elementType=self.elementType,
                                     elementRes=self.elementRes,
                                     minCoord=minCoord,
                                     maxCoord=maxCoord,
                                     periodic=self.periodic)

        self.mesh_fields = {}
        self.submesh_fields = {}
        self.swarm_fields = {}

        # Add common mesh variables
        # Temperature Field is initialised to None
        self.temperature = None
        self.add_submesh_field("pressureField", nodeDofCount=1)
        self.add_mesh_field("velocityField", nodeDofCount=self.mesh.dim)
        self.add_mesh_field("tractionField", nodeDofCount=self.mesh.dim)
        self.add_submesh_field("_strainRateField", nodeDofCount=1)

        # symmetric component of the gradient of the flow velocityField.
        self.strainRate = fn.tensor.symmetric(self.velocityField.fn_gradient)
        self._strainRate_2ndInvariant = fn.tensor.second_invariant(
            self.strainRate
        )

        # Create the material swarm
        self.swarm = Swarm(mesh=self.mesh, particleEscape=True)
        self._swarmLayout = uw.swarm.layouts.PerCellSpaceFillerLayout(
            swarm=self.swarm,
            particlesPerCell=rcParams["swarm.particles.per.cell"])

        self.swarm.populate_using_layout(layout=self._swarmLayout)

        # timing and checkpointing
        self.checkpointID = 0
        self.time = 0.0 * u.megayears
        self.step = 0
        self._dt = None

        # viscosity limiter
        self.minViscosity = rcParams["minimum.viscosity"]
        self.maxViscosity = rcParams["maximum.viscosity"]

        self._defaultMaterial = self.index
        self.materials = materials
        self.materials.append(self)

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
        self.velocityBCs = velocityBCs
        self.temperatureBCs = temperatureBCs
        self.frictionalBCs = frictionalBCs
        self._isostasy = isostasy
        self.surfaceProcesses = surfaceProcesses

        self.pressSmoother = PressureSmoother(self.mesh, self.pressureField)

        # Passive Tracers
        self.passive_tracers = {}

        # Visualisation
        self.plot = Plots(self)
        self._visugrid = visugrid

        # Mesh advector
        self._advector = None

        # Initialise remaining attributes
        self._default_strain_rate = 1e-15 / u.second
        self._solution_exist = fn.misc.constant(False)
        self._isYielding = None
        self._temperatureDot = None
        self._temperature = None
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

        self.population_control = uw.swarm.PopulationControl(
            self.swarm,
            aggressive=rcParams["popcontrol.aggressive"],
            splitThreshold=rcParams["popcontrol.split.threshold"],
            maxSplits=rcParams["popcontrol.max.splits"],
            particlesPerCell=rcParams["popcontrol.particles.per.cell"])

        # Add Common Swarm Variables
        self.add_swarm_field("materialField", dataType="int", count=1,
                             init_value=self.index)
        self.add_swarm_field("plasticStrain", dataType="double", count=1)
        self.add_swarm_field("_viscosityField", dataType="double", count=1)
        self.add_swarm_field("_densityField", dataType="double", count=1)
        self.add_swarm_field("meltField", dataType="double", count=1)
        self.add_swarm_field("timeField", dataType="double", count=1)
        self.timeField.data[...] = 0.0

        if self.mesh.dim == 3:
            stress_dim = 6
        else:
            stress_dim = 3

        self.add_swarm_field("_previousStressField", dataType="double",
                             count=stress_dim)
        self.add_swarm_field("_stressField", dataType="double",
                             count=stress_dim)

        #self._add_surface_tracers()

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
        self._outputDir = value

    def restart(self, step=None, restartDir=None, badlands_prefix="outbdls",
                badlands_step=None):
        """ Restart a Model

        parameters:
        -----------
            restartDir: (string)
                Directory that contains the outputs of the model
                you want to restart
            step: (int)
                Step from which you want to restart the model.

        """
        if step and not isinstance(step, int):
            raise ValueError("step must be an int")

        if restartDir and not isinstance(restartDir, str):
            raise ValueError("restartDir must be a path to a folder")

        if not restartDir:
            restartDir = self.outputDir

        if not os.path.exists(restartDir):
            return

        if step is None:
            indices = [int(os.path.splitext(filename)[0].split("-")[-1])
                    for filename in os.listdir(restartDir) if "-" in
                    filename]

            if indices:
                step = max(indices) - 1

        if not step or step < 1:
            return

        # Get time from swarm-%.h5 file
        import h5py
        with h5py.File(os.path.join(restartDir, "swarm-%s.h5" % step), "r") as h5f:
            self.time = u.Quantity(h5f.attrs.get("time"))

        if uw.rank() == 0:
            print(80*"="+"\n")
            print("Restarting Model from Step {0} at Time = {1}\n".format(step, self.time))
            print(80*"="+"\n")

        self.checkpointID = step
        self.mesh.load(os.path.join(restartDir, "mesh.h5"))
        self.swarm = Swarm(mesh=self.mesh, particleEscape=True)
        self.swarm.load(os.path.join(restartDir, 'swarm-%s.h5' % step))
        self._initialize()

        # Reload all the restart fields
        for field in rcParams["restart.fields"]:
            "Temporary !!!!"
            if field == "temperature":
                continue
            obj = getattr(self, field)
            path = os.path.join(restartDir, field + "-%s.h5" % step)
            if uw.rank() == 0:
                print("Reloading field {0} from {1}".format(field, path))
            obj.load(str(path))

        # Temperature is a special case...
        path = os.path.join(restartDir, "temperature" + "-%s.h5" % step)
        if os.path.exists(path):
            if not self.temperature:
                self.temperature = MeshVariable(mesh=self.mesh,
                                                nodeDofCount=1)
                self._temperatureDot = MeshVariable(mesh=self.mesh,
                                                    nodeDofCount=1)
                self._heatFlux = MeshVariable(mesh=self.mesh,
                                              nodeDofCount=1)
                self._temperatureDot.data[...] = 0.
                self._heatFlux.data[...] = 0.
            obj = getattr(self, "temperature")
            if uw.rank() == 0:
                print("Reloading field {0} from {1}".format("temperature", path))
            obj.load(str(path))

        # Reload Passive Tracers
        for (key, tracer) in iteritems(self.passive_tracers):

            if uw.rank() == 0:
                print("Reloading {0} passive tracers".format(tracer.name))

            fname = tracer.name + '-%s.h5' % step
            fpath = os.path.join(restartDir, fname)
            with h5py.File(fpath, "r") as h5f:
                vertices = h5f["data"].value * u.Quantity(h5f.attrs["units"])
                vertices = [vertices[:, dim] for dim in range(self.mesh.dim)]
                obj = PassiveTracers(self.mesh,
                                     self.velocityField,
                                     tracer.name,
                                     vertices=vertices,
                                     particleEscape=tracer.particleEscape)

            attr_name = tracer.name.lower() + "_tracers"
            setattr(self, attr_name, obj)
            self.passive_tracers[key] = obj

        # Restart Badlands if we are running a coupled model
        if isinstance(self.surfaceProcesses, surfaceProcesses.Badlands):
            badlands_model = self.surfaceProcesses
            restartFolder = badlands_model.restartFolder
            restartStep = badlands_model.restartStep

            # Parse xmf for the last timestep time
            import xml.etree.ElementTree as etree
            xmf = restartFolder+"/xmf/tin.time"+str(restartStep)+".xmf"
            tree = etree.parse(xmf)
            root = tree.getroot()
            badlands_time = float(root[0][0][0].attrib["Value"])
            uw_time = self.time.to(u.years).magnitude

            if np.abs(badlands_time - uw_time) > 1:
                raise ValueError("""Time in Underworld and Badlands outputs
                                 differs:\n
                                 Badlands: {0}\n
                                 Underworld: {1}""".format(badlands_time,
                                                           uw_time))

            airIndex = badlands_model.airIndex
            sedimentIndex = badlands_model.sedimentIndex
            XML = badlands_model.XML
            resolution = badlands_model.resolution
            checkpoint_interval = badlands_model.checkpoint_interval

            self.surfaceProcesses = surfaceProcesses.Badlands(
                airIndex, sedimentIndex,
                XML, resolution,
                checkpoint_interval,
                restartFolder=restartFolder,
                restartStep=restartStep
                )

        return

    @property
    def projMaterialField(self):
        """ Material field projected on the mesh """
        self._materialFieldProjector.solve()
        self._projMaterialField.data[:] = np.rint(
            self._projMaterialField.data[:]
        )
        return self._projMaterialField

    @property
    def projPlasticStrain(self):
        """ Plastic Strain Field projected on the mesh """
        self._plasticStrainProjector.solve()
        return self._projPlasticStrain

    @property
    def projTimeField(self):
        """ Plastic Strain Field projected on the mesh """
        self._timeFieldProjector.solve()
        return self._projTimeField

    @property
    def strainRateField(self):
        """ Strain Rate Field """
        self._strainRateField.data[:] = self._strainRate_2ndInvariant.evaluate(
            self.mesh.subMesh)
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
    def projStressField(self):
        """ Stress Tensor on mesh """
        self._stressField.data[...] = self._stressFn.evaluate(self.swarm)
        self._stressFieldProjector.solve()
        return self._projStressField

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
        if value:
            self._surfaceProcesses.timeField = self.timeField
        if isinstance(value, surfaceProcesses.Badlands):
            self._surfaceProcesses.Model = self

    def set_temperatureBCs(self, left=None, right=None,
                           top=None, bottom=None,
                           front=None, back=None,
                           indexSets=None, materials=None,
                           bottom_material=None, top_material=None,
                           left_material=None, right_material=None,
                           back_material=None, front_material=None):

        """ Set Model thermal conditions

        Parameters:
            left:
                Temperature or flux along the left wall.
                Flux must be a vector (Fx, Fy, [Fz])
                Default is 'None'
            right:
                Temperature or flux along the right wall.
                Flux must be a vector (Fx, Fy, [Fz])
                Default is 'None'
            top:
                Temperature or flux along the top wall.
                Flux must be a vector (Fx, Fy, [Fz])
                Default is 'None'
            bottom:
                Temperature or flux along the bottom wall.
                Flux must be a vector (Fx, Fy, [Fz])
                Default is 'None'
            front:
                Temperature or flux along the front wall.
                Flux must be a vector (Fx, Fy, [Fz])
                Default is 'None'
            back:
                Temperature or flux along the back wall.
                Flux must be a vector (Fx, Fy, [Fz])
                Default is 'None'
            indexSets: (set, temperature)
                underworld mesh index set and associate temperature
            materials:
                list of materials for which temperature need to be
                fixed. [(material, temperature)]
            bottom_material:
                Specify Material that lays at the bottom of the model
                Default is 'None' (required if a flux is defined at the
                base).
            top_material:
                Specify Material that lays at the top of the model
                Default is 'None' (required if a flux is defined at the top).
            left_material:
                Specify Material that lays at the left of the model
                Default is 'None' (required if a flux is defined at the left).
            right_material:
                Specify Material that lays at the right of the model
                Default is 'None' (required if a flux is defined at the right).
            back_material:
                Specify Material that lays at the back of the model
                Default is 'None' (required if a flux is defined at the back).
            front_material:
                Specify Material that lays at the front of the model
                Default is 'None' (required if a flux is defined at the front).

        """

        if not self.temperature:
            self.temperature = MeshVariable(mesh=self.mesh,
                                            nodeDofCount=1)
            self._temperatureDot = MeshVariable(mesh=self.mesh,
                                                nodeDofCount=1)
            self._heatFlux = MeshVariable(mesh=self.mesh,
                                          nodeDofCount=1)
            self._temperatureDot.data[...] = 0.
            self._heatFlux.data[...] = 0.

        self.temperatureBCs = TemperatureBCs(self, left=left, right=right,
                                             top=top, bottom=bottom,
                                             back=back, front=front,
                                             indexSets=indexSets,
                                             materials=materials,
                                             bottom_material=bottom_material,
                                             top_material=top_material,
                                             left_material=left_material,
                                             right_material=right_material,
                                             back_material=back_material,
                                             front_material=front_material)
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
            if material.diffusivity:
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
                    (self._densityFn * nd(material.capacity)))
            else:
                HeatProdMap[material.index] = 0.

        self.HeatProdFn = fn.branching.map(fn_key=self.materialField,
                                           mapping=HeatProdMap)

        obj = uw.systems.AdvectionDiffusion(
            self.temperature,
            self._temperatureDot,
            velocityField=self.velocityField,
            fn_diffusivity=self.DiffusivityFn,
            fn_sourceTerm=self.HeatProdFn,
            conditions=self._temperatureBCs
        )
        return obj

    def _stokes(self):
        """ Stokes solver """

        gravity = tuple([nd(val) for val in self.gravity])
        self._buoyancyFn = self._densityFn * gravity

        if any([material.viscosity for material in self.materials]):

            stokes_object = uw.systems.Stokes(
                velocityField=self.velocityField,
                pressureField=self.pressureField,
                conditions=self._velocityBCs,
                fn_viscosity=self._viscosityFn,
                fn_bodyforce=self._buoyancyFn,
                fn_stresshistory=self._elastic_stressFn,
                fn_one_on_lambda=self._lambdaFn)
                #useEquationResidual=rcParams["useEquationResidual"])

            solver = uw.systems.Solver(stokes_object)
            solver.set_inner_method(rcParams["solver"])

            if rcParams["penalty"]:
                solver.set_penalty(rcParams["penalty"])

            if rcParams["mg.levels"]:
                solver.options.mg.levels = rcParams["mg.levels"]

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

    set_mechanicalBCs = set_velocityBCs

    @property
    def _velocityBCs(self):
        """ Retrieve kinematic boundary conditions """
        if not self.velocityBCs:
            raise ValueError("Set Boundary Conditions")
        return self.velocityBCs.get_conditions()

    def add_material(self, material=None, shape=None, name="unknown",
                     reset=False, fill=True):
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

        if fill:
            self._fill_model()

        return mat

    def _fill_model(self):
        """ Initialize the Material Field from the list of materials"""

        conditions = [(obj.shape.fn, obj.index)
                      for obj in self.materials if obj.shape is not None]

        conditions.append((True, self._defaultMaterial))
        vals = fn.branching.conditional(conditions).evaluate(self.swarm)
        self.materialField.data[:] = vals

    def add_swarm_field(self, name, dataType="double", count=1,
                        init_value=0., **kwargs):
        newField = self.swarm.add_variable(dataType, count, **kwargs)
        setattr(self, name, newField)
        newField.data[...] = init_value
        self.swarm_fields[name] = newField

        # Create mesh variable for projection
        if name.startswith("_"):
            proj_name = "_proj" + name[1].upper() + name[2:]
        else:
            proj_name = "_proj" + name[0].upper() + name[1:]
        projected = self.add_mesh_field(proj_name, nodeDofCount=count,
                                        dataType="double")
        # Create a projector
        if name.startswith("_"):
            projector_name = name + "Projector"
        else:
            projector_name = "_" + name + "Projector"
        projector = uw.utils.MeshVariable_Projection(projected,
                                                     newField, type=0)
        setattr(self, projector_name, projector)

        return newField

    def add_mesh_field(self, name, nodeDofCount,
                       dataType="double", init_value=0., **kwargs):
        newField = self.mesh.add_variable(nodeDofCount, dataType, **kwargs)
        setattr(self, name, newField)
        newField.data[...] = init_value
        self.mesh_fields[name] = newField
        return newField

    def add_submesh_field(self, name, nodeDofCount,
                          dataType="double", init_value=0., **kwargs):
        newField = MeshVariable(self.mesh.subMesh, nodeDofCount,
                                dataType, **kwargs)
        setattr(self, name, newField)
        newField.data[...] = init_value
        self.submesh_fields[name] = newField

    @property
    def _densityFn(self):
        densityMap = {}
        for material in self.materials:

            if self.temperature:
                dens_handler = material.density
                dens_handler.temperatureField = self.temperature
                dens_handler.pressureField = self.pressureField
                densityMap[material.index] = dens_handler.effective_density()
            else:
                dens_handler = material.density
                densityMap[material.index] = nd(dens_handler.reference_density)

            if material.meltExpansion:
                fact = material.meltExpansion * self.meltField
                densityMap[material.index] = (densityMap[material.index] *
                                              (1.0 - fact))

        return fn.branching.map(fn_key=self.materialField, mapping=densityMap)

    def set_frictional_boundary(self, right=None, left=None,
                                top=None, bottom=None, front=None,
                                back=None, thickness=2):
        """ Set Frictional Boundary conditions

        Frictional boundaries are implemented as a thin layer of frictional
        material along the chosen boundaries.

        Parameters:
        -----------
        Returns:
        --------
            Underworld mesh variable that maps the boundaries.
            (Values are set to 1 if the mesh element applies
            a frictional condition, 0 otherwise).
        """

        self.frictionalBCs = FrictionBoundaries(self, rightFriction=right,
                                                leftFriction=left,
                                                topFriction=top,
                                                bottomFriction=bottom,
                                                frontFriction=front,
                                                backFriction=back,
                                                thickness=thickness)

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

                if not material.minViscosity:
                    minViscosity = self.minViscosity
                else:
                    minViscosity = material.minViscosity

                if not material.maxViscosity:
                    maxViscosity = self.maxViscosity
                else:
                    maxViscosity = material.maxViscosity

                ViscosityHandler._viscosity_limiter = ViscosityLimiter(minViscosity,
                                                                       maxViscosity)

                ViscosityMap[material.index] = ViscosityHandler.muEff

        # Elasticity
        for material in self.materials:
            if material.elasticity:

                # We are dealing with viscous material
                # Raise an error is no viscosity has been defined.
                if not material.viscosity:
                    raise ValueError("""Viscosity undefined for
                                     {0}""".format(material.name))
                ElasticityHandler = material.elasticity
                ElasticityHandler.viscosity = ViscosityMap[material.index]
                if not material.minViscosity:
                    minViscosity = self.minViscosity
                else:
                    minViscosity = material.minViscosity

                if not material.maxViscosity:
                    maxViscosity = self.maxViscosity
                else:
                    maxViscosity = material.maxViscosity

                ElasticityHandler._viscosity_limiter = ViscosityLimiter(minViscosity,
                                                                        maxViscosity)
                ViscosityMap[material.index] = ElasticityHandler.muEff

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

                if material.stressLimiter:
                    stressLimiter = StressLimiter(material.stressLimiter)
                    yieldStress = stressLimiter.apply(yieldStress)
                elif self.stressLimiter:
                    stressLimiter = StressLimiter(self.stressLimiter)
                    yieldStress = stressLimiter.apply(yieldStress)

                if material.elasticity:
                    ElasticityHandler = material.elasticity
                    mu = nd(ElasticityHandler.shear_modulus)
                    dt_e = nd(ElasticityHandler.observation_time)
                    strainRate = fn.tensor.symmetric(self.velocityField.fn_gradient)
                    D_eff = strainRate + 0.5 * self._previousStressField / (mu * dt_e)
                    SRInv = fn.tensor.second_invariant(D_eff)
                else:
                    SRInv = self._strainRate_2ndInvariant

                eij = fn.branching.conditional(
                    [(SRInv < 1e-20, 1e-20),
                     (True, SRInv)])

                muEff = 0.5 * yieldStress / eij
                if not material.minViscosity:
                    minViscosity = self.minViscosity
                else:
                    minViscosity = material.minViscosity

                if not material.maxViscosity:
                    maxViscosity = self.maxViscosity
                else:
                    maxViscosity = material.maxViscosity

                viscosity_limiter = ViscosityLimiter(minViscosity,
                                                     maxViscosity)
                muEff = viscosity_limiter.apply(muEff)
                PlasticityMap[material.index] = muEff

            if self.frictionalBCs is not None:
                from copy import copy

                # Only affect plastic materials
                if material.plasticity:

                    YieldHandler = copy(material.plasticity)
                    YieldHandler.frictionCoefficient = self.frictionalBCs.friction
                    YieldHandler.frictionAfterSoftening = self.frictionalBCs.friction
                    YieldHandler.pressureField = self.pressureField
                    YieldHandler.plasticStrain = self.plasticStrain

                    if self.mesh.dim == 2:
                        yieldStress = YieldHandler._get_yieldStress2D()

                    if self.mesh.dim == 3:
                        yieldStress = YieldHandler._get_yieldStress3D()

                    eij = fn.branching.conditional(
                        [(self._strainRate_2ndInvariant <= 1e-20, 1e-20),
                         (True, self._strainRate_2ndInvariant)])

                    muEff = 0.5 * yieldStress / eij

                    if not material.minViscosity:
                        minViscosity = self.minViscosity
                    else:
                        minViscosity = material.minViscosity

                    if not material.maxViscosity:
                        maxViscosity = self.maxViscosity
                    else:
                        maxViscosity = material.maxViscosity

                    viscosity_limiter = ViscosityLimiter(minViscosity,
                                                         maxViscosity)
                    muEff = viscosity_limiter.apply(muEff)

                    conditions = [(self.frictionalBCs._mask > 0.0, muEff),
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
        if self._solution_exist.value:
            self._isYielding = (fn.branching.conditional(yieldConditions) *
                                self._strainRate_2ndInvariant)
        else:
            self._isYielding = fn.misc.constant(0.0)

        return viscosityFn

    @property
    def _stressFn(self):
        """ Calculate Stress """
        stressMap = {}
        for material in self.materials:
            # Calculate viscous stress
            viscousStressFn = self._viscous_stressFn()
            elasticStressFn = self._elastic_stressFn
            stressMap[material.index] = viscousStressFn + elasticStressFn

        return fn.branching.map(fn_key=self.materialField,
                                mapping=stressMap)

    def _viscous_stressFn(self):
        return 2. * self._viscosityFn * self.strainRate

    @property
    def _elastic_stressFn(self):
        if any([material.elasticity for material in self.materials]):
            stressMap = {}
            for material in self.materials:
                if material.elasticity:
                    ElasticityHandler = material.elasticity
                    ElasticityHandler.viscosity = self._viscosityFn
                    ElasticityHandler.previousStress = self._previousStressField
                    elasticStressFn = ElasticityHandler.elastic_stress
                else:
                    # Set elastic stress to zero if material
                    # has no elasticity
                    elasticStressFn = [0.0] * 3 if self.mesh.dim == 2 else [0.0] * 6
                stressMap[material.index] = elasticStressFn
            return fn.branching.map(fn_key=self.materialField,
                                    mapping=stressMap)
        else:

            elasticStressFn = [0.0] * 3 if self.mesh.dim == 2 else [0.0] * 6
            return elasticStressFn

    def _update_stress_history(self, dt):
        dt_e = []
        for material in self.materials:
            if material.elasticity:
                dt_e.append(nd(material.elasticity.observation_time))
        dt_e = np.array(dt_e).min()
        phi = dt / dt_e
        veStressFn_data = self._stressFn.evaluate(self.swarm)
        self._previousStressField.data[:] *= (1. - phi)
        self._previousStressField.data[:] += phi * veStressFn_data[:]

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
                if material.diffusivity:
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
                        (self._densityFn *
                         nd(material.capacity))
                    )

                else:
                    HeatProdMap[material.index] = 0.

                # Melt heating
                if material.latentHeatFusion and self._dt:
                    dynamicHeating = self._get_dynamic_heating(material)
                    HeatProdMap[material.index] += dynamicHeating

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
        area, = surfaceArea.evaluate()
        p0, = surfacepressureFieldIntegral.evaluate()
        offset = p0 / area
        self.pressureField.data[:] -= offset
        self.pressSmoother.smooth()

        for material in self.materials:
            if material.viscosity:
                material.viscosity.firstIter.value = False

        self._solution_exist.value = True

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
        mat = self.projMaterialField.evaluate(self.mesh.data[0:self.mesh.nodesLocal])
        # Round to closest integer
        mat = np.rint(mat)
        # convert array as integer
        mat = mat.astype("int")[:, 0]
        # Get nodes corresponding to material
        nodes = np.arange(0, self.mesh.nodesLocal)[mat == material.index]
        return uw.mesh.FeMesh_IndexSet(self.mesh, topologicalIndex=0,
                                       size=self.mesh.nodesGlobal,
                                       fromObject=nodes)

    def solve(self):
        """ Solve Stokes """

        if self.step == 0:
            tol = rcParams["initial.nonlinear.tolerance"]
            minIterations = rcParams["initial.nonlinear.min.iterations"]
            maxIterations = rcParams["initial.nonlinear.max.iterations"]
        else:
            tol = rcParams["nonlinear.tolerance"]
            minIterations = rcParams["nonlinear.min.iterations"]
            maxIterations = rcParams["nonlinear.max.iterations"]

        if uw.__version__.split(".")[1] > 5:
            # The minimum number of iteration only works with version 2.6
            # 2.6 is still in development...
            self._stokes().solve(
                nonLinearIterate=True,
                nonLinearMinIterations=minIterations,
                nonLinearMaxIterations=maxIterations,
                callback_post_solve=self._calibrate_pressureField,
                nonLinearTolerance=tol)
        else:
            self._stokes().solve(
                nonLinearIterate=True,
                nonLinearMaxIterations=maxIterations,
                callback_post_solve=self._calibrate_pressureField,
                nonLinearTolerance=tol)

        self._solution_exist.value = True

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

    def run_for(self, duration=None, checkpoint_interval=None,
                timeCheckpoints=None, swarm_checkpoint=None, dt=None,
                glucifer_outputs=False):
        """ Run the Model

        Parameters:
        -----------

            duration: duration of the Model in Time units or nb of steps
            checkpoint_interval: checkpoint interval
            timeCheckpoints: Additional checkpoints
            dt: force time interval.
            glucifer_outputs: output glucifer figures [False]
        """

        if uw.rank()==0 and not os.path.exists(self.outputDir):
            os.mkdir(self.outputDir)
        uw.barrier()

        step = self.step
        time = nd(self.time)
        units = duration.units
        duration = time + nd(duration)
        user_dt = nd(dt)

        next_checkpoint = None

        if timeCheckpoints:
            timeCheckpoints = [nd(val) for val in timeCheckpoints]

        if checkpoint_interval:
            next_checkpoint = time + nd(checkpoint_interval)

        if checkpoint_interval or timeCheckpoints:
            # Check point initial set up
            self.checkpoint()
            if glucifer_outputs:
                self.output_glucifer_figures(self.checkpointID)

        # Save model to json
        # Comment as it does not work in parallel
        # self.save()

        while time < duration:

            self.preSolveHook()

            self.solve()

            # Whats the longest we can run before reaching the end
            # of the model or a checkpoint?
            # Need to generalize that
            dt = rcParams["CFL"] * self.swarm_advector.get_max_dt()

            if self.temperature:
                dt = min(dt, self._advdiffSystem.get_max_dt())
                dt *= rcParams["CFL"]

            if checkpoint_interval:
                dt = min(dt, next_checkpoint - time)

            self._dt = min(dt, duration - time)

            if user_dt:
                self._dt = min(self._dt, user_dt)

            dte = []
            for material in self.materials:
                if material.elasticity:
                    dte.append(nd(material.elasticity.observation_time))

            if dte:
                dte = np.array(dte).min()
                # Cap dt for observation time, dte / 3.
                if dte and self._dt > (dte / 3.):
                    self._dt = dte / 3.

            uw.barrier()

            self._update()

            self.step += 1
            self.time += Dimensionalize(self._dt, units)
            time += self._dt

            if time == next_checkpoint:
                self.checkpointID += 1
                self.checkpoint()
                if glucifer_outputs:
                    self.output_glucifer_figures(self.checkpointID)
                next_checkpoint += nd(checkpoint_interval)

            if checkpoint_interval or step % 1 == 0:
                if uw.rank() == 0:
                    print("Time: ", str(self.time.to(units)),
                          'dt:', str(Dimensionalize(self._dt, units)))

            self.postSolveHook()

        return 1

    def preSolveHook(self):
        return

    def postSolveHook(self):
        return

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
            self.update_melt_fraction()

        # Solve for temperature
        if self.temperature:
            self._advdiffSystem.integrate(dt)

        if self._advector:
            self.swarm_advector.integrate(dt)
            self._advector.advect_mesh(dt)
        else:
            # Integrate Swarms in time
            self.swarm_advector.integrate(dt, update_owners=True)

        # Update stress
        if any([material.elasticity for material in self.materials]):
            self._update_stress_history(dt)

        if self.passive_tracers:
            for (key, tracers) in iteritems(self.passive_tracers):
                tracers.integrate(dt)

        # Do pop control
        self.population_control.repopulate()
        self.swarm.update_particle_owners()

        if self.surfaceProcesses:
            self.surfaceProcesses.solve(dt)

        # Update Time Field
        self.timeField.data[...] += dt

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

    def add_passive_tracers(self, name, vertices=None,
                            particleEscape=False):
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
        if name in self.passive_tracers.keys():
            print("{0} tracers exists already".format(name))
            return

        tracers = PassiveTracers(self.mesh,
                                 self.velocityField,
                                 name=name,
                                 vertices=vertices,
                                 particleEscape=particleEscape)

        self.passive_tracers[name] = tracers

        setattr(self, name.lower() + "_tracers", tracers)

        return tracers

    def _add_surface_tracers(self):
        """ Add tracers at the surface """

        if self.mesh.dim < 3:
            xcoords = np.linspace(nd(self.minCoord[0]), nd(self.maxCoord[0]), 1000)
            ycoords = 0.
            self.add_passive_tracers(name="Surface", vertices=[xcoords, ycoords])
        else:
            xcoords = np.linspace(nd(self.minCoord[0]), nd(self.maxCoord[0]), 100)
            ycoords = np.linspace(nd(self.minCoord[1]), nd(self.maxCoord[1]), 100)
            xcoords, ycoords = np.meshgrid(xcoords, ycoords)
            xcoords = xcoords.ravel()
            ycoords = ycoords.ravel()
            zcoords = 0.
            self.add_passive_tracers(name="Surface",
                                     vertices=[xcoords, ycoords, zcoords])

        return

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

    def update_melt_fraction(self):
        # Calculate New meltField
        meltFraction = self._get_melt_fraction()
        self.meltField.data[:] = meltFraction.evaluate(self.swarm)

    def _get_dynamic_heating(self, material):
        """ Calculate additional heating source due to melt

        Returns:
        --------
            Underworld function

        """

        ratio = material.latentHeatFusion / material.capacity

        if not ratio.dimensionless:
            raise ValueError("""Unit Error in either Latent Heat Fusion or
                             Capacity (Material: """ + material.name)
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
                    materialMap[material.index] = nd(material.compressibility)

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

        self._save_fields(variables, checkpointID)
        self._save_swarms(variables, checkpointID)

        # Checkpoint passive tracers and associated tracked fields
        if self.passive_tracers:
            for (key, tracers) in iteritems(self.passive_tracers):
                tracers.save(self.outputDir, checkpointID, self.time)

    def profile(self,
                start_point,
                end_point,
                field,
                field_units=u.dimensionless,
                distance_units=u.dimensionless,
                npoints=1000):

        if (field not in rcParams["mesh.variables"] and
           field not in rcParams["swarm.variables"]):
            raise ValueError("""{0} is not a valid field, \n
                                Valid fields are {1} or {2}""".format(
                                    field, rcParams["mesh.variables"],
                                    rcParams["swarm.variables"]))

        start_point = np.array([nd(val) for val in start_point])
        end_point = np.array([nd(val) for val in end_point])

        points = [start_point, end_point]
        x_coords, y_coords = zip(*points)

        x = np.linspace(x_coords[0], x_coords[1], npoints)

        if x_coords[0] == x_coords[1]:
            y = np.linspace(y_coords[0], y_coords[1], npoints)
        else:
            f = interpolate.interp1d(x_coords, y_coords)
            y = f(x)

        dx = np.diff(x)
        dy = np.diff(y)
        distances = np.zeros(x.shape)
        distances[1:] = np.sqrt(dx**2 + dy**2)
        distances = np.cumsum(distances)

        obj = getattr(self, field)

        if isinstance(obj, SwarmVariable):
            obj = getattr(self, "proj" + field[0].upper() + field[1:])

        pts = np.array(zip(x, y))

        values = obj.evaluate(pts)

        distance_fact = Dimensionalize(1.0, distance_units)
        field_fact = Dimensionalize(1.0, field_units)

        distances *= distance_fact.magnitude
        values *= field_fact.magnitude

        return distances, values

    def output_glucifer_figures(self, step):
        """ Output glucifer Figures to the gldb store """

        import glucifer
        GluciferStore = glucifer.Store(os.path.join(self.outputDir,
                                       "glucifer"))
        GluciferStore.step = step

        for field in rcParams["glucifer.outputs"]:
            if getattr(self, field):
                func = eval("self.plot." + field)
                fig = func(store=GluciferStore, show=False)
                fig.save()

    def add_visugrid(self, elementRes, minCoord=None, maxCoord=None):
        """ Add a tracking grid to the Model

        This is essentially a lagrangian grid that deforms with the materials.

        Parameters:
        -----------
            elementRes:
                Grid resolution in number of elements
                along each axis (x, y, z).
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
                try:
                    units = rcParams[field + ".SIunits"]
                except KeyError:
                    units = None

            if self._advector:
                mesh_name = 'mesh-%s' % checkpointID
                mesh_prefix = os.path.join(self.outputDir, mesh_name)
            else:
                mesh_name = 'mesh'
                mesh_prefix = os.path.join(self.outputDir, mesh_name)

            mH = self.mesh.save('%s.h5' % mesh_prefix, units=u.kilometers,
                                time=self.time)
            file_prefix = os.path.join(self.outputDir, field + '-%s' % checkpointID)
            obj = getattr(self, field)
            handle = obj.save('%s.h5' % file_prefix, units=units, time=self.time)
            obj.xdmf('%s.xdmf' % file_prefix, handle, field, mH, mesh_name,
                     modeltime=self.time.magnitude)

        elif field in rcParams["swarm.variables"]:
            if not units:
                try:
                    units = rcParams[field + ".SIunits"]
                except KeyError:
                    units = None

            sH = self.swarm.save(os.path.join(self.outputDir,
                                 'swarm-%s.h5' % checkpointID),
                                 units=u.kilometers, time=self.time)
            file_prefix = os.path.join(self.outputDir,
                                       field + '-%s' % checkpointID)
            obj = getattr(self, field)
            handle = obj.save('%s.h5' % file_prefix, units=units, time=self.time)
            obj.xdmf('%s.xdmf' % file_prefix, handle, field, sH, 'swarm',
                     modeltime=self.time.magnitude)
        else:
            raise ValueError(field, ' is not a valid variable name \n')

    def _save_fields(self, fields, checkpointID):

        if self._advector:
            mesh_name = 'mesh-%s' % checkpointID
            mesh_prefix = os.path.join(self.outputDir, mesh_name)
        else:
            mesh_name = 'mesh'
            mesh_prefix = os.path.join(self.outputDir, mesh_name)

        mH = self.mesh.save('%s.h5' % mesh_prefix, units=u.kilometers,
                            time=self.time)

        filename = "XDMF.fields."+str(checkpointID).zfill(5)+".xmf"
        filename = os.path.join(self.outputDir, filename)

        # First write the XDMF header
        string = uw.utils._xdmfheader()
        string += uw.utils._spacetimeschema(mH, mesh_name, self.time)

        for field in fields:
            if field == "temperature" and not self.temperature:
                continue
            if field in rcParams["mesh.variables"]:
                field = str(field)
                try:
                    units = rcParams[field + ".SIunits"]
                except KeyError:
                    units = None

                # Save the h5 file and write the field schema for
                # each one of the field variables
                obj = getattr(self, field)
                file_prefix = os.path.join(self.outputDir, field + '-%s' % checkpointID)
                handle = obj.save('%s.h5' % file_prefix, units=units,
                                  time=self.time)
                string += uw.utils._fieldschema(handle, field)

        # Write the footer to the xmf
        string += uw.utils._xdmffooter()

        # Write the string to file - only proc 0
        xdmfFH = open(filename, "w")
        xdmfFH.write(string)
        xdmfFH.close()

    def _save_swarms(self, fields, checkpointID):

        swarm_name = 'swarm-%s.h5' % checkpointID

        sH = self.swarm.save(os.path.join(self.outputDir,
                             swarm_name),
                             units=u.kilometers,
                             time=self.time)

        filename = "XDMF.swarms."+str(checkpointID).zfill(5)+".xmf"
        filename = os.path.join(self.outputDir, filename)

        # First write the XDMF header
        string = uw.utils._xdmfheader()
        string += uw.utils._swarmspacetimeschema(sH, swarm_name, self.time)

        for field in fields:
            if field in rcParams["swarm.variables"]:
                field = str(field)
                try:
                    units = rcParams[field + ".SIunits"]
                except KeyError:
                    units = None

                # Save the h5 file and write the field schema for
                # each one of the field variables
                obj = getattr(self, field)
                file_prefix = os.path.join(self.outputDir, field + '-%s' % checkpointID)
                handle = obj.save('%s.h5' % file_prefix, units=units, time=self.time)
                string += uw.utils._swarmvarschema(handle, field)

        # Write the footer to the xmf
        string += uw.utils._xdmffooter()

        # Write the string to file - only proc 0
        xdmfFH = open(filename, "w")
        xdmfFH.write(string)
        xdmfFH.close()

    def save(self, filename=None):
        save_model(self, filename)

    def geometry_from_shapefile(self, filename, units=None):
        from ._utils import MoveImporter
        Importer = MoveImporter(filename, units=units)

        shape_dict = {name: []  for name in Importer.names}

        for polygon in Importer.generator:
            name = polygon["properties"]["Name"]
            vertices = polygon["coordinates"]
            shape = shapes.Polygon(vertices)
            shape_dict[name].append(shape)

        for name, shape in shape_dict.iteritems():
            mshape = shapes.MultiShape(shape)
            self.add_material(name=name,
                              shape=mshape,
                              fill=False)

        self._fill_model()


    def velocity_rms(self):
        vdotv_fn = uw.function.math.dot(self.velocityField, self.velocityField)
        fn_2_integrate = (1., vdotv_fn)
        (v2, vol) = self.mesh.integrate(fn=fn_2_integrate)
        import math
        vrms = math.sqrt(v2/vol)
        os.write(1, "Velocity rms (vrms): {0}".format(vrms))
        return vrms


def save_model(model, filename):
    if not filename:
        filename = model.name + ".json"

    path = os.path.join(model.outputDir, filename)
    with open(path, "w") as f:
        json.dump(model, f, sort_keys=True, indent=4, cls=json_encoder.ObjectEncoder)

def load_model(filename, step=None):
    """ Reload Model from json file """
    global scaling_coefficients

    with open(filename, "r") as f:
        model = json.load(f)

    # Set rcParams
    rcParams = model.pop("rcParams")

    # Set scaling
    scaling = model.pop("scaling")
    for elem in scaling_coefficients:
        decoder = json_encoder.ObjectDecoder()
        scaling_coefficients[elem] = decoder.json_to_model(scaling[elem])

    decoder = json_encoder.ObjectDecoder()
    Model = decoder.json_to_model(model)
    return Model


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
