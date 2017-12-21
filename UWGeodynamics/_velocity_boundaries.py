import UWGeodynamics.scaling as sca
import underworld as uw
import underworld.function as fn
from .LecodeIsostasy import LecodeIsostasy
from .scaling import nonDimensionalize as nd

u = UnitRegistry = sca.UnitRegistry

class VelocityBCs(object):
    """ Class to define the mechanical boundary conditions """


    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, indexSets=None):
        """ Defines mechanical boundary conditions

        The type of conditions is determined through the units used do define
        the parameters:
            * Units of velocity ([length] / [time]) represent a kinematic
            condition (Dirichlet)
            * Units of stress / pressure ([force] / [area]) are set as
            stress condition (Neumann).

        parameters
        ----------

        Model: (UWGeodynamics.Model)
            An UWGeodynamics Model (See UWGeodynamics.Model)
        left:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define mechanical conditions on the left side of the Model.
            Conditions are defined for each Model direction (x, y, [z])
        right:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define mechanical conditions on the right side of the Model.
            Conditions are defined for each Model direction (x, y, [z])
        top:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define mechanical conditions on the top side of the Model.
            Conditions are defined for each Model direction (x, y, [z])
        bottom:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define mechanical conditions on the bottom side of the Model.
            Conditions are defined for each Model direction (x, y, [z])
        indexSets: (list)
            List of node where to apply predefined velocities.

        Only valid for 3D Models:

        front:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define mechanical conditions on the front side of the Model.
            Conditions are defined for each Model direction (x, y, [z])
        back:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define mechanical conditions on the front side of the Model.
            Conditions are defined for each Model direction (x, y, [z])


        examples:
        ---------

        The following example defines a (2x1) meter Underworld model with
        freeslip conditions on all the sides.

        >>> import UWGeodynamics as GEO
        >>> u = GEO.u
        >>> Model = GEO.Model(elementRes=(64, 64),
                              minCoord=(-1. * u.meter, -50. * u.centimeter),
                              maxCoord=(1. * u.meter, 50. * u.centimeter))
        >>> Model.set_velocityBCs(left=[0, None], right=[0,None], top=[None,0],
                                  bottom=[None, 0])

        """

        self.Model = Model
        self.left = left
        self.right = right
        self.top = top
        self.bottom = bottom
        self.back = back
        self.front = front
        self.indexSets = indexSets

    def _isNeumann(self, x):
        """ Returns true if x as units of stress """

        if not isinstance(x, u.Quantity):
            return False
        val = x.to_base_units()
        return val.units == u.kilogram / (u.meter * u.second**2)

    def get_conditions(self):
        """ Get the mechanical boundary conditions

        Returns
        -------

        List of conditions as:
            [<underworld.conditions._conditions.DirichletCondition,
             <underworld.conditions._conditions.NeumannCondition]
        or
            [<underworld.conditions._conditions.DirichletCondition]

        """

        # A lot of that code is redundant...
        # Could be cleaned up a bit...

        Model = self.Model
        dirichletIndices = []
        neumannIndices = []

        for dim in range(Model.mesh.dim):
            dirichletIndices.append(Model.mesh.specialSets["Empty"])
            neumannIndices.append(Model.mesh.specialSets["Empty"])

        if self.left is not None:
            for dim in range(Model.mesh.dim):
                if isinstance(self.left[dim], (list, tuple)):
                    func = fn.branching.conditional(self.left[dim])
                    Model.velocityField.data[Model._left_wall.data, dim] = func.evaluate(Model.mesh.data[Model._left_wall.data])[:,dim]
                    dirichletIndices[dim] += Model._left_wall
                    continue
                if self.left[dim] is not None:
                    if self._isNeumann(self.left[dim]):
                        Model.tractionField.data[Model._left_wall.data, dim] = nd(self.left[dim])
                        neumannIndices[dim] += Model._left_wall
                    else:
                        Model.velocityField.data[Model._left_wall.data, dim] = nd(self.left[dim])
                        dirichletIndices[dim] += Model._left_wall


        if self.right is not None:
            for dim in range(Model.mesh.dim):
                if isinstance(self.right[dim], (list, tuple)):
                    func = fn.branching.conditional(self.right[dim])
                    Model.velocityField.data[Model._right_wall.data, dim] = func.evaluate(Model.mesh.data[Model._right_wall.data])[:,dim]
                    dirichletIndices[dim] += Model._right_wall
                    continue
                if self.right[dim] is not None:
                    if self._isNeumann(self.right[dim]):
                        Model.tractionField.data[Model._right_wall.data, dim] = nd(self.right[dim])
                        neumannIndices[dim] += Model._right_wall
                    else:
                        Model.velocityField.data[Model._right_wall.data, dim] = nd(self.right[dim])
                        dirichletIndices[dim] += Model._right_wall

        if self.top is not None:
            for dim in range(Model.mesh.dim):
                if isinstance(self.top[dim], (list, tuple)):
                    func = fn.branching.conditional(self.top[dim])
                    Model.velocityField.data[Model._top_wall.data, dim] = func.evaluate(Model.mesh.data[Model._top_wall.data])[:,dim]
                    dirichletIndices[dim] += Model._top_wall
                    continue
                if self.top[dim] is not None:
                    if self._isNeumann(self.top[dim]):
                        Model.tractionField.data[Model._top_wall.data, dim] = nd(self.top[dim])
                        neumannIndices[dim] += Model._top_wall
                    else:
                        Model.velocityField.data[Model._top_wall.data, dim] = nd(self.top[dim])
                        dirichletIndices[dim] += Model._top_wall

        if self.bottom is not None:
            if isinstance(self.bottom, LecodeIsostasy):
                Model._isostasy = self.bottom
                Model._isostasy.mesh = Model.mesh
                Model._isostasy.swarm = Model.swarm
                Model._isostasy.velocityField = Model.velocityField
                Model._isostasy.materialIndexField = Model.materialField
                Model._isostasy._densityFn = Model._densityFn
                dirichletIndices[-1] += Model._bottom_wall
            else:
                Model._isostasy = None
                for dim in range(Model.mesh.dim):
                    if isinstance(self.bottom[dim], (list, tuple)):
                        func = fn.branching.conditional(self.bottom[dim])
                        Model.velocityField.data[Model._bottom_wall.data, dim] = func.evaluate(Model.mesh.data[Model._bottom_wall.data])[:,dim]
                        dirichletIndices[dim] += Model._bottom_wall
                        continue
                    if self.bottom[dim] is not None:
                        if self._isNeumann(self.bottom[dim]):
                            Model.tractionField.data[Model._bottom_wall.data, dim] = nd(self.bottom[dim])
                            neumannIndices[dim] += Model._bottom_wall
                        else:
                            Model.velocityField.data[Model._bottom_wall.data, dim] = nd(self.bottom[dim])
                            dirichletIndices[dim] += Model._bottom_wall

        if self.front is not None and Model.mesh.dim > 2:
            for dim in range(Model.mesh.dim):
                if isinstance(self.front[dim], (list, tuple)):
                    func = fn.branching.conditional(self.front[dim])
                    Model.velocityField.data[Model._front_wall.data, dim] = func.evaluate(Model.mesh.data[Model._front_wall.data])[:,dim]
                    dirichletIndices[dim] += Model._front_wall
                    continue
                if self.front[dim] is not None:
                    if self._isNeumann(self.front[dim]):
                        Model.tractionField.data[Model._front_wall.data, dim] = nd(self.front[dim])
                        neumannIndices[dim] += Model._front_wall
                    else:
                        Model.velocityField.data[Model._front_wall.data, dim] = nd(self.front[dim])
                        dirichletIndices[dim] += Model._front_wall

        if self.back is not None and Model.mesh.dim > 2:
            for dim in range(Model.mesh.dim):
                if isinstance(self.back[dim], (list, tuple)):
                    func = fn.branching.conditional(self.back[dim])
                    Model.velocityField.data[Model._back_wall.data, dim] = func.evaluate(Model.mesh.data[Model._back_wall.data])[:,dim]
                    dirichletIndices[dim] += Model._back_wall
                    continue
                if self.back[dim] is not None:
                    if self._isNeumann(self.back[dim]):
                        Model.tractionField.data[Model._back_wall.data, dim] = nd(self.back[dim])
                        neumannIndices[dim] += Model._back_wall
                    else:
                        Model.velocityField.data[Model._back_wall.data, dim] = nd(self.back[dim])
                        dirichletIndices[dim] += Model._back_wall
        
        if self.indexSets:
            for indexSet in self.indexSets:
                for dim in range(Model.mesh.dim):
                    if indexSet[dim] is not None:
                        Model.velocityField.data[indexSet.data, dim] = nd(indexSet[dim])
                        dirichletIndices[dim] += indexSet

        conditions = []

        conditions.append(uw.conditions.DirichletCondition(variable=Model.velocityField,
                                                           indexSetsPerDof=dirichletIndices))

        neumannIndices = tuple([val if val.data.size > 0 else None for val in neumannIndices])

        if neumannIndices != tuple([None for val in range(Model.mesh.dim)]):
            conditions.append(uw.conditions.NeumannCondition(fn_flux=Model.tractionField,
                                                             variable=Model.velocityField,
                                                             indexSetsPerDof=neumannIndices))

        if not conditions:
            raise ValueError("Undefined conditions, please check your condition")

        return conditions
