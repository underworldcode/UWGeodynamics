from __future__ import print_function,  absolute_import
import numpy as np
import underworld as uw
import underworld.function as fn
from .scaling import nonDimensionalize as nd
from .scaling import UnitRegistry as u
from .LecodeIsostasy import LecodeIsostasy
from ._material import Material
from ._utils import Balanced_InflowOutflow
from ._utils import MovingWall
from .shapes import Shape


class BoundaryConditions(object):
    """Base Class for Boundary Condition objects"""

    def __init__(self, Model, field, varfield, left=None, right=None,
                 top=None, bottom=None, front=None, back=None,
                 nodeSets=tuple(), materials=None, order_wall_conditions=None,
                 condition_type="Dirichlet"):

        """ Define boundaries condition

        A condition can be a temperature (float, int, Pint Quantity) or
        an Underworld function which evaluates as a temperature.

        parameters
        ----------

            Model: (UWGeodynamics.Model)
                An UWGeodynamics Model (See UWGeodynamics.Model)

            field: Underworld Mesh Variable
                Field on which the condition applies

            varfield: Underworld Mesh Variable
                Related Field, needed for neumann conditions

            left:
                Define conditions on the left side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            right:
                Define conditions on the right side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            top:
                Define conditions on the top side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            bottom:
                Define conditions on the bottom side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            nodeSets: list of tuples: [(nodes, condition)]
                List of node where to apply predefined condition.

                The nodes can be a list or a numpy array containing the
                local index of the nodes. It can also be an Underworld
                IndexSet. You can also pass an UWGeodynamics shape.

                The condition can be an Underworld Function, A Pint
                Quantity of a scalar.

            materials: list of tuples: [(Material, condition)]
                List of material on which to apply a condition.
                The materials must be UWGeodynamics Material objects.

                The condition can be an Underworld Function, A Pint
                Quantity of a scalar.

            order_wall_conditions: list of str, [left, right, top, bottom,
                front, back]
                Order in which the boundaries are processed.

            condition_type: str
                can be either "Dirichlet" or "Neumann"

            Only valid for 3D Models:

            front:
                Define conditions on the front side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            back:
                Define conditions on the front side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

        """

        self.Model = Model
        self.left = left
        self.right = right
        self.top = top
        self.bottom = bottom
        self.back = back
        self.front = front
        self.nodeSets = nodeSets
        self.materials = materials
        self.field = field
        self.varfield = varfield

        self.condition_type = condition_type

        if self.condition_type not in ["Dirichlet", "Neumann"]:
            raise ValueError("""Condition type can only be 'Dirichlet' or
                             'Neumann'""")

        self._indices = []
        for _ in range(Model.mesh.dim):
            self._indices.append(Model.mesh.specialSets["Empty"])

        self._wall_indexSets = {"bottom": (self.bottom,
                                           self.Model.bottom_wall),
                                "top": (self.top,
                                        self.Model.top_wall),
                                "left": (self.left,
                                         self.Model.left_wall),
                                "right": (self.right,
                                          self.Model.right_wall),
                                "front": (self.front,
                                          self.Model.front_wall),
                                "back": (self.back,
                                         self.Model.back_wall)}
        self.order_wall_conditions = ["bottom", "top", "front", "back",
                                      "left", "right"]

    def __getitem__(self, name):
        return self.__dict__[name]

    def _apply_conditions_nodes(self, condition, nodes):
        """ Apply condition to a set of nodes

        Parameters:
        -----------
            condition:
                condition
            nodes: np.array, list, IndexSet, Shape, Function
                set of nodes
        """

        # Expect a list or tuple of dimension equal to the dof of the
        # field on which the condition is used

        if isinstance(condition, (float, int, u.Quantity)):
            condition = [condition]

        # If nodes is a shape get the corresponding local
        # nodes
        if isinstance(nodes, (Shape, fn.Function)):
            if isinstance(nodes, Shape):
                mask = nodes.fn.evaluate(self.Model.mesh)
            else:
                mask = nodes.evaluate(self.Model.mesh)
            nodes = np.arange(self.Model.mesh.nodesLocal)
            nodes = nodes.astype("int")[mask.ravel()]

        # Convert nodes numpy arrays / list to UW IndexSet
        if isinstance(nodes, (list, np.ndarray)):
            nodes = uw.mesh.FeMesh_IndexSet(
                self.Model.mesh, topologicalIndex=0,
                size=self.Model.mesh.nodesGlobal,
                fromObject=nodes)

        # Check that the domain actually contains some boundary nodes
        if isinstance(condition, (list, tuple)) and nodes.data.size > 0:
            for dim in range(self.field.data.shape[1]):

                # Scalar condition
                if isinstance(condition[dim], (u.Quantity, float, int)):
                    self.field.data[nodes.data, dim] = nd(condition[dim])
                    self._indices[dim] += nodes

                # If it's an Underworld function
                elif isinstance(condition[dim], fn.Function):
                    func = condition[dim]
                    self.field.data[nodes.data, dim] = (
                        func.evaluate(
                            self.Model.mesh.data[nodes.data])[:, dim])
                    self._indices[dim] += nodes
                elif condition[dim] is not None:
                    raise ValueError("""Wrong condition""")

        return

    def get_conditions(self):

        self._indices = []

        for _ in range(self.field.data.shape[1]):
            self._indices.append(self.Model.mesh.specialSets["Empty"])

        for set_ in self.order_wall_conditions:
            (condition, nodes) = self._wall_indexSets[set_]
            if nodes is not None:
                self._apply_conditions_nodes(condition, nodes)

        if self.nodeSets:
            for (nodes, condition) in self.nodeSets:
                self._apply_conditions_nodes(condition, nodes)

        if self.materials:
            for (material, condition) in self.materials:
                nodes = self.Model._get_material_indices(material)
                self._apply_conditions_nodes(condition, nodes)

        if self.condition_type is "Dirichlet":
            return uw.conditions.DirichletCondition(
                variable=self.field, indexSetsPerDof=self._indices)
        elif self.condition_type is "Neumann":
            return uw.conditions.NeumannCondition(
                fn_flux=self.field,
                variable=self.varfield,
                indexSetsPerDof=self._indices)


class VelocityBCs(BoundaryConditions):
    """ Class to define the mechanical boundary conditions """

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, nodeSets=None, materials=None,
                 order_wall_conditions=None):
        """ Defines kinematic boundary conditions

        A condition can be a velocity (float, int, Pint Quantity) or
        an Underworld function which evaluates as a velocity.

        parameters
        ----------

        Model: (UWGeodynamics.Model)
            An UWGeodynamics Model (See UWGeodynamics.Model)

        left:(tuple)
            Define kinematic conditions on the left side of the Model.
            Conditions are defined for each Model direction (x, y, [z])

        right:(tuple)
            Define kinematic conditions on the right side of the Model.
            Conditions are defined for each Model direction (x, y, [z])

        top:(tuple)
            Define kinematic conditions on the top side of the Model.
            Conditions are defined for each Model direction (x, y, [z])

        bottom:(tuple)
            Define kinematic conditions on the bottom side of the Model.
            Conditions are defined for each Model direction (x, y, [z])

        nodeSets: list of tuples: [(nodes, condition)]
            List of node where to apply predefined condition.

            The nodes can be a list or a numpy array containing the
            local index of the nodes. It can also be an Underworld
            IndexSet. You can also pass an UWGeodynamics shape.

            The condition can be an Underworld Function, A Pint
            Quantity of a scalar.

        materials: list of tuples: [(Material, condition)]
            List of material on which to apply a condition.
            The materials must be UWGeodynamics Material objects.

            The condition can be an Underworld Function, A Pint
            Quantity of a scalar.

        order_wall_conditions: list of str, [left, right, top, bottom, front, back]
            Order in which the boundaries are processed.

        Only valid for 3D Models:

        front:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define kinematic conditions on the front side of the Model.
            Conditions are defined for each Model direction (x, y, [z])

        back:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define kinematic conditions on the front side of the Model.
            Conditions are defined for each Model direction (x, y, [z])

        examples:
        ---------

        The following example defines a (2x1) meter Underworld model with
        freeslip conditions on all the sides.

        >>> import UWGeodynamics as GEO
        >>> u = GEO.u
        >>> Model = GEO.Model(elementRes=(64, 64),
        ...                   minCoord=(-1. * u.meter, -50. * u.centimeter),
        ...                   maxCoord=(1. * u.meter, 50. * u.centimeter))
        >>> velocityBCs = Model.set_velocityBCs(left=[0, None],
        ...                                     right=[0,None],
        ...                                     top=[None,0],
        ...                                     bottom=[None, 0])

        """

        super(VelocityBCs, self).__init__(Model, Model.velocityField, None,
                                          left, right, top, bottom,
                                          front, back, nodeSets, materials,
                                          order_wall_conditions, "Dirichlet")

        # Link Moving Walls
        for arg in [self.left, self.right, self.top, self.bottom, self.front,
                    self.back]:
            if isinstance(arg, MovingWall):
                arg.Model = self.Model

    def _apply_conditions_nodes(self, condition, nodes):
        """ Apply condition to a set of nodes

        Parameters:
        -----------
            condition:
                condition
            nodes: np.array, list, IndexSet, Shape, Function
                set of nodes
        """

        # Special case (Bottom LecodeIsostasy)
        if isinstance(condition, LecodeIsostasy):

            # Apply support condition
            self.Model._isostasy = self.bottom
            self.Model._isostasy.mesh = self.Model.mesh
            self.Model._isostasy.swarm = self.Model.swarm
            self.Model._isostasy._mesh_advector = self.Model._advector
            self.Model._isostasy.velocityField = self.Model.velocityField
            self.Model._isostasy.materialIndexField = self.Model.materialField
            self.Model._isostasy._densityFn = self.Model._densityFn
            vertical_walls_conditions = {
                "left": self.left,
                "right": self.right,
                "front": self.front,
                "back": self.back
            }
            self.Model._isostasy.vertical_walls_conditions = (
                vertical_walls_conditions)
            self._indices[-1] += self.Model.bottom_wall
            return

        if isinstance(condition, MovingWall):
            condition.wall = nodes
            set_, axis = condition.get_wall_indices()
            func = condition.velocityFn

            # Intersect of wall nodes and current local domain
            intersect = np.intersect1d(self.Model.mesh.data_nodegId, set_.data)
            ISet = uw.mesh.FeMesh_IndexSet(
                self.Model.mesh, topologicalIndex=0,
                size=self.Model.mesh.nodesGlobal,
                fromObject=intersect)

            for dim in range(self.Model.mesh.dim):
                if ISet.data.size > 0:
                    if (dim == axis):
                        self.field.data[ISet.data, dim] = func.evaluate(ISet)[:, 0]
                        self._indices[dim] += ISet
                    else:
                        self.field.data[ISet.data, dim] = 0.
                        self._indices[dim] += ISet

            return

        # Expect a list or tuple of dimension mesh.dim.
        # Check that the domain actually contains some boundary nodes
        # (nodes is not None)
        if isinstance(condition, (list, tuple)) and nodes.data.size > 0:
            for dim in range(self.Model.mesh.dim):

                # Inflow Outflow
                if isinstance(condition[dim], Balanced_InflowOutflow):
                    obj = condition[dim]
                    obj.ynodes = self.Model.mesh.data[nodes.data, 1]
                    obj._get_side_flow()
                    self.Model.velocityField.data[nodes.data, dim] = (
                        obj._get_side_flow())
                    self._indices[dim] += nodes

        super(VelocityBCs, self)._apply_conditions_nodes(condition, nodes)

        return


class StressBCs(BoundaryConditions):
    """ Class to define the stress boundary conditions """

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, nodeSets=tuple(), materials=None,
                 order_wall_conditions=None):
        """ Defines stress boundary conditions

        A condition can be a stress (float, int, Pint value) or an
        Underworld function which evaluates as a stress.

        parameters
        ----------

            Model: (UWGeodynamics.Model)
                An UWGeodynamics Model (See UWGeodynamics.Model)

            left:(tuple)
                Define conditions on the left side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            right:(tuple)
                Define conditions on the right side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            top:(tuple)
                Define conditions on the top side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            bottom:(tuple)
                Define conditions on the bottom side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            nodeSets: list of tuples: [(nodes, condition)]
                List of node where to apply predefined condition.

                The nodes can be a list or a numpy array containing the
                local index of the nodes. It can also be an Underworld
                IndexSet. You can also pass an UWGeodynamics shape.

                The condition can be an Underworld Function, A Pint
                Quantity of a scalar.

            materials: list of tuples: [(Material, condition)]
                List of material on which to apply a condition.
                The materials must be UWGeodynamics Material objects.

                The condition can be an Underworld Function, A Pint
                Quantity of a scalar.

            order_wall_conditions: list of str, [left, right, top, bottom,
                front, back]
                Order in which the boundaries are processed.

            Only valid for 3D Models:

            front:(tuple) with length 2 in 2D, and length 3 in 3D.
                Define mechanical conditions on the front side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            back:(tuple) with length 2 in 2D, and length 3 in 3D.
                Define mechanical conditions on the front side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

        """
        super(StressBCs, self).__init__(Model, Model.tractionField,
                                        Model.velocityField,
                                        left, right, top, bottom,
                                        front, back, nodeSets, materials,
                                        order_wall_conditions, "Neumann")


_dim_temp = {"[temperature]": 1.0}
_dim_heat_flux = {"[mass]": 1.0, "[time]": -3.0}


class TemperatureBCs(BoundaryConditions):

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, nodeSets=None, materials=None,
                 order_wall_conditions=None):
        """ Define temperature boundaries condition

        A condition can be a temperature (float, int, Pint Quantity) or
        an Underworld function which evaluates as a temperature.

        parameters
        ----------

            Model: (UWGeodynamics.Model)
                An UWGeodynamics Model (See UWGeodynamics.Model)

            left:
                Define conditions on the left side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            right:
                Define conditions on the right side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            top:
                Define conditions on the top side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            bottom:
                Define conditions on the bottom side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            nodeSets: list of tuples: [(nodes, condition)]
                List of node where to apply predefined condition.

                The nodes can be a list or a numpy array containing the
                local index of the nodes. It can also be an Underworld
                IndexSet. You can also pass an UWGeodynamics shape.

                The condition can be an Underworld Function, A Pint
                Quantity of a scalar.

            materials: list of tuples: [(Material, condition)]
                List of material on which to apply a condition.
                The materials must be UWGeodynamics Material objects.

                The condition can be an Underworld Function, A Pint
                Quantity of a scalar.

            order_wall_conditions: list of str, [left, right, top, bottom,
                front, back]
                Order in which the boundaries are processed.

            Only valid for 3D Models:

            front:
                Define conditions on the front side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            back:
                Define conditions on the front side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

        """

        super(TemperatureBCs, self).__init__(
            Model, Model.temperature, None, left, right, top, bottom,
            front, back, nodeSets, materials, order_wall_conditions,
            "Dirichlet")


class HeatFlowBCs(BoundaryConditions):

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, nodeSets=None, materials=None,
                 order_wall_conditions=None):
        """ Define heat flow boundaries condition

        A condition can be a heat flow (float, int, Pint Quantity) or
        an Underworld function which evaluates as a heat flow.

        parameters
        ----------

            Model: (UWGeodynamics.Model)
                An UWGeodynamics Model (See UWGeodynamics.Model)

            left:
                Define conditions on the left side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            right:
                Define conditions on the right side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            top:
                Define conditions on the top side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            bottom:
                Define conditions on the bottom side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            nodeSets: list of tuples: [(nodes, condition)]
                List of node where to apply predefined condition.

                The nodes can be a list or a numpy array containing the
                local index of the nodes. It can also be an Underworld
                IndexSet. You can also pass an UWGeodynamics shape.

                The condition can be an Underworld Function, A Pint
                Quantity of a scalar.

            materials: list of tuples: [(Material, condition)]
                List of material on which to apply a condition.
                The materials must be UWGeodynamics Material objects.

                The condition can be an Underworld Function, A Pint
                Quantity of a scalar.

            order_wall_conditions: list of str, [left, right, top, bottom,
                front, back]
                Order in which the boundaries are processed.

            Only valid for 3D Models:

            front:
                Define conditions on the front side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

            back:
                Define conditions on the front side of the Model.
                Conditions are defined for each Model direction (x, y, [z])

        """

        args = [left, right, top, bottom, front, back]

        message = """Heat Flux condition on the wall must be entered as
                     (Flux Value, Material) where Material is the material
                     located on the external side of the wall"""

        for arg in args:
            if arg:
                if not isinstance(arg, tuple):
                    raise ValueError(message)
                if not isinstance(arg[1], Material):
                    raise ValueError(message)
                if isinstance(arg[0], u.Quantity):
                    if not (arg[0].dimensionality == _dim_heat_flux):
                        raise ValueError(message)

        if left:
            left = self._get_heat_flux(left[0], left[1])

        if right:
            right = self._get_heat_flux(right[0], right[1])

        if top:
            top = self._get_heat_flux(top[0], right[1])

        if bottom:
            bottom = self._get_heat_flux(bottom[0], bottom[1])

        if front:
            front = self._get_heat_flux(front[0], front[1])

        if back:
            back = self._get_heat_flux(back[0], back[1])

        super(HeatFlowBCs, self).__init__(
            Model, Model._heatFlux, Model.temperature,
            left, right, top, bottom,
            front, back, nodeSets, materials,
            order_wall_conditions, "Neumann")

    @staticmethod
    def _get_heat_flux(heat_flow, material):

        if not (heat_flow and material):
            return

        if not material.capacity:
            raise ValueError("""Material {0} has no capacity
                             defined""".format(material.name))

        if not material.density.reference_density:
            raise ValueError("""Material {0} has no density
                             defined""".format(material.name))

        cp = material.capacity
        rho = material.density.reference_density
        return heat_flow / (rho * cp)
