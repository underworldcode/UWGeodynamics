from __future__ import print_function,  absolute_import
import underworld as uw
import underworld.function as fn
import numpy as np
from mpi4py import MPI
from .LecodeIsostasy import LecodeIsostasy
from .scaling import nonDimensionalize as nd
from .scaling import UnitRegistry as u
from ._utils import Balanced_InflowOutflow
from ._utils import MovingWall
from ._boundary_conditions import BoundaryConditions

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


class VelocityBCs(BoundaryConditions):
    """ Class to define the mechanical boundary conditions """

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, nodeSets=tuple(),
                 order_wall_conditions=None):
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
        indexSets: [(condition, IndexSet)]
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

        super(VelocityBCs, self).__init__(Model, left, right, top, bottom,
                                          front, back, nodeSets,
                                          order_wall_conditions)
        # Link Moving Walls
        for arg in [self.left, self.right, self.top, self.bottom, self.front,
                    self.back]:
            if isinstance(arg, MovingWall):
                arg.Model = self.Model

    def __getitem__(self, name):
        return self.__dict__[name]

    def _apply_conditions_nodes(self, condition, nodes):
        """ Apply condition to a set of nodes

        Parameters:
        -----------
            condition:
                velocity condition
            nodes:
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
            self._dirichlet_indices[-1] += self.Model.bottom_wall
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
                        self.Model.velocityField.data[ISet.data, dim] = (
                            func.evaluate(ISet)[:, 0])
                        self._dirichlet_indices[dim] += ISet
                    else:
                        self.Model.velocityField.data[ISet.data, dim] = 0.
                        self._dirichlet_indices[dim] += ISet

            return

        # Expect a list or tuple of dimension mesh.dim.
        # Check that the domain actually contains some boundary nodes
        # (nodes is not None)
        if isinstance(condition, (list, tuple)) and nodes.data.size > 0:
            for dim in range(self.Model.mesh.dim):

                if isinstance(condition[dim], fn.Function):
                    func = condition[dim]
                    self.Model.velocityField.data[nodes.data, dim] = (
                        func.evaluate(
                            self.Model.mesh.data[nodes.data])[:, dim])
                    self._dirichlet_indices[dim] += nodes

                # User defined function
                if isinstance(condition[dim], (list, tuple)):
                    func = fn.branching.conditional(condition[dim])
                    self.Model.velocityField.data[nodes.data, dim] = (
                        func.evaluate(
                            self.Model.mesh.data[nodes.data])[:, 0])
                    self._dirichlet_indices[dim] += nodes

                # Scalar condition
                if isinstance(condition[dim], (u.Quantity, float, int)):
                    self.Model.velocityField.data[nodes.data, dim] = (
                        nd(condition[dim]))
                    self._dirichlet_indices[dim] += nodes

                # Inflow Outflow
                if isinstance(condition[dim], Balanced_InflowOutflow):
                    obj = condition[dim]
                    obj.ynodes = self.Model.mesh.data[nodes.data, 1]
                    obj._get_side_flow()
                    self.Model.velocityField.data[nodes.data, dim] = (
                        obj._get_side_flow())
                    self._dirichlet_indices[dim] += nodes

                if isinstance(condition[dim], LecodeIsostasy):
                    # Apply support condition
                    if self.Model.mesh.dim - 1 != dim:
                        raise ValueError("""Can not apply LecodeIsostasy on that
                                         dimension""")

                    self.Model._isostasy = condition[dim]
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
                    self._dirichlet_indices[dim] += nodes

        return

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

        Model = self.Model

        # Reinitialise neumnann and dirichlet condition
        self._dirichlet_indices = []

        for _ in range(Model.mesh.dim):
            self._dirichlet_indices.append(Model.mesh.specialSets["Empty"])

        for set_ in self.order_wall_conditions:
            (condition, nodes) = self._wall_indexSets[set_]
            if nodes is not None:
                self._apply_conditions_nodes(condition, nodes)

        if self.nodeSets:
            for (condition, nodes) in self.nodeSets:
                self._apply_conditions_nodes(condition, nodes)

        dirichlet_conditions = None

        dirichlet_conditions = uw.conditions.DirichletCondition(
            variable=Model.velocityField,
            indexSetsPerDof=self._dirichlet_indices)

        return dirichlet_conditions


class StressBCs(BoundaryConditions):
    """ Class to define the stress boundary conditions """

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, nodeSets=tuple(),
                 order_wall_conditions=None):
        """ Defines stress boundary conditions

        parameters
        ----------

        Model: (UWGeodynamics.Model)
            An UWGeodynamics Model (See UWGeodynamics.Model)
        left:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define conditions on the left side of the Model.
            Conditions are defined for each Model direction (x, y, [z])
        right:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define conditions on the right side of the Model.
            Conditions are defined for each Model direction (x, y, [z])
        top:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define conditions on the top side of the Model.
            Conditions are defined for each Model direction (x, y, [z])
        bottom:(tuple) with length 2 in 2D, and length 3 in 3D.
            Define conditions on the bottom side of the Model.
            Conditions are defined for each Model direction (x, y, [z])
        indexSets: [(condition, IndexSet)]
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

        """
        super(StressBCs, self).__init__(Model, left, right, top, bottom,
                         front, back, nodeSets, order_wall_conditions)

    def __getitem__(self, name):
        return self.__dict__[name]

    def _apply_conditions_nodes(self, condition, nodes):
        """ Apply condition to a set of nodes

        Parameters:
        -----------
            condition:
                velocity condition
            nodes:
                set of nodes

        """

        if not nodes:
            return

        # Expect a list or tuple of dimension mesh.dim.
        # Check that the domain actually contains some boundary nodes
        # (nodes is not None)
        if isinstance(condition, (list, tuple)) and nodes.data.size > 0:
            for dim in range(self.Model.mesh.dim):

                # Scalar condition
                if isinstance(condition[dim], (u.Quantity, float, int)):
                    self.Model.tractionField.data[nodes.data, dim] = (
                        nd(condition[dim]))
                    self._neumann_indices[dim] += nodes

        return

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

        Model = self.Model

        # Reinitialise neumnann condition
        self._neumann_indices = []

        for _ in range(Model.mesh.dim):
            self._neumann_indices.append(Model.mesh.specialSets["Empty"])

        for set_ in self.order_wall_conditions:
            (condition, nodes) = self._wall_indexSets[set_]
            if nodes is not None:
                self._apply_conditions_nodes(condition, nodes)

        if self.nodeSets:
            for (condition, nodes) in self.nodeSets:
                self._apply_conditions_nodes(condition, nodes)

        self.neumann_conditions = None
        _neumann_indices = []

        # Remove empty Sets
        for val in self._neumann_indices:
            if val.data.size > 0:
                _neumann_indices.append(val)
            else:
                _neumann_indices.append(None)
        self._neumann_indices = tuple(_neumann_indices)

        # Now we only create a Neumann condition if we have a stress condition
        # somewhere, on any of the procs.
        local_procs_has_neumann = np.zeros((uw.nProcs()))
        global_procs_has_neumann = np.zeros((uw.nProcs()))

        if self._neumann_indices != tuple([None] * Model.mesh.dim):
            local_procs_has_neumann[uw.rank()] = 1

        comm.Allreduce(local_procs_has_neumann, global_procs_has_neumann)
        comm.Barrier()

        if any(global_procs_has_neumann):
            self.neumann_conditions = uw.conditions.NeumannCondition(
                fn_flux=Model.tractionField,
                variable=Model.velocityField,
                indexSetsPerDof=self._neumann_indices)

        return self.neumann_conditions

