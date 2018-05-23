import underworld as uw
import underworld.function as fn
import numpy as np
from mpi4py import MPI
from .LecodeIsostasy import LecodeIsostasy
from .scaling import nonDimensionalize as nd
from .scaling import UnitRegistry as u
from ._utils import Balanced_InflowOutflow

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def _is_neumann(val):
    """ Returns true if x as units of stress """

    if not isinstance(val, u.Quantity):
        return False
    val = val.to_base_units()
    return val.units == u.kilogram / (u.meter * u.second**2)


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

        self.dirichlet_indices = []
        self.neumann_indices = []

    def __getitem__(self, name):
        return self.__dict__[name]

    def apply_condition_nodes(self, condition, nodes):
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

                # User defined function
                if isinstance(condition[dim], (list, tuple)):
                    func = fn.branching.conditional(condition[dim])
                    self.Model.velocityField.data[nodes.data, dim] = (
                        func.evaluate(
                            self.Model.mesh.data[nodes.data])[:, dim])
                    self.dirichlet_indices[dim] += nodes
                    continue

                # Scalar condition
                if isinstance(condition[dim], (u.Quantity, float, int)):

                    # Process dirichlet condition
                    if not _is_neumann(condition[dim]):
                        self.Model.velocityField.data[nodes.data, dim] = (
                            nd(condition[dim]))
                        self.dirichlet_indices[dim] += nodes
                    # Process neumann condition
                    else:
                        self.Model.tractionField.data[nodes.data, dim] = (
                            nd(condition[dim]))
                        self.neumann_indices[dim] += nodes

                # Inflow Outflow
                if isinstance(condition[dim], Balanced_InflowOutflow):
                    obj = condition[dim]
                    obj.ynodes = self.Model.mesh.data[nodes.data, 1]
                    obj._get_side_flow()
                    self.Model.velocityField.data[nodes.data, dim] = (
                        obj._get_side_flow())
                    self.dirichlet_indices[dim] += nodes

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
                    self.dirichlet_indices[dim] += nodes

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
        self.dirichlet_indices = []
        self.neumann_indices = []

        for dim in range(Model.mesh.dim):
            self.dirichlet_indices.append(Model.mesh.specialSets["Empty"])
            self.neumann_indices.append(Model.mesh.specialSets["Empty"])

        # Apply condition to the bottom first
        if isinstance(self.bottom, LecodeIsostasy):
            # Apply support condition
            Model._isostasy = self.bottom
            Model._isostasy.mesh = Model.mesh
            Model._isostasy.swarm = Model.swarm
            Model._isostasy._mesh_advector = Model._advector
            Model._isostasy.velocityField = Model.velocityField
            Model._isostasy.materialIndexField = Model.materialField
            Model._isostasy._densityFn = Model._densityFn
            vertical_walls_conditions = {
                "left": self.left,
                "right": self.right,
                "front": self.front,
                "back": self.back
            }
            Model._isostasy.vertical_walls_conditions = (
                vertical_walls_conditions)
            self.dirichlet_indices[-1] += Model._bottom_wall
        else:
            Model._isostasy = None
            self.apply_condition_nodes(self.bottom, Model._bottom_wall)

        # Apply condition at the top.
        self.apply_condition_nodes(self.top, Model._top_wall)

        # Apply conditions to each wall.
        self.apply_condition_nodes(self.left, Model._left_wall)
        self.apply_condition_nodes(self.right, Model._right_wall)
        self.apply_condition_nodes(self.indexSets, self.indexSets)

        if Model.mesh.dim > 2:
            self.apply_condition_nodes(self.front, Model._front_wall)
            self.apply_condition_nodes(self.back, Model._back_wall)

        conditions = []

        conditions.append(uw.conditions.DirichletCondition(
            variable=Model.velocityField,
            indexSetsPerDof=self.dirichlet_indices))

        neumann_indices = []

        # Remove empty Sets
        for val in self.neumann_indices:
            if val.data.size > 0:
                neumann_indices.append(val)
            else:
                neumann_indices.append(None)
        self.neumann_indices = tuple(neumann_indices)

        # Now we only create a Neumann condition if we have a stress condition
        # somewhere, on any of the procs.
        local_procs_has_neumann = np.zeros((uw.nProcs()))
        global_procs_has_neumann = np.zeros((uw.nProcs()))
        if self.neumann_indices != tuple([None for val in range(Model.mesh.dim)]):
            local_procs_has_neumann[uw.rank()] = 1

        comm.Allreduce(local_procs_has_neumann, global_procs_has_neumann)
        comm.Barrier()

        if any(global_procs_has_neumann):
            conditions.append(uw.conditions.NeumannCondition(
                fn_flux=Model.tractionField,
                variable=Model.velocityField,
                indexSetsPerDof=self.neumann_indices))

        if not conditions:
            raise ValueError("Undefined conditions")

        return conditions

