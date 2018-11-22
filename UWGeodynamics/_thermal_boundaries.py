from __future__ import print_function,  absolute_import
import underworld as uw
import numpy as np
from mpi4py import MPI
from .scaling import nonDimensionalize as nd
from .scaling import UnitRegistry as u
from ._boundary_conditions import BoundaryConditions
from ._material import Material

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

_dim_temp = {"[temperature]": 1.0}
_dim_heat_flux = {"[mass]": 1.0, "[time]": -3.0}


class TemperatureBCs(BoundaryConditions):

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, nodeSets=None, materials=None,
                 order_wall_conditions=None):

        self.materials = materials

        super(TemperatureBCs, self).__init__(
            Model, left, right, top, bottom, front, back,
            nodeSets, order_wall_conditions)

    def __getitem__(self, name):
        return self.__dict__[name]

    def get_conditions(self):

        Model = self.Model

        self.dirichlet_indices = []
        self.dirichlet_indices.append(Model.mesh.specialSets["Empty"])

        for wall in self.order_wall_conditions:
            condition, indexSet = self._wall_indexSets[wall]
            if condition and indexSet is not None:
                Model.temperature.data[indexSet.data] = nd(condition)
                self.dirichlet_indices[0] += indexSet

        if self.nodeSets:
            for elem in self.nodeSets:
                nodeSet, temp = elem
                isThere = np.in1d(Model.mesh.data_nodegId, list(nodeSet))
                if isThere.any():
                    local_indices = np.arange(Model.mesh.nodesDomain)[isThere]
                    indexSet = uw.mesh.FeMesh_IndexSet(
                        self.Model.mesh, topologicalIndex=0,
                        size=self.Model.mesh.nodesGlobal,
                        fromObject=local_indices)
                    Model.temperature.data[local_indices] = nd(temp)
                    self.dirichlet_indices[0] += indexSet

        if self.materials:
            for (material, temp) in self.materials:
                if material and nd(temp):
                    indexSet = Model._get_material_indices(material)
                    Model.temperature.data[indexSet.data] = nd(temp)
                    self.dirichlet_indices[0] += indexSet

        return uw.conditions.DirichletCondition(
            variable=Model.temperature,
            indexSetsPerDof=self.dirichlet_indices)


class HeatFlowBCs(BoundaryConditions):

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, order_wall_conditions=None):

        super(HeatFlowBCs, self).__init__(
            Model, left, right, top, bottom,
            front, back, nodeSets=None,
            order_wall_conditions=order_wall_conditions)

        args = [left, right, top, bottom, front, back]

        message = """Heat Flux condition on the wall must be entered as
                     (Flux Value, Material) where Material is the material
                     located on the external side of the wall"""

        for arg in args:
            if arg:
                if not isinstance(arg, tuple):
                    raise ValueError(message)
                if not isinstance(arg[0], (float, int, u.Quantity)):
                    raise ValueError
                if not isinstance(arg[1], Material):
                    raise ValueError(message)
                if isinstance(arg[0], u.Quantity):
                    if not (arg[0].dimensionality == _dim_heat_flux):
                        raise ValueError(message)

    def __getitem__(self, name):
        return self.__dict__[name]

    @staticmethod
    def _get_heat_flux(heat_flow, material):

        if not material.capacity:
            raise ValueError("""Material {0} has no capacity
                             defined""".format(material.name))

        if not material.density.reference_density:
            raise ValueError("""Material {0} has no density
                             defined""".format(material.name))

        cp = material.capacity
        rho = material.density.reference_density
        return heat_flow / (rho * cp)

    def get_conditions(self):

        Model = self.Model
        # neumnann and dirichlet condition
        self._neumann_indices = []

        self._neumann_indices.append(Model.mesh.specialSets["Empty"])

        for wall in self.order_wall_conditions:
            if wall:
                condition, indexSet = self._wall_indexSets[wall]
                if condition and indexSet is not None:
                    (condition, material) = condition
                    self._neumann_indices[0] += indexSet
                    heat_flow = self._get_heat_flux(condition, material)
                    Model._heatFlux.data[indexSet.data] = nd(heat_flow)

        self.neumann_conditions = None
        _neumann_indices = []
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
                fn_flux=Model._heatFlux,
                variable=Model.temperature,
                indexSetsPerDof=self._neumann_indices)

        return self.neumann_conditions
