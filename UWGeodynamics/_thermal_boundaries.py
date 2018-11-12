import underworld as uw
import numpy as np
from .scaling import nonDimensionalize as nd
from .scaling import UnitRegistry as u
from ._boundary_conditions import BoundaryConditions

_dim_temp = {"[temperature]":1.0}
_dim_heat_flux = {"[mass]":1.0, "[time]":-3.0}

class TemperatureBCs(BoundaryConditions):

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, nodeSets=None, materials=None,
                 bottom_material=None, top_material=None, left_material=None,
                 right_material=None, back_material=None, front_material=None,
                 order_wall_conditions=None):

        super().__init__(Model, left, right, top, bottom, front, back,
                         nodeSets, order_wall_conditions)

        self.materials = materials
        self.bottom_material = bottom_material
        self.top_material = top_material
        self.left_material = left_material
        self.right_material = right_material
        self.back_material = back_material
        self.front_material = front_material

        self._material_boundaries = {"left": left_material,
                                     "right": right_material,
                                     "top": top_material,
                                     "bottom": bottom_material,
                                     "back": back_material,
                                     "front": front_material}

        for key in self._material_boundaries.keys():
            wall = self.__dict__[key]
            if wall:
                wall_material = self._material_boundaries[key]
                if (wall.dimensionality == _dim_heat_flux) and not(wall and wall_material):
                    raise ValueError(
                        """left wall is {0} and left_material is {1}""".format(
                            wall, wall_material))

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
        # Reinitialise neumnann and dirichlet condition
        self.dirichlet_indices = []
        self.neumann_indices = []

        self.dirichlet_indices.append(Model.mesh.specialSets["Empty"])
        self.neumann_indices.append(Model.mesh.specialSets["Empty"])

        for wall in self.order_wall_conditions:
            condition, indexSet = self._wall_indexSets[wall]
            if condition:
                if condition.dimensionality == _dim_temp:
                    Model.temperature.data[indexSet.data] = nd(condition)
                    self.dirichlet_indices[0] += indexSet
                elif condition.dimensionality == _dim_heat_flux:
                    material = self._material_boundaries[wall]
                    self.neumann_indices[0] += indexSet
                    heat_flow = self._get_heat_flux(condition, material)
                    Model._heatFlux.data[indexSet.data] = nd(heat_flow)

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

        conditions = []

        conditions.append(uw.conditions.DirichletCondition(
            variable=Model.temperature,
            indexSetsPerDof=self.dirichlet_indices))

        neumann_indices = []
        for val in self.neumann_indices:
            if val.data.size > 0:
                neumann_indices.append(val)
            else:
                neumann_indices.append(None)
        neumann_indices = tuple(neumann_indices)

        if neumann_indices != tuple([None]):
            conditions.append(uw.conditions.NeumannCondition(
                fn_flux=Model._heatFlux,
                variable=Model.temperature,
                indexSetsPerDof=self.neumann_indices))

        if not conditions:
            raise ValueError("Undefined conditions")

        return conditions

