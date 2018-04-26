import underworld as uw
import numpy as np
from .scaling import nonDimensionalize as nd
from .scaling import UnitRegistry as u


class TemperatureBCs(object):

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, indexSets=None, materials=None,
                 bottom_material=None, top_material=None, left_material=None,
                 right_material=None, back_material=None, front_material=None):

        self.Model = Model
        self.left = left
        self.right = right
        self.top = top
        self.bottom = bottom
        self.front = front
        self.back = back
        self.indexSets = indexSets
        self.materials = materials
        self.bottom_material = bottom_material
        self.top_material = top_material
        self.left_material = left_material
        self.right_material = right_material
        self.back_material = back_material
        self.front_material = front_material

    def __getitem__(self, name):
        return self.__dict__[name]

    def _get_heat_flux(self, side, heat_flow, material):

        if not material.capacity:
            raise ValueError("""Material {0} has no capacity
                             defined""".format(material.name))

        if not material.density.reference_density:
            raise ValueError("""Material {0} has no density
                             defined""".format(material.name))

        cp = material.capacity
        rho = material.density.reference_density
        return heat_flow / (rho * cp)

    def _check_temp(self, val):
        if isinstance(val, u.Quantity):
            check = (val.units == u.degK) or (val.units == u.degC)
        return check

    def _check_flux(self, heat_flow):
        if isinstance(heat_flow, u.Quantity):
            val = heat_flow.to_base_units()
            check = (val.units == u.kilogram * u.second**-3)
        return check

    def get_conditions(self):

        Model = self.Model
        # Reinitialise neumnann and dirichlet condition
        self.dirichlet_indices = []
        self.neumann_indices = []

        self.dirichlet_indices.append(Model.mesh.specialSets["Empty"])
        self.neumann_indices.append(Model.mesh.specialSets["Empty"])

        # A lot of the following code is redundant and could be simplified.
        # It is however rather easy to follow

        self.dirichlet_indices = [Model.mesh.specialSets["Empty"]]
        if self.left is not None:
            if self._check_temp(self.left):
                Model.temperature.data[Model._left_wall.data] = nd(self.left)
                self.dirichlet_indices[0] += Model._left_wall
            elif self._check_flux(self.left):
                if not self.left_material:
                    raise ValueError("left_material is missing")
                material = self.left_material
                self.neumann_indices[0] += Model._left_wall
                heat_flow = self._get_heat_flux("left", self.left, material)
                Model._heatFlux.data[Model._left_wall.data] = (
                    nd(heat_flow))

        if self.right is not None:
            if self._check_temp(self.right):
                Model.temperature.data[Model._right_wall.data] = nd(self.right)
                self.dirichlet_indices[0] += Model._right_wall
            elif self._check_flux(self.right):
                if not self.right_material:
                    raise ValueError("right_material is missing")
                material = self.right_material
                self.neumann_indices[0] += Model._right_wall
                heat_flow = self._get_heat_flux("right", self.right, material)
                Model._heatFlux.data[Model._right_wall.data] = (
                    nd(heat_flow))

        if self.top is not None:
            if self._check_temp(self.top):
                Model.temperature.data[Model._top_wall.data] = nd(self.top)
                self.dirichlet_indices[0] += Model._top_wall
            elif self._check_flux(self.top):
                if not self.top_material:
                    raise ValueError("top_material is missing")
                material = self.top_material
                self.neumann_indices[0] += Model._top_wall
                heat_flow = self._get_heat_flux("top", self.top, material)
                Model._heatFlux.data[Model._top_wall.data] = (
                    nd(heat_flow))

        if self.bottom is not None:
            if self._check_temp(self.bottom):
                Model.temperature.data[Model._bottom_wall.data] = (
                    nd(self.bottom))
                self.dirichlet_indices[0] += Model._bottom_wall
            elif self._check_flux(self.bottom):
                if not self.bottom_material:
                    raise ValueError("bottom_material is missing")
                material = self.bottom_material
                self.neumann_indices[0] += Model._bottom_wall
                heat_flow = self._get_heat_flux(
                    "bottom", self.bottom, material)
                Model._heatFlux.data[Model._bottom_wall.data] = (
                    nd(heat_flow))

        if self.back is not None and Model.mesh.dim > 2:
            if self._check_temp(self.back):
                Model.temperature.data[Model._back_wall.data] = nd(self.back)
                self.dirichlet_indices[0] += Model._back_wall
            elif self._check_flux(self.back):
                if not self.back_material:
                    raise ValueError("back_material is missing")
                material = self.back_material
                self.neumann_indices[0] += Model._back_wall
                heat_flow = self._get_heat_flux("back", self.back, material)
                Model._heatFlux.data[Model._back_wall.data] = (
                    nd(heat_flow))

        if self.front is not None and Model.mesh.dim > 2:
            if self._check_temp(self.front):
                Model.temperature.data[Model._front_wall.data] = nd(self.front)
                self.dirichlet_indices[0] += Model._front_wall
            elif self._check_flux(self.front):
                if not self.front_material:
                    raise ValueError("front_material is missing")
                material = self.front_material
                self.neumann_indices[0] += Model._front_wall
                heat_flow = self._get_heat_flux("front", self.front, material)
                Model._heatFlux.data[Model._front_wall.data] = (
                    nd(heat_flow))

        if self.indexSets:
            for elem in self.indexSets:
                indexSet, temp = elem
                isThere = np.in1d(Model.mesh.data_nodegId, list(indexSet))
                if isThere.any():
                    local_indices = np.arange(Model.mesh.nodesDomain)[isThere]
                    indexSet = uw.mesh.FeMesh_IndexSet(
                        self.Model.mesh, topologicalIndex=0,
                        size=self.Model.mesh.nodesGlobal,
                        fromObject=local_indices)
                    self._check_temp(temp)
                    Model.temperature.data[local_indices] = nd(temp)
                    self.dirichlet_indices[0] += indexSet

        if self.materials:
            for (material, temp) in self.materials:
                if material and nd(temp):
                    self._check_temp(temp)
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

