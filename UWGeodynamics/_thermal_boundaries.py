import underworld as uw
from .scaling import nonDimensionalize as nd

class TemperatureBCs(object):

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None,
                 indexSets=[], materials=None):


        self.Model = Model
        self.left = left
        self.right = right
        self.top = top
        self.bottom = bottom
        self.front = front
        self.back = back
        self.indexSets = indexSets
        self.materials = materials

    def get_conditions(self):

        Model = self.Model

        indices = [Model.mesh.specialSets["Empty"]]
        if self.left is not None:
            Model.temperature.data[Model._left_wall.data] = nd(self.left)
            indices[0] += Model._left_wall
        if self.right is not None:
            Model.temperature.data[Model._right_wall.data] = nd(self.right)
            indices[0] += Model._right_wall
        if self.top is not None:
            Model.temperature.data[Model._top_wall.data] = nd(self.top)
            indices[0] += Model._top_wall
        if self.bottom is not None:
            Model.temperature.data[Model._bottom_wall.data] = nd(self.bottom)
            indices[0] += Model._bottom_wall
        if self.back is not None and Model.mesh.dim > 2:
            Model.temperature.data[Model._back_wall.data] = nd(self.back)
            indices[0] += Model._back_wall
        if self.front is not None and Model.mesh.dim > 2:
            Model.temperature.data[Model._front_wall.data] = nd(self.front)
            indices[0] += Model._front_wall
        if self.indexSets:
            for indexSet, temp in self.indexSets:
                Model.temperature.data[indexSet.data] = nd(temp)
                indices[0] += indexSet
        if self.materials:
            for (material, temp) in self.materials:
                if material and nd(temp):
                    indexSet = Model._get_material_indices(material)
                    Model.temperature.data[indexSet.data] = nd(temp)
                    indices[0] += indexSet

        return uw.conditions.DirichletCondition(variable=Model.temperature,
                                                 indexSetsPerDof=indices)
