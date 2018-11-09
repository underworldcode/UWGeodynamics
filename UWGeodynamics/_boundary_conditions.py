
class BoundaryConditions(object):

    def __init__(self, Model, left=None, right=None, top=None, bottom=None,
                 front=None, back=None, nodeSets=tuple(),
                 order_wall_conditions=None):

        self.Model = Model
        self.left = left
        self.right = right
        self.top = top
        self.bottom = bottom
        self.back = back
        self.front = front
        self.nodeSets = nodeSets

        if self.Model.mesh.dim == 2:
            self._wall_indexSets = {"bottom": (self.bottom,
                                               self.Model.bottom_wall),
                                    "top": (self.top,
                                            self.Model.top_wall),
                                    "left": (self.left,
                                             self.Model.left_wall),
                                    "right": (self.right,
                                              self.Model.right_wall)}
            if order_wall_conditions:
                if len(order_wall_conditions) <= 5:
                    self.order_wall_conditions = order_wall_conditions
            else:
                self.order_wall_conditions = ["bottom", "top", "left", "right"]

        if self.Model.mesh.dim == 3:
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
            if order_wall_conditions:
                if len(order_wall_conditions) <= 7:
                    self.order_wall_conditions = order_wall_conditions
            else:
                self.order_wall_conditions = ["bottom", "top", "front", "back",
                                              "left", "right"]

