from __future__ import print_function,  absolute_import


class BoundaryConditions(object):
    """Base Class for BOundary Condition objects"""

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
