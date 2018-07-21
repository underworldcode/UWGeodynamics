import numpy as np
import underworld as uw
import underworld.function as fn
from .scaling import nonDimensionalize as nd


class Shape(object):

    def __init__(self):
        self._fn = None
        return

    @property
    def fn(self):
        return self._fn

    @fn.setter
    def fn(self, value):
        self._fn = value

    def __and__(self, B):
        newShape = Shape()
        newShape.fn = self._fn & B._fn
        return newShape

    def __add__(self, B):
        newShape = MultiShape([self, B])
        return newShape

    def __or__(self, B):
        newShape = Shape()
        newShape.fn = self._fn | B._fn
        return newShape

class Polygon(Shape):

    def __init__(self, vertices):
        self.vertices = vertices
        vertices = [(nd(x), nd(y)) for x, y in self.vertices]
        self._fn = uw.function.shape.Polygon(np.array(vertices))


class HalfSpace(Shape):
    """ Class to define a HalfSpace

        A plan defined by a normal vector is used to split the space
        in two half spaces. By default the origin of the coordinate
        system is at (0., 0.).

        Particles tested against this class are assigned a boolean value.
    """

    def __init__(self, normal, origin=None, reverse=False):
        """ HalfSpace

        arguments:

            normal: A vector defining the normal to the plan.
            origin: Origin
            reverse: by default, particles tested against this class are
                     assigned "True" if they lay on or below the plan.
                     You can reverse than behavior by setting reverse=True.
        """

        if isinstance(normal, (tuple, list)):
            self.normal = fn.misc.constant([float(nd(val)) for val in normal])
        else:
            raise ValueError("{0} must be a list or tuple".format(normal))

        if isinstance(origin, (tuple, list)):
            self.origin = fn.misc.constant([float(nd(val)) for val in origin])
        else:
            self.origin = fn.misc.constant([0.] * len(normal))

        self.reverse = reverse

    @property
    def _fn(self):
        coords = fn.input()
        new_coords = coords - self.origin
        func = fn.math.dot(self.normal, new_coords)

        # True if below, False if above
        if not self.reverse:
            conditions = [(func <= 0., True), (func > 0., False)]
        else:
            conditions = [(func >= 0., True), (func < 0., False)]

        return fn.branching.conditional(conditions)


class MultiShape(Shape):

    def __init__(self, shapes):
        self.shapes = shapes

    @property
    def _fn(self):
        import operator
        import functools
        self._fnlist = []
        for shape in self.shapes:
            self._fnlist.append(shape._fn)
        func = functools.reduce(
            operator.or_,
            self._fnlist,
            fn.misc.constant(False))
        return func


class CombinedShape(Shape):

    def __init__(self, shapes):
        self.shapes = shapes

    @property
    def _fn(self):
        import operator
        import functools
        self._fnlist = []
        for shape in self.shapes:
            self._fnlist.append(shape._fn)
        func = functools.reduce(
            operator.and_,
            self._fnlist,
            fn.misc.constant(True))
        return func


class Layer(Shape):
    def __init__(self, top, bottom):
        self.top = top
        self.bottom = bottom


class Layer2D(Shape):

    def __init__(self, top, bottom):
        self.top = top
        self.bottom = bottom

    @property
    def _fn(self):
        coord = fn.input()
        func = ((coord[1] <= nd(self.top)) &
                (coord[1] >= nd(self.bottom)))
        return func

    @property
    def top(self):
        return self._top

    @top.setter
    def top(self, value):
        self._top = value

    @property
    def bottom(self):
        return self._bottom

    @bottom.setter
    def bottom(self, value):
        self._bottom = value


class Layer3D(Shape):

    def __init__(self, top, bottom):
        self.top = top
        self.bottom = bottom

    @property
    def _fn(self):
        coord = fn.input()
        func = ((coord[2] <= nd(self.top)) &
                (coord[2] >= nd(self.bottom)))
        return func

    @property
    def top(self):
        return self._top

    @top.setter
    def top(self, value):
        self._top = value

    @property
    def bottom(self):
        return self._bottom

    @bottom.setter
    def bottom(self, value):
        self._bottom = value


class Box(Shape):

    def __init__(self, top, bottom, minX=0., maxX=0., minY=None, maxY=None):
        self.top = top
        self.bottom = bottom
        self.minX = minX
        self.maxX = maxX
        self.minY = minY
        self.maxY = maxY

    @property
    def _fn(self):
        coord = fn.input()
        if (self.minY is not None) and (self.maxY is not None):
            func = ((coord[1] <= nd(self.maxY)) &
                    (coord[1] >= nd(self.minY)) &
                    (coord[0] <= nd(self.maxX)) &
                    (coord[0] >= nd(self.minX)) &
                    (coord[2] <= nd(self.top)) &
                    (coord[2] >= nd(self.bottom)))
        else:
            func = ((coord[1] <= nd(self.top)) &
                    (coord[1] >= nd(self.bottom)) &
                    (coord[0] <= nd(self.maxX)) &
                    (coord[0] >= nd(self.minX)))
        return func

    @property
    def minX(self):
        return self._minX

    @minX.setter
    def minX(self, value):
        self._minX = value

    @property
    def maxX(self):
        return self._maxX

    @maxX.setter
    def maxX(self, value):
        self._maxX = value

    @property
    def top(self):
        return self._top

    @top.setter
    def top(self, value):
        self._top = value

    @property
    def bottom(self):
        return self._bottom

    @bottom.setter
    def bottom(self, value):
        self._bottom = value


class Disk(Shape):

    def __init__(self, center, radius):
        self.center = center
        self.radius = radius

    @property
    def _fn(self):
        center = tuple(nd(x) for x in list(self.center))
        radius = nd(self.radius)
        coord = fn.input() - center
        return fn.math.dot(coord, coord) < radius**2


Sphere = Disk


class Annulus(Shape):

    def __init__(self, center, r1, r2):
        self.center = center
        self.r1 = r1
        self.r2 = r2

    @property
    def _fn(self):
        center = tuple(nd(x) for x in list(self.center))
        r1 = nd(self.r1)
        r2 = nd(self.r2)
        coord = fn.input() - center
        return (fn.math.dot(coord, coord) < r2**2) & (fn.math.dot(coord, coord) > r1**2)

