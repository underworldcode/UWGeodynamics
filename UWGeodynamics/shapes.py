from __future__ import print_function,  absolute_import
import numpy as np
import underworld as uw
import underworld.function as fn
from .scaling import nonDimensionalize as nd


class Shape(object):
    """Base Class for Shape objects"""

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
        newShape = CombinedShape([self, B])
        return newShape

    def __or__(self, B):
        newShape = Shape()
        newShape.fn = self._fn | B._fn
        return newShape


class Polygon(Shape):
    """Polygon Shape Class"""

    def __init__(self, vertices):
        """Create a polygon shape

        Parameters
        ----------

        vertices : vertices of the polygon as (x,y) pairs.

        Returns
        -------

        Polygon Shape Class
        """
        self.vertices = vertices
        self.top = min([y for x, y in self.vertices])
        self.bottom = max([y for x, y in self.vertices])
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

        Parameters:
        -----------

        normal: A vector defining the normal to the plan.
        origin: Origin
        reverse: by default, particles tested against this class are
                 assigned "True" if they lay on or below the plan.
                 You can reverse than behavior by setting reverse=True.

        Returns:
        --------

        A UWGeodynamics Shape object.

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


class CombinedShape(Shape):
    """CombinedShape class"""

    def __init__(self, shapes):
        """Combine a list of shapes into one

        Parameters
        ----------

        shapes : list of Shape objects

        Returns
        -------

        An UWGeodynamics Shape object
        """
        self.shapes = shapes
        self.top = max([shape.top for shape in shapes])
        self.bottom = min([shape.top for shape in shapes])

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


class MultiShape(CombinedShape):

    def __init__(self, shapes):

        Warning("MultiShape is now deprecated, use CombinedShape instead")
        super(MultiShape, self).__init__(shapes)


class IntersectShape(Shape):
    """IntersectShape class"""

    def __init__(self, shapes):
        """ Take the intersection of some shapes

        Parameters
        ----------

        shapes : A list of Shapes

        Returns
        -------

        An UWGeodynamics Shape object
        """
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
    """Layer 2D"""

    def __init__(self, top, bottom):
        """Create a 2D Layer object

        Parameters
        ----------

        top : top of the layer
        bottom : bottom of the layer

        Returns
        -------

        AN UWGeodynamics Shape object
        """
        self.top = top
        self.bottom = bottom

    @property
    def _fn(self):
        coord = fn.input()
        func = ((coord[1] <= nd(self.top)) &
                (coord[1] >= nd(self.bottom)))
        return func


class Layer3D(Layer):
    """Layer3D"""

    def __init__(self, top, bottom):
        """Create a 3D Layer object

        Parameters
        ----------

        top : top of the layer
        bottom : bottom of the layer

        Returns
        -------

        AN UWGeodynamics Shape object
        """
        self.top = top
        self.bottom = bottom

    @property
    def _fn(self):
        coord = fn.input()
        func = ((coord[2] <= nd(self.top)) &
                (coord[2] >= nd(self.bottom)))
        return func

class Layer2D(Layer):

    def __init__(self, top, bottom):

        Warning("Layer2D is now deprecated, use Layer instead")
        super(Layer2D, self).__init__(top, bottom)

class Box(Shape):
    """Box"""

    def __init__(self, top, bottom, minX=0., maxX=0., minY=None, maxY=None):
        """Create a Box Shape

        Parameters
        ----------

        top : Top of the Box
        bottom : Bottom of the Box
        minX : Minimum extent of the Box along the x-axis
        maxX : Maximum extent of the Box along the x-axis

        Only in 3D:

        minY : Minimum extent of the Box along the y-axis
        maxY : Maximum extent of the Box along the y-axis

        Returns
        -------
        """
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


class Disk(Shape):
    """Disk"""

    def __init__(self, center, radius):
        """Create a Disk shape

        Parameters
        ----------

        center : center of the disk
        radius : radius of the disk

        Returns
        -------

        An UWGeodynamics Shape
        """
        self.center = center
        self.radius = radius
        self.top = center[1] + self.radius
        self.bottom = center[1] - self.radius

    @property
    def _fn(self):
        center = tuple(nd(x) for x in list(self.center))
        radius = nd(self.radius)
        coord = fn.input() - center
        return fn.math.dot(coord, coord) < radius**2


Sphere = Disk


class Annulus(Shape):
    """Annulus"""

    def __init__(self, center, r1, r2):
        """Create an Annulus shape

        Parameters
        ----------

        center : center of the annulus
        r1 : Internal radius
        r2 : External radius

        Returns
        -------

        An UWGeodynamics Shape object
        """
        self.center = center
        self.r1 = r1
        self.r2 = r2
        self.bottom = center[1] - self.r2
        self.top = center[1] + self.r2

    @property
    def _fn(self):
        center = tuple(nd(x) for x in list(self.center))
        r1 = nd(self.r1)
        r2 = nd(self.r2)
        coord = fn.input() - center
        return (fn.math.dot(coord, coord) < r2**2) & (fn.math.dot(coord, coord) > r1**2)
