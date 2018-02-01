import glucifer
import underworld as uw
from mpi4py import MPI
from .scaling import nonDimensionalize as nd
from .scaling import Dimensionalize
from .scaling import UnitRegistry as u
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def visugrid_drawing_object2(Model):

    visugrid = Model.visugrid.mesh

    m_minx, m_maxx = _get_minmax_coordinates_mesh(visugrid, 0)
    m_miny, m_maxy = _get_minmax_coordinates_mesh(visugrid, 1)

    v_minx, v_maxx = _get_minmax_coordinates_mesh(Model.mesh, 0)
    v_miny, v_maxy = _get_minmax_coordinates_mesh(Model.mesh, 1)

    minx = min(m_minx, v_minx)
    miny = min(m_miny, v_miny)
    maxx = max(m_maxx, v_maxx)
    maxy = max(m_maxy, v_maxy)

    norm_x1 = (m_minx - minx) / (maxx - minx)
    norm_x2 = (m_maxx - minx) / (maxx - minx)

    norm_y1 = (m_miny - miny) / (maxy - miny)
    norm_y2 = (m_maxy - miny) / (maxy - miny)

    return norm_x1, norm_x2, norm_y1, norm_y2


def _get_minmax_coordinates_mesh(mesh, axis=0):
    """ Return the minimum and maximum coordinates along axis

    parameter:
    ----------
        axis:
            axis

    returns:
    -------
        tuple: minV, maxV

    """
    maxVal = np.zeros((1))
    minVal = np.zeros((1))
    maxVal[0] = mesh.data[:, axis].max()
    minVal[0] = mesh.data[:, axis].min()

    uw.barrier
    comm.Allreduce(MPI.IN_PLACE, maxVal, op=MPI.MAX)
    comm.Allreduce(MPI.IN_PLACE, minVal, op=MPI.MIN)
    uw.barrier

    return minVal, maxVal


class Plots(object):

    def __init__(self, Model):

        self.Model = Model

    @property
    def _boundingBox(self):
        minCoord = tuple([nd(val) for val in self.Model.minCoord])
        maxCoord = tuple([nd(val) for val in self.Model.maxCoord])
        boundingBox = (minCoord, maxCoord)
        return boundingBox

    def materialField(self, figsize=None, title="Material Field",
                      colours=None, script=None, cullface=False,
                      mask=None, visugrid=None, onMesh=False,
                      tracers=[], show=True, store=None, **kwargs):

        Fig = glucifer.Figure(store=store, figsize=figsize)
        Fig["title"] = title
        Fig["boundingBox"] = self._boundingBox

        if onMesh:
            Fig.Surface(self.Model.mesh, self.Model.projMaterialField,
                        **kwargs)
        else:
            pts = Fig.Points(self.Model.swarm,
                             fn_colour=self.Model.materialField,
                             colours=colours,
                             cullface=cullface, **kwargs)
            # pts.colourBar["binlabels"] = True
            # pts.colourBar["align"] = "right"
            # pts.colourBar["vertical"] = True
            pts.colourBar.colourMap["discrete"] = True

        if visugrid:
            Fig.Mesh(visugrid.mesh)

        Fig.script(script)
        if show:
            #Fig.viewer().window()
            Fig.show()

        return Fig

    material = materialField

    def viscosityField(self, figsize=None, title="Viscosity Field",
                       units=u.pascal * u.second, logScale=True,
                       projected=False, cullface=False,
                       script=None, show=True,
                       store=None, **kwargs):

        Fig = glucifer.Figure(store=store, figsize=figsize)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = Dimensionalize(1.0, units).magnitude
        Fig.Points(self.Model.swarm, self.Model._viscosityFn * fact,
                   logScale=logScale, cullface=cullface, **kwargs)

        Fig.script(script)
        if show:
            Fig.show()

        return Fig

    viscosity = viscosityField

    def strainRate(self, figsize=None, title="Strain Rate Field",
                   units=1.0 / u.second,
                   cullface=False, onMesh=True,
                   logScale=True, colours="coolwarm",
                   script=None, show=True,
                   store=None, **kwargs):

        Fig = glucifer.Figure(store=store, figsize=figsize)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh.subMesh,
                    self.Model._strainRate_2ndInvariant * fact,
                    cullface=cullface, logScale=logScale,
                    onMesh=onMesh, colours=colours,
                    **kwargs)

        Fig.script(script)
        if show:
            Fig.show()

        return Fig

    def density(self, figsize=None, title="Density Field",
                units=u.kilogram / u.metre**3,
                script=None, cullface=False, show=True,
                store=None, **kwargs):

        Fig = glucifer.Figure(store=store, figsize=figsize)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = Dimensionalize(1.0, units).magnitude
        Fig.Points(self.Model.swarm,
                   fn_colour=self.Model._densityFn * fact,
                   cullface=cullface, **kwargs)

        Fig.script(script)
        if show:
            Fig.show()

        return Fig

    def temperature(self, figsize=None, title="Temperature Field",
                    units=u.degK, script=None, cullface=False,
                    colours="coolwarm", show=True,
                    store=None, **kwargs):

        Fig = glucifer.Figure(store=store, figsize=figsize)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh, self.Model.temperature * fact,
                    colours=colours, cullface=cullface,
                    **kwargs)

        Fig.script(script)
        if show:
            Fig.show()

        return Fig

    def pressureField(self, figsize=None, title="Pressure Field",
                      units=u.pascal, cullface=False,
                      script=None, show=True,
                      store=None, **kwargs):

        Fig = glucifer.Figure(store=store, figsize=figsize)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh, self.Model.pressureField * fact,
                    cullface=cullface, **kwargs)

        Fig.script(script)
        if show:
            Fig.show()

        return Fig

    def velocityField(self, figsize=None, title="Velocity Field",
                      units=u.centimeter / u.year, cullface=False,
                      script=None, show=True, store=None, **kwargs):

        Fig = glucifer.Figure(store=store, figsize=figsize)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = Dimensionalize(1.0, units).magnitude
        velmagfield = uw.function.math.sqrt(
            uw.function.math.dot(self.Model.velocityField,
                                 self.Model.velocityField))
        Fig.Surface(self.Model.mesh, velmagfield * fact,
                    cullface=cullface, **kwargs)
        Fig.VectorArrows(self.Model.mesh, self.Model.velocityField,
                         **kwargs)

        Fig.script(script)
        if show:
            Fig.show()

        return Fig

    def plasticStrain(self, figsize=None, title="Plastic Strain",
                      cullface=False, script=None, show=True,
                      store=None, **kwargs):

        Fig = glucifer.Figure(store=store, figsize=figsize)
        Fig["title"] = title
        Fig["boundingBox"] = self._boundingBox

        Fig.Points(self.Model.swarm, fn_colour=self.Model.plasticStrain,
                   cullface=cullface, **kwargs)

        Fig.script(script)
        if show:
            Fig.show()

        return Fig

    plastic_strain = plasticStrain

    def melt_fraction(self, figsize=None, title="Melt fraction",
                      cullface=False, script=None, show=True,
                      store=None, **kwargs):

        Fig = glucifer.Figure(store=store, figsize=figsize)
        Fig["title"] = title
        Fig["boundingBox"] = self._boundingBox

        Fig.Points(self.Model.swarm, fn_colour=self.Model.meltField,
                   cullface=cullface, **kwargs)

        Fig.script(script)

        if show:
            Fig.show()

        return Fig
