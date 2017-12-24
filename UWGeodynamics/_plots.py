import os
import glucifer
import scaling as sca
from scaling import nonDimensionalize as nd
import matplotlib.pyplot as plt
import numpy as np

u = UnitRegistry = sca.UnitRegistry

def _clean_local(locals_):
    del(locals_["self"])
    kwargs = locals_.pop("kwargs")
    for key, val in kwargs.iteritems():
        locals_[key] = val
    return locals_


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
        boundingBox =(minCoord, maxCoord)
        return boundingBox

    def material(self, figsize=None, title="Material Field",
                 colours=None, script=None, cullface=False,
                 mask=None, visugrid=None, onMesh=False,
                 tracers=[], **kwargs):

        saved_args = _clean_local(locals())
        Fig = glucifer.Figure(**saved_args)
        Fig["title"] = title
        Fig["boundingBox"] = self._boundingBox

        if onMesh:
            Fig.Surface(self.Model.mesh, self.Model.projMaterialField)
        else:
            pts = Fig.Points(self.Model.swarm,
                             fn_colour=self.Model.materialField,
                             **saved_args)
            #pts.colourBar["binlabels"] = True
            #pts.colourBar["align"] = "right"
            #pts.colourBar["vertical"] = True
            pts.colourBar.colourMap["discrete"] = True

        Fig.script(script)
        Fig.show()

        return Fig

    def viscosity(self, figsize=None, title="Viscosity Field",
                  units=u.pascal*u.second, logScale=True,
                  projected=False, cullface=False,
                  script=None, **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Points(self.Model.swarm, self.Model._viscosityFn*fact,
                   **saved_args)

        Fig.script(script)
        Fig.show()

        return Fig

    def strainRate(self, figsize=None, title="Strain Rate Field",
                   units=1.0/u.second,
                   cullface=False, onMesh=True,
                   logScale=True, colours="coolwarm",
                   script=None, **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh,
                    self.Model._strainRate_2ndInvariant*fact,
                    **saved_args)

        Fig.script(script)
        Fig.show()

        return Fig

    def density(self, figsize=None, title="Density Field",
                units=u.kilogram/u.metre**3,
                script=None, cullface=False, **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Points(self.Model.swarm,
                   fn_colour=self.Model._densityFn*fact,
                   **saved_args)

        Fig.script(script)
        Fig.show()

        return Fig

    def temperature(self, figsize=None, title="Temperature Field",
                    units=u.degK, script=None, cullface=False,
                    colourScale="coolwarm", **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh, self.Model.temperature*fact,
                    **saved_args)

        Fig.script(script)
        Fig.show()

        return Fig

    def pressureField(self, figsize=None, title="Pressure Field",
                      units=u.pascal, cullface=False,
                      script=None, **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh, self.Model.pressureField*fact,
                    **saved_args)

        Fig.script(script)
        Fig.show()

        return Fig

    def velocityField(self, figsize=None, title="Velocity Field",
                      units=u.centimeter/u.year, cullface=False,
                      script=None, **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        Fig["title"] = title + " " + str(units)
        Fig["boundingBox"] = self._boundingBox

        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh,self.Model.velocityField[0]*fact,
                   **saved_args)
        Fig.VectorArrows(self.Model.mesh, self.Model.velocityField,
                         **saved_args)

        Fig.script(script)
        Fig.show()

        return Fig

    def plastic_strain(self, figsize=None, title="Plastic Strain",
                       cullface=False, script=None, **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        Fig["title"] = title
        Fig["boundingBox"] = self._boundingBox

        Fig.Points(self.Model.swarm, fn_colour=self.Model.plasticStrain,
                   **saved_args)

        Fig.script(script)
        Fig.show()

        return Fig

    def melt_fraction(self, figsize=None, title="Melt fraction",
                      cullface=False, script=None, **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        Fig["title"] = title
        Fig["boundingBox"] = self._boundingBox

        Fig.Points(self.Model.swarm, fn_colour=self.Model.meltField,
                   **saved_args)

        Fig.script(script)
        Fig.show()

        return Fig
