import os
import glucifer
import scaling as sca
from scaling import nonDimensionalize as nd
import matplotlib.pyplot as plt
import numpy as np

u = UnitRegistry = sca.UnitRegistry

default_parameters = {
        "figsize": (1200,400),
        "title": "",
        "store": None,
        "axis": True,
        "quality": 1
        }

def _clean_local(locals_):
    del(locals_["self"])
    kwargs = locals_.pop("kwargs")
    for key, val in kwargs.iteritems():
        locals_[key] = val
    return locals_

def visugrid_drawing_object(Model): 
    # TODO Use a class based on glucifer.objects
    visugrid = Model.visugrid.mesh

    xrange = visugrid.maxCoord[0] - visugrid.minCoord[0]
    yrange = visugrid.maxCoord[1] - visugrid.minCoord[1]
    dx = np.abs(nd(Model.minCoord[0]) - visugrid.minCoord[0])
    dy = np.abs(nd(Model.minCoord[1]) - visugrid.minCoord[1])
    xmin = (Model.mesh.minCoord[0] - Model.mesh.minCoord[0] + dx) / xrange
    xmax = (Model.mesh.maxCoord[0] - Model.mesh.minCoord[0] + dx) / xrange
    ymin = (Model.mesh.minCoord[1] - Model.mesh.minCoord[1] + dy) / yrange
    ymax = (Model.mesh.maxCoord[1] - Model.mesh.minCoord[1] + dy) / yrange
    xmin += 0.007
    xmax -= 0.007
    
    return glucifer.objects.Mesh(visugrid.mesh, xmin=xmin,
                                 xmax=xmax, ymin=ymin, ymax=ymax)

class Plots(object):

    def __init__(self, Model):

        self.Model = Model

    @property
    def _get_bounding_box(self):
        minCoord = tuple([nd(val) for val in self.Model.minCoord])
        maxCoord = tuple([nd(val) for val in self.Model.maxCoord])
        boundingBox =(minCoord, maxCoord)
        return boundingBox

    def material(self, figsize=None,
                 colours=None, script=None, cullface=False,
                 mask=None, visugrid=None, projected=False,
                 tracers=[], **kwargs):
        
        saved_args = _clean_local(locals())
        Fig = glucifer.Figure(**saved_args)

        if projected:
            pass
        else:
            Fig.Points(self.Model.swarm, fn_colour=self.Model.materialField,
                       **saved_args)
        
        Fig.script(script)
        Fig.show()

        return Fig

    def viscosity(self, figsize=None, units=u.pascal*u.second,
                       projected=False, cullface=False,
                       script=None, **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
       
        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Points(self.Model.swarm, self.Model._viscosityFn*fact,
                   **saved_args)
        
        Fig.script(script)
        Fig.show()

        return Fig
    
    def strainRate(self, figsize=None, units=1.0/u.second,
                   cullface=False, onMesh=True,
                   logScale=True, colours="coolwarm",
                   script=None, **kwargs):
        
        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        
        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh, self.Model._strainRate_2ndInvariant*fact,
                    **saved_args)
        Fig.script(script)
        Fig.show()
        return Fig

    def density(self, figsize=None, units=u.kilogram/u.metre**3,
                script=None, cullface=False, **kwargs):
        
        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)

        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Points(self.Model.swarm, fn_colour=self.Model._densityFn*fact,
                   **saved_args)
        Fig.script(script)
        Fig.show()
        return Fig

    def temperature(self, figsize=None, units=u.degK, script=None,
                    cullface=False, colourScale="coolwarm", **kwargs):
        
        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
       
        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh, self.Model.temperature*fact,
                    **saved_args)

        Fig.script(script)
        Fig.show()
        return Fig
    
    def pressureField(self, figsize=None, units=u.pascal,
                      cullface=False, script=None, **kwargs):
        
        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh, self.Model.pressureField*fact,
                    **saved_args)
        
        Fig.script()
        Fig.show()
        return Fig

    def velocityField(self, figsize=None, units=u.centimeter/u.year,
                      cullface=False, script=None, **kwargs):
        
        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)

        fact = sca.Dimensionalize(1.0, units).magnitude
        Fig.Surface(self.Model.mesh,self.Model.velocityField[0]*fact,
                   **saved_args)
        Fig.VectorArrows(self.Model.mesh, self.Model.velocityField,
                         **saved_args)
        
        Fig.script(script)
        Fig.show()
        return Fig
    
    def strain(self, figsize=None, cullface=False,
               script=None, **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        
        Fig.Points(self.Model.swarm, fn_colour=self.Model.plasticStrain,
                   **saved_args)
        
        Fig.script(script)
        Fig.show()
        return Fig
    
    def melt(self, figsize=None, cullface=False,
             script=None, **kwargs):

        saved_args = _clean_local(locals())
        units = saved_args.pop("units")
        Fig = glucifer.Figure(**saved_args)
        
        Fig.Points(self.Model.swarm, fn_colour=self.Model.meltField,
                   **saved_args)

        Fig.script(script)
        Fig.show()
        return Fig
