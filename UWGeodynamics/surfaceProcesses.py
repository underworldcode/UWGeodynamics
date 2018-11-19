from __future__ import print_function,  absolute_import
import abc
import underworld.function as fn
from .scaling import nonDimensionalize as nd


ABC = abc.ABCMeta('ABC', (object,), {})


class SurfaceProcesses(ABC):

    def __init__(self, Model=None):

        self.Model = Model

    @property
    def Model(self):
        return self._Model

    @Model.setter
    def Model(self, value):
        self._Model = value
        if value:
            self._init_model()

    @abc.abstractmethod
    def _init_model(self):
        pass

    @abc.abstractmethod
    def solve(self, dt):
        pass


class Badlands(SurfaceProcesses):
    """ A wrapper class for Badlands Linkage"""

    def __init__(self, airIndex,
                 sedimentIndex, XML, resolution, checkpoint_interval,
                 surfElevation=0., verbose=True, Model=None, restartFolder=None,
                 restartStep=None, timeField=None, surfaceTracers=None):
        try:

            import pyBadlands

        except ImportError :
            raise ImportError("""pyBadlands import as failed. Please check your
                              installation, PYTHONPATH and PATH environment
                              variables""")

        self.airIndex = airIndex
        self.sedimentIndex = sedimentIndex
        self.XML = XML
        self.resolution = resolution
        self.checkpoint_interval = checkpoint_interval
        self.surfElevation = surfElevation
        self.verbose = verbose
        self.restartFolder = restartFolder
        self.restartStep = restartStep
        self.surfaceTracers = surfaceTracers

        self.timeField = timeField
        self.Model = Model

    def _init_model(self):

        from .linkage import SPM

        self.mesh = self._Model.mesh
        self.velocityField = self._Model.velocityField
        self.swarm = self._Model.swarm
        self.materialField = self._Model.materialField

        self._BadlandsModel = SPM(self.mesh, self.velocityField, self.swarm,
                                  self.materialField, self.airIndex,
                                  self.sedimentIndex,
                                  self.XML, nd(self.resolution),
                                  nd(self.checkpoint_interval),
                                  nd(self.surfElevation), self.verbose,
                                  self.restartFolder, self.restartStep)
        return

    def solve(self, dt):
        self._BadlandsModel.solve(dt)
        return


class ErosionThreshold(SurfaceProcesses):

    def __init__(self, air=None, threshold=None, surfaceTracers=None,
                 Model=None, **kwargs):

        super(ErosionThreshold, self).__init__(Model=Model)

        self.Model = Model
        self.threshold = threshold
        self.air = air
        self.surfaceTracers = surfaceTracers
        self.Model = Model

    def _init_model(self):

        materialField = self.Model.materialField

        materialMap = {}
        for material in self.air:
            materialMap[material.index] = 1.0

        isAirMaterial = fn.branching.map(fn_key=materialField,
                                         mapping=materialMap,
                                         fn_default=0.0)

        belowthreshold = [(((isAirMaterial < 0.5) & (fn.input()[1] > nd(self.threshold))), self.air[0].index),
                          (True, materialField)]

        self._fn = fn.branching.conditional(belowthreshold)

    def solve(self, dt):

        if not self.Model:
            raise ValueError("Model is not defined")

        self.Model.materialField.data[:] = self._fn.evaluate(self.Model.swarm)
        if self.surfaceTracers:
            if self.surfaceTracers.swarm.particleCoordinates.data.size > 0:
                coords = self.surfaceTracers.swarm.particleCoordinates
                coords.data[coords.data[:, -1] > nd(self.threshold), -1] = nd(self.threshold)
        return


class SedimentationThreshold(SurfaceProcesses):

    def __init__(self, air=None, sediment=None,
                 threshold=None, timeField=None, Model=None,
                 surfaceTracers=None, **kwargs):

        super(SedimentationThreshold, self).__init__(Model=Model)

        self.timeField = timeField
        self.air = air
        self.sediment = sediment
        self.threshold = threshold
        self.surfaceTracers = surfaceTracers
        self.Model = Model

    def _init_model(self):

        materialField = self.Model.materialField

        materialMap = {}
        for material in self.air:
            materialMap[material.index] = 1.0

        isAirMaterial = fn.branching.map(fn_key=materialField,
                                         mapping=materialMap,
                                         fn_default=0.0)

        belowthreshold = [(((isAirMaterial > 0.5) & (fn.input()[1] < nd(self.threshold))), 0.),
                          (True, 1.)]

        self._change_material = fn.branching.conditional(belowthreshold)

        conditions = [(self._change_material < 0.5, self.sediment[0].index),
                      (True, materialField)]

        self._fn = fn.branching.conditional(conditions)

    def solve(self, dt):

        if not self.Model:
            raise ValueError("Model is not defined")

        if self.timeField:
            fn = self._change_material * self.timeField
            self.timeField.data[...] = fn.evaluate(self.Model.swarm)

        self.Model.materialField.data[:] = self._fn.evaluate(self.Model.swarm)

        if self.surfaceTracers:
            if self.surfaceTracers.swarm.particleCoordinates.data.size > 0:
                coords = self.surfaceTracers.swarm.particleCoordinates
                coords.data[coords.data[:, -1] < nd(self.threshold), -1] = nd(self.threshold)


class ErosionAndSedimentationThreshold(SedimentationThreshold, ErosionThreshold):

    def __init__(self, air=None, sediment=None,
                 threshold=None, timeField=None,
                 surfaceTracers=None, Model=None, **kwargs):

        super(ErosionAndSedimentationThreshold, self).__init__(Model=Model)

        self.timeField = timeField
        self.air = air
        self.sediment = sediment
        self.threshold = threshold
        self.surfaceTracers = surfaceTracers
        self.Model = Model

    def _init_model(self):

        ErosionThreshold._init_model(self)
        SedimentationThreshold._init_model(self)

    def solve(self, dt):

        ErosionThreshold.solve(self, dt)
        SedimentationThreshold.solve(self, dt)
