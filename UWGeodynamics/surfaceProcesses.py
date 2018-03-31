try:
    from .linkage import SPM
except:
    pass

import underworld.function as fn
from scaling import nonDimensionalize as nd


class Badlands(object):
    """ A wrapper class for Badlands Linkage"""

    def __init__(self, airIndex,
                 sedimentIndex, XML, resolution, checkpoint_interval,
                 surfElevation=0., verbose=True, Model=None, restartFolder=None,
                 restartStep=None, timeField=None):

        self.airIndex = airIndex
        self.sedimentIndex = sedimentIndex
        self.XML = XML
        self.resolution = resolution
        self.checkpoint_interval = checkpoint_interval
        self.surfElevation = surfElevation
        self.verbose = verbose
        self.restartFolder = restartFolder
        self.restartStep = restartStep

        self.timeField = timeField
        self.Model = Model
        return

    @property
    def Model(self):
        return self._Model

    @Model.setter
    def Model(self, value):
        self._Model = value
        if value:
            self._init_Badlands()

    def _init_Badlands(self):
        self.mesh = self._Model.mesh
        self.velocityField = self._Model.velocityField
        self.swarm = self._Model.swarm
        self.materialField = self._Model.materialField

        self._BadlandsModel = SPM(self.mesh, self.velocityField, self.swarm,
                                 self.materialField, self.airIndex, self.sedimentIndex,
                                 self.XML, nd(self.resolution),
                                 nd(self.checkpoint_interval),
                                 nd(self.surfElevation), self.verbose,
                                 self.restartFolder, self.restartStep)
        return

    def solve(self, dt):
        self._BadlandsModel.solve(dt)
        return


class ErosionThreshold(object):

    def __init__(self, swarm=None, materialIndexField=None,
                 air=None, sediment=None, threshold=None):

        self.materialIndexField = materialIndexField
        self.swarm = swarm
        self.threshold = threshold
        self.air = air
        self.sediment = sediment


        materialMap = {}
        for material in air:
            materialMap[material.index] = 1.0

        isAirMaterial = fn.branching.map(fn_key=materialIndexField, mapping=materialMap, fn_default=0.0)

        belowthreshold = [(((isAirMaterial < 0.5) & (fn.input()[1] > nd(threshold))), air[0].index),
                         (True, materialIndexField)]

        self._fn = fn.branching.conditional(belowthreshold)

    def solve(self, dt):
        self.materialIndexField.data[:] = self._fn.evaluate(self.swarm)
        return


class SedimentationThreshold(object):

    def __init__(self, swarm=None, materialIndexField=None,
                 air=None, sediment=None, threshold=None, timeField=None):

        self.materialIndexField = materialIndexField
        self.timeField = timeField
        self.swarm = swarm
        self.air = air
        self.sediment = sediment
        self.threshold = threshold

        materialMap = {}
        for material in air:
            materialMap[material.index] = 1.0

        isAirMaterial = fn.branching.map(fn_key=materialIndexField, mapping=materialMap, fn_default=0.0)

        belowthreshold = [(((isAirMaterial > 0.5) & (fn.input()[1] < nd(threshold))), 0.),
                         (True, 1.)]

        self._change_material = fn.branching.conditional(belowthreshold)

        conditions = [(self._change_material < 0.5, sediment[0].index),
                      (True, materialIndexField)]

        self._fn = fn.branching.conditional(conditions)

    def solve(self, dt):

        if self.timeField:
            fn = self._change_material * self.timeField
            self.timeField.data[...] = fn.evaluate(self.swarm)

        self.materialIndexField.data[:] = self._fn.evaluate(self.swarm)

        return


class ErosionAndSedimentationThreshold(object):

    def __init__(self, swarm=None, materialIndexField=None,
                 air=None, sediment=None, threshold=None, timeField=None):

        self.materialIndexField = materialIndexField
        self.timeField = timeField
        self.swarm = swarm
        self.air = air
        self.sediment = sediment
        self.threshold = threshold

        materialMap = {}
        for material in air:
            materialMap[material.index] = 1.0

        isAirMaterial = fn.branching.map(fn_key=materialIndexField, mapping=materialMap, fn_default=0.0)

        sedimentation = [(((isAirMaterial > 0.5) & (fn.input()[1] < nd(threshold))), 0.),
                         (True, 1.0)]

        self._sedimented = fn.branching.conditional(sedimentation)

        conditions = [(self._sedimented < 0.5, sediment[0].index),
                      (True, materialIndexField)]

        erosion = [(((isAirMaterial < 0.5) & (fn.input()[1] > nd(threshold))), sediment[0].index),
                         (True, materialIndexField)]

        self._fn1 = fn.branching.conditional(conditions)
        self._fn2 = fn.branching.conditional(erosion)

    def solve(self, dt):

        if self.timeField:
            fn = self._sedimented * self.timeField
            self.timeField.data[...] = fn.evaluate(self.swarm)

        self.materialIndexField.data[:] = self._fn1.evaluate(self.swarm)
        self.materialIndexField.data[:] = self._fn2.evaluate(self.swarm)

        return
