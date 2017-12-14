

class strainWeakening(object):

    def __init__(self, swarm, plasticStrain,
                 healingRate, 
                 initialDamageFraction,
                 initialDamageWaveNumber,
                 initialDamageWaveNumberSinI,
                 initialDamageWaveNumberCosI,
                 initialDamageWaveNumberSinJ,
                 initialDamageWaveNumberCosJ,
                 initialDamageWaveNumberSinK,
                 initialDamageWaveNumberCosK,
                 initialDamageFactor,
                 initialStrainShape,
                 ):

        self._swarm = swarm
        self._plasticStrain = plasticStrain
        self.healingRate = healingRate
        self.initialDamageFraction = initialDamageFraction
        self.initialDamageWaveNumber = initialDamageWaveNumber
        self.initialDamageWaveNumberSinI = initialDamageWaveNumberSinI
        self.initialDamageWaveNumberCosI = initialDamageWaveNumberCosI
        self.initialDamageWaveNumberSinJ = initialDamageWaveNumberSinJ
        self.initialDamageWaveNumberCosJ = initialDamageWaveNumberCosJ
        self.initialDamageWaveNumberSinK = initialDamageWaveNumberSinK
        self.initialDamageWaveNumberCosK = initialDamageWaveNumberCosK
        self.initialDamageFactor = initialDamageFactor
        self.initialStrainShape = initialStrainShape


    def strainWeakeningInitialize(self):
        pass

