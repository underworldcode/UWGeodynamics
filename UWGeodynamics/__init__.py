try:
    import underworld as uw
except ImportError:
    raise ImportError("Can not find Underworld, please check your installation")


import scaling
import shapes
import utils
import surfaceProcesses as SPM
from scaling import UnitRegistry
from scaling import nonDimensionalize
from scaling import Dimensionalize
from LecodeIsostasy import LecodeIsostasy
from lithopress import LithostaticPressure
from rheology import Rheology, ConstantViscosity, ViscousCreep, DruckerPrager
from rheology import VonMises
from rheology import ViscousCreepRegistry, PlasticityRegistry
from Material import Material
from Model import Model
from Melt import Solidus, Liquidus, SolidusRegistry, LiquidusRegistry
from rcParams import rcParams as config

nd = nonDimensionalize
sca = scaling
u = UnitRegistry

rheologies = ViscousCreepRegistry()
yieldCriteria = PlasticityRegistry()

