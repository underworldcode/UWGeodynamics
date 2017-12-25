try:
    import underworld
except ImportError:
    raise ImportError("Can not find Underworld, please check your installation")

import shapes
import surfaceProcesses
from .scaling import COEFFICIENTS as scaling_coefficients
from .scaling import UnitRegistry
from .scaling import nonDimensionalize
from .scaling import Dimensionalize
from .LecodeIsostasy import LecodeIsostasy
from .lithopress import LithostaticPressure
from ._rheology import Rheology, ConstantViscosity, ViscousCreep, DruckerPrager
from ._rheology import VonMises
from ._rheology import ViscousCreepRegistry, PlasticityRegistry
from ._material import Material, MaterialRegistry
from ._model import Model
from ._melt import Solidus, Liquidus, SolidusRegistry, LiquidusRegistry
from ._rcParams import rcParams as config

nd = nonDimensionalize
u = UnitRegistry

rheologies = ViscousCreepRegistry()
yieldCriteria = PlasticityRegistry()
materials = MaterialRegistry()
scaling = scaling_coefficients

