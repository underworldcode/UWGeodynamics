from pint import UnitRegistry
from _utils import TransformedDict

u = UnitRegistry = UnitRegistry()

COEFFICIENTS = TransformedDict()
COEFFICIENTS["[length]"] = 1.0 * u.meter
COEFFICIENTS["[mass]"] = 1.0 * u.kilogram 
COEFFICIENTS["[time]"] = 1.0 * u.second
COEFFICIENTS["[temperature]"] = 1.0 * u.degK
COEFFICIENTS["[substance]"] = 1.0 * u.mole
