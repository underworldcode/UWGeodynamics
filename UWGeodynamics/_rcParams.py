# Configuration file

from scaling import u
import six
try:
    import collections.abc as abc
except ImportError:
    # python 2
    import collections as abc

def validate_gravity(s):
    return u.check('[length]')(s)

def _listify_validator(scalar_validator, allow_stringlist=False):
    def f(s):
        if isinstance(s, six.string_types):
            try:
                return [scalar_validator(v.strip()) for v in s.split(',')
                        if v.strip()]
            except Exception:
                if allow_stringlist:
                    # Sometimes, a list of colors might be a single string
                    # of single-letter colornames. So give that a shot.
                    return [scalar_validator(v.strip()) for v in s if v.strip()]
                else:
                    raise
        # We should allow any generic sequence type, including generators,
        # Numpy ndarrays, and pandas data structures.  However, unordered
        # sequences, such as sets, should be allowed but discouraged unless the
        # user desires pseudorandom behavior.
        elif isinstance(s, abc.Iterable) and not isinstance(s, abc.Mapping):
            # The condition on this list comprehension will preserve the
            # behavior of filtering out any empty strings (behavior was
            # from the original validate_stringlist()), while allowing
            # any non-string/text scalar values such as numbers and arrays.
            return [scalar_validator(v) for v in s
                    if not isinstance(v, six.string_types) or v]
        else:
            msg = "{0!r} must be of type: string or non-dictionary iterable.".format(s)
            raise ValueError(msg)
    f.__doc__ = scalar_validator.__doc__
    return f

def validate_quantity(s):
    # Convert to quantity
    s = u.Quantity(s)
    if s.dimensionless:
        return s.magnitude
    return s

def validate_float(s):
    try:
        return float(s)
    except:
        raise ValueError("Could not convert value to float") 

def validate_solver(s):
    if s in ["mg"]:
        return s
    else:
        raise ValueError("Wrong solver option")

def validate_int(s):
    try:
        return int(s)
    except:
        raise ValueError("Could not convert value to int") 

def validate_path(s):
    return s

def validate_bool(b):
    """Convert b to a boolean or raise"""
    if isinstance(b, six.string_types):
        b = b.lower()
    if b in ('t', 'y', 'yes', 'on', 'true', '1', 1, True):
        return True
    elif b in ('f', 'n', 'no', 'off', 'false', '0', 0, False):
        return False
    else:
        raise ValueError('Could not convert "%s" to boolean' % b)

validate_stringlist = _listify_validator(six.text_type)
validate_stringlist.__doc__ = 'return a list'

rcParams = {
"output.directory": ["outputs", validate_path],

"minimum.viscosity": [1e19 * u.pascal * u.second, validate_quantity],
"maximum.viscosity": [1e25 * u.pascal * u.second, validate_quantity],

"swarm.variables" : [["materialField", "plasticStrain", "viscosityField", "densityField"], validate_stringlist],
"mesh.variables" :  [["velocityField", "temperature", "pressureField",
                     "strainRateField", "projMaterialField",
                     "projViscosityField", 
                     "projPlasticStrain",
                     "projDensityField"], validate_stringlist],
"default.outputs" : [["materialField", "temperature", "pressureField", "plasticStrain", "velocityField"], validate_stringlist], 

"gravity": [9.81 * u.meter / u.second**2, validate_quantity],
"swarm.particles.per.cell": [50, validate_int],

"popcontrol.aggressive" : [True, validate_bool],
"popcontrol.split.threshold" : [0.15, validate_float],
"popcontrol.max.splits" : [10, validate_int],
"popcontrol.particles.per.cell" : [50, validate_int],

"solver" : ["mg", validate_solver],
"penalty" : [0.0, validate_float],
"nonlinear.tolerance": [1e-2, validate_float],
"maximum.timestep" : [200000, validate_int],
"nonlinear.min.iterations": [3, validate_int],
"nonlinear.max.iterations": [500, validate_int],

"scaling.length": [1.0 * u.meter, validate_quantity],
"scaling.mass": [1.0 * u.kilogram, validate_quantity],
"scaling.time": [1.0 * u.second, validate_quantity],
"scaling.temperature": [1.0 * u.degK, validate_quantity],
"scaling.substance": [1.0 * u.mole, validate_quantity],

"velocity.units" : [u.centimeter / u.year, validate_quantity],
"temperature.units"   : [u.degK, validate_quantity],
"pressure.units" : [u.pascal, validate_quantity],
"strain.rate" : [1.0/u.second, validate_quantity],
"viscosity.units": [u.pascal * u.second, validate_quantity],
"density.units": [u.kilogram / u.meter**3, validate_quantity],


"viscosityField" : [u.pascal * u.second, validate_quantity],
"densityField" : [u.kilogram / u.metre**3, validate_quantity],
"velocityField" : [u.centimeter / u.year, validate_quantity],
"temperature" : [u.degK, validate_quantity],
"pressureField" : [u.pascal , validate_quantity],
"strainRateField" : [1.0 / u.second, validate_quantity],
"projViscosityField"  : [u.pascal * u.second, validate_quantity],
"projDensityField" : [u.kilogram / u.metre**3, validate_quantity]
        }

