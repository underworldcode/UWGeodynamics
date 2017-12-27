# Configuration file

from scaling import u


def validate_gravity(s):
    return u.check('[length]')(s)



rcParams = {

"gravity": 9.81 * u.meter / u.second**2,
"particles.per.cell": 20,
"courant.factor": "unused",
"maximum.time.step": 30000 * u.years,

"default.upper.crust.viscosity": None,
"default.lower.crust.viscosity": None,


"output.directory": "outputs",
"default.solver" : "mg",
"linear.tolerance": 1e-5,
"nonlinear.tolerance": 1e-2,
"nonlinear.min.iterations": 3,
"nonlinear.max.iterations": 500,
"maximum.timestep" : 200000,


"default.outputs" : ["materialField", "temperature", "pressureField", "plasticStrain", "velocityField"], 
"swarm.variables" : ["materialField", "plasticStrain", "viscosityField", "densityField"],
"mesh.variables" :  ["velocityField", "temperature", "pressureField",
                     "strainRateField", "projMaterialField",
                     "projViscosityField", 
                     "projPlasticStrain",
                     "projDensityField"],




"velocityField.SIunits" : u.centimeter / u.year,
"temperature.SIunits"   : u.degK,
"pressureField.SIunits" : u.pascal,
"strainRateField.SIunits" : 1.0/u.second,
"materialField.SIunits" : None,
"plasticStrain.SIunits" : None,
"viscosityField.SIunits": u.pascal * u.second,
"densityField.SIunits": u.kilogram / u.meter**3,

"projMaterialField.SIunits" : None,
"projViscosityField.SIunits" : u.pascal * u.second,
"projPlasticStrain.SIunits" : None,
"projDensityField.SIunits" : u.kilogram / u.meter**3
        }

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

rcParamsNew = {

"gravity": [9.81 * u.meter / u.second**2, validate_quantity],
"nonlinear.tolerance": [1e-2, validate_float],
"linear.tolerance": [1e-5, validate_float],

"scaling.length": [1.0 * u.meter, validate_quantity],
"scaling.mass": [1.0 * u.kilogram, validate_quantity],
"scaling.time": [1.0 * u.second, validate_quantity],
"scaling.temperature": [1.0 * u.degK, validate_quantity],
"scaling.substance": [1.0 * u.mole, validate_quantity]
        }

