# Configuration file
from scaling import u
from _validate import *


rcParams = {

"model.name": ["Model", validate_string],
"output.directory": ["outputs", validate_path],
"element.type" : ["Q1/dQ0", validate_string],

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
"swarm.particles.per.cell": [20, validate_int],

"popcontrol.aggressive" : [True, validate_bool],
"popcontrol.split.threshold" : [0.15, validate_float],
"popcontrol.max.splits" : [10, validate_int],
"popcontrol.particles.per.cell" : [20, validate_int],

"solver" : ["mg", validate_solver],
"penalty" : [0.0, validate_float],
"nonlinear.tolerance": [1e-2, validate_float],
"maximum.timestep" : [200000, validate_int],
"nonlinear.min.iterations": [3, validate_int],
"nonlinear.max.iterations": [500, validate_int],

"rheology.default.uppercrust": ["Patterson et al., 1990", validate_viscosity],
"rheology.default.midcrust": ["Patterson et al., 1990", validate_viscosity],
"rheology.default.lowercrust": ["Wang et al., 2012", validate_viscosity],
"rheology.default.mantlelithosphere": ["Hirth et al., 2003", validate_viscosity],
"rheology.default.mantle": ["Karato and Wu, 1990", validate_viscosity],

"scaling.length": [1.0 * u.meter, validate_quantity],
"scaling.mass": [1.0 * u.kilogram, validate_quantity],
"scaling.time": [1.0 * u.second, validate_quantity],
"scaling.temperature": [1.0 * u.degK, validate_quantity],
"scaling.substance": [1.0 * u.mole, validate_quantity],

# "velocity.SIunits" : [u.centimeter / u.year, validate_quantity],
# "temperature.SIunits"   : [u.degK, validate_quantity],
# "pressure.SIunits" : [u.pascal, validate_quantity],
# "strainrate.SIunits" : [1.0/u.second, validate_quantity],
# "viscosity.SIunits": [u.pascal * u.second, validate_quantity],
# "density.SIunits": [u.kilogram / u.meter**3, validate_quantity],


"viscosityField.SIunits" : [u.pascal * u.second, validate_quantity],
"densityField.SIunits" : [u.kilogram / u.metre**3, validate_quantity],
"velocityField.SIunits" : [u.centimeter / u.year, validate_quantity],
"temperature.SIunits" : [u.degK, validate_quantity],
"pressureField.SIunits" : [u.pascal , validate_quantity],
"strainRateField.SIunits" : [1.0 / u.second, validate_quantity],
"projViscosityField.SIunits"  : [u.pascal * u.second, validate_quantity],
"projDensityField.Siunits" : [u.kilogram / u.metre**3, validate_quantity]
        }

