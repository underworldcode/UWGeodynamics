# Configuration file

from scaling import u

rcParams = {

"gravity": 9.81 * u.meter / u.second**2,
"particles.per.cell": 20,
"courant.factor": "unused",
"maximum.time.step": 30000 * u.years,



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

"default.length.scaling": 1.0 * u.meter,
"default.mass.scaling": 1.0 * u.kilogram,
"default.time.scaling": 1.0 * u.second,
"default.temperature.scaling": 1.0 * u.degK,
"default.substance.scaling": 1.0 * u.mole,



"velocityField.SIunits" : u.centimeter / u.year,
"temperature.SIunits"   : u.degK,
"pressureField.SIunits" : u.pascal,
"strainRateField.SIunits"    : 1.0/u.second,
"materialField.SIunits" : None,
"plasticStrain.SIunits" : None,
"viscosityField.SIunits": u.pascal * u.second,
"densityField.SIunits": u.kilogram / u.meter**3,

"projMaterialField.SIunits" : None,
"projViscosityField.SIunits" : u.pascal * u.second,
"projPlasticStrain.SIunits" : None,
"projDensityField.SIunits" : u.kilogram / u.meter**3
        }
