# Configuration file

from scaling import u

rcParams = {

"output.directory": "outputs",
"default.solver" : "mg",

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
