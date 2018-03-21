from itertools import count
from copy import copy
from collections import OrderedDict
import json
import pkg_resources
from .scaling import u
from ._rheology import ConstantViscosity
from ._density import ConstantDensity


class Material(object):
    _ids = count(0)

    def __init__(self, name="Undefined", density=0.0,
                 diffusivity=None, capacity=None,
                 radiogenicHeatProd=0.0, shape=None, viscosity=None,
                 plasticity=None, elasticity=None, solidus=None, liquidus=None,
                 minViscosity=None, maxViscosity=None, stressLimiter=None,
                 latentHeatFusion=0.0, meltExpansion=0.0, meltFraction=0.0,
                 meltFractionLimit=1.0, viscosityChangeX1=0.15,
                 viscosityChangeX2=0.3, viscosityChange=1.0):

        self.index = self._ids.next()

        self.name = name
        self.top = None
        self.bottom = None

        self.shape = shape
        self._density = None
        self.density = density
        self.referenceDensity = self.density
        self.diffusivity = diffusivity
        self.capacity = capacity
        self._thermalExpansivity = None
        self.radiogenicHeatProd = radiogenicHeatProd

        self.minViscosity = minViscosity
        self.maxViscosity = maxViscosity
        self.stressLimiter = stressLimiter

        self.compressibility = None

        self.solidus = solidus
        self.liquidus = liquidus
        self.latentHeatFusion = latentHeatFusion
        self.meltFraction = meltFraction
        self.meltFractionLimit = meltFractionLimit
        self.meltExpansion = meltExpansion
        self.viscosityChangeX1 = viscosityChangeX1
        self.viscosityChangeX2 = viscosityChangeX2
        self.viscosityChange = viscosityChange
        self.melt = False

        self._viscosity = None
        self.viscosity = viscosity
        self.plasticity = plasticity

        self.elasticity = elasticity

    def _repr_html_(self):
        return _material_html_repr(self)

    def __getitem__(self, name):
        return self.__dict__[name]

    @property
    def viscosity(self):
        return self._viscosity

    @viscosity.setter
    def viscosity(self, value):
        if isinstance(value, (u.Quantity, float)):
            self._viscosity = ConstantViscosity(value)
        else:
            self._viscosity = value

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        if isinstance(value, (u.Quantity, float)):
            self._density = ConstantDensity(value)
        else:
            self._density = value
            self._thermalExpansivity = value.thermalExpansivity

    @property
    def thermalExpansivity(self):
        return self._thermalExpansivity

    def add_melt_modifier(self, solidus, liquidus, latentHeatFusion,
                          meltExpansion,
                          meltFraction=0.0,
                          meltFractionLimit=1.0,
                          viscosityChangeX1=0.15,
                          viscosityChangeX2=0.3,
                          viscosityChange=1e3):

        self.solidus = solidus
        self.liquidus = liquidus
        self.latentHeatFusion = latentHeatFusion
        self.meltFraction = meltFraction
        self.meltFractionLimit = meltFractionLimit
        self.meltExpansion = meltExpansion
        self.viscosityChangeX1 = viscosityChangeX1
        self.viscosityChangeX2 = viscosityChangeX2
        self.viscosityChange = viscosityChange
        self.melt = True


_default = OrderedDict()
_default["Radiogenic Heat Production"] = "radiogenicHeatProd"
_default["Diffusivity"] = "diffusivity"
_default["Capacity"] = "capacity"
_default["Min Viscosity Limit"] = "minViscosity"
_default["Max Viscosity Limit"] = "maxViscosity"

_melt = OrderedDict()
_melt["Solidus"] = ""
_melt["Liquidus"] = ""
_melt["Latent Heat Fusion"] = "latentHeatFusion"
_melt["Melt Expansion"] = "meltExpansion"
_melt["Melt Fraction"] = "meltFraction"
_melt["Melt Fraction Limit"] = "meltFractionLimit"
_melt["Viscosity Change"] = "viscosityChange"
_melt["Viscosity Change X1"] = "viscosityChangeX1"
_melt["Viscosity Change X2"] = "viscosityChangeX2"


def _material_html_repr(Material):
    header = "<table>"
    footer = "</table>"
    html = """
    <tr>
      <th colspan="2">{0}</th>
    </tr>""".format(Material.name)

    if Material.density:
        key = "Density"
        value = str(Material.density.name)
        html += "<tr><td>{0}</td><td>{1}</td></tr>".format(key, value)

    for key, val in _default.iteritems():
        value = Material.__dict__.get(val)
        if value is None:
            value = Material.__dict__.get("_"+val)
        html += "<tr><td>{0}</td><td>{1}</td></tr>".format(key, value)

    filename = None

    if Material.viscosity and Material.plasticity:
        type_ = "(Visco-plastic)"
    elif Material.viscosity:
        type_ = "(Viscous)"
    elif Material.plasticity:
        type_ = "(Plastic)"
    else:
        type_ = ""

    html += """
    <tr>
      <th colspan="2">Rheology {}</th>
    </tr>""".format(type_)

    if Material.viscosity:
        html += "<tr><td>{0}</td><td>{1}</td></tr>".format(
            "Viscosity", Material.viscosity.name)
    if Material.plasticity:
        html += "<tr><td>{0}</td><td>{1}</td></tr>".format(
            "Plasticity", Material.plasticity.name)

    return header + html + footer

class MaterialRegistry(object):

    def __init__(self, filename=None):

        if not filename:
            filename = pkg_resources.resource_filename(__name__, "ressources/Materials.json")

        with open(filename, "r") as infile:
            _materials = json.load(infile)

        self._dir = {}
        for material, parameters in _materials.iteritems():
            name = material.replace(" ", "_").replace(",", "").replace(".", "")
            name = name.replace(")", "").replace("(", "")

            for key in parameters.keys():
                if parameters[key] is not None:
                    value = parameters[key]["value"]
                    units = parameters[key]["units"]
                    if units != "None":
                        parameters[key] = u.Quantity(value, units)
                    else:
                        parameters[key] = value

            self._dir[name] = Material(name=material, **parameters)

    def __dir__(self):
        # Make all the rheology available through autocompletion
        return list(self._dir.keys())

    def __getattr__(self, item):
        # Make sure to return a new instance of ViscousCreep
        return copy(self._dir[item])
