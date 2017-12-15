from itertools import count
from scaling import u
from rheology import ConstantViscosity

class Material(object):
    _ids = count(0)

    def __init__(self, name="Undefined", vertices=None, density=None,
                 diffusivity=None, capacity=None, thermalExpansivity=None,
                 radiogenicHeatProd=0.0, shape=None, viscosity=None,
                 plasticity=None, solidus=None, liquidus=None,
                 latentHeatFusion=0.0, meltExpansion=0.0, meltFraction=0.0,
                 meltFractionLimit=1.0, viscosityChangeX1=0.15, 
                 viscosityChangeX2=0.3, viscosityChange=1.0):

        self.index = self._ids.next()
        
        self.name = name
        self.top = None
        self.bottom = None

        self.shape = shape
        self.density = density
        self.diffusivity = diffusivity
        self.capacity = capacity
        self.thermalExpansivity = thermalExpansivity
        self.radiogenicHeatProd = radiogenicHeatProd
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

        self.rheology = None  # For backward compatibility

        self._viscosity = viscosity
        self._plasticity = plasticity

    def __repr__(self):
        return self.name

    @property
    def shape(self):
        return self._shape

    @shape.setter
    def shape(self, value=None):
        self._shape = value

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        self._density = value

    @property
    def diffusivity(self):
        return self._diffusivity

    @diffusivity.setter
    def diffusivity(self, value):
        self._diffusivity = value

    @property
    def capacity(self):
        return self._capacity

    @capacity.setter
    def capacity(self, value):
        self._capacity = value

    @property
    def radiogenicHeatProd(self):
        return self._radiogenicHeatProd

    @radiogenicHeatProd.setter
    def radiogenicHeatProd(self, value):
        self._radiogenicHeatProd = value
    
    @property
    def thermalExpansivity(self):
        return self._thermalExpansivity

    @thermalExpansivity.setter
    def thermalExpansivity(self, value):
        self._thermalExpansivity = value
    
    @property
    def viscosity(self):
        return self._viscosity

    @viscosity.setter
    def viscosity(self, value):
        if isinstance(value, u.Quantity) or isinstance(value, float):
            self._viscosity = ConstantViscosity(value)
        else:
            self._viscosity = value

    @property
    def plasticity(self):
        return self._plasticity

    @plasticity.setter
    def plasticity(self, value):
        self._plasticity = value
    
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
