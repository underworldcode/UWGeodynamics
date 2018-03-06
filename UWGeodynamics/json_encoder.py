import json
import inspect
from collections import OrderedDict
import UWGeodynamics as GEO
from .scaling import UnitRegistry as u
from ._frictional_boundary import FrictionBoundaries
from ._thermal_boundaries import TemperatureBCs
from ._velocity_boundaries import VelocityBCs
from ._visugrid import Visugrid
from json import JSONDecoder


class ObjectEncoder(json.JSONEncoder):

    def default(self, obj):

        if isinstance(obj, GEO.Model):
            
            model = OrderedDict()

            model["Type"] = "Model"
            model["name"] = obj["name"]

            # Save rcparams (use file)
            with open(GEO.uwgeodynamics_fname(), "r") as f:
                rcParams = f.read()

            model["rcParams"] = rcParams

            model["elementRes"] = obj["elementRes"]
            model["minCoord"] = obj["minCoord"]
            model["maxCoord"] = obj["maxCoord"]
            model["gravity"] = obj["gravity"]
            model["periodic"] = obj["periodic"]
            model["elementType"] = obj["elementType"]

            # Encode Scaling
            scaling = {}
            for key, val in GEO.scaling_coefficients.iteritems():
                scaling[key] = val
            model["scaling"] = scaling

            # Encode velocity boundary conditions
            model["velocityBCs"] = obj.velocityBCs

            # Encode temperature boundary conditions
            model["temperatureBCs"] = obj.temperatureBCs

            model["materials"] = []
            # Encode materials
            for material in reversed(obj.materials):
                if material is obj:
                    #val = super(Model, obj).to_json()
                    #model["materials"].append(val)
                    continue
                else:
                    model["materials"].append(material)

            return model

        # Material
        if isinstance(obj, GEO.Material):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        # Constant Density
        if isinstance(obj, GEO.ConstantDensity):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        # Linear Density
        if isinstance(obj, GEO.LinearDensity):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        # Boundary conditions
        # Frictional Boundary
        if isinstance(obj, FrictionBoundaries):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        # Solidus, Liquidus
        if isinstance(obj, (GEO.Solidus, GEO.Liquidus)):
            d = OrderedDict()
            return d

        # Constant Viscosity
        if isinstance(obj, GEO.ConstantViscosity):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        # Drucker Prager
        if isinstance(obj, GEO.DruckerPrager):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        # Von Mises
        if isinstance(obj, GEO.VonMises):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        # Viscous Creep
        if isinstance(obj, GEO.ViscousCreep):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, TemperatureBCs):
            d = OrderedDict()
            d["Type"] = "Temperature Boundary Conditions"

            attributes = [
                "left",
                "right",
                "top",
                "bottom",
                "back",
                "front",
                "indexSets"]

            for attribute in attributes:
                d[attribute] = obj[attribute]

            return d

        if isinstance(obj, GEO.Balanced_InflowOutflow):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, VelocityBCs):
            d = OrderedDict()
            d["Type"] = "Velocity Boundary Conditions"

            attributes = [
                "left",
                "right",
                "top",
                "bottom",
                "back",
                "front",
                "indexSets"]

            for attribute in attributes:
                d[attribute] = obj[attribute]

            return d

        if isinstance(obj, Visugrid):
            d = OrderedDict()
            d["Type"] = "Visugrid"
            d["minCoord"] = obj.minCoord
            d["maxCoord"] = obj.maxCoord
            d["elementRes"] = obj.elementRes
            return d

        # Shapes
        if isinstance(obj, GEO.shapes.Polygon):
            d = OrderedDict()
            d["Type"] = "Polygon"
            d["vertices"] = obj.vertices
            return d

        if isinstance(obj, GEO.shapes.MultiShape):
            d = OrderedDict()
            d["Type"] = "Multishape"
            d["shapes"] = obj.shapes
            return d

        if isinstance(obj, GEO.shapes.Layer):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, GEO.shapes.Box):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, GEO.shapes.Disk):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, GEO.shapes.Annulus):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, GEO.surfaceProcesses.Badlands):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, GEO.surfaceProcesses.ErosionThreshold):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, GEO.surfaceProcesses.SedimentationThreshold):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, GEO.surfaceProcesses.ErosionAndSedimentationThreshold):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, GEO.LecodeIsostasy):
            args, _, _, _ = inspect.getargspec(obj.__init__)
            d = OrderedDict()
            d["Type"] = type(obj).__name__
            for arg in args[1:]:
                try:
                    d[arg] = getattr(obj, arg)
                except AttributeError:
                    d[arg] = getattr(obj, "_" + arg)
            return d

        if isinstance(obj, GEO._mesh_advector._mesh_advector):
            d = OrderedDict()
            d["Type"] = "_mesh_advector"
            d["axis"] = obj.axis
            return d

        if isinstance(obj, GEO.u.Quantity):
            d = OrderedDict()
            d["value"] = obj.magnitude
            d["units"] = obj.units
            return d

        if isinstance(obj, GEO.u.Unit):
            return str(obj)

        return super(ObjectEncoder, self).default(obj)


class ObjectDecoder(json.JSONDecoder):

    def __init__(self, *args, **kwargs):
        JSONDecoder.__init__(self, object_hook=self.json_to_model,
                             *args, **kwargs)

    def json_to_model(self, d):

        if "Type" in d and d["Type"] == "Model":
            model = GEO.Model
            args, _, _, _ = inspect.getargspec(model.__init__)
            kwargs = {}
            for arg in args:
                if arg != "self":
                    kwargs[arg] = self.json_to_model(d[arg])

            return model(**kwargs)

        if "units" in d:
            value = d["value"]
            units = d["units"]
            return GEO.u.Quantity(value, units)




