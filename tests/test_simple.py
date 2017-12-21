
def test_module_loading():
    import UWGeodynamics as GEO

def test_unit_registry():
    import UWGeodynamics as GEO

    u = GEO.u

def test_scaling():
    import UWGeodynamics as GEO

    u = GEO.u
    velocity = 1.0 * u.centimeter / u.hour
    model_length = 2. * u.meter
    model_height = 1. * u.meter
    refViscosity = 1e6 * u.pascal * u.second
    bodyforce = 200 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

    KL = model_height
    Kt = KL / velocity
    KM = bodyforce * KL**2 * Kt**2

    GEO.scaling.scaling["[length]"] = KL
    GEO.scaling.scaling["[time]"] = Kt
    GEO.scaling.scaling["[mass]"]= KM

def test_model_creation():
    import UWGeodynamics as GEO

    u = GEO.u
    velocity = 1.0 * u.centimeter / u.hour
    model_length = 2. * u.meter
    model_height = 1. * u.meter
    refViscosity = 1e6 * u.pascal * u.second
    bodyforce = 200 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

    KL = model_height
    Kt = KL / velocity
    KM = bodyforce * KL**2 * Kt**2

    GEO.scaling.scaling["[length]"] = KL
    GEO.scaling.scaling["[time]"] = Kt
    GEO.scaling.scaling["[mass]"]= KM

    Model = GEO.Model(elementRes=(64, 64), 
                      minCoord=(-1. * u.meter, -50. * u.centimeter), 
                      maxCoord=(1. * u.meter, 50. *
                      u.centimeter))
    return Model
    
def test_adding_materials():
    import UWGeodynamics as GEO
    Model = test_model_creation()
    air = Model.add_material(
        name="Air", 
        shape=GEO.shapes.Layer(top=Model.top,
                               bottom=Model.bottom)
    )
    return air
