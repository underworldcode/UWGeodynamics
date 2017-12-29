
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

    GEO.scaling_coefficients["[length]"] = KL
    GEO.scaling_coefficients["[time]"] = Kt
    GEO.scaling_coefficients["[mass]"]= KM

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

    GEO.scaling_coefficients["[length]"] = KL
    GEO.scaling_coefficients["[time]"] = Kt
    GEO.scaling_coefficients["[mass]"]= KM

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

def test_set_velocity_boundary_conditions():
    import UWGeodynamics as GEO
    u = GEO.u
    Model = GEO.Model(elementRes=(10, 10),
                      minCoord=(0. * u.kilometer, 0. * u.kilometer),
                      maxCoord=(10. * u.kilometer, 10. * u.kilometer))
    conditions = Model.set_velocityBCs(left=[1.0*u.centimetre/u.year, None],
                                       right=[-1.0*u.centimetre/u.year, None],
                                       bottom=[None, 0.],
                                       top=[None,0.])
    assert(len(conditions) == 1)

def test_set_velocity_boundary_conditions_in_3D():
    import UWGeodynamics as GEO
    u = GEO.u
    Model = GEO.Model(elementRes=(10, 10, 10),
                      minCoord=(0. * u.kilometer,
                                0. * u.kilometer,
                                0. *u.kilometer),
                      maxCoord=(10. * u.kilometer,
                                10. * u.kilometer,
                                10. * u.kilometer))

    conditions = Model.set_velocityBCs(
        left=[1.0*u.centimetre/u.year, None, 0.],
        right=[-1.0*u.centimetre/u.year, None, 0.],
        bottom=[None, None, 0.],
        top=[None, None, 0.],
        front=[None, 0., None],
        back=[None, 0., None])
    assert len(conditions) == 1
