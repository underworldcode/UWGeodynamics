from .utils import _notebook_run

def test_steady_state_example():
    _notebook_run("examples/1_01_Steady_State_Heat.ipynb")


def test_convection_example():
    _notebook_run("examples/1_02_Convection_Example.ipynb")


def test_convection_blankenbach():
    _notebook_run("examples/1_03_BlankenbachBenchmark.ipynb")


def test_stokes_sinker():
    _notebook_run("examples/1_05_StokesSinker.ipynb")


def test_hypnic_jerk():
    _notebook_run("examples/1_06_HypnicJerk.ipynb")

# ### Not sure why this is not working. This triggers non-linear solver.
# ### def test_slab_subduction():
# ###    _notebook_run("examples/1_07_SlabSubduction.ipynb")


def test_viscoelastic_halfspace():
    _notebook_run("examples/1_08_ViscoElasticHalfSpace.ipynb")


def test_viscoelastic_shear():
    _notebook_run("examples/1_09_ViscoElasticShear.ipynb")


def test_viscoplasticity_simple_shear():
    _notebook_run("examples/1_10_Viscoelastoplasticity-in-simple-shear.ipynb")


def test_stokes_sinker_3D():
   _notebook_run("examples/1_11_StokesSinker3D.ipynb")


def test_columns_traction_bottom():
    _notebook_run("examples/1_20_ColumnsTractionBottom.ipynb")


def test_columns_traction_bottom_3D():
   _notebook_run("examples/1_21_3D_ColumnsTractionBottom.ipynb")


def test_freesurface_simple():
    _notebook_run("examples/1_23_01_FreeSurface_Simple_Example.ipynb")

    
def test_freesurface_Kaus2010_RT_instability():
    _notebook_run("examples/1_23_02_FreeSurface_Kaus2010_Rayleigh-Taylor_Instability.ipynb")
    
    
def test_freesurface_Crameri2012_Case1():
    _notebook_run("examples/1_23_03_FreeSurface_Crameri2012Case1_Relaxation.ipynb")
    
    
def test_freesurface_Crameri2012_Case1():
    _notebook_run("examples/1_23_04_FreeSurface_Crameri2012Case2_Rising_Plume.ipynb")


def test_define_3D_volume():
   _notebook_run("examples/1_24_Define_3D_volumes.ipynb")


def test_hot_canon_ball():
    _notebook_run("examples/1_25_Hot_Canon_Ball.ipynb")


def test_numerical_sandbox_moving_wall():
    _notebook_run("examples/1_26_NumericalSandboxCompression-MovingWall.ipynb")


def test_column_pure_thermal_advection():
    _notebook_run("examples/1_27_ColumnPureThermalAdvection.ipynb")


def test_poiseuille_under_pressure():
    _notebook_run("examples/1_30_Poiseuille_Under_Pressure.ipynb")


def test_user_defined_geotherm():
    _notebook_run("examples/1_31_User_defined_geotherm_and_TP_dependent_densities.ipynb")


def test_passive_tracers():
    _notebook_run("examples/1_32_Passive_Tracers_tests.ipynb")


def test_shear_band_pure_shear():
   _notebook_run("examples/2_09_ShearBandsPureShear.ipynb")


def test_rayleigh_taylor_kekken():
   _notebook_run("examples/2_15_Rayleigh-Taylor_van_Keken_et_al_1997.ipynb")

