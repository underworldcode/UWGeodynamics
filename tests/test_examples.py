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


#def test_slab_subduction():
#    _notebook_run("examples/1_07_SlabSubduction.ipynb")


def test_viscoelastic_halfspace():
    _notebook_run("examples/1_08_ViscoElasticHalfSpace.ipynb")


def test_viscoelastic_shear():
    _notebook_run("examples/1_09_ViscoElasticShear.ipynb")


def test_viscoplasticity_simple_shear():
    _notebook_run("examples/1_10_Viscoelastoplasticity-in-simple-shear.ipynb")


#def test_stokes_sinker_3D():
#    _notebook_run("examples/1_11_StokesSinker3D.ipynb")


def test_columns_traction_bottom():
    _notebook_run("examples/1_20_ColumnsTractionBottom.ipynb")


#def test_columns_traction_bottom_3D():
#    _notebook_run("examples/1_21_3D_ColumnsTractionBottom.ipynb")


def test_freesurface_simple():
    _notebook_run("examples/1_23_FreeSurface_Simple_Example.ipynb")


#def test_define_3D_volume():
#    _notebook_run("1_24_Define_3D_volumes.ipynb")


def test_hot_canon_ball():
    _notebook_run("examples/1_25_Hot_Canon_Ball.ipynb")


def test_shear_band_pure_shear():
   _notebook_run("examples/2_09_ShearBandsPureShear.ipynb")


def test_rayleigh_taylor_kekken():
   _notebook_run("examples/2_15_Rayleigh-Taylor_van_Keken_et_al_1997.ipynb")

