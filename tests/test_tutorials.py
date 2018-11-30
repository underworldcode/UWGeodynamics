from .utils import _notebook_run

def test_tutorial1():
    _notebook_run("tutorials/Tutorial_1_ThermoMechanical_Model.ipynb")


def test_tutorial2():
    _notebook_run("tutorials/Tutorial_2_Melt.ipynb")


def test_tutorial3():
    _notebook_run("tutorials/Tutorial_3_SandboxExtension_static_mesh.ipynb")


def test_tutorial3B():
    _notebook_run("tutorials/Tutorial_3B_SandboxExtension_deform_mesh.ipynb")


def test_tutorial4():
    _notebook_run("tutorials/Tutorial_4_NumericalSandboxCompression.ipynb")


def test_tutorial5():
    _notebook_run("tutorials/Tutorial_5_Convergence_Model.ipynb")


def test_tutorial6():
    _notebook_run("tutorials/Tutorial_6_Simple_Surface_Processes.ipynb")


def test_tutorial7():
    _notebook_run("tutorials/Tutorial_7_3D_Lithospheric_Model.ipynb")


def test_tutorial8():
    _notebook_run("tutorials/Tutorial_8_Subduction_ViscoElastic.ipynb")


def test_tutorial9():
    _notebook_run("tutorials/Tutorial_9_passive_margins.ipynb")


def test_tutorial10():
    _notebook_run("tutorials/Tutorial_10_Thrust_Wedges.ipynb")
