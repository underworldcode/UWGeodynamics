import os
import subprocess
import tempfile
import sys

import nbformat

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(TEST_DIR, os.pardir))
sys.path.insert(0, PROJECT_DIR)


def _notebook_run(path):
    """Execute a notebook via nbconvert and collect output.
       :returns (parsed nb object, execution errors)
    """
    dirname, __ = os.path.split(path)
    os.chdir(dirname)
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        args = ["jupyter", "nbconvert",
                "--to", "notebook", "--execute",
                "--output", fout.name, path]
        subprocess.check_call(args)

        fout.seek(0)
        nb = nbformat.read(fout, nbformat.current_nbformat)

    errors = [output for cell in nb.cells if "outputs" in cell
              for output in cell["outputs"]
              if output.output_type == "error"]

    return nb, errors


def test_tutorial1():
    path = os.path.join(PROJECT_DIR,
                        "tutorials/Tutorial_1_ThermoMechanical_Model.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_tutorial2():
    path = os.path.join(PROJECT_DIR, "tutorials/Tutorial_2_Melt.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_tutorial3():
    path = os.path.join(
        PROJECT_DIR,
        "tutorials/Tutorial_3_SandboxExtension_static_mesh.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_tutorial4():
    path = os.path.join(
        PROJECT_DIR,
        "tutorials/Tutorial_4_NumericalSandboxCompression.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_tutorial5():
    path = os.path.join(PROJECT_DIR,
                        "tutorials/Tutorial_5_Convergence_Model.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_tutorial6():
    path = os.path.join(PROJECT_DIR,
                        "tutorials/Tutorial_6_Simple_Surface_Processes.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_tutorial7():
    path = os.path.join(PROJECT_DIR,
                        "tutorials/Tutorial_7_3D_Lithospheric_Model.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_tutorial8():
    path = os.path.join(PROJECT_DIR,
                        "tutorials/Tutorial_8_Subduction_ViscoElastic.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []
