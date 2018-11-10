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


def test_steady_state_example():
    path = os.path.join(
        PROJECT_DIR,
        "examples/1_01_Steady_State_Heat.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_convection_example():
    path = os.path.join(
        PROJECT_DIR,
        "examples/1_02_Convection_Example.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_convection_blankenbach():
    path = os.path.join(
        PROJECT_DIR,
        "examples/1_03_BlankenbachBenchmark.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_stokes_sinker():
    path = os.path.join(
        PROJECT_DIR,
        "examples/1_05_StokesSinker.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_slab_subduction():
    path = os.path.join(
        PROJECT_DIR,
        "examples/1_07_SlabSubduction.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_column_traction_bottom():
    path = os.path.join(
        PROJECT_DIR,
        "examples/1_20_ColumnsTractionBottom.ipynb")
    _, errors = _notebook_run(path)
    assert errors == []


def test_shear_bands_pure_shear():
    path = os.path.join(
        PROJECT_DIR,
        "examples/2_09_ShearBandsPureShear")
    _, errors = _notebook_run(path)
    assert errors == []
