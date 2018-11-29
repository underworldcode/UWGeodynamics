import os
import subprocess
import sys
from shutil import copyfile
import re

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(TEST_DIR, os.pardir))
RESULT_DIR = os.path.join(PROJECT_DIR, "tests/test_results")
sys.path.insert(0, PROJECT_DIR)


def _notebook_run(example):
    path = os.path.join(PROJECT_DIR, example)
    outpath = os.path.join(RESULT_DIR, example)
    copyfile(path, outpath)

    s = open(outpath).read()
    s = re.sub(r'Model\.run\_for\(.*\)', 'Model.run_for(nstep=2)', s)
    f = open(outpath, 'w')
    f.write(s)
    f.close()

    args = ["jupyter", "nbconvert",
            "--to", "notebook", "--execute",
            "--ExecutePreprocessor.timeout=0",
            "--output", outpath, outpath]
    subprocess.check_call(args)


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
