import os
import subprocess
import tempfile
import sys

import nbformat

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(TEST_DIR, os.pardir))
sys.path.insert(0, PROJECT_DIR)


def _notebook_run(example):
    path = os.path.join(PROJECT_DIR, example)
    args = ["jupyter", "nbconvert",
            "--to", "notebook", "--execute",
            "--output", path, path]
    subprocess.check_call(args)


def test_steady_state_example():
    _notebook_run("examples/1_01_Steady_State_Heat.ipynb")

#def test_convection_example():
#    _notebook_run("examples/1_02_Convection_Example.ipynb")

#def test_convection_blankenbach():
#    _notebook_run("examples/1_03_BlankenbachBenchmark.ipynb")

def test_stokes_sinker():
    _notebook_run("examples/1_05_StokesSinker.ipynb")

def test_hypnic_jerk():
    _notebook_run("examples/1_06_HypnicJerk.ipynb")

def test_slab_subduction():
    _notebook_run("examples/1_07_SlabSubduction.ipynb")

def test_viscoelastic_halfspace():
    _notebook_run("examples/1_08_ViscoElasticHalfSpace.ipynb")

def test_viscoelastic_shear():
    _notebook_run("examples/1_09_ViscoElasticShear.ipynb")

def test_viscoplasticity_simple_shear():
    _notebook_run("examples/1_10_Viscoelastoplasticity-in-simple-shear.ipynb")

def test_indentor_benchmark():
    _notebook_run("examples/1_22_Indentor_Benchmark.ipynb")

def test_freesurface_simple():
    _notebook_run("examples/1_23_FreeSurface_Simple_Example.ipynb")

def test_canon_ball():
    _notebook_run("examples/1_25_Hot_Canon_Ball.ipynb")


