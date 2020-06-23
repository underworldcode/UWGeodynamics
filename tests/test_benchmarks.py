from .utils import _notebook_run


def test_slab_detachment_benchmark():
    _notebook_run("docs/benchmarks/1_12_Slab_Detachment_Benchmark.ipynb")

def test_indentor_benchmark():
    _notebook_run("docs/benchmarks/1_22_Indentor_Benchmark.ipynb")

def test_self_subduction_case1():
    _notebook_run("docs/benchmarks/2D_Self_Subduction_Case1.ipynb")

def test_self_subduction_case2():
    _notebook_run("docs/benchmarks/2D_Self_Subduction_Case2.ipynb")

def test_brick_compression():
    _notebook_run("docs/benchmarks/Kaus_BrickBenchmark-Compression.ipynb")

def test_brick_extension():
    _notebook_run("docs/benchmarks/Kaus_BrickBenchmark_Extension.ipynb")
