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
