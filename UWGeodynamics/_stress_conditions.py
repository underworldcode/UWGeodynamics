import underworld as uw
import underworld.function as fn
import numpy as np
from mpi4py import MPI
from .LecodeIsostasy import LecodeIsostasy
from .scaling import nonDimensionalize as nd
from .scaling import UnitRegistry as u
from ._utils import Balanced_InflowOutflow
from ._utils import MovingWall

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

