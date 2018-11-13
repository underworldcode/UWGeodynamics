from __future__ import print_function,  absolute_import
import underworld as uw
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def _visugrid_drawing_object(Model, visugrid):

    m_minx, m_maxx = _get_minmax_coordinates_mesh(Model.mesh, 0)
    m_miny, m_maxy = _get_minmax_coordinates_mesh(Model.mesh, 1)
    model_domain = np.array([(m_minx, m_miny), (m_maxx, m_maxy)])

    v_minx, v_maxx = _get_minmax_coordinates_mesh(visugrid.mesh, 0)
    v_miny, v_maxy = _get_minmax_coordinates_mesh(visugrid.mesh, 1)
    # visugrid_domain = np.array([(v_minx, v_miny), (v_maxx, v_maxy)])

    minx = min(m_minx, v_minx)
    miny = min(m_miny, v_miny)
    maxx = max(m_maxx, v_maxx)
    maxy = max(m_maxy, v_maxy)

    full_domain = np.array([(minx, miny), (maxx, maxy)])
    bounds = full_domain[1] - full_domain[0]

    minx = (model_domain[0][0] - full_domain[0][0]) / bounds[0]
    maxx = 1 + (model_domain[1][0] - full_domain[1][0]) / bounds[0]
    miny = (model_domain[0][1] - full_domain[0][1]) / bounds[1]
    maxy = 1 + (model_domain[1][1] - full_domain[1][1]) / bounds[1]

    return (minx[0], maxx[0]), (miny[0], maxy[0])


def _get_minmax_coordinates_mesh(mesh, axis=0):
    """ Return the minimum and maximum coordinates along axis

    parameter:
    ----------
        axis:
            axis

    returns:
    -------
        tuple: minV, maxV

    """
    maxVal = np.zeros((1))
    minVal = np.zeros((1))
    maxVal[0] = mesh.data[:, axis].max()
    minVal[0] = mesh.data[:, axis].min()

    uw.barrier()
    comm.Allreduce(MPI.IN_PLACE, maxVal, op=MPI.MAX)
    comm.Allreduce(MPI.IN_PLACE, minVal, op=MPI.MIN)
    uw.barrier()

    return minVal, maxVal

