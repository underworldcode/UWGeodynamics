import numpy as np
from .scaling import nonDimensionalize as nd


def remesh(mesh, x=None, y=None, z=None):
    nx, ny = mesh.elementRes

    def new_points(intervals, elements):
        pts = []
        for idx in range(len(intervals) - 1):
            pts.append(np.linspace(intervals[idx], intervals[idx + 1], elements[idx] + 1))
        pts = np.unique(np.hstack(pts))
        pts.sort()
        return pts

    def check_vals(intervals, elements, axis):
        if ((intervals[0] != mesh.minCoord[axis]) or
            (intervals[-1] != mesh.maxCoord[axis])):
            raise ValueError("""Intervals do not match mesh extent""")

        if np.sum(np.array(elements)) != mesh.elementRes[axis]:
            raise ValueError("""Total nb of elements do not match the nb of elements in the mesh""")

    if x:
        intervals, elements = x
        intervals = [nd(val) for val in intervals]
        elements = [nd(val) for val in elements]
        check_vals(intervals, elements, 0)
        pts = new_points(intervals, elements)

        with mesh.deform_mesh():
            new_vals = np.tile(pts, ny + 1)
            mesh.data[:, 0] = new_vals[mesh.data_nodegId.flatten()]

    if y:
        intervals, elements = y
        intervals = [nd(val) for val in intervals]
        elements = [nd(val) for val in elements]
        check_vals(intervals, elements, 1)
        pts = new_points(intervals, elements)

        with mesh.deform_mesh():
            vals = np.repeat(pts, nx + 1)
            mesh.data[:, 1] = vals[mesh.data_nodegId.flatten()]
    return mesh

