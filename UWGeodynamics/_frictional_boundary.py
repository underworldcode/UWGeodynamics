import underworld as uw
from .scaling import nonDimensionalize as nd
import numpy as np

class FrictionBoundaries(object):
    """ This class flags elements at the boundaries

    """
    def __init__(self, Model, right=False, left=False,
                 top=False, bottom=False, front=False,
                 back=False, friction=0.0, thickness=2):

        self.Model = Model
        self.friction = friction
        self.thickness = thickness

        # Build borders
        globalIndices = np.arange(np.prod(Model.mesh.elementRes))
        globalIndices = globalIndices.reshape((Model.mesh.elementRes[::-1]))

        self.bottom = bottom
        self.right = right
        self.left = left
        self.top = top
        self.front = front
        self.back = back

        if self.bottom:
            self.bottom = globalIndices[:thickness].ravel()

        if self.top:
            self.top = globalIndices[-thickness:].ravel()

        if self.right and Model.mesh.dim < 3:
            self.right = globalIndices[:, -thickness:].ravel()

        if self.left and Model.mesh.dim < 3:
            self.left = globalIndices[:, :thickness].ravel()

        if self.right and Model.mesh.dim > 2:
            self.right = globalIndices[:, :, -thickness:].ravel()

        if self.left and Model.mesh.dim > 2:
            self.left = globalIndices[:, :, :thickness].ravel()

        if self.front and Model.mesh.dim > 2:
            self.front = globalIndices[:, :thickness, :].ravel()

        if self.back and Model.mesh.dim > 2:
            self.back = globalIndices[:, -thickness:, :].ravel()

        boundaries = [self.bottom, self.right, self.left, self.top, self.back,
                      self.front]
        boundaries = [boundary for boundary in boundaries if boundary is not False]

        self.boundaries = np.concatenate(boundaries)

        subMesh = Model.mesh.subMesh
        self._mask = uw.mesh.MeshVariable(mesh=subMesh, nodeDofCount=1)
        self._mask.data[:] = 0
        self._mask.data[np.intersect1d(subMesh.data_nodegId.ravel(), self.boundaries)] = 1

