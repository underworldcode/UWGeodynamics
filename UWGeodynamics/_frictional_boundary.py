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

        if self.right:
            if Model.mesh.dim < 3:
                self.right = globalIndices[:, -thickness:].ravel()
            else:
                self.right = globalIndices[:, :, -thickness:].ravel()

        if self.left:
            if Model.mesh.dim < 3:
                self.left = globalIndices[:, :thickness].ravel()
            else:
                self.left = globalIndices[:, :, :thickness].ravel()

        if self.front:
            if Model.mesh.dim > 2:
                self.front = globalIndices[:, :thickness, :].ravel()
            else:
                raise ValueError("Mesh is 2D")

        if self.back:
            if Model.mesh.dim > 2:
                self.back = globalIndices[:, -thickness:, :].ravel()
            else:
                raise ValueError("Mesh is 2D")

        boundaries = [self.bottom, self.right, self.left, self.top, self.back,
                      self.front]
        boundaries = [boundary for boundary in boundaries if boundary is not False]

        self.boundaries = np.concatenate(boundaries)

        subMesh = Model.mesh.subMesh
        self._mask = uw.mesh.MeshVariable(mesh=subMesh, nodeDofCount=1)
        self._mask.data[:] = 0
        # Take the intersection between the globalID and the boundaries where
        # friction is to be applied
        intersect = np.intersect1d(subMesh.data_nodegId.ravel(), self.boundaries)
        # Create a mask to highlight those elements in the local domain
        self._mask.data[:,0] = np.in1d(subMesh.data_nodegId.ravel(), intersect)

