import underworld as uw
from .scaling import nonDimensionalize as nd
import numpy as np

class FrictionBoundaries(object):
    """ This class flags elements at the boundaries

        Works in 2D only for now
    """
    def __init__(self, Model, right=False, left=False,
                 top=False, bottom=False, friction=0.0, thickness=2):

        if Model.mesh.dim > 2:
            raise NotImplemented("Frictional boundaries are not yet implemented in 3D")

        self.Model = Model
        self.friction = friction
        self.thickness = thickness

        # Build borders
        globalIndices = np.arange(reduce(lambda x, y: x*y, Model.mesh.elementRes))
        globalIndices = globalIndices.reshape((Model.mesh.elementRes[::-1]))

        self.bottom = bottom
        self.right  = right
        self.left   = left
        self.top    = top

        if self.bottom:
            self.bottom = globalIndices[:thickness].ravel()

        if self.top:
            self.top = globalIndices[-thickness:].ravel()

        if self.right:
            self.right = globalIndices[:,-thickness:].ravel()

        if self.left:
            self.left = globalIndices[:,:thickness].ravel()

        boundaries = [self.bottom, self.right, self.left, self.top]
        boundaries = [boundary for boundary in boundaries if boundary is not False]

        self.boundaries = np.concatenate(boundaries)

        subMesh = Model.mesh.subMesh
        self._mask = uw.mesh.MeshVariable(mesh=subMesh, nodeDofCount=1)
        self._mask.data[:] = 0
        self._mask.data[np.intersect1d(subMesh.data_nodegId.ravel(), self.boundaries)] = 1





