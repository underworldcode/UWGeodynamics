from copy import copy
import underworld as uw
import numpy as np
import sys


class _mesh_advector(object):

    def __init__(self, Model, axis):

        self._mesh2nd = Model.mesh
        self.Model = Model

    def _advect_along_axis(self, dt, axis=0):
        if axis != 0:
            raise ValueError("Axis not supported yet")

        # Get minimum and maximum coordinates for the current mesh
        minx, maxx = _get_minmax_coordinates_mesh(axis)

        minVxLeftWall, maxVxLeftWall   = self._get_minmax_velocity_wall(self.Model.leftWall, axis)
        minVxrightWall, maxVxrightWall = self._get_minmax_velocity_wall(self.Model.rightWall, axis)

        vxright = maxVxrightWall if (np.abs(maxVxrightWall) > np.abs(minVxrightWall)) else minVxrightWall
        vxleft  = maxVxleftWall  if (np.abs(maxVxleftWall)  > np.abs(minVxleftWall))  else minVxleftWall

        newMinx = minx + vxleft * dt
        newMaxX = maxx + vxright * dt
        newLength = np.abs(newMinx - newMaxx)

        newValues = np.linspace(newMinx, newMaxx, self.Model.mesh.elementRes[axis]+1)

        with self._mesh2nd.deform_mesh():
            indices = self._mesh2nd.data_nodegId % self.Model.mesh.elementRes[axis]+1
            self._mesh2nd.data[:, axis] = newValues[indices.flatten()]
        
        uw.barrier

        with self.Model.mesh.deform_mesh():
            self.Model.mesh.data[:, axis] = self._mesh2nd.data[:, axis]

        self.Model.velocityField.data[...] = np.copy(self.Model.velocityField.evaluate(self.Model.mesh))
        self.Model.pressureField.data[...] = np.copy(self.Model.pressureField.evaluate(self.Model.subMesh))
       
        if self.Model.rightWall:
            self.Model.velocityField.data[self.Model.rightWall.data,:] = vxright
        
        if self.Model.leftWall:
            self.Model.velocityField.data[self.Model.leftWall.data,:]  = vxleft

    def advect_mesh(self, dt):
        self.advect_along_axis(dt)

    def _get_minmax_velocity_wall(self, wall, axis=0):
        """ Return the minimum and maximum velocity component on the wall 
        
        parameters:
        -----------
            wall: (indexSet)
                The wall.
            axis: 
                axis (velocity component). 
        """

        # Initialise value to max and min sys values
        maxV = np.ones((1)) * sys.float_info.min
        minV = np.ones((1)) * sys.float_info.max


        # if local domain has wall, get velocities
        if wall:
            velocities  = self.velocityField.data[wall.data, axis]
       
        # get local min and max
        maxV[0] = velocities.max()
        minV[0] = velocities.min()

        # reduce operation
        uw.barrier
        comm.Allreduce(MPI.IN_PLACE, maxV, op=MPI.MAX)
        comm.Allreduce(MPI.IN_PLACE, minV, op=MPI.MIN)
        uw.barrier

        return minV, maxV

    def _get_minmax_coordinates_mesh(self, axis=0):
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
        maxVal[0] = self.Model.mesh.data[:, axis].max()
        minVal[0] = self.Model.mesh.data[:, axis].min()

        uw.barrier
        comm.Allreduce(MPI.IN_PLACE, maxVal, op=MPI.MAX)
        comm.Allreduce(MPI.IN_PLACE, minVal, op=MPI.MIN)
        uw.barrier
        
        return minx, maxx

