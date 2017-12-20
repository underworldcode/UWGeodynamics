from copy import copy
import underworld as uw
import numpy as np
import sys
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


class _mesh_advector(object):

    def __init__(self, Model, axis):

        self._mesh2nd = Model.mesh
        self.Model = Model

    def _advect_along_axis(self, dt, axis=0):
        if axis != 0:
            raise ValueError("Axis not supported yet")

        # Get minimum and maximum coordinates for the current mesh
        minX, maxX = self._get_minmax_coordinates_mesh(axis)

        minvxLeftWall, maxvxLeftWall   = self._get_minmax_velocity_wall(self.Model.leftWall, axis)
        minvxRightWall, maxvxRightWall = self._get_minmax_velocity_wall(self.Model.rightWall, axis)

        if np.abs(maxvxRightWall) > np.abs(minvxRightWall):
            vxRight = maxvxRightWall   
        else:
            vxRight = minvxRightWall
        
        if (np.abs(maxvxLeftWall)  > np.abs(minvxLeftWall)):
            vxLeft = maxvxLeftWall  
        else:
            vxLeft = minvxLeftWall

        minX += vxLeft * dt
        maxX += vxRight * dt
        length = np.abs(minX - maxX)

        newValues = np.linspace(minX, maxX, self.Model.mesh.elementRes[axis]+1)
        newValues = np.repeat(newValues[np.newaxis,:], self.Model.mesh.elementRes[1] + 1, axis)

        with self._mesh2nd.deform_mesh():
            self._mesh2nd.data[:, axis] = newValues.flatten()
        
        uw.barrier

        with self.Model.mesh.deform_mesh():
            self.Model.mesh.data[:, axis] = self._mesh2nd.data[:, axis]

        self.Model.velocityField.data[...] = np.copy(self.Model.velocityField.evaluate(self.Model.mesh))
        self.Model.pressureField.data[...] = np.copy(self.Model.pressureField.evaluate(self.Model.mesh.subMesh))
       
        if self.Model.rightWall:
            self.Model.velocityField.data[self.Model.rightWall.data,:] = vxRight
        
        if self.Model.leftWall:
            self.Model.velocityField.data[self.Model.leftWall.data,:]  = vxLeft

    def advect_mesh(self, dt):
        self._advect_along_axis(dt)

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
            velocities  = self.Model.velocityField.data[wall.data, axis]
       
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
        
        return minVal, maxVal

