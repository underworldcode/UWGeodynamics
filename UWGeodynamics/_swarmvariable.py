import underworld as uw
import h5py
import numpy as np
from mpi4py import MPI
from .scaling import Dimensionalize
from .scaling import nonDimensionalize
from .scaling import UnitRegistry as u


class SwarmVariable(uw.swarm.SwarmVariable):

    def __init__(self, swarm, dataType, count, writeable=True, **kwargs):
        super(SwarmVariable, self).__init__(swarm, dataType, count,
                                            writeable=True, **kwargs)

    def save(self, filename, units=None, time=None):
        """
        Save the swarm variable to disk.

        Parameters
        ----------
        filename : str
            The filename for the saved file. Relative or absolute paths may be
            used, but all directories must exist.
        swarmHandle :uw.utils.SavedFileData , optional
            The saved swarm file handle. If provided, a reference to the swarm file
            is made. Currently this doesn't provide any extra functionality.

        Returns
        -------
        underworld.utils.SavedFileData
            Data object relating to saved file. This only needs to be retained
            if you wish to create XDMF files and can be ignored otherwise.

        Notes
        -----
        This method must be called collectively by all processes.

        Example
        -------
        First create the swarm, populate, then add a variable:

        >>> mesh = uw.mesh.FeMesh_Cartesian( elementType='Q1/dQ0', elementRes=(16,16), minCoord=(0.,0.), maxCoord=(1.,1.) )
        >>> swarm = uw.swarm.Swarm(mesh)
        >>> swarm.populate_using_layout(uw.swarm.layouts.PerCellGaussLayout(swarm,2))
        >>> svar = swarm.add_variable("int",1)

        Write something to variable

        >>> import numpy as np
        >>> svar.data[:,0] = np.arange(swarm.particleLocalCount)

        Save to a file:

        >>> ignoreMe = swarm.save("saved_swarm.h5")
        >>> ignoreMe = svar.save("saved_swarm_variable.h5")

        Now let's try and reload. First create a new swarm and swarm variable,
        and then load both:

        >>> clone_swarm = uw.swarm.Swarm(mesh)
        >>> clone_svar = clone_swarm.add_variable("int",1)
        >>> clone_swarm.load("saved_swarm.h5")
        >>> clone_svar.load("saved_swarm_variable.h5")

        Now check for equality:

        >>> import numpy as np
        >>> np.allclose(svar.data,clone_svar.data)
        True

        >>> # clean up:
        >>> if uw.rank() == 0:
        ...     import os;
        ...     os.remove( "saved_swarm.h5" )
        ...     os.remove( "saved_swarm_variable.h5" )

        """

        if not isinstance(filename, str):
            raise TypeError("'filename' parameter must be of type 'str'")

        # setup mpi basic vars
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nProcs = comm.Get_size()

        # allgather the number of particles each proc has
        swarm = self.swarm
        procCount = comm.allgather(swarm.particleLocalCount)
        particleGlobalCount = np.sum(procCount) #swarm.particleGlobalCount

        # calculate the hdf5 file offset
        offset=0
        for i in xrange(rank):
            offset += procCount[i]

        # open parallel hdf5 file
        h5f = h5py.File(name=filename, mode="w", driver='mpio', comm=MPI.COMM_WORLD)

        # attribute of the proc offsets - used for loading from checkpoint
        h5f.attrs["proc_offset"] = procCount

        # write the entire local swarm to the appropriate offset position
        globalShape = (particleGlobalCount, self.data.shape[1])
        dset = h5f.create_dataset("data",
                                   shape=globalShape,
                                   dtype=self.data.dtype)
        fact = 1.0
        if units:
            fact = Dimensionalize(1.0, units=units).magnitude
            h5f.attrs['units'] = str(units)

        if time is not None:
            h5f.attrs['time'] = str(time)

        if swarm.particleLocalCount > 0: # only add if there are local particles
            dset[offset:offset+swarm.particleLocalCount] = self.data[:] * fact

        h5f.close()

        return uw.utils.SavedFileData( self, filename )

    def load(self, filename):
        """
        Load the swarm variable from disk. This must be called *after* the swarm.load().

        Parameters
        ----------
        filename : str
            The filename for the saved file. Relative or absolute paths may be
            used, but all directories must exist.

        Notes
        -----
        This method must be called collectively by all processes.


        Example
        -------
        Refer to example provided for 'save' method.

        """

        if not isinstance(filename, str):
            raise TypeError("'filename' parameter must be of type 'str'")

        if self.swarm._checkpointMapsToState != self.swarm.stateId:
            raise RuntimeError("'Swarm' associate with this 'SwarmVariable' does not appear to be in the correct state.\n" \
                               "Please ensure that you have loaded the swarm prior to loading any swarm variables.")
        gIds = self.swarm._local2globalMap

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nProcs = comm.Get_size()

        # open hdf5 file
        h5f = h5py.File(name=filename, mode="r", driver='mpio', comm=MPI.COMM_WORLD)

        # get units
        try:
            units = h5f.attrs["units"]
        except:
            units = None

        if units and units != "None":
            units = u.parse_expression(units)
        else:
            units = None

        dset = h5f.get('data')
        if dset == None:
            raise RuntimeError("Can't find 'data' in file '{}'.\n".format(filename))

        if dset.shape[1] != self.data.shape[1]:
            raise RuntimeError("Cannot load file data on current swarm. Data in file '{0}', " \
                               "has {1} components -the particlesCoords has {2} components".format(filename, dset.shape[1], self.particleCoordinates.data.shape[1]))

        particleGobalCount = self.swarm.particleGlobalCount

        if dset.shape[0] != particleGobalCount:
            if rank == 0:
                import warnings
                warnings.warn("Warning, it appears {} particles were loaded, but this h5 variable has {} data points". format(particleGobalCount, dset.shape[0]), RuntimeWarning)

        size = len(gIds) # number of local2global mapped indices
        if size > 0:     # only if there is a non-zero local2global do we load
            if units:
                vals = dset[gIds,:]
                self.data[:] = nonDimensionalize(vals * units)
            else:
                self.data[:] = dset[gIds,:]

        h5f.close();

