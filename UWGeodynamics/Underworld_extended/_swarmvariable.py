from __future__ import print_function,  absolute_import
import underworld as uw
import h5py
import numpy as np
from mpi4py import MPI
from UWGeodynamics.scaling import Dimensionalize
from UWGeodynamics.scaling import nonDimensionalize
from UWGeodynamics.scaling import UnitRegistry as u
from UWGeodynamics.version import git_revision as __git_revision__


class SwarmVariable(uw.swarm.SwarmVariable):

    def __init__(self, swarm, dataType, count, writeable=True, **kwargs):
        super(SwarmVariable, self).__init__(swarm, dataType, count,
                                            writeable=True, **kwargs)

    def load(self, filename, collective=False):
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
        rank = comm.rank

        # open hdf5 file
        h5f = h5py.File(name=filename, mode="r", driver='mpio', comm=MPI.COMM_WORLD)


        dset = h5f.get('data')
        if dset == None:
            raise RuntimeError("Can't find 'data' in file '{}'.\n".format(filename))

        if dset.shape[1] != self.data.shape[1]:
            raise RuntimeError("Cannot load file data on current swarm. Data in file '{0}', " \
                               "has {1} components -the particlesCoords has {2} components".format(filename, dset.shape[1], self.particleCoordinates.data.shape[1]))

        if dset.shape[0] != self.swarm.particleGlobalCount:
            raise RuntimeError("It appears that the swarm has {} particles, but provided h5 file has {} data points. Please check that " \
                               "both the Swarm and the SwarmVariable were saved at the same time, and that you have reloaded using " \
                               "the correct files.".format(particleGobalCount, dset.shape[0]))

        # for efficiency, we want to load swarmvariable data in the largest stride chunks possible.
        # we need to determine where required data is contiguous.
        # first construct an array of gradients. the required data is contiguous
        # where the indices into the array are increasing by 1, ie have a gradient of 1.
        gradIds = np.zeros_like(gIds)            # creates array of zeros of same size & type
        if len(gIds) > 1:
            gradIds[:-1] = gIds[1:] - gIds[:-1]  # forward difference type gradient

        # note that we do only the first read into dset collective. this call usually
        # does the entire read, but if it doesn't we won't know how many calls will
        # be necessary, hence only collective calling the first.
        done_collective = False
        guy = 0
        while guy < len(gIds):
            # do contiguous
            start_guy = guy
            while gradIds[guy]==1:  # count run of contiguous. note bounds check not required as last element of gradIds is always zero.
                guy += 1
            # copy contiguous chunk if found.. note that we are copying 'plus 1' items
            if guy > start_guy:
                if collective and not done_collective:
                    with dset.collective:
                        self.data[start_guy:guy+1] = dset[gIds[start_guy]:gIds[guy]+1]
                        done_collective = True
                else:
                    self.data[start_guy:guy+1] = dset[gIds[start_guy]:gIds[guy]+1]
                guy += 1

            # do non-contiguous
            start_guy = guy
            while guy<len(gIds) and gradIds[guy]!=1:  # count run of non-contiguous
                guy += 1
            # copy non-contiguous items (if found) using index array slice
            if guy > start_guy:
                if collective and not done_collective:
                    with dset.collective:
                        self.data[start_guy:guy,:] = dset[gIds[start_guy:guy],:]
                        done_collective = True
                else:
                    self.data[start_guy:guy,:] = dset[gIds[start_guy:guy],:]

        # if we haven't entered a collective call, do so now to
        # avoid deadlock. we just do an empty read/write.
        if collective and not done_collective:
            with dset.collective:
                self.data[0:0, :] = dset[0:0, :]

        # get units
        try:
            units = h5f.attrs["units"]
        except KeyError:
            units = None

        if units and units != "None":
            units = u.parse_expression(units)
        else:
            units = None

        h5f.close()

        if units:
            self.data[:] = nonDimensionalize(self.data * units)

    def save( self, filename, collective=False, units=None, time=None):
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
        rank = comm.rank

        # allgather the number of particles each proc has
        swarm = self.swarm
        procCount = comm.allgather(swarm.particleLocalCount)
        particleGlobalCount = np.sum(procCount)

        # calculate the hdf5 file offset
        offset=0
        for i in range(comm.rank):
            offset += procCount[i]

        # open parallel hdf5 file
        with h5py.File(name=filename, mode="w", driver='mpio', comm=MPI.COMM_WORLD) as h5f:
            # write the entire local swarm to the appropriate offset position
            globalShape = (particleGlobalCount, self.data.shape[1])
            dset = h5f.create_dataset("data",
                                       shape=globalShape,
                                       dtype=self.data.dtype)
            fact = 1.0
            if units:
                fact = Dimensionalize(1.0, units=units).magnitude

            if collective:
                with dset.collective:
                    dset[offset:offset+swarm.particleLocalCount] = self.data[:] * fact
            else:
                dset[offset:offset + swarm.particleLocalCount] = self.data[:] * fact

        # let's reopen in serial to write the attrib.
        # not sure if this really is necessary.
        comm.barrier()
        if comm.rank == 0:
            with h5py.File(name=filename, mode="a") as h5f:
                # attribute of the proc offsets - used for loading from checkpoint
                h5f.attrs["proc_offset"] = procCount
                h5f.attrs['units'] = str(units)
                h5f.attrs['time'] = str(time)
                h5f.attrs["git commit"] = __git_revision__

        return uw.utils.SavedFileData( self, filename )

    def copy(self, deepcopy=False):
        """
        This method returns a copy of the swarmvariable.

        Parameters
        ----------
        deepcopy: bool
            If True, the variable's data is also copied into
            new variable.

        Returns
        -------
        underworld.swarm.SwarmVariable
            The swarm variable copy.

        Example
        -------
        >>> import math
        >>> mesh = uw.mesh.FeMesh_Cartesian()
        >>> swarm = uw.swarm.Swarm()
        >>> swarm.populate.using_layout(uw.swarm.layouts.PerCellGaussLayout(swarm, 2)
        >>> svar = swarm.add_variable("double", 1)
        >>> svar.data[:] = math.exp(1.)
        >>> svarCopy = svar.copy()
        >>> svarCopy.swarm == var.swarm
        True
        >>> svarCopy.dataType == svar.dataType
        True
        >>> import numpy as np
        >>> np.allclose(svar.data,svarCopy.data)
        False
        >>> svarCopy2 = svar.copy(deepcopy=True)
        >>> np.allclose(svar.data,svarCopy2.data)
        True

        """

        if not isinstance(deepcopy, bool):
            raise TypeError("'deepcopy' parameter is expected to be of type 'bool'.")

        newSv = SwarmVariable(self.swarm, self.dataType, self.count)

        if deepcopy:
            newSv.data[:] = self.data[:]

        return newSv
