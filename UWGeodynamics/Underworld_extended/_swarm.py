from __future__ import print_function,  absolute_import
import underworld as uw
import h5py
import numpy as np
from mpi4py import MPI
from UWGeodynamics.scaling import nonDimensionalize
from UWGeodynamics.scaling import UnitRegistry as u
from . import _swarmvariable as svar


class Swarm(uw.swarm.Swarm):
    def __init__(self, mesh, particleEscape=False, **kwargs):
        super(Swarm, self).__init__(mesh, particleEscape, **kwargs)

    def _setup(self):
        if self._cself.particleCoordVariable:
            self._particleCoordinates = svar.SwarmVariable(
                self, "double", self.mesh.dim,
                _cself=self._cself.particleCoordVariable,
                writeable=False)

    def add_variable(self, dataType, count):
        """
        Add a variable to each particle in this swarm. Variables can be added
        at any point. Removal of variables is however not currently supported.
        See help(SwarmVariable) for further information.

        Parameters
        ----------
        dataType: str
            The data type for the variable. Available types are  "char",
            "short", "int", "float" or "double".
        count: unsigned
            The number of values to be stored for each particle.

        Returns
        -------
        underworld.swarm.SwarmVariable
            The newly created swarm variable.

        Example
        -------
        >>> # first we need a mesh
        >>> mesh = uw.mesh.FeMesh_Cartesian( elementType='Q1/dQ0', elementRes=(16,16), minCoord=(0.,0.), maxCoord=(1.,1.) )
        >>> # create swarm
        >>> swarm = uw.swarm.Swarm(mesh)
        >>> # add a variable
        >>> svar = swarm.add_variable("char",1)
        >>> # add another
        >>> svar = swarm.add_variable("double",3)
        >>> # add some particles
        >>> swarm.populate_using_layout(uw.swarm.layouts.PerCellGaussLayout(swarm,2))
        >>> # add another variable
        >>> svar = swarm.add_variable("double",5)

        """
        return svar.SwarmVariable( self, dataType, count )

    def save(self, filename, collective=False, units=None, time=None):
        """
        Save the swarm to disk.

        Parameters
        ----------
        filename : str
            The filename for the saved file. Relative or absolute paths may be
            used, but all directories must exist.

        Returns
        -------
        underworld.utils.SavedFileData
            Data object relating to saved file. This only needs to be retained
            if you wish to create XDMF files and can be ignored otherwise.

        Notes
        -----
        This method must be called collectively by all processes

        Example
        -------
        First create the swarm, and populate with layout:

        >>> mesh = uw.mesh.FeMesh_Cartesian( elementType='Q1/dQ0', elementRes=(16,16), minCoord=(0.,0.), maxCoord=(1.,1.) )
        >>> swarm = uw.swarm.Swarm(mesh)
        >>> swarm.populate_using_layout(uw.swarm.layouts.PerCellGaussLayout(swarm,2))

        Save to a file:

        >>> ignoreMe = swarm.save("saved_swarm.h5")

        Now let's try and reload. First create an empty swarm, and then load:

        >>> clone_swarm = uw.swarm.Swarm(mesh)
        >>> clone_swarm.load( "saved_swarm.h5" )

        Now check for equality:

        >>> import numpy as np
        >>> np.allclose(swarm.particleCoordinates.data,clone_swarm.particleCoordinates.data)
        True

        >>> # clean up:
        >>> if uw.rank() == 0:
        ...     import os;
        ...     os.remove( "saved_swarm.h5" )

        """

        if not isinstance(filename, str):
            raise TypeError("Expected filename to be provided as a string")

        # just save the particle coordinates SwarmVariable
        self.particleCoordinates.save(filename, collective, units=units, time=time)

        return uw.utils.SavedFileData( self, filename )

    def load( self, filename, collective=False, try_optimise=True, verbose=False ):
        """
        Load a swarm from disk. Note that this must be called before any SwarmVariable
        members are loaded.

        Parameters
        ----------
        filename : str
            The filename for the saved file. Relative or absolute paths may be
            used.
        try_optimise : bool, Default=True
            Will speed up the swarm load time but warning - this algorithm assumes the
            previously saved swarm data was made on an identical mesh and mesh partitioning
            (number of processors) with respect to the current mesh. If this isn't the case then
            the reloaded particle ordering will be broken, leading to an invalid swarms.
            One can disable this optimisation and revert to a brute force algorithm, much slower,
            by setting this option to False.
        verbose : bool
            Prints a swarm load progress bar.

        Notes
        -----
        This method must be called collectively by all processes.

        Example
        -------
        Refer to example provided for 'save' method.

        """

        if not isinstance(filename, str):
            raise TypeError("Expected 'filename' to be provided as a string")

        # open hdf5 file
        h5f = h5py.File(name=filename, mode="r", driver='mpio', comm=MPI.COMM_WORLD)

        # get units
        try:
            units = h5f.attrs["units"]
        except KeyError:
            units = None

        if units and units != "None":
            units = u.parse_expression(units)
        else:
            units = None

        dset = h5f.get('data')
        if dset == None:
            raise RuntimeError("Can't find 'data' in file '{0}'.\n".format(filename))
        if dset.shape[1] != self.particleCoordinates.data.shape[1]:
            raise RuntimeError("Cannot load file data on current swarm. Data in file '{0}', " \
                               "has {1} components -the particlesCoords has {2} components".format(filename, dset.shape[1], self.particleCoordinates.data.shape[1]))
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nProcs = comm.Get_size()

        if rank == 0 and verbose:
            bar = uw.utils._ProgressBar( start=0, end=dset.shape[0]-1, title="loading "+filename)

        # try and read the procCount attribute & assume that if nProcs in .h5 file
        # is equal to the current no. procs then the particles will be distributed the
        # same across the processors. (Danger if different discretisations are used... i think)
        # else try and load the whole .h5 file.
        # we set the 'offset' & 'size' variables to achieve the above

        offset = 0
        size = dset.shape[0] # number of particles in h5 file

        if try_optimise:
            procCount = h5f.attrs.get('proc_offset')
            if procCount is not None and nProcs == len(procCount):
                for p_i in range(rank):
                    offset += procCount[p_i]
                size = procCount[rank]

        valid = np.zeros(0, dtype='i') # array for read in
        chunk=int(2e7) # read in this many points at a time

        firstChunk = True
        (multiples, remainder) = divmod( size, chunk )
        for ii in range(multiples+1):
            # setup the points to begin and end reading in
            chunkStart = offset + ii*chunk
            if ii == multiples:
                chunkEnd = chunkStart + remainder
                if remainder == 0: # in the case remainder is 0
                    break
            else:
                chunkEnd = chunkStart + chunk

            # Add particles to swarm, ztmp is the corresponding local array
            # non-local particles are not added and their ztmp index is -1.
            # Note that for the first chunk, we do collective read, as this
            # is the only time that we can guaranteed that all procs will
            # take part, and usually most (if not all) particles are loaded
            # in this step.
            if firstChunk and collective:
                with dset.collective:
                    if units:
                        vals = nonDimensionalize(dset[ chunkStart : chunkEnd] * units)
                        ztmp = self.add_particles_with_coordinates(vals)
                    else:
                        ztmp = self.add_particles_with_coordinates(dset[ chunkStart : chunkEnd ])
                    firstChunk = False
            else:
                if units:
                    vals = nonDimensionalize(dset[ chunkStart : chunkEnd] * units)
                    ztmp = self.add_particles_with_coordinates(vals)
                else:
                    ztmp = self.add_particles_with_coordinates(dset[ chunkStart : chunkEnd ])
            tmp = np.copy(ztmp) # copy because ztmp is 'readonly'

            # slice out -neg bits and make the local indices global
            it = np.nditer(tmp, op_flags=['readwrite'], flags=['f_index'])
            while not it.finished:
                if it[0] >= 0:
                    it[0] = chunkStart+it.index # local to global
                it.iternext()

            # slice out -neg bits
            tmp = tmp[tmp[:]>=0]
            # append to valid
            valid = np.append(valid, tmp)

            if rank == 0 and verbose:
                bar.update(chunkEnd)

        h5f.close()
        self._local2globalMap = valid
        # record which swarm state this corresponds to
        self._checkpointMapsToState = self.stateId
