from __future__ import print_function,  absolute_import
import underworld as uw
import h5py
from mpi4py import MPI
from UWGeodynamics.scaling import Dimensionalize
from UWGeodynamics.scaling import nonDimensionalize
from UWGeodynamics.scaling import UnitRegistry as u
from UWGeodynamics.version import git_revision as __git_revision__
from . import _meshvariable as var

class FeMesh_Cartesian(uw.mesh.FeMesh_Cartesian):

    def __init__(self, elementType="Q1/dQ0",
                 elementRes=(4,4), minCoord=(0.,0.),
                 maxCoord=(1.,1.), periodic=None,
                 partitioned=True, **kwargs):

        super(FeMesh_Cartesian, self).__init__(elementType,
                                               elementRes,
                                               minCoord,
                                               maxCoord,
                                               periodic,
                                               partitioned,
                                               **kwargs)

    def add_variable(self, nodeDofCount, dataType='double', **kwargs):
        """
        Creates and returns a mesh variable using the discretisation of the given mesh.

        To set / read nodal values, use the numpy interface via the 'data' property.

        Parameters
        ----------
        dataType : string
            The data type for the variable.
            Note that only 'double' type variables are currently
            supported.
        nodeDofCount : int
            Number of degrees of freedom per node the variable will have

        Returns
        -------
        underworld.mesh.MeshVariable
            The newly created mesh variable.

        Example
        -------
        >>> linearMesh  = uw.mesh.FeMesh_Cartesian( elementType='Q1/dQ0', elementRes=(16,16), minCoord=(0.,0.), maxCoord=(1.,1.) )
        >>> scalarFeVar = linearMesh.add_variable( nodeDofCount=1, dataType="double" )
        >>> q0field     = linearMesh.subMesh.add_variable( 1 )  # adds variable to secondary elementType discretisation
        """

        return  var.MeshVariable(self, nodeDofCount, dataType, **kwargs)

    def save(self, filename, units=None, time=None):
        """
        Save the mesh to disk

        Parameters
        ----------
        filename : string
            The name of the output file.

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
        First create the mesh:

        >>> mesh = uw.mesh.FeMesh_Cartesian( elementType='Q1/dQ0', elementRes=(16,16), minCoord=(0.,0.), maxCoord=(1.,1.) )

        Save to a file (note that the 'ignoreMe' object isn't really required):

        >>> ignoreMe = mesh.save("saved_mesh.h5")

        Now let's try and reload. First create new mesh (note the different spatial size):

        >>> clone_mesh = uw.mesh.FeMesh_Cartesian( elementType='Q1/dQ0', elementRes=(16,16), minCoord=(0.,0.), maxCoord=(1.5,1.5) )

        Confirm clone mesh is different from original mesh:

        >>> import numpy as np
        >>> np.allclose(mesh.data,clone_mesh.data)
        False

        Now reload using saved file:

        >>> clone_mesh.load("saved_mesh.h5")

        Now check for equality:

        >>> np.allclose(mesh.data,clone_mesh.data)
        True

        >>> # clean up:
        >>> if uw.rank() == 0:
        ...     import os;
        ...     os.remove( "saved_mesh.h5" )

        """

        if hasattr(self.generator, 'geometryMesh'):
            raise RuntimeError("Cannot save this mesh as it's a subMesh. "
                                + "Most likely you only need to save its geometryMesh")
        if not isinstance(filename, str):
            raise TypeError("'filename', must be of type 'str'")

        h5f = h5py.File(name=filename, mode="w", driver='mpio', comm=MPI.COMM_WORLD)

        fact = 1.0
        if units:
            fact = Dimensionalize(1.0, units=units).magnitude
            h5f.attrs['units'] = str(units)

        # save attributes and simple data - MUST be parallel as driver is mpio
        h5f.attrs['dimensions'] = self.dim
        h5f.attrs['mesh resolution'] = self.elementRes
        h5f.attrs['max'] = tuple([fact*x for x in self.maxCoord])
        h5f.attrs['min'] = tuple([fact*x for x in self.minCoord])
        h5f.attrs['regular'] = self._cself.isRegular
        h5f.attrs['elementType'] = self.elementType
        h5f.attrs['time'] = str(time)
        h5f.attrs["git commit"] = __git_revision__

        # write the vertices
        globalShape = ( self.nodesGlobal, self.data.shape[1] )
        dset = h5f.create_dataset("vertices",
                                  shape=globalShape,
                                  dtype=self.data.dtype)

        local = self.nodesLocal
        # write to the dset using the local set of global node ids
        with dset.collective:
            dset[self.data_nodegId[0:local],:] = self.data[0:local] * fact

        # write the element node connectivity
        globalShape = ( self.elementsGlobal, self.data_elementNodes.shape[1] )
        dset = h5f.create_dataset("en_map",
                                  shape=globalShape,
                                  dtype=self.data_elementNodes.dtype)

        local = self.elementsLocal
        # write to the dset using the local set of global node ids
        with dset.collective:
            dset[self.data_elgId[0:local], :] = self.data_elementNodes[0:local]

        h5f.close()

        # return our file handle
        return uw.utils.SavedFileData(self, filename)

    def load(self, filename):
        """
        Load the mesh from disk.

        Parameters
        ----------
        filename: str
            The filename for the saved file. Relative or absolute paths may be
            used, but all directories must exist.

        Notes
        -----
        This method must be called collectively by all processes.

        If the file data array is the same length as the current mesh
        global size, it is assumed the file contains compatible data. Note that
        this may not be the case where for example where you have saved using a
        2*8 resolution mesh, then loaded using an 8*2 resolution mesh.

        Provided files must be in hdf5 format, and use the correct schema.

        Example
        -------
        Refer to example provided for 'save' method.

        """
        self.reset()
        if not isinstance(filename, str):
            raise TypeError("Expected filename to be provided as a string")

        # get field and mesh information
        h5f = h5py.File( filename, "r", driver='mpio', comm=MPI.COMM_WORLD );

        # get resolution of old mesh
        res = h5f.attrs['mesh resolution']
        if res is None:
            raise RuntimeError("Can't read the 'mesh resolution' for the field hdf5 file,"+
                   " was it created correctly?")

        if (res == self.elementRes).all() == False:
            raise RuntimeError("Provided file mesh resolution does not appear to correspond to\n"\
                               "resolution of mesh object.")

        # get units
        try:
            units = h5f.attrs["units"]
        except KeyError:
            units = None

        if units and units != "None":
            units = u.parse_expression(units)
        else:
            units = None

        dset = h5f.get('vertices')
        if dset == None:
            raise RuntimeError("Can't find the 'vertices' dataset in hdf5 file '{0}'".format(filename) )

        dof = dset.shape[1]
        if dof != self.data.shape[1]:
            raise RuntimeError("Can't load hdf5 '{0}', incompatible data shape".format(filename))

        if len(dset) != self.nodesGlobal:
            raise RuntimeError("Provided data file appears to be for a different resolution mesh.")

        with self.deform_mesh(isRegular=h5f.attrs['regular']):
            with dset.collective:
                if units:
                    self.data[0:self.nodesLocal] = nonDimensionalize(
                        dset[self.data_nodegId[0:self.nodesLocal], :] * units)
                else:
                    self.data[0:self.nodesLocal] = dset[self.data_nodegId[0:self.nodesLocal], :]

        h5f.close()
