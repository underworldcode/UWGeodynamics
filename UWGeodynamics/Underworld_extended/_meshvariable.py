from __future__ import print_function,  absolute_import
import underworld as uw
import h5py
import numpy as np
import os
from mpi4py import MPI
from UWGeodynamics.scaling import Dimensionalize
from UWGeodynamics.scaling import nonDimensionalize
from UWGeodynamics.scaling import UnitRegistry as u
from UWGeodynamics.version import git_revision as __git_revision__


class MeshVariable(uw.mesh.MeshVariable):

    def __init__(self, mesh, nodeDofCount, dataType="double", **kwargs):
        super(MeshVariable, self).__init__(mesh, nodeDofCount, dataType, **kwargs)

    def load(self, filename, interpolate=False):
        """
        Load the MeshVariable from disk.

        Parameters
        ----------
        filename: str
            The filename for the saved file. Relative or absolute paths may be
            used, but all directories must exist.
        interpolate: bool
            Set to True to interpolate a file containing different resolution data.
            Note that a temporary MeshVariable with the file data will be build
            on **each** processor. Also note that the temporary MeshVariable
            can only be built if its corresponding mesh file is available.
            Also note that the supporting mesh mush be regular.

        Notes
        -----
        This method must be called collectively by all processes.

        If the file data array is the same length as the current variable
        global size, it is assumed the file contains compatible data. Note that
        this may not be the case where for example where you have saved using a
        2*8 resolution mesh, then loaded using an 8*2 resolution mesh.

        Provided files must be in hdf5 format, and use the correct schema.

        Example
        -------
        Refer to example provided for 'save' method.

        """
        if not isinstance(filename, str):
            raise TypeError("Expected filename to be provided as a string")

        # get field and mesh information
        h5f = h5py.File( filename, "r", driver='mpio', comm=MPI.COMM_WORLD );
        dset = h5f.get('data')

        # get units
        try:
            units = h5f.attrs["units"]
        except KeyError:
            units = None

        if units and units != "None":
            units = u.parse_expression(units)
        else:
            units = None

        if dset == None:
            raise RuntimeError("Can't find the 'data' in hdf5 file '{0}'".format(filename) )

        dof = dset.shape[1]
        if dof != self.data.shape[1]:
            raise RuntimeError("Can't load hdf5 '{0}', incompatible data shape".format(filename))

        if len(dset) == self.mesh.nodesGlobal:

            # assume dset matches field exactly
            mesh = self.mesh
            local = mesh.nodesLocal

            with dset.collective:
                self.data[0:local] = dset[mesh.data_nodegId[0:local],:]

        else:
            if not interpolate:
                raise RuntimeError("Provided data file appears to be for a different resolution MeshVariable.\n"\
                                   "If you would like to interpolate the data to the current variable, please set\n" \
                                   "the 'interpolate' parameter. Check docstring for important caveats of interpolation method.")

            # if here then we build a local version of the entire file field and interpolate it's values

            # first get file field's mesh
            if h5f.get('mesh') == None:
                raise RuntimeError("The hdf5 field to be loaded with interpolation must have an associated "+
                        "'mesh' hdf5 file. Resave the field with its associated mesh."+
                        "i.e. myField.save(\"filename.h5\", meshFilename)" )
            # get resolution of old mesh
            res = h5f['mesh'].attrs.get('mesh resolution')
            if res is None:
                raise RuntimeError("Can't read the 'mesh resolution' for the field hdf5 file,"+
                       " was it created correctly?")

            # get max of old mesh
            inputMax = h5f['mesh'].attrs.get('max')
            if inputMax is None:
                raise RuntimeError("Can't read the 'max' for the field hdf5 file,"+
                       " was it created correctly?")

            inputMin = h5f['mesh'].attrs.get('min')
            if inputMin is None:
                raise RuntimeError("Can't read the 'min' for the field hdf5 file,"+
                       " was it created correctly?")
            regular = h5f['mesh'].attrs.get('regular')
            if regular and regular!=True:
                raise RuntimeError("Saved mesh file appears to correspond to a irregular mesh.\n"\
                                   "Interpolating from irregular mesh not currently supported." )

            elType = h5f['mesh'].attrs.get('elementType')
            # for backwards compatiblity, the 'elementType' attribute was added Feb2017
            if elType == None:
                elType = 'Q1'

            # build the NON-PARALLEL field and mesh
            inputMesh = uw.mesh.FeMesh_Cartesian( elementType = (elType+"/DQ0"), # only geometryMesh can be saved
                                          elementRes  = tuple(res),
                                          minCoord    = tuple(inputMin),
                                          maxCoord    = tuple(inputMax),
                                          partitioned=False)

            # load data onto MeshVariable
            if len(dset) == inputMesh.nodesGlobal:
                inputField = uw.mesh.MeshVariable( mesh=inputMesh, nodeDofCount=dof )
            elif  dset.shape[0] == inputMesh.subMesh.nodesGlobal:
                # load as a subMesh
                # assume the dset field belongs to the subMesh
                inputField = uw.mesh.MeshVariable( mesh=inputMesh.subMesh, nodeDofCount=dof )
            else:
                # raise error
                raise RuntimeError("The saved mesh file can't be read onto the interpolation grid.\n" \
                                   "Note: only subMesh variable with elementType 'DQ0' can be used presently used")

            # copy hdf5 numpy array onto serial inputField
            inputField.data[:] = dset[:]

            # interpolate 'inputField' onto the self nodes
            self.data[:] = inputField.evaluate(self.mesh.data)

        if units:
            if units.units == "degC":
                units = u.degK
                self.data[:] = nonDimensionalize((self.data[:]+273.15)*units)
            else:
                self.data[:] = nonDimensionalize(self.data[:]*units)

        uw.libUnderworld.StgFEM._FeVariable_SyncShadowValues( self._cself )
        h5f.close()

    def save(self, filename, meshHandle=None, units=None, time=None):
        """
        Save the MeshVariable to disk.

        Parameters
        ----------
        filename : string
            The name of the output file. Relative or absolute paths may be
            used, but all directories must exist.
        meshHandle :uw.utils.SavedFileData , optional
            The saved mesh file handle. If provided, a link is created within the
            mesh variable file to this saved mesh file. Important for checkpoint when
            the mesh deforms.

        Notes
        -----
        This method must be called collectively by all processes.

        Returns
        -------
        underworld.utils.SavedFileData
            Data object relating to saved file. This only needs to be retained
            if you wish to create XDMF files and can be ignored otherwise.

        Example
        -------
        First create the mesh add a variable:

        >>> mesh = uw.mesh.FeMesh_Cartesian( elementType='Q1/dQ0', elementRes=(16,16), minCoord=(0.,0.), maxCoord=(1.,1.) )
        >>> var = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1, dataType="double" )

        Write something to variable

        >>> import numpy as np
        >>> var.data[:,0] = np.arange(var.data.shape[0])

        Save to a file (note that the 'ignoreMe' object isn't really required):

        >>> ignoreMe = var.save("saved_mesh_variable.h5")

        Now let's try and reload.

        >>> clone_var = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1, dataType="double" )
        >>> clone_var.load("saved_mesh_variable.h5")

        Now check for equality:

        >>> np.allclose(var.data,clone_var.data)
        True

        >>> # clean up:
        >>> if uw.rank() == 0:
        ...     import os;
        ...     os.remove( "saved_mesh_variable.h5" )

        """
        if not isinstance(filename, str):
            raise TypeError("Expected 'filename' to be provided as a string")

        mesh = self.mesh
        h5f = h5py.File(name=filename, mode="w", driver='mpio', comm=MPI.COMM_WORLD)

        # ugly global shape def
        globalShape = ( mesh.nodesGlobal, self.data.shape[1] )
        # create dataset
        dset = h5f.create_dataset("data",
                                  shape=globalShape,
                                  dtype=self.data.dtype)
        fact = 1.0
        if units:
            fact = Dimensionalize(1.0, units=units).magnitude
            if units == "degC":
                fact = Dimensionalize(1.0, units=u.degK).magnitude
            # Save unit type as attribute
            h5f.attrs['units'] = str(units)

        if time:
            h5f.attrs['time'] = str(time)

        h5f.attrs["git commit"] = __git_revision__

        # write to the dset using the global node ids
        local = mesh.nodesLocal

        with dset.collective:
            if units == "degC":
                dset[mesh.data_nodegId[0:local],:] = self.data[0:local] * fact - 273.15
            else:
                dset[mesh.data_nodegId[0:local],:] = self.data[0:local] * fact

        # save a hdf5 attribute to the elementType used for this field - maybe useful
        h5f.attrs["elementType"] = np.string_(mesh.elementType)

        if hasattr( mesh.generator, "geometryMesh"):
            mesh = mesh.generator.geometryMesh

        if meshHandle:
            if not isinstance(meshHandle, (str, uw.utils.SavedFileData)):
                raise TypeError("Expected 'meshHandle' to be of type 'uw.utils.SavedFileData'")

            if isinstance(meshHandle, str):
                # DEPRECATION check
                import warnings
                warnings.warn("'meshHandle' paramater should be of type uw.utils.SaveFileData. Please update your models. "+
                              "Accepting 'meshHandle' as a string parameter will be removed in the next release.")
                meshFilename = meshHandle
            else:
                meshFilename = meshHandle.filename

            if not os.path.exists(meshFilename):
                raise ValueError("You are trying to link against the mesh file '{}'\n\
                                  that does not appear to exist. If you need to link \n\
                                  against a mesh file, please make sure it is created first.".format(meshFilename))
            # set reference to mesh (all procs must call following)
            h5f["mesh"] = h5py.ExternalLink(meshFilename, "./")

        h5f.close()

        # return our file handle
        return uw.utils.SavedFileData(self, filename)
