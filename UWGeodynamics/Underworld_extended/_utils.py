import os
import h5py
from underworld.utils._utils import _xdmfAttributeschema


_dtypes_to_xdmf = {
    '<f8': ("Float", "8"),
    '<f4': ("Float", "4"),
    '<f1': ("Float", "1"),
    '<i8': ("Int", "8"),
    '<i4': ("Int", "4"),
    '<i2': ("Int", "2"),
    '<i1': ("Int", "1"),
    '<u8': ("UInt", "8"),
    '<u4': ("UInt", "4"),
    '<u2': ("UInt", "2"),
    '<u1': ("UInt", "1"),
}


def _swarmvarschema(varSavedData, varname):
    """"
    Writes the attribute schema for a swarm variable xdmf file

    Parameters:
    ----------
    varSavedData : SavedFileData
        The SavedFileData handle to the saved SwarmVariable
    varname : str
        The name, in xdmf, to give the SwarmVariable

    Returns
    -------
    out : str
        string containing the xdmf schema
    """

    # retrieve bits from varSavedData
    var = varSavedData.pyobj
    varfilename = varSavedData.filename
    refName = os.path.basename(varSavedData.filename)

    # set parameters - serially open the varfilename
    h5f = h5py.File(name=varfilename, mode="r")
    dset = h5f.get('data')
    if dset is None:
        raise RuntimeError("Can't find 'data' in file '{}'.\n".format(varfilename))
    globalCount = len(dset)
    h5f.close()

    dof_count = var.data.shape[1]
    vartype, precision = _dtypes_to_xdmf[var.data.dtype.str]
    variableType = "NumberType=\"{0}\" Precision=\"{1}\"".format(vartype, precision)

    out = _xdmfAttributeschema(varname, variableType, "Node", globalCount, dof_count, refName )

    return out
