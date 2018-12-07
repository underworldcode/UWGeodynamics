
if __name__ == "__main__":
    import doctest
    import UWGeodynamics
    doctest.testmod(UWGeodynamics._model)
    doctest.testmod(UWGeodynamics._boundary_conditions)
    doctest.testmod(UWGeodynamics._freesurface)
    doctest.testmod(UWGeodynamics._frictional_boundary)
    doctest.testmod(UWGeodynamics._material)
    doctest.testmod(UWGeodynamics._melt)
    doctest.testmod(UWGeodynamics._mesh_advector)
    doctest.testmod(UWGeodynamics._rheology)
    doctest.testmod(UWGeodynamics.shapes)
    doctest.testmod(UWGeodynamics.surfaceProcesses)
    doctest.testmod(UWGeodynamics._utils)

