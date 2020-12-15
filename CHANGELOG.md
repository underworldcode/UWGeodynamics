# Version 2.11

- Bug fixes
- Add option to pass UW function (or mesh variable) to the density material properties.
- Add option to initiliase temperature and pressure fields with an UW function or a mesh variable.


# Version 2.10


- Layer2D which has now been deprecated for a while has been removed


# Version 2.9

### User related changes:
- The 'visualisation' python module replaces the 'glucifer' python module (as per UW 2.9) 
To update in your model, replace the line 
     `import glucifer`
  with
     `from UWGeodynamics import visualisation as vis`
  'vis' can be used as 'glucifer' was.

- Compatible with UW 2.9
- New Numpy Pressure Smoother

### Enhancements:
- Upgrade to Dockerfile configuration, using "muli-stage" builds.

### Bug Fixes:
- Passive tracer initialisation in 3D.
- Restarting with track fields and swarm global indices.
- Passive tracer interface changes. vertices must be passed as a numpy array
