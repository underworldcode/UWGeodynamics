UWGeodynamics User Guide
========================

Docker
------

We strongly encourage users to run UWGeodynamics using the docker images
we provide on `Docker
Hub <https://hub.docker.com/r/underworldcode/uwgeodynamics>`__.

Docker containers provide and easy-way to set up and distribute
applications. They also provide a safe and consistent environment which
facilitate debugging and reproducibility of models. The image we provide
contains all the dependencies and configuration files required to run
Underworld models. Users can start developping model as soon as they
have downloaded the image, independently of the operating system running
on their machine.

Docker container using `Kitematic <https://kitematic.com/>`__
-------------------------------------------------------------

Getting UWGeodynamics docker image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download and Install `Kitematic <https://kitematic.com/>`__ The
   software is available for Windows, MacOsx and Linux. Be aware that on
   linux the installation may differ depending on the distribution you
   are running.

2. Open Kitematic and search for the **uwgeodynamics** image.
3. Create a container by clicking on the create button.

Starting a container:
~~~~~~~~~~~~~~~~~~~~~

You should now have a container appearing on the left side of your
kitematic window. The first thing to do now is to create a link between
a local directory (A directory on your physical hard drive) and a volume
directory inside the docker container. A volume is a special directory
that can be accessed from outside the container. It is the location you
will use to save your results.

Local Installation
------------------

Requirements
------------

-  Python >= 3.5
-  A Working version of Underworld2 >=2.7.0 (Please refer to the
   Underworld documentation)
-  pint >= 0.8

Install UWGeodynamics
---------------------

Pip install (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~

The UWGeodynamics module can be installed directly from the Python
Package Index:

::

       pip3 install UWGeodynamics

Install from sources
~~~~~~~~~~~~~~~~~~~~

The module source files are available through
`github <https://github.com/rbeucher/UWGeodynamics.git>`__

::

       git clone https://github.com/rbeucher/UWGeodynamics.git

It can then be installed globally on your system using

::

       pip3 install -e UWGeodynamics

The Jupyter notebook
--------------------

The [[Jupyter|https://jupyter.org]] notebook provides a powerfull
environment for the development and analysis of Underworld model.

Importing UWGeodynamics
-----------------------

*UWGeodynamics* can be imported as follow:

.. code:: python

   import UWGeodynamics as GEO

Working with units
------------------

*UWGeodynamics* uses *[[Pint|https://pint.readthedocs.io/en/latest]]*, a
Python package to define, operate and manipulate **physical quantities**
(A numerical value with unit of measurement). Pint is a very powerful
package that handles conversion and operation between units.

We recommend using SI units but other systems are also available.

*[[Pint|https://pint.readthedocs.io/en/latest]]* **Unit Registry** can
be used as follow:

.. code:: python

   import UWGeodynamics as GEO
   u = GEO.UnitRegistry

or simply

.. code:: python

   u = GEO.u

You can have a quick overview of all the units available by hitting tab
after the “.” of the u object.

[[img/tabtab.gif]]

Quantities can then be defined as follow:

.. code:: python

   length = 100. * u.kilometer
   width = 50. * u.kilometer
   gravity = 9.81 * u.meter / u.second**2

*[[Pint|https://pint.readthedocs.io/en/latest]]* offers the possibility
to append a prefix to the units. 1 million year can thus be defined as
follow:

.. code:: python

   1.0 * u.megayear

Model Scaling
-------------

Model can be scaled using a series of scaling coefficients

.. code:: python

   import UWGeodynamics as GEO

The default scaling coefficients are defined as follow:

+---------------+--------------+
| Dimension     | value        |
+===============+==============+
| [mass]        | 1.0 kilogram |
+---------------+--------------+
| [length]      | 1.0 meter    |
+---------------+--------------+
| [temperature] | 1.0 kelvin   |
+---------------+--------------+
| [time]        | 1.0 second   |
+---------------+--------------+
| [substance]   | 1.0 mole     |
+---------------+--------------+

The scaling value can be changed by accessing each scaling coefficient
as follow

.. code:: python

   u = GEO.u

   GEO.scaling_coefficients["[length]"] = 3.0 * u.kilometer
   GEO.scaling_coefficients["[mass]"] = 4.0 * u.kilogram
   GEO.scaling_coefficients["[temperature]"] = 273.15 * u.degK
   GEO.scaling_coefficients["[time]"] = 300 * u.years

The unit entered are checked internally and an error is raised if the
units are incompatible. The value is automatically converted to the base
units (meter, second, degree, etc).

Scaling a Model
---------------

To scale a model, the user must define a serie of characteristic
physical values and assign them to the scaling object.

Arguments with units will be scaled by the UWGeodynamics functions.

.. code:: python

   KL = 100 * u.kilometer  # Characteristic length
   Kt = 1. * u.year        # Characteristic time
   KM = 3000. * u.kilogram # Characteristic mass
   KT = 1200. * u.degK     # Characteristic temperature

   GEO.scaling_coefficients["[length]"] = KL
   GEO.scaling_coefficients["[time]"] = Kt
   GEO.scaling_coefficients["[mass]"]= KM
   GEO.scaling_coefficients["[temperature]"] = KT

Tools
~~~~~

It is sometime necessary to scale or convert values back to units.

We provide 2 function to process the conversion:

.. code:: python

   GEO.nonDimensionalize
   GEO.Dimensionalize

The nonDimensionalize function is also available as:

.. code:: python

   GEO.nd

Example
~~~~~~~

1. define a length of 300 kilometers.
2. use the GEO.nd function to scale it.
3. convert the value back to SI units.

.. code:: python

   length = 300. * u.kilometers
   scaled_length = GEO.nd(length)
   length_meters = GEO.Dimensionalize(scaled_length, u.meters)

The Model object
----------------

The central element or “object” of the UWGeodynamics module is the
**Model** object.

It has several uses: - It defines the extent and the outside geometry of
your problem. - It works as a container for the field variables.

It basically defines the universe on which you are going to apply
physical rules (Gravity field, boundary condition, composition,
temperature etc.) It is the equivalent of the box in which you would put
the sand and silicon if you were to build an analog experiment in a lab.
One important difference is that the “box” his not empty, it is
populated with particles that have already some properties. The
properties are changed by defining new materials.

.. code:: python

   import UWGeodynamics as GEO
   u = GEO.u
   Model = GEO.Model(elementRes=(64, 64), 
                     minCoord=(0. * u.kilometer, 0. * u.kilometer),
                     maxCoord=(64. * u.kilometer, 64. * u.kilometer))

The Material object
-------------------

The *UWGeodynamics* module is designed around the idea of materials,
which are essentially a way to define physical properties across the
Model domain.

Materials are defined using the **Material** object as follow:

.. code:: python

   import UWGeodynamics as GEO

   crust = GEO.Material(name="Crust")

Typing the name of the material in an empty cell will return a table
which summarizes the property of the material:

[[img/Material1.png]]

As you can see, most of the property are undefined.

They are several ways to define the physical parameters of our Material.

-  The first one is to add them directly when creating the object
   itself:

.. code:: python

   import UWGeodynamics as GEO

   u = GEO.u
   crust = GEO.Material(name="Crust", density=3000*u.kilogram/u.metre**3)

-  The second option is to change the property after creating the
   **Material**:

.. code:: python

   import UWGeodynamics as GEO

   u = GEO.u
   crust = GEO.Material(name="Crust")
   crust.density = 3000. * u.kilogram / u.metre **3

The second option is often easier to read.

**UWGeodynamics contains some basic dimensionality checks. Entering
wrong units will raise an error**

Material can be added to a model as follow:

.. code:: python

   import UWGeodynamics as GEO
   u = GEO.u
   Model = GEO.Model()
   crust = Model.add_material(name="Crust")

Although optional, tt is a good idea to give a **name** to the material.
The **Model.add_material** method will return a Material object. That
object is a python object that will then be used to define the property
of the material.

Material shape
~~~~~~~~~~~~~~

A material (or a phase) is first defined by the space it takes in the
box (its shape).

There is a range of shapes available

2D:

-  `Layer <#layer>`__ (2D/3D)
-  `Polygon <#polygon>`__ (2D)
-  `Box <#box>`__ (2D)
-  `Disk <#disk>`__ (2D) / `Spheres <#spheres>`__ (3D)
-  `Annulus <#annulus>`__ (2D)
-  `MultiShape <#multishape>`__ (Combination of any of the above) (2D)
-  `HalfSpace <#halfspace>`__ (3D)

Layer
^^^^^

.. code:: python

   import UWGeodynamics as GEO

   u = GEO.u
   Model = GEO.Model()
   shape = GEO.shapes.Layer(top=30.*u.kilometer, bottom=0.*u.kilometer)
   material = Model.add_material(name="Material", shape=shape)

   Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   Fig.save("layers.png")

[[/img/layers.png]]

Polygon
^^^^^^^

.. code:: python

   import UWGeodynamics as GEO

   u = GEO.u
   Model = GEO.Model()
   polygon = GEO.shapes.Polygon(vertices=[(10.* u.kilometer, 10.*u.kilometer),
                                          (20.* u.kilometer, 35.*u.kilometer),
                                          (35.* u.kilometer, 5.*u.kilometer)])
   material = Model.add_material(name="Material", shape=polygon)

   Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   Fig.save("polygon.png")

[[/img/polygon.png]]

Box
^^^

.. code:: python

   import UWGeodynamics as GEO

   u = GEO.u
   Model = GEO.Model()
   box = GEO.shapes.Box(top=10.* u.kilometer, bottom=5*u.kilometer,
                        minX=10.*u.kilometer, maxX=15*u.kilometer)
   material = Model.add_material(name="Material", shape=box)

   Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   Fig.save("box.png")

[[/img/box.png]]

Disk
^^^^

.. code:: python

   import UWGeodynamics as GEO

   u = GEO.u
   Model = GEO.Model()
   disk = GEO.shapes.Disk(center=(32. * u.kilometer, 32. * u.kilometer), radius=10.*u.kilometer)
   material = Model.add_material(name="Material", shape=disk)

   Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   Fig.save("disk.png")

[[/img/disk.png]]

Annulus
^^^^^^^

.. code:: python

   import UWGeodynamics as GEO

   u = GEO.u
   Model = GEO.Model()
   annulus = GEO.shapes.Annulus(center=(35.*u.kilometer, 50.*u.kilometer),
                                r1=5.*u.kilometer, 
                                r2=10.*u.kilometer)
   material = Model.add_material(name="Material", shape=annulus)

   Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   Fig.save("annulus.png")

[[/img/annulus.png]]

MultiShape
^^^^^^^^^^

Several shapes can be combined to form a material shape:

.. code:: python

   import UWGeodynamics as GEO

   u = GEO.u
   Model = GEO.Model()
   disk1 = GEO.shapes.Disk(center=(10. * u.kilometer, 10. * u.kilometer), radius=10.*u.kilometer)
   disk2 = GEO.shapes.Disk(center=(20. * u.kilometer, 20. * u.kilometer), radius=5.*u.kilometer)

   shape = GEO.shapes.MultiShape([disk1, disk2])
   material = Model.add_material(name="Material", shape=shape)
   Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   Fig.save("multishape.png")

[[/img/multishape.png]]

Material Attributes
~~~~~~~~~~~~~~~~~~~

.. code:: python

   Model.density = 200. * u.kg / u.m**3
   material.density = 3000 * u.kilogram / u.meter**3
   Fig = Model.plot.density(figsize=(400, 400))

[[/img/density.png]]

Material Rheologies
-------------------

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u

Viscous Rheology
~~~~~~~~~~~~~~~~

Registry / Database
^^^^^^^^^^^^^^^^^^^

[[/img/ViscousCreepRegistry.gif]]

User Defined
^^^^^^^^^^^^

Constant viscosity
''''''''''''''''''

.. code:: python

   >>> viscosity = 1e21 * u.pascal * u.second
   >>> viscosity = GEO.ConstantViscosity(1e21 * u.pascal * u.second)

Viscous Creep
'''''''''''''

.. code:: python

   >>> viscosity = GEO.ViscousCreep(preExponentialFactor=1.0, 
                                    stressExponent=1.0, 
                                    activationVolume=0., 
                                    activationEnergy=200 * u.kilojoules,
                                    waterFugacity=0.0,
                                    grainSize=0.0,
                                    meltFraction=0.,
                                    grainSizeExponent=0., 
                                    waterFugacityExponent=0., 
                                    meltFractionFactor=0.0,
                                    f=1.0)

Single parameters can then be modified

.. code:: python

   >>> viscosity.activationEnergy = 300. * u.kilojoule

Plastic Rheology (Yield)
~~~~~~~~~~~~~~~~~~~~~~~~

.. _registry-database-1:

Registry / Database
^^^^^^^^^^^^^^^^^^^

[[/img/PlasticityRegistry.gif]]

.. _user-defined-1:

User Defined
^^^^^^^^^^^^

.. code:: python

   >>> plasticity = GEO.DruckerPrager(cohesion=10. * u.megapascal, 
                                      cohesionAfterSoftening=10. * u.megapascal, 
                                      frictionCoefficient = 0.3,
                                      frictionAfterSoftening = 0.2,
                                      epsilon1=0.5,
                                      epsilon2=1.5)

Mechanical Boundary Conditions
------------------------------

Kinematic or mechanichal boundary conditions are a critical part of any
geodynamic model design. In the following, we quickly detail the options
available to define boundary conditions in Underworld using the
UWGeodynamics module.

How to define boundary conditions and how to make sure those are
consistent are questions beyond the scope of this manual.

We will define a simple model for the sake of the example.

.. code:: python

   >>> import UWGeodynamics as GEO

   >>> u = GEO.u

   >>> Model = GEO.Model(elementRes=(64, 64),
                         minCoord=(0. * u.kilometer, 0. * u.kilometer),
                         maxCoord=(64. * u.kilometer, 64. * u.kilometer))

Kinematic or mechanical boundary conditions are set using the
**set_velocityBCs** method. Conditions are defined for each wall (left,
right, bottom, top, back and front (3D only)). For each wall, the user
must define the condition for each degree of freedom (2 in 2D (x,y), 3
in 3D (x,y,z).

if :math:`V` is a vector (:math:`V_x`, :math:`V_y`, :math:`V_z`) that we
want to apply on the left wall, the *left* parameter must be defined as
follow:

::

   left=[Vx, Vy, Vz]

In the following example we set the boundary condition to be:

-  left wall: \\V_x = -1.0\\ :raw-latex:`\text{cm / yr}`$,
   :math:`Vy=None`
-  right wall: :math:`V_x = 1.0 \text{cm / yr}`, :math:`Vy=None`
-  bottom wall: :math:`V_x = None`, :math:`V_y= 0.` (free slip)

It is an extension model with a total rate of extension equal to 2.0
centimetre / year. No :math:`V_x` is prescribed at the bottom, while
:math:`V_y` is set to :math:`0.` no material will be able to enter or
leave the model domain from that side. The material is free to move
vertically along the side walls.

.. code:: python

   >>> Model.set_velocityBCs(left=[1.0*u.centimetre/u.year, None],
                             right=[-1.0*u.centimetre/u.year, None],
                             bottom=[None, 0.],
                             top=[None,0.])

   >>> Fig = Model.plot.velocityField()

[[/img/mechanicalBCs1.png]]

Support Conditions
------------------

Isostasy is an important concept in geodynamics. It is essentially a
consequence of the redistribution of mass within a deforming Earth. One
important limitation of our geodynamic model is that we model special
cases inside rectangular boxes while earth is actually a sphere. One may
however need to provide a way to maintain the volume / mass inside the
domain in order to mimic isostasy. There is no ideal way to model
isostasy in a boxed model, it is however possible to approach isostasy
using a support condition.

Options are to:

-  Balance flows using a kinematic condition at the base of the model.
-  Balance flows using a stress condition at the base of the model.
-  Balance flows along the sides.

Lecode Isostasy
~~~~~~~~~~~~~~~

The Lecode Isostasy submodule provides a way to model isostatic support
at the base of the model. It calculates the velocity to apply at the
base of each elemental column. It applies the principles of Airy
isostatic model by approximating the weight of each column. The
calculation is done dynamically and velocities will change from one step
to the next. It is a good option to use in most cases.

The option can be used by creating a LecodeIsostasy object using the
``GEO.LecodeIsostasy`` class. The object requires the index of the
material of reference (the material number). One can apply an average
velocity (calculated across each column base) using the ``average``
parameter (default to False).

.. code:: python

   >>> Model.set_velocityBCs(left=[1.0*u.centimetre/u.year, None],
                             right=[-1.0*u.centimetre/u.year, None],
                             bottom=[None, GEO.LecodeIsostasy(reference_mat=Model.index)],
                             top=[None,0.])

Traction Condition (Stress)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another approach to model isostasy is to defined a stress at the base of
the model. This is done using units of stress (derived SI units =
pascal). The model will then maintain the stress by adjusting the flow
(and thus velocities) across the border. Units are important: units of
stress defines a stress while units of velocity define a velocity.

.. code:: python

   >>> Model.set_velocityBCs(left=[1.0*u.centimetre/u.year, None],
                             right=[-1.0*u.centimetre/u.year, None],
                             bottom=[None, 10.*u.megapascal],
                             top=[None,0.])

Velocity varying along a wall
-----------------------------

It is sometime necessary to define a velocity only for a section of a
wall. That can be done using a **condition**. A condition is a set of
rule to apply on a wall.

As an example, we will apply a velocity of :math:`5.0\text{cm/yr}` for
the part of the left wall below 32 kilometre. Velocity is set to be
:math:`1.0\text{cm/yr}` above.

.. code:: python

   >>> conditions = [(Model.y < GEO.nd(32 * u.kilometer), GEO.nd(5.0 * u.centimeter/u.year)),
                     (True, GEO.nd(1.0*u.centimeter/u.year))]
     
   >>> Model.set_velocityBCs(left=[conditions, None],
                             right=[-1.0*u.centimetre/u.year, None],
                             bottom=[None, 10.*u.megapascal],
                             top=[None,0.])
   >>> Fig = Model.plot.velocityField()

[[/img/mechanicalBCs2.png]]

nodeSets
~~~~~~~~

(to be implemented)

3D Model
--------

Defining boundary conditions for a 3D model is no different than above.
The user must define the velocity components with 3 degree of freedom
instead of 2.

.. code:: python

   >>> Model2 = GEO.Model(elementRes=(16, 16, 16),
                          minCoord=(0. * u.kilometer, 0. * u.kilometer, 0. * u.kilometer),
                          maxCoord=(64. * u.kilometer, 64. * u.kilometer, 64. * u.kilometer))

.. code:: python

   >>> Model2.set_velocityBCs(left=[1.0*u.centimetre/u.year, None, 0.],
                              right=[-1.0*u.centimetre/u.year, None, 0.],
                              bottom=[None, None, 0.],
                              top=[None, None, 0.],
                              front=[None, 0., None],
                              back=[None, 0., None])

Thermal Boundary Conditions
---------------------------

Setting absolute temperatures at the boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setting the temperature at the top of a model to be
:math:`500 \text{kelvin}` at the top and :math:`1600 \text{kelvin}` at
the bottom is done as follow.

.. code:: python

   >>> Model.set_temperatureBCs(top=500. * u.degK, bottom=1600. * u.degK)

[[/img/thermalBCs1.png]]

You can of course define temperatures on the sidewalls:

.. code:: python

   >>> Model.set_temperatureBCs(right=500. * u.degK, left=1600. * u.degK)

[[/img/thermalBCs2.png]]

Setting Heat flux at the boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   >>> Model.set_temperatureBCs(top=500. * u.degK, bottom=-0.22 * u.milliwatt / u.metre**2, bottom_material=Model)

Fix the temperature of a Material
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   >>> Model.set_temperatureBCs(top=500. * u.degK, bottom=-0.22 * u.milliwatt / u.metre**2, bottom_material=Model,
                                materials=[(air, 273. * u.Kelvin)])

Fix the temperature of some nodes using an IndexSet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can assign a temperature to a list of nodes by passing a list of
node indices (global).

.. code:: python

   >>> nodes = [0, 1, 2]
   >>> Model.set_temperatureBCs(top=500. * u.degK, bottom=-0.22 * u.milliwatt / u.metre**2, bottom_material=Model,
                                nodeSets=[(nodes, 273. * u.Kelvin)])

Frictional Boundaries
---------------------

Frictional Boundaries can be set as follow:

.. code:: python

   Model.set_frictional_boundary(left=True, 
                                 right=True, 
                                 bottom=True, 
                                 top=False, 
                                 friction=19.0, 
                                 thickness=3)

Where *left*, *right*, *top*, *bottom*, parameters are the side you want
to apply a frictional boundary condition on. *friction* is the angle of
friction (in degrees). *thickness* is the thickness of the boundary.

Surface Processes
-----------------

A range of basic surface processes function are available from the
*surfaceProcesses* submodule. Surface processes are turned on once you
have passed a valid surface processes function to the
``surfaceProcesses`` method of the ``Model`` object.

Example:

.. code:: python

     import UWGeodynamics as GEO

     Model.surfaceProcesses = GEO.surfaceProcesses.SedimentationThreshold(air=[air], sediment=[sediment], threshold=0. * u.meter)

Three simple function are available:

1. Total Erosion Above Threshold (``ErosionThreshold``).
2. Total Sedimentation Below Threshold (``SedimentationThreshold``)
3. Combination of the 2 above. (``ErosionAndSedimentationThreshold``)

Coupling with Badlands
~~~~~~~~~~~~~~~~~~~~~~

UWGeodynamics provide a way to couple an Underworld model to Badlands.
**More documentation needed**

.. code:: python

   import UWGeodynamics as GEO

   Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(
       airIndex=[air.index], sedimentIndex=sediment.index,
       XML="ressources/badlands.xml", resolution=1. * u.kilometer,
       checkpoint_interval=0.01 * u.megayears)

Passive Tracers and Grids
-------------------------

Passive tracers
~~~~~~~~~~~~~~~

.. code:: python

   import UWGeodynamics as GEO

   u = GEO.u

   Model = GEO.Model(elementRes=(64,64),
                     minCoord=(0.*u.kilometer, 0.* u.kilometer),
                     maxCoord=(64.* u.kilometer, 64 * u.kilometer))

   x = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), 1000)
   y = 32. * u.kilometer

   P = Model.add_passive_tracers(vertices=[x,y])

Tracer patterns / Finite Strain Ellipsis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Visugrid
~~~~~~~~

Running the Model
-----------------

Model initialization
~~~~~~~~~~~~~~~~~~~~

Once your model is set up and initialize. You can run it using the
*Model.run_for* method.

You have 2 options:

1. Run the model for some given number of steps:

.. code:: python

      Model.run_for(nstep=10)

1. Specify an endTime

.. code:: python

      Model.run_for(endTime=1.0* u.megayears)
      # which is equivalent to
      Model.run_for(1.0*u.megayears)

Specify a timestep
^^^^^^^^^^^^^^^^^^

UWGeodynamics calculates the time step automatically based on some
numerical stability criteria. You can force a specific time step or
force the time step to be constant throughou

Saving data
^^^^^^^^^^^

As your model is running you will need to save the results to files.

The *Model.run_for* command provides a series of arguments to help you
save the results at some regular. You can define:

1. A *checkpoint_interval*

.. code:: python

      Model.run_for(endTime=1.0*u.megayears, 
                    checkpoint_interval=0.1* u.megayears)

**The value passed to the checkpoint_interval must have units of time**
1. A list of checkpoint times:

.. code:: python

      Model.run_for(endTime=1.0*u.megayears, 
                    checkpoint_interval=0.1* u.megayears,
                    checkpoint_times=[0.85 * u.megayears, 
                                      0.21 * u.megayears])

**This can be used together or without the checkpoint_interval**

UWGeodynamics will save all the fields defined in the
GEO.rcParams[“default.outputs”] list. You can change that list before
running the model.

Checkpointing
~~~~~~~~~~~~~

By checkpointing we mean saving the data required to restart the Model.
This includes the *mesh*, the *swarm* and all the associated variables.

However, as the swarm and the swarm variables can be very large and can
take a lot of space on disk, the user can decide to save them only every
second, third, fourth etc. checkpoint step.

This is done passing the *restart_checkpoint* argument to the
*Model.run_for* function:

.. code:: python

   Model.run_for(endTime=1.0*u.megayears, 
                 checkpoint_interval=0.1* u.megayears,
             restart_checkpoint=2

By default, the swarm and the swarm variables are saved every time the
model reaches a checkpoint time (``restart_checkpoint=1``).

Restarting the Model
--------------------

When checkpointing a model, the model state is not explicitely saved,
only the mesh, swarms and explicitely saved variables are… We thus need
to recreate the **Model** object before restarting it.

In practice, that means the user must run all commands preceding the
**Model.run_for** command.

When running the **Model.run_for** command *UWGeodynamics* will first
check if an output already exists in the output folder. If it does, the
program will attempt to reload the last available step.

The user can alter this behavior using the **restartStep** and
**restartFolder** arguments:

-  **restartStep** is None by default. The default behaviour is to look
   into **restartFolder** for an existing output and attempt a restart
   from there. Setting it to False will overwrite any existing outputs
   in the *output* folder. If its value is an integer, this corresponds
   to the step number you want to restart from.

-  **restartFolder** is the folder where the program should look for
   previously saved data. It is set to **Model.outputs** by default.

.. code:: python

   import UWGeodynamics as GEO

   u = GEO.u

   Model = GEO.Model(elementRes=(64, 64),
                     minCoord=(0. * u.kilometer, 0. * u.kilometer),
                     maxCoord=(64. * u.kilometer, 64. * u.kilometer))

   # normal model setup

   # Default (restart, restartFolder are optional in this case)
   Model.run_for(2.0 * u.megayears, restartStep=True, restartFolder="your_restart_folder") 

   # Restart from step 10
   Model.run_for(2.0 * u.megayears, restartStep=10, restartFolder="your_restart_folder") 

   # Overwrite existing outputs
   Model.run_for(2.0 * u.megayears, restartStep=False) 
