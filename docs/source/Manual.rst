UWGeodynamics User Guide
========================

Docker_
-------

Docker containers provide and easy-way to set up and distribute
applications. They also provide a safe and consistent environment which
facilitate debugging and reproducibility of models. The image we provide
contains all the dependencies and configuration files required to run
Underworld models. Users can start developping model as soon as they
have downloaded the image, independently of the operating system running
on their machine.

We strongly encourage users to run UWGeodynamics using the docker images
we provide on `Docker Hub`_

Different version of the `underworldcode/uwgeodynamics` image can be
pulled using a tag:

1. The *latest* tag points to the github master branch and uses the latest
   *underworld* release.
2. The *dev* tag points to the github development and uses the development
   branch of *underworld*.
3. release tags such as *0.9.8* points to the specified version.

Command line
~~~~~~~~~~~~

Once you have installed docker on your system you can *pull* the
*UWGeodynamics* official image as follow:

.. code:: bash

  docker pull underworldcode/uwgeodynamics

You can then start a docker container. (An instance of
an image).

.. code:: bash

  docker run -d \
     --name my_container \
     --port 8888:8888 \
     --mount source=myvol,target=/workspace/user-data \
     underworldcode/uwgeodynamics

You can access the container via your browser at the following
address: http://localhost:8888


It is also possible to ssh into the container as follow:

.. code:: bash

  docker exec -it my_container /bin/bash


Kitematic_
~~~~~~~~~~

Kitematic_ is a program that provides a graphical user interface to
the *docker* daemon and to Docker Hub.
It is available on Linux, MacOSX and Windows.

1. Download and Install Kitematic_
   The software is available for Windows, MacOsx and Linux. Be aware that on
   linux the installation may differ depending on the distribution you
   are running.

2. Open Kitematic and search for the **uwgeodynamics** image.
3. Create a container by clicking on the create button.

You should now have a container appearing on the left side of your
kitematic window. The first thing to do now is to create a link between
a local directory (A directory on your physical hard drive) and a volume
directory inside the docker container. A volume is a special directory
that can be accessed from outside the container. It is the location you
will use to save your results.

Local Installation
------------------

If you really need to install the software natively on your system...
This is not recommended and involves installing *Underworld* and all
its dependencies. Docker is highly recommended!!!

Requirements
~~~~~~~~~~~~

-  Python >= 2.7
-  A Working version of Underworld2 >=2.6.0 (Please refer to the
   Underworld documentation)
-  pint >= 0.8

**Note on Python 3 compatibility**:
The bleeding edge version of *Underworld* (development branch)
is now python 3 compatible only.
*UWGeodynamics* is python 3 ready and can thus be used with it.

Install
~~~~~~~

**from Pip**

The UWGeodynamics module can be installed directly from the Python
Package Index:

.. code:: bash

  pip3 install UWGeodynamics

**from sources**

The module source files are available through github_

.. code:: bash

  git clone https://github.com/underworldcode/UWGeodynamics.git

It can then be installed globally on your system using

.. code:: bash

  pip3 install -e UWGeodynamics

The Jupyter notebook
--------------------

The Jupyter_ notebook provides a powerfull
environment for the development and analysis of Underworld model.
*Underworld* and *UWGeodynamics* recommend using Jupyter notebooks
for the development of geodynamic models.

If you are not familiar with Jupyter notebooks, we suggest you follow
a quick introduction `here <https://mybinder.org/v2/gh/ipython/ipython-in-depth/master?filepath=binder/Index.ipynb>`_.


import UWGeodynamics
--------------------

*UWGeodynamics* can be imported as follow:

.. code:: python

   >>> import UWGeodynamics as GEO

Working with units
------------------

*UWGeodynamics* uses Pint_, a
Python package to define, operate and manipulate **physical quantities**
(A numerical value with unit of measurement). Pint is a very powerful
package that handles conversion and operation between units.

We recommend using SI units but other systems are also available.

Pint_ **Unit Registry** can be used as follow:

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.UnitRegistry

or simply

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u

You can have a quick overview of all the units available by hitting tab
after the “.” of the u object.

.. image:: img/tabtab.gif

Quantities can then be defined as follow:

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u
   >>> length = 100. * u.kilometer
   >>> width = 50. * u.kilometer
   >>> gravity = 9.81 * u.meter / u.second**2

Pint_ offers the possibility to append a prefix to the units.
1 million year can thus be defined as follow:

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u
   >>> 1.0 * u.megayear

Model Scaling
-------------

Model can be scaled using a series of scaling coefficients

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> GEO.scaling_coefficients

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

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u

   >>> GEO.scaling_coefficients["[length]"] = 3.0 * u.kilometer
   >>> GEO.scaling_coefficients["[mass]"] = 4.0 * u.kilogram
   >>> GEO.scaling_coefficients["[temperature]"] = 273.15 * u.degK
   >>> GEO.scaling_coefficients["[time]"] = 300 * u.years

The unit entered are checked internally and an error is raised if the
units are incompatible. The value is automatically converted to the base
units (meter, second, degree, etc).

To scale a model, the user must define a serie of characteristic
physical values and assign them to the scaling object.

Arguments with units will be scaled by the UWGeodynamics functions.

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u

   >>> KL = 100 * u.kilometer
   >>> Kt = 1. * u.year
   >>> KM = 3000. * u.kilogram
   >>> KT = 1200. * u.degK

   >>> GEO.scaling_coefficients["[length]"] = KL
   >>> GEO.scaling_coefficients["[time]"] = Kt
   >>> GEO.scaling_coefficients["[mass]"]= KM
   >>> GEO.scaling_coefficients["[temperature]"] = KT

Dimensionalize / non-Dimensionalize
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We provide 2 functions :code:`GEO.nonDimensionalize` and :code:`GEO.Dimensionalize`
to convert between non-dimensional and dimensional values.
The function are also available respectively as :code:`GEO.nd` and
:code:`GEO.dim`.

**Example:**

1. define a length of 300 kilometers.
2. use the GEO.nd function to scale it.
3. convert the value back to SI units.

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u

   >>> GEO.scaling_coefficients["[length]"] = 300. * u.kilometer

   >>> length = 300. * u.kilometers
   >>> scaled_length = GEO.nd(length)
   >>> print(scaled_length)
   1.0
   >>> length_meters = GEO.Dimensionalize(scaled_length, u.meters)
   >>> print(length_meters)
   300.0 kilometer


Building a Model
----------------

Design principles
~~~~~~~~~~~~~~~~~



The Model object
~~~~~~~~~~~~~~~~

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

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u
   >>> Model = GEO.Model(elementRes=(64, 64),
                         minCoord=(0. * u.kilometer, 0. * u.kilometer),
                         maxCoord=(64. * u.kilometer, 64. * u.kilometer))

The Material object
~~~~~~~~~~~~~~~~~~~

The *UWGeodynamics* module is designed around the idea of materials,
which are essentially a way to define physical properties across the
Model domain.

Materials are defined using the **Material** object as follow:

.. code:: python

   >>> import UWGeodynamics as GEO

   >>> crust = GEO.Material(name="Crust")

Typing the name of the material in an empty cell will return a table
which summarizes the property of the material:

.. image:: img/Material1.png

As you can see, most of the property are undefined.

They are several ways to define the physical parameters of our Material.

-  The first one is to add them directly when creating the object
   itself:

.. code:: python

   >>> import UWGeodynamics as GEO

   >>> u = GEO.u
   >>> crust = GEO.Material(name="Crust", density=3000*u.kilogram/u.metre**3)

-  The second option is to change the property after creating the
   **Material**:

.. code:: python

   >>> import UWGeodynamics as GEO

   >>> u = GEO.u
   >>> crust = GEO.Material(name="Crust")
   >>> crust.density = 3000. * u.kilogram / u.metre **3

The second option is often easier to read.

**UWGeodynamics contains some basic dimensionality checks. Entering
wrong units will raise an error**

Material can be added to a model as follow:

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> crust = Model.add_material(name="Crust")

Although optional, tt is a good idea to give a **name** to the material.
The **Model.add_material** method will return a Material object. That
object is a python object that will then be used to define the property
of the material.

Material Attributes
~~~~~~~~~~~~~~~~~~~

The Material object comes with a series of attribute that can
be used to define its physical behavior.

.. table:: Materials attributes
  :widths: auto

  =================== ==================
  Name                    Description
  =================== ==================
  shape
  density             Density
  diffusivity         Thermal Diffusivity
  capacity            Thermal Capacity
  radiogenicHeatProd  Radiogenic Heat Production
  viscosity           Viscous behavior
  plasticity          Plastic behavior
  elasticity          Elastic behavior
  minViscosity        Minimum Viscosity allowed
  maxViscosity        Maximum Viscosity allowed
  stressLimiter       Maximum sustainable stress
  healingRate         Plastic Strain Healing Rate
  solidus             Solidus
  liquidus            Liquidus
  latentHeatFusion    Latent Heat Fusion (Enthalpy of Fusion)
  meltExpansion       Melt Expansion
  meltFraction        Initial Melt Fraction
  meltFractionLimit   Maximum Fraction of Melt
  viscosityChange     Change in Viscosity over Melt Fraction range
  viscosityChangeX1   Melt Fraction Range begin
  viscosityChangeX2   Melt Fraction Range end
  =================== ==================

**Examples**

.. code:: python

   >>> Model.density = 200. * u.kg / u.m**3
   >>> myMaterial = GEO.Material(name="My Material")
   >>> myMaterial.density = 3000 * u.kilogram / u.meter**3
   >>> myMaterial.viscosity = 1e19 * u.pascal * u.second
   >>> myMaterial.radiogenicHeatProd = 0.7 * u.microwatt / u.meter**3
   >>> myMaterial.diffusivity = 1.0e-6 * u.metre**2 / u.second

Global properties
^^^^^^^^^^^^^^^^^

The user can define attributes on the *Model* itself.
The values will be used as global values for materials with undefined
attributes

**Example**

.. code:: python

   >>> Model.density = 200. * u.kg / u.m**3
   >>> myMaterial = GEO.Material(name="My Material")

The density of myMaterial will default to 200. kilogram / cubic metre unless
its *density* attribute is specified.


Material shape
^^^^^^^^^^^^^^

The *shape* attribute essentially describe the initial
location of a material.
It is used to build a initial geometry of the model.

There is a range of shapes available

-  Layer (2D/3D)
-  Polygon (2D)
-  Box (2D)
-  Disk (2D) / Spheres (3D)
-  Annulus (2D)
-  MultiShape (Combination of any of the above) (2D)
-  HalfSpace (3D)

**Layer**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> shape = GEO.shapes.Layer(top=30.*u.kilometer, bottom=0.*u.kilometer)
   >>> material = Model.add_material(name="Material", shape=shape)

   >>> Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   >>> Fig.save("layers.png")

.. image:: /img/layers.png

**Polygon**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> polygon = GEO.shapes.Polygon(vertices=[(10.* u.kilometer, 10.*u.kilometer),
                                              (20.* u.kilometer, 35.*u.kilometer),
                                              (35.* u.kilometer, 5.*u.kilometer)])
   >>> material = Model.add_material(name="Material", shape=polygon)

   >>> Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   >>> Fig.save("polygon.png")

.. image:: /img/polygon.png

**Box**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> box = GEO.shapes.Box(top=10.* u.kilometer, bottom=5*u.kilometer,
                            minX=10.*u.kilometer, maxX=15*u.kilometer)
   >>> material = Model.add_material(name="Material", shape=box)

   >>> Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   >>> Fig.save("box.png")

.. image:: /img/box.png

**Disk**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> disk = GEO.shapes.Disk(center=(32. * u.kilometer, 32. * u.kilometer), radius=10.*u.kilometer)
   >>> material = Model.add_material(name="Material", shape=disk)

   >>> Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   >>> Fig.save("disk.png")

.. image:: /img/disk.png

**Annulus**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> annulus = GEO.shapes.Annulus(center=(35.*u.kilometer, 50.*u.kilometer),
                                    r1=5.*u.kilometer,
                                    r2=10.*u.kilometer)
   >>> material = Model.add_material(name="Material", shape=annulus)

   >>> Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   >>> Fig.save("annulus.png")

.. image:: /img/annulus.png

**MultiShape**

Several shapes can be combined to form a material shape:

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> disk1 = GEO.shapes.Disk(center=(10. * u.kilometer, 10. * u.kilometer), radius=10.*u.kilometer)
   >>> disk2 = GEO.shapes.Disk(center=(20. * u.kilometer, 20. * u.kilometer), radius=5.*u.kilometer)

   >>> shape = GEO.shapes.MultiShape([disk1, disk2])
   >>> material = Model.add_material(name="Material", shape=shape)

   >>>Fig = Model.plot.material(figsize=(400, 400), fn_size=3.0)
   >>>Fig.save("multishape.png")

.. image:: /img/multishape.png

the following is equivalent:

.. code:: python

  >>> import UWGeodynamics as GEO
  >>> import glucifer

  >>> u = GEO.u
  >>> Model = GEO.Model()
  >>> disk1 = GEO.shapes.Disk(center=(32. * u.kilometer, 32. * u.kilometer), radius=10.*u.kilometer)
  >>> disk2 = GEO.shapes.Disk(center=(32. * u.kilometer, 22. * u.kilometer), radius=10.*u.kilometer)

  >>> shape = disk1 + disk2
  >>> material = Model.add_material(name="Material", shape=shape)

  >>> Fig = glucifer.Figure(figsize=(400,400))
  >>> Fig.Points(Model.swarm, Model.materialField, fn_size=3.)
  >>> Fig.show()
  >>> Fig.save("multishape-2.png")


You can also take the intersection of some shapes:

.. code:: python

  >>> import UWGeodynamics as GEO
  >>> import glucifer

  >>> u = GEO.u
  >>> Model = GEO.Model()
  >>> disk1 = GEO.shapes.Disk(center=(32. * u.kilometer, 32. * u.kilometer), radius=10.*u.kilometer)
  >>> disk2 = GEO.shapes.Disk(center=(32. * u.kilometer, 22. * u.kilometer), radius=10.*u.kilometer)

  >>> shape = disk1 & disk2
  >>> material = Model.add_material(name="Material", shape=shape)

  >>> Fig = glucifer.Figure(figsize=(400,400))
  >>> Fig.Points(Model.swarm, Model.materialField, fn_size=3.)
  >>> Fig.show()
  >>> Fig.save("multishape-3.png")

**Multiple materials**

You can add as many materials as needed:

.. code:: python

  >>> import UWGeodynamics as GEO
  >>> import glucifer

  >>> u = GEO.u
  >>> Model = GEO.Model()
  >>> shape = GEO.shapes.Layer(top=30.*u.kilometer, bottom=0.*u.kilometer)
  >>> material1 = Model.add_material(name="Material 1", shape=shape)

  >>> polygon = GEO.shapes.Polygon(vertices=[(10.* u.kilometer, 10.*u.kilometer),
  >>>                                        (20.* u.kilometer, 35.*u.kilometer),
  >>>                                        (35.* u.kilometer, 5.*u.kilometer)])

  >>> material2 = Model.add_material(name="Material 2", shape=polygon)

  >>> Fig = glucifer.Figure(figsize=(400,400))
  >>> Fig.Points(Model.swarm, Model.materialField, fn_size=3.)
  >>> Fig.show()
  >>> Fig.save("multiple_materials.png")


Rheologies
~~~~~~~~~~

Newtonian Rheology
^^^^^^^^^^^^^^^^^^

A newtonian rheology can be applied by assigning a viscosity

.. code:: python

  >>> import UWGeodynamics as GEO

  >>> myMaterial = GEO.Material(name="Newtonian Material")
  >>> myMaterial.viscosity = 1e19 * u.pascal * u.second

Non-Newtonian Rheology
^^^^^^^^^^^^^^^^^^^^^^

*UWGeodynamics* provides a library of commonly used Viscous Creep Flow Laws.
They can be accessed using the `GEO.ViscousCreepRegistry` registry:

.. image:: /img/ViscousCreepRegistry.gif


**Example:**

.. code:: python

  >>> import UWGeodynamics as GEO
  >>> material = GEO.Material(name="Material")

  >>> rh = GEO.ViscousCreepRegistry()
  >>> material.viscosity = rh.Gleason_and_Tullis_1995

You can scale viscosity by using a multiplier.
For example to make the **Gleason and Tullis, 1995** rheology
30 times stronger you can do:

.. code:: python

  >>> import UWGeodynamics as GEO
  >>> material = GEO.Material(name="Material")

  >>> rh = GEO.ViscousCreepRegistry()
  >>> material.viscosity = 30 * rh.Gleason_and_Tullis_1995

The user can of course define its own rheology.

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

Plastic Behavior (Yield)
~~~~~~~~~~~~~~~~~~~~~~~~

As for Viscous Creep, we provide a registry of commmonly used
plastic behaviors.
They can be accessed using the `GEO.PlasticityRegistry` registry.

.. image:: /img/PlasticityRegistry.gif

The user can define its own parameters:

.. code:: python

   >>> plasticity = GEO.DruckerPrager(cohesion=10. * u.megapascal,
                                      cohesionAfterSoftening=10. * u.megapascal,
                                      frictionCoefficient = 0.3,
                                      frictionAfterSoftening = 0.2,
                                      epsilon1=0.5,
                                      epsilon2=1.5)

   >>> plasticity = GEO.VonMises(cohesion=10. * u.megapascal)

Mechanical Boundary ConditionVs
-------------------------------

Mechanical boundary conditions are a critical part of any
geodynamic model design. In the following, we quickly detail the options
available to define mechanical boundary conditions in Underworld using the
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

Kinematic boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kinematic boundary conditions are set using the **set_velocityBCs** method.
Conditions are defined for each wall (left, right, bottom, top, back and front (3D only)).
For each wall, the user must define the condition for each degree of freedom
(2 in 2D (x,y), 3 in 3D (x,y,z).

if :math:`V` is a vector :math:`(V_x, V_y, V_z)` that we
want to apply on the left wall, the *left* parameter must be defined as
:code:`left=[Vx, Vy, Vz]`.

In the following example we set the boundary condition to be:

-  left wall: :math:`V_x = -1.0 \text{cm / yr}`,
   :math:`Vy = None`
-  right wall: :math:`V_x = 1.0 \text{cm / yr}`, :math:`Vy=None`
-  bottom wall: :math:`V_x = None`, :math:`V_y= 0.` (free slip)

It is an extension model with a total rate of extension equal to 2.0
centimetre / year. No :math:`V_x` is prescribed at the bottom, while
:math:`V_y` is set to :math:`0.` no material will be able to enter or
leave the model domain from that side. The material is free to move
vertically along the side walls.

.. code:: python

   >>> Model.set_velocityBCs(left=[1.0*u.centimetre/u.year, None],
   >>>                       right=[-1.0*u.centimetre/u.year, None],
   >>>                       bottom=[None, 0.],
   >>>                       top=[None,0.])

.. image:: /img/mechanicalBCs1.png

3D
^^

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

Velocity varying along a wall
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

.. image:: /img/mechanicalBCs2.png

Assign Viscosity to Internal nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Stress Conditions
~~~~~~~~~~~~~~~~~

Stress conditions can be applied to the boundaries using the
**set_stressBCs** method:

In the following example we apply a stress of 200.0 megapascal to the
bottom of our model:

.. code:: python

   >>> Model.set_stressBCs(bottom=[None, 200. * u.megapascal])

Note that you will have to make sure that kinematic and stress conditions
are compatible.


Isostasy
~~~~~~~~

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

Lecode Isostasy (kinematic)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Traction Condition (stress)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another approach to model isostasy is to defined a stress at the base of
the model. This is done using units of stress (derived SI units =
pascal). The model will then maintain the stress by adjusting the flow
across the border.

.. code:: python

   >>> Model.set_stressBCs(bottom=[None, 10.*u.megapascal])


Thermal Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Absolute temperatures
^^^^^^^^^^^^^^^^^^^^^

Setting the temperature at the top of a model to be
:math:`500 \text{kelvin}` at the top and :math:`1600 \text{kelvin}` at
the bottom is done as follow.

.. code:: python

   >>> Model.set_temperatureBCs(top=500. * u.degK, bottom=1600. * u.degK)

You can of course define temperatures on the sidewalls:

.. code:: python

   >>> Model.set_temperatureBCs(right=500. * u.degK, left=1600. * u.degK)

**Fix the temperature of a Material**

.. code:: python

   >>> Model.set_temperatureBCs(top=500. * u.degK, bottom=-0.22 * u.milliwatt / u.metre**2, bottom_material=Model,
                                materials=[(air, 273. * u.Kelvin)])

**Fix the temperature of internal nodes**

You can assign a temperature to a list of nodes by passing a list of
node indices (global).

.. code:: python

   >>> nodes = [0, 1, 2]
   >>> Model.set_temperatureBCs(top=500. * u.degK, bottom=-0.22 * u.milliwatt / u.metre**2, bottom_material=Model,
                                nodeSets=[(nodes, 273. * u.Kelvin)])

Heat flux
^^^^^^^^^

.. code:: python

   >>> Model.set_temperatureBCs(top=500. * u.degK, bottom=-0.22 * u.milliwatt / u.metre**2, bottom_material=Model)


Frictional Boundaries
~~~~~~~~~~~~~~~~~~~~~

Frictional Boundaries can be set as follow:

.. code:: python

   >>> Model.set_frictional_boundary(left=True,
   >>>                               right=True,
   >>>                               bottom=True,
   >>>                               top=False,
   >>>                               friction=19.0,
   >>>                               thickness=3)

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

   >>> import UWGeodynamics as GEO

   >>> Model.surfaceProcesses = GEO.surfaceProcesses.SedimentationThreshold(air=[air], sediment=[sediment], threshold=0. * u.meter)

Three simple function are available:

1. Total Erosion Above Threshold (``ErosionThreshold``).
2. Total Sedimentation Below Threshold (``SedimentationThreshold``)
3. Combination of the 2 above. (``ErosionAndSedimentationThreshold``)

Coupling with Badlands
~~~~~~~~~~~~~~~~~~~~~~

UWGeodynamics provide a way to couple an Underworld model to Badlands.
**More documentation needed**

.. code:: python

   >>> import UWGeodynamics as GEO

   >>> Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(
   >>>     airIndex=[air.index], sedimentIndex=sediment.index,
   >>>     XML="ressources/badlands.xml", resolution=1. * u.kilometer,
   >>>     checkpoint_interval=0.01 * u.megayears)

Passive Tracers and Grids
-------------------------

Passive tracers
~~~~~~~~~~~~~~~

.. code:: python

   >>> import UWGeodynamics as GEO

   >>> u = GEO.u

   >>> Model = GEO.Model(elementRes=(64,64),
   >>>                   minCoord=(0.*u.kilometer, 0.* u.kilometer),
   >>>                   maxCoord=(64.* u.kilometer, 64 * u.kilometer))

   >>> x = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), 1000)
   >>> y = 32. * u.kilometer

   >>> P = Model.add_passive_tracers(vertices=[x,y])

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

   >>> Model.run_for(nstep=10)

1. Specify an endTime

.. code:: python

   >>> Model.run_for(endTime=1.0* u.megayears)

which is equivalent to

.. code:: python

   >>> Model.run_for(1.0*u.megayears)

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

   >>> Model.run_for(endTime=1.0*u.megayears,
                     checkpoint_interval=0.1* u.megayears)

**The value passed to the checkpoint_interval must have units of time**
1. A list of checkpoint times:

.. code:: python

   >>> Model.run_for(endTime=1.0*u.megayears,
   >>>                  checkpoint_interval=0.1* u.megayears,
   >>>                  checkpoint_times=[0.85 * u.megayears,
   >>>                                    0.21 * u.megayears])

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

   >>> Model.run_for(endTime=1.0*u.megayears,
   >>>               checkpoint_interval=0.1* u.megayears,
   >>>               restart_checkpoint=2

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

-  **restartStep** is *-1* by default. The default behaviour is to look
   into **restartFolder** for an existing output and attempt a restart
   from the last output available.
   Setting it to False will overwrite any existing outputs
   in the *output* folder. If its value is an integer, this corresponds
   to the step number you want to restart from.

-  **restartFolder** is the folder where the program should look for
   previously saved data. It is set to **Model.outputs** by default.

.. code:: python

   >>> import UWGeodynamics as GEO

   >>> u = GEO.u

   >>> Model = GEO.Model(elementRes=(64, 64),
   >>>                   minCoord=(0. * u.kilometer, 0. * u.kilometer),
   >>>                   maxCoord=(64. * u.kilometer, 64. * u.kilometer))

   >>> # Default (restart, restartFolder are optional in this case)
   >>> Model.run_for(2.0 * u.megayears, restartStep=-1, restartFolder="your_restart_folder")

   >>> # Restart from step 10
   >>> Model.run_for(2.0 * u.megayears, restartStep=10, restartFolder="your_restart_folder")

   >>> # Overwrite existing outputs
   >>> Model.run_for(2.0 * u.megayears, restartStep=False)


Dynamic rc settings
-------------------

You can dynamically change the default rc settings in a python script or
interactively from the python shell. All of the rc settings are stored
in a dictionary-like variable called `UWGeodynamics.rcParams`, which
is global to the UWGeodynamics package. rcParams can be modified
directly, for example:

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> GEO.rcParams['solver'] = "mumps"
   >>> GEO.rcParams['penalty'] = 1e6


The ``UWGeodynamics.rcdefaults`` command will restore the standard
UWGeodynamics default settings.

There is some degree of validation when setting the values of rcParams,
see ``UWGeodynamics.rcsetup`` for details.

The ``uwgeodynamicsrc`` file
----------------------------

UWGeodynamics uses ``uwgeodynamicsrc`` configuration files to customize
all kinds of properties, which we call ``rc settings`` or
``rc parameters``. For now, you can control the defaults of a limited
set of property in matplotlib UWGeodynamics looks for
``uwgeodynamicsrc`` in four locations, in the following order:

1. ``uwgeodynamicsrc`` in the current working directory, usually used
   for specific customizations that you do not want to apply elsewhere.

2. ``$UWGEODYNAMICSRC`` if it is a file, else
   ``$UWGEODYNAMICSRC/uwgeodynamicsrc``.

3. It next looks in a user-specific place, depending on your platform:

   -  On Linux, it looks in ``.config/uwgeodynamics/uwgeodynamicsrc``
      (or ``$XDG_CONFIG_HOME/uwgeodynamics/uwgeodynamicsrc``) if you’ve
      customized your environment.

   -  On other platforms, it looks in
      ``.uwgeodynamics/uwgeodynamicsrc``.

4. ``{INSTALL}/UWGeodynamics/uwgeo-data/uwgeodynamicsrc``, where
   ``{INSTALL}`` is something like ``/usr/lib/python2.5/site-packages``
   on Linux, and maybe ``C:\\Python25\\Lib\\site-packages`` on Windows.
   Every time you install matplotlib, this file will be overwritten, so
   if you want your customizations to be saved, please move this file to
   your user-specific matplotlib directory.

To display where the currently active ``uwgeodynamicsrc`` file was
loaded from, one can do the following:

.. code:: python

     >>> import UWGeodynamics as GEO
     >>> GEO.uwgeodynamics_fname()
     '/home/foo/.config/uwgeodynamics/uwgeodynamicsrc'

See below for a sample.

\_uwgeodynamicsrc-sample:
-------------------------

.. _Jupyter: http://jupyter.org/
.. _Docker Hub: https://hub.docker.com/r/underworldcode/uwgeodynamics
.. _Kitematic: https://kitematic.com/
.. _github: https://github.com/underworldcode/UWGeodynamics.git
.. _Pint: https://pint.readthedocs.io/en/latest
