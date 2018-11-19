User Guide
==========

The Jupyter notebook
--------------------

The Jupyter_ notebook provides a powerfull
environment for the development and analysis of Underworld model.
*Underworld* and *UWGeodynamics* recommend using Jupyter notebooks
for the development of geodynamic models.

If you are not familiar with Jupyter notebooks, we suggest you follow
a quick introduction `here <https://mybinder.org/v2/gh/ipython/ipython-in-depth/master?filepath=binder/Index.ipynb>`_.


Design principles
-----------------

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
   >>> length = 100. * u.kilometre
   >>> width = 50. * u.kilometre
   >>> gravity = 9.81 * u.metre / u.second**2

Pint_ offers the possibility to append a prefix to the units.
1 million year can thus be defined as follow:

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u
   >>> 1.0 * u.megayear

.. note::

   Unit abbreviation is also possible :code:`u.km` is equivalent to :code:`u.kilometer`.
   You can refer to the Pint_ documentation for all abbreviations available.


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
| [length]      | 1.0 metre    |
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

   >>> GEO.scaling_coefficients["[length]"] = 3. * u.kilometre
   >>> GEO.scaling_coefficients["[mass]"] = 4. * u.kilogram
   >>> GEO.scaling_coefficients["[temperature]"] = 273.15 * u.degK
   >>> GEO.scaling_coefficients["[time]"] = 300. * u.years

The unit entered are checked internally and an error is raised if the
units are incompatible. The value is automatically converted to the base
units (metre, second, degree, etc).

To scale a model, the user must define a serie of characteristic
physical values and assign them to the scaling object.

Arguments with units will be scaled by the UWGeodynamics functions.

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u

   >>> KL = 100 * u.kilometre
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

1. define a length of 300 kilometres.
2. use the GEO.nd function to scale it.
3. convert the value back to SI units.

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u

   >>> GEO.scaling_coefficients["[length]"] = 300. * u.kilometre

   >>> length = 300. * u.kilometre
   >>> scaled_length = GEO.nd(length)
   >>> print(scaled_length)
   1.0
   >>> length_metres = GEO.Dimensionalize(scaled_length, u.metre)
   >>> print(length_metres)
   300.0 kilometre


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

   >>> import UWGeodynamics as GEO
   >>> u = GEO.u
   >>> Model = GEO.Model(elementRes=(64, 64),
   ...                   minCoord=(0. * u.kilometre, 0. * u.kilometre),
   ...                   maxCoord=(64. * u.kilometre, 64. * u.kilometre))

The Material object
-------------------

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

They are several ways to define the physical parametres of our Material.

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

.. warning::

   UWGeodynamics contains some basic dimensionality checks. Entering
   wrong units will raise an error

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
   >>> myMaterial.density = 3000 * u.kilogram / u.metre**3
   >>> myMaterial.viscosity = 1e19 * u.pascal * u.second
   >>> myMaterial.radiogenicHeatProd = 0.7 * u.microwatt / u.metre**3
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
-  Disk (2D)
-  Spheres (3D)
-  Annulus (2D)
-  MultiShape (Combination of any of the above) (2D)
-  HalfSpace (3D)

**Layer**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> shape = GEO.shapes.Layer(top=30.*u.kilometre, bottom=0.*u.kilometre)
   >>> material = Model.add_material(name="Material", shape=shape)

   >>> Fig = glucifer.Figure(figsize(1200,400))
   >>> Fig.Points(Model.swarm, Model.materialField)
   >>> Fig.show()

.. image:: /img/layers.png

**Polygon**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> polygon = GEO.shapes.Polygon(vertices=[(10.* u.kilometre, 10.*u.kilometre),
                                              (20.* u.kilometre, 35.*u.kilometre),
                                              (35.* u.kilometre, 5.*u.kilometre)])
   >>> material = Model.add_material(name="Material", shape=polygon)

   >>> Fig = glucifer.Figure(figsize(1200,400))
   >>> Fig.Points(Model.swarm, Model.materialField)
   >>> Fig.show()

.. image:: /img/polygon.png

**Box**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> box = GEO.shapes.Box(top=10.* u.kilometre, bottom=5*u.kilometre,
                            minX=10.*u.kilometre, maxX=15*u.kilometre)
   >>> material = Model.add_material(name="Material", shape=box)

   >>> Fig = glucifer.Figure(figsize(1200,400))
   >>> Fig.Points(Model.swarm, Model.materialField)
   >>> Fig.show()

.. image:: /img/box.png

**Disk**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> disk = GEO.shapes.Disk(center=(32. * u.kilometre, 32. * u.kilometre),
   ...                        radius=10.*u.kilometre)
   >>> material = Model.add_material(name="Material", shape=disk)

   >>> Fig = glucifer.Figure(figsize(1200,400))
   >>> Fig.Points(Model.swarm, Model.materialField)
   >>> Fig.show()

.. image:: /img/disk.png


**Sphere (3D)**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model(elementRes=(16, 16, 16),
                         minCoord=(-1. * u.m, -1. * u.m, -50. * u.cm),
                         maxCoord=(1. * u.m, 1. * u.m, 50. * u.cm))

   >>> sphereShape = GEO.shapes.Sphere(center=(0., 0., 20.*u.centimetre),
                                       radius=20. * u.centimetre))

**Annulus**

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> annulus = GEO.shapes.Annulus(center=(35.*u.kilometre, 50.*u.kilometre),
   ...                              r1=5.*u.kilometre,
   ...                              r2=10.*u.kilometre)
   >>> material = Model.add_material(name="Material", shape=annulus)

   >>> Fig = glucifer.Figure(figsize(400,400))
   >>> Fig.Points(Model.swarm, Model.materialField)
   >>> Fig.show()

.. image:: /img/annulus.png


**MultiShape**

Several shapes can be combined to form a material shape:

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.u
   >>> Model = GEO.Model()
   >>> disk1 = GEO.shapes.Disk(center=(10. * u.kilometre, 10. * u.kilometre),
   ...                         radius=10.*u.kilometre)
   >>> disk2 = GEO.shapes.Disk(center=(20. * u.kilometre, 20. * u.kilometre),
   ...                         radius=5.*u.kilometre)

   >>> shape = GEO.shapes.MultiShape([disk1, disk2])
   >>> material = Model.add_material(name="Material", shape=shape)

   >>> Fig = glucifer.Figure(figsize(400,400))
   >>> Fig.Points(Model.swarm, Model.materialField)
   >>> Fig.show()

.. image:: /img/multishape.png

The following is equivalent:

.. code:: python

  >>> import UWGeodynamics as GEO
  >>> import glucifer

  >>> u = GEO.u
  >>> Model = GEO.Model()
  >>> disk1 = GEO.shapes.Disk(center=(32. * u.kilometre, 32. * u.kilometre),
  ...                         radius=10.*u.kilometre)
  >>> disk2 = GEO.shapes.Disk(center=(32. * u.kilometre, 22. * u.kilometre),
  ...                         radius=10.*u.kilometre)

  >>> shape = disk1 + disk2
  >>> material = Model.add_material(name="Material", shape=shape)

  >>> Fig = glucifer.Figure(figsize(400,400))
  >>> Fig.Points(Model.swarm, Model.materialField)
  >>> Fig.show()


You can also take the intersection of some shapes:

.. code:: python

  >>> import UWGeodynamics as GEO
  >>> import glucifer

  >>> u = GEO.u
  >>> Model = GEO.Model()
  >>> disk1 = GEO.shapes.Disk(center=(32. * u.kilometre, 32. * u.kilometre),
  ...                         radius=10.*u.kilometre)
  >>> disk2 = GEO.shapes.Disk(center=(32. * u.kilometre, 22. * u.kilometre),
  ...                         radius=10.*u.kilometre)

  >>> shape = disk1 & disk2
  >>> material = Model.add_material(name="Material", shape=shape)

  >>> Fig = glucifer.Figure(figsize(400,400))
  >>> Fig.Points(Model.swarm, Model.materialField)
  >>> Fig.show()


**HalfSpace**

HalfSpaces can be used to divide space in 2 domains. The divide is a plan defined
by its normal vector. The convention is to keep the domain opposite to direction
defined by the normal vector.

.. note::

   HalfSpaces can be combined to define 3D shapes / volumes.

.. image:: /img/3D_halfspaces.png

.. code:: python

   >>> import UWGeodynamics as GEO
   >>> import glucifer

   >>> u = GEO.UnitRegistry

   >>> Model = GEO.Model(elementRes=(34, 34, 12),
   ...                   minCoord=(0. * u.km, 0. * u.km, -2880. * u.km),
   ...                   maxCoord=(9000. * u.km, 2000. * u.km, 20. * u.km))

   >>> halfspace1 = GEO.shapes.HalfSpace(normal=(-1.,0.,1.), origin=(4000. * u.km, 0. * u.km, -1000. * u.km))
   >>> halfspace2 = GEO.shapes.HalfSpace(normal=(0.,0.,1.), origin=(7000. * u.km, 1000. * u.km, 0. * u.km))
   >>> halfspace3 = GEO.shapes.HalfSpace(normal=(1.,0.,0.), origin=(9000. * u.km, 1000. * u.km, -500. * u.km))
   >>> halfspace4 = GEO.shapes.HalfSpace(normal=(0.,0.,-1.), origin=(6500. * u.km, 1000. * u.km, -1000. * u.km))

   >>> Fig = glucifer.Figure()
   >>> Fig.Points(Model.swarm, Model.materialField, cullface=False, opacity=1.)
   >>> Fig.Mesh(Model.mesh)
   >>> viewer = Fig.viewer(resolution=(1200,600))
   >>> viewer = Fig.viewer(axis=True)
   >>> viewer.rotatex(-70)
   >>> viewer.rotatey(-10)
   >>> viewer.window()


.. image:: /img/3D_halfspaces2.png

**Multiple materials**

You can add as many materials as needed:

.. code:: python

  >>> import UWGeodynamics as GEO
  >>> import glucifer

  >>> u = GEO.u
  >>> Model = GEO.Model()
  >>> shape = GEO.shapes.Layer(top=30.*u.kilometre, bottom=0.*u.kilometre)
  >>> material1 = Model.add_material(name="Material 1", shape=shape)

  >>> polygon = GEO.shapes.Polygon(vertices=[(10.* u.kilometre, 10.*u.kilometre),
  ...                                        (20.* u.kilometre, 35.*u.kilometre),
  ...                                        (35.* u.kilometre, 5.*u.kilometre)])

  >>> material2 = Model.add_material(name="Material 2", shape=polygon)

  >>> Fig = glucifer.Figure(figsize=(400,400))
  >>> Fig.Points(Model.swarm, Model.materialField, fn_size=3.)
  >>> Fig.show()
  >>> Fig.save("multiple_materials.png")


Rheologies
----------

Newtonian Rheology
~~~~~~~~~~~~~~~~~~

A newtonian rheology can be applied by assigning a viscosity

.. code:: python

  >>> import UWGeodynamics as GEO

  >>> myMaterial = GEO.Material(name="Newtonian Material")
  >>> myMaterial.viscosity = 1e19 * u.pascal * u.second

Non-Newtonian Rheology
~~~~~~~~~~~~~~~~~~~~~~

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
   ...                              stressExponent=1.0,
   ...                              activationVolume=0.,
   ...                              activationEnergy=200 * u.kilojoules,
   ...                              waterFugacity=0.0,
   ...                              grainSize=0.0,
   ...                              meltFraction=0.,
   ...                              grainSizeExponent=0.,
   ...                              waterFugacityExponent=0.,
   ...                              meltFractionFactor=0.0,
   ...                              f=1.0)

Single parametres can then be modified

.. code:: python

   >>> viscosity.activationEnergy = 300. * u.kilojoule

Plastic Behavior (Yield)
~~~~~~~~~~~~~~~~~~~~~~~~

As for Viscous Creep, we provide a registry of commmonly used
plastic behaviors.
They can be accessed using the `GEO.PlasticityRegistry` registry.

.. image:: /img/PlasticityRegistry.gif

The user can define its own parametres:

.. code:: python

   >>> material.plasticity = GEO.DruckerPrager(
   ...     cohesion=10. * u.megapascal,
   ...     cohesionAfterSoftening=10. * u.megapascal,
   ...     frictionCoefficient = 0.3,
   ...     frictionAfterSoftening = 0.2,
   ...     epsilon1=0.5,
   ...     epsilon2=1.5)

   >>> material.plasticity = GEO.VonMises(cohesion=10. * u.megapascal)

Elasticity
~~~~~~~~~~

Elastic behavior can be added to a material:

.. code:: python

   >>> material.elasticity(shear_modulus=10e9 * u.pascal,
                           observation_time=10000 * u.year)


Mechanical Boundary Conditions
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
   ...                   minCoord=(0. * u.kilometre, 0. * u.kilometre),
   ...                   maxCoord=(64. * u.kilometre, 64. * u.kilometre))

Kinematic boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kinematic boundary conditions are set using the **set_velocityBCs** method.
Conditions are defined for each wall (left, right, bottom, top, back and front (3D only)).
For each wall, the user must define the condition for each degree of freedom
(2 in 2D (x,y), 3 in 3D (x,y,z).

if :math:`V` is a vector :math:`(V_x, V_y, V_z)` that we
want to apply on the left wall, the *left* parametre must be defined as
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
   ...                       right=[-1.0*u.centimetre/u.year, None],
   ...                       bottom=[None, 0.],
   ...                       top=[None,0.])

.. image:: /img/mechanicalBCs1.png

3D
^^

Defining boundary conditions for a 3D model is no different than above.
The user must define the velocity components with 3 degree of freedom
instead of 2.

.. code:: python

   >>> Model2 = GEO.Model(elementRes=(16, 16, 16),
   ...                    minCoord=(0. * u.kilometre, 0. * u.kilometre, 0. * u.kilometre),
   ...                    maxCoord=(64. * u.kilometre, 64. * u.kilometre, 64. * u.kilometre))

.. code:: python

   >>> Model2.set_velocityBCs(left=[1.0*u.centimetre/u.year, None, 0.],
   ...                        right=[-1.0*u.centimetre/u.year, None, 0.],
   ...                        bottom=[None, None, 0.],
   ...                        top=[None, None, 0.],
   ...                        front=[None, 0., None],
   ...                        back=[None, 0., None])

Velocity varying along a wall
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is sometime necessary to define a velocity only for a section of a
wall. That can be done using a **condition**. A condition is a set of
rule to apply on a wall.

As an example, we will apply a velocity of :math:`5.0\text{cm/yr}` for
the part of the left wall below 32 kilometre. Velocity is set to be
:math:`1.0\text{cm/yr}` above.

.. code:: python

   >>> conditions = [(Model.y < GEO.nd(32 * u.kilometre), GEO.nd(5.0 * u.centimetre/u.year)),
                     (True, GEO.nd(1.0*u.centimetre/u.year))]

   >>> Model.set_velocityBCs(left=[conditions, None],
   ...                       right=[-1.0*u.centimetre/u.year, None],
   ...                       bottom=[None, 10.*u.megapascal],
   ...                       top=[None,0.])
   >>> Fig = Model.plot.velocityField()

.. image:: /img/mechanicalBCs2.png


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

Frictional Boundaries
~~~~~~~~~~~~~~~~~~~~~

Frictional Boundaries can be set as follow:

.. code:: python

   >>> Model.set_frictional_boundary(left=True,
   ...                               right=True,
   ...                               bottom=True,
   ...                               top=False,
   ...                               friction=19.0,
   ...                               thickness=3)

Where *left*, *right*, *top*, *bottom*, parametres are the side you want
to apply a frictional boundary condition on. *friction* is the angle of
friction (in degrees). *thickness* is the thickness of the boundary.

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
parametre (default to False).

.. code:: python

   >>> Model.set_velocityBCs(left=[1.0*u.centimetre/u.year, None],
   ...                       right=[-1.0*u.centimetre/u.year, None],
   ...                       bottom=[None, GEO.LecodeIsostasy(reference_mat=Model.index)],
   ...                       top=[None,0.])

Traction Condition (stress)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another approach to model isostasy is to defined a stress at the base of
the model. This is done using units of stress (derived SI units =
pascal). The model will then maintain the stress by adjusting the flow
across the border.

.. code:: python

   >>> Model.set_stressBCs(bottom=[None, 10.*u.megapascal])


Thermal Boundary Conditions
---------------------------

Absolute temperatures
~~~~~~~~~~~~~~~~~~~~~

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

   >>> Model.set_temperatureBCs(top=500. * u.degK,
   ...                          bottom=-0.022 * u.milliwatt / u.metre**2,
   ...                          bottom_material=Model,
   ...                          materials=[(air, 273. * u.Kelvin)])

.. Note::

   Model inflow is negative, outflow is positive.


**Fix the temperature of internal nodes**

You can assign a temperature to a list of nodes by passing a list of
node indices (global).

.. code:: python

   >>> nodes = [0, 1, 2]
   >>> Model.set_temperatureBCs(top=500. * u.degK,
   ...                          bottom=-0.022 * u.milliwatt / u.metre**2,
   ...                          bottom_material=Model,
   ...                          nodeSets=[(nodes, 273. * u.Kelvin)])

Heat flux
~~~~~~~~~

.. code:: python

   >>> Model.set_temperatureBCs(top=500. * u.degK, bottom=-0.22 * u.milliwatt / u.metre**2, bottom_material=Model)


Model initialization
--------------------

Initialization of the pressure and temperature field is done using the
::code::python`Model.init_model` method.

The default behavior is to initialise the temperature field to a steady-state
while the pressure field is initialize to the lithostatic pressure.

You can deactivate pressure or temperature initialization by setting the
corresponding argument to `False` (`Model.init_model(temperature=False)`)

.. warning::

   The lithostatic pressure calculation relies on a regular quadratic mesh.
   Most of the time this is fine for model initialization as models often
   starts on a regular mesh. However, this will not work on a deformed mesh

Running the Model
-----------------

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
~~~~~~~~~~~~~~~~~~

UWGeodynamics calculates the time step automatically based on some
numerical stability criteria. You can force a specific time step or
force the time step to be constant throughou

Saving data
~~~~~~~~~~~

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
   ...                  checkpoint_interval=0.1* u.megayears,
   ...                  checkpoint_times=[0.85 * u.megayears,
   ...                                    0.21 * u.megayears])

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
   ...               checkpoint_interval=0.1* u.megayears,
   ...               restart_checkpoint=2

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
   ...                   minCoord=(0. * u.kilometre, 0. * u.kilometre),
   ...                   maxCoord=(64. * u.kilometre, 64. * u.kilometre))

   >>> # Default (restart, restartFolder are optional in this case)
   >>> Model.run_for(2.0 * u.megayears, restartStep=-1, restartFolder="your_restart_folder")

   >>> # Restart from step 10
   >>> Model.run_for(2.0 * u.megayears, restartStep=10, restartFolder="your_restart_folder")

   >>> # Overwrite existing outputs
   >>> Model.run_for(2.0 * u.megayears, restartStep=False)


Parallel run
------------

Model can be run on multiple processors:

You first need to convert your jupyter notebook to a python script:

.. code:: bash

  jupyter nbconvert --to python my_script.ipynb


You can then run the python script as follow:

.. code:: bash

  mpirun -np 4 python my_script.py


.. warning::

   Underworld and UWGeodynamics functions are parallel safe and
   can be run on multiple CPUs. This might not be the case of other
   python library you might be interested in using in your Model.
   For example, matplotlib plots will not work in parallel and must
   be processed in serial.
   *Tutorial 1* has examples of matplotlib plots which are only done
   on the rank 0 CPU.


Passive Tracers
---------------

.. code:: python

   >>> import UWGeodynamics as GEO

   >>> u = GEO.u

   >>> Model = GEO.Model(elementRes=(64,64),
   ...                   minCoord=(0.*u.kilometre, 0.* u.kilometre),
   ...                   maxCoord=(64.* u.kilometre, 64 * u.kilometre))

   >>> x = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), 1000)
   >>> y = 32. * u.kilometre

   >>> P = Model.add_passive_tracers(vertices=[x,y])

.. note::

   You can pass a list of centroids to the `Model.add_passive_tracers` method.
   In that case, the coordinates of the passive tracers are relative to the
   position of the centroids. The pattern is repeated around each centroid.

Surface Processes
-----------------

A range of basic surface processes function are available from the
*surfaceProcesses* submodule. Surface processes are turned on once you
have passed a valid surface processes function to the
``surfaceProcesses`` method of the ``Model`` object.

Example:

.. code:: python

   >>> import UWGeodynamics as GEO

   >>> Model.surfaceProcesses = GEO.surfaceProcesses.SedimentationThreshold(air=[air], sediment=[sediment], threshold=0. * u.metre)

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
   ...     airIndex=[air.index], sedimentIndex=sediment.index,
   ...     XML="ressources/badlands.xml", resolution=1. * u.kilometre,
   ...     checkpoint_interval=0.01 * u.megayears)


Deforming Mesh
--------------

Uniaxial deformation can be turned on using the ``Model.mesh_advector()``
method. The method takes an ``axis`` argument which defines the direction
of deformation (x=0, y=1, z=2)

.. code:: python

   >>> Model.mesh_advector(axis=0)

Element are stretched or compressed uniformly across the model.
This will results in a change in resolution with time.

Top Free surface
----------------

Free surface can be turned on using the ``Model.freesurface`` switch.

.. code:: python

   >>> Model.freesurface = True

.. warning::

   No stabilization algorithm has been implemented yet.


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


.. table::

   ======== =========== ============
   name      function    default val.
   ======== =========== ============


The ``uwgeodynamicsrc`` file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

UWGeodynamics uses ``uwgeodynamicsrc`` configuration files to customize
all kinds of properties, which we call ``rc settings`` or
``rc parametres``. For now, you can control the defaults of a limited
set of property in UWGeodynamics looks for
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
   ``{INSTALL}`` is something like ``/usr/lib/python2.7/site-packages``
   on Linux, and maybe ``C:\\Python27\\Lib\\site-packages`` on Windows.
   Every time you install UWgeodynamics, this file will be overwritten, so
   if you want your customizations to be saved, please move this file to
   your user-specific directory.

To display where the currently active ``uwgeodynamicsrc`` file was
loaded from, one can do the following:

.. code:: python

     >>> import UWGeodynamics as GEO
     >>> GEO.uwgeodynamics_fname()
     '/home/foo/.config/uwgeodynamics/uwgeodynamicsrc'

See below for a sample.

\_uwgeodynamicsrc-sample:
~~~~~~~~~~~~~~~~~~~~~~~~~

.. _Jupyter: http://jupyter.org/
.. _Docker Hub: https://hub.docker.com/r/underworldcode/uwgeodynamics
.. _Kitematic: https://kitematic.com/
.. _github: https://github.com/underworldcode/UWGeodynamics.git
.. _Pint: https://pint.readthedocs.io/en/latest
