Welcome to UWGeodynamics's documentation!
=========================================

.. image:: ./img/Tutorial1.gif

The UWGeodynamics module intents to facilitate rapid prototyping of geodynamics
models using Underworld. It can be seen as a set of high-level functions within
the Underworld ecosystem. It is a means to quickly get the user into Underworld
modelling and assumes very little knowledge in coding. The module make some
assumptions based on how the user defines the boundary conditions and the
properties of the materials (rocks, phases). Its simplicity comes with a
relatively more rigid workflow (compared to the classic Underworld functions).
However, the user can easily break the high level objects and get back to core
Underworld function at any step of model design.

The UWGeodynamics is inspired by the Lithospheric Modelling Recipe (LMR)
originally developed by Luke Mondy, Guillaume Duclaux and Patrice Rey
for Underworld 1. Some of the naming conventions have been reused to facilitate
the transition from LMR. The Rheological libraries are also taken from LMR.

As we think the low-level interface is more flexible, and in so allows for more
complex models, we strongly encourage users to explore and break
the High Level functions.

We hope that the user will naturally move to the low-level functionalities as
he or her gets more confident, and by doing so will access the wide range
of possibilities offered by Underworld.

.. image:: ./img/SandboxCompression.gif

.. toctree::
   :maxdepth: 2

   Installation
   User Guide
   Troubleshoot

   Examples
   Benchmarks
   Tutorials
