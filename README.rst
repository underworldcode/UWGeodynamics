
.. image:: /docs/source/img/logos.png
    :target: https://www.earthbyte.org/the-basin-genesis-hub

Underworld Geodynamics project
==============================

.. image:: http://joss.theoj.org/papers/10.21105/joss.01136/status.svg
   :target: https://doi.org/10.21105/joss.01136

.. image:: https://zenodo.org/badge/114189389.svg
    :target: https://zenodo.org/badge/latestdoi/114189389
    :alt: DOI      

.. image:: https://api.codacy.com/project/badge/Grade/85b5f7736d03441db786549d6e357c9e
    :target: https://www.codacy.com/app/romainbeucher/UWGeodynamics?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=rbeucher/UWGeodynamics&amp;utm_campaign=Badge_Grade
    :alt: Codacy

.. image:: https://img.shields.io/pypi/v/uwgeodynamics.svg
    :target: https://pypi.python.org/pypi/uwgeodynamics
    :alt: Pip

.. image:: https://readthedocs.org/projects/uwgeodynamics/badge
    :target: http://uwgeodynamics.readthedocs.org/
    :alt: Docs

.. image:: https://travis-ci.org/underworldcode/UWGeodynamics.svg?branch=development
    :target: https://travis-ci.org/underworldcode/UWGeodynamics

.. image:: https://github.com/underworldcode/UWGeodynamics/blob/development/tutorials/images/Tutorial1.gif

.. image:: https://github.com/underworldcode/UWGeodynamics/blob/development/docs/source/img/collision_wedge.gif

The UWGeodynamics module facilitates prototyping of geodynamics models using Underworld. 
It can be seen as a set of high-level functions within the Underworld ecosystem.
It is a means to quickly get the user into Underworld modelling and assumes very
little knowledge in coding. The module make some assumptions based on how the user
defines the boundary conditions and the properties of the materials (rocks, phases).
Its simplicity comes with a relatively more rigid workflow (compared to the classic Underworld functions).
However, the user can easily break the high level objects and get back to core
Underworld function at any step of model design.

The UWGeodynamics is inspired by the [Lithospheric Modelling Recipe (LMR)](https://github.com/LukeMondy/lithospheric_modelling_recipe) originally developed by
Luke Mondy, Guillaume Duclaux and Patrice Rey for Underworld 1. 
Some of the naming conventions have been reused to facilitate the transition from LMR.
The Rheological libraries is also taken from LMR.

As we think the low-level interface is more flexible, and in so allows for more complex models,
we strongly encourage users to explore and break the High Level functions.

We hope that the user will naturally move to the low-level functionalities as he
or her gets more confident, and by doing so will access the wide range of 
possibilities offered by Underworld.

.. image:: /docs/source/img/SandboxCompression.gif

UWGeodynamics and Underworld
----------------------------

*UWGeodynamics* uses the Underworld_  Application Programming Interface (API).
Both projects are supported by The Underworld development team led by Louis Moresi and based in Melbourne, Australia
at the University of Melbourne and at Monash University.

*Underworld* and *UWGeodynamics* both provide powerful tools to develop numerical geodynamic models.
But their approaches are different: *UWGeodynamics* largely guides users into a way of doing things.
The Underworld API provides a series of tools and components (Mesh, Mesh variables, system of equations, functions)
and leaves the responsibility to arrange those components to the user. The main advantage of the Underworld API is its flexibility.
The main inconvenient resides in a somewhat steeper learning curve. *UWGeodynamics* components are
designed to be more natural to non-experimented numerical modellers or people with little knowledge in programming.
It is a way to quickly get started and design numerical models. Developing complex models can also be facilitated
by the *UWGeodynamics* high-level interface as it requires less time and less involvement
with the details of the Underworld API.

The two approaches are complementary and mixing the two approaches is possible and highly encouraged.

Note on versioning
------------------

Since version 1.0 The Underworld development team has decided to match the *UWGeodynamics* version number with
the latest supported version of Underworld. 
UWGeodynamics v2.7 is then supporing Underworld up to version 2.7.

The third number is used for *UWGeodynamics* only (v2.7.1, v2.7.2 etc.)

The development branch is based on the current *Underworld* development branch.

The Current release (**DOI** citable): 

`DOI <https://zenodo.org/badge/114189389.svg)](https://zenodo.org/badge/latestdoi/114189389>`_

Quick Start / Testing
----------------------

We provide a docker container via binder_.
This is a quick solution to get you started and run the examples and tutorials
without installing anything on your machine. That is a good way to see if the
software can actually be useful to you. 
The ressource are however limited and you should not try to run model with high resolution.
3D models can not be run in the binder.

Where to find documentation?
----------------------------

The full documentation is available on `ReadTheDocs <http://uwgeodynamics.readthedocs.org/>`_

Additional documentation and function specific documentation can be find in the python doctrings.
You can acces them in the Jupyter_ notebook by prepending or appending the method, variable or function with ``?``.

Installation
-------------

Docker_ installation
~~~~~~~~~~~~~~~~~~~~

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
3. release tags such as *v2.7.1* points to a specific version.

**Command line**

Once you have installed docker on your system you can *pull* the
*UWGeodynamics* official image as follow:

.. code:: bash

  docker pull underworldcode/uwgeodynamics

You can list all the images available on your system as follow:

.. code:: bash

  docker images

An image can be deleted as follow:

.. code:: bash

  docker rmi underworldcode/uwgeodynamics

You can then start a docker container. (An instance of
an image).

.. code:: bash

  docker run -d \
     --name my_container \
     -p 8888:8888 \
     --mount source=myvol,target=/workspace/user-data \
     underworldcode/uwgeodynamics

You can access the container via your browser at the following
address: http://localhost:8888

It is also possible to ssh into the container as follow:

.. code:: bash

  docker exec -it my_container /bin/bash

You can list the containers currently existing on your machine by running:

.. code:: bash

  docker ps -a

The "a" means "all container". The :code:`docker ps` command only list
running containers.

Docker containers can be stop (so that they do not use CPU or RAM ressource):

.. code:: bash

  docker stop my_container

They can also be deleted:

.. code:: bash

  docker rm my_container

.. warning::

  It's a good idea to keep track of how many containers have been created as
  they can rapidly take a lot of space on your machine.

Kitematic_
~~~~~~~~~~

Kitematic_ is a program that provides a graphical user interface to
the *docker* daemon and to Docker Hub.
The software is available for Windows, MacOsx and Linux. Be aware that on
linux the installation may differ depending on the distribution you
are running.

1. Download and Install Kitematic_
2. Open Kitematic and search for the **uwgeodynamics** image.
3. Create a container by clicking on the create button.

You should now have a container appearing on the left side of your
kitematic window. The first thing to do now is to create a link between
a local directory (A directory on your physical hard drive) and a volume
directory inside the docker container. A volume is a special directory
that can be accessed from outside the container. It is the location you
will use to save your results.

Local Installation
~~~~~~~~~~~~~~~~~~~~

This is not recommended and involves installing *Underworld* and all
its dependencies. Docker is highly recommended!!!

**Requirements**

-  Python >= 2.7
-  A Working version of Underworld2 >=2.6.0 (Please refer to the
   Underworld documentation)
-  pint >= 0.8

.. note::
  The bleeding edge version of *Underworld* (development branch)
  is now python 3 compatible only.
  *UWGeodynamics* is python 3 ready and can thus be used with it.

**Install**

**from Pip**

The UWGeodynamics module can be installed directly from the Python
Package Index:

.. code:: bash

  pip install UWGeodynamics

**from sources**

The module source files are available through github_

.. code:: bash

  git clone https://github.com/underworldcode/UWGeodynamics.git

It can then be installed globally on your system using

.. code:: bash

  pip install UWGeodynamics/


Seeking Support?
----------------

Error messages are useful to understand the source of a problem.

If you cannot solve the problem by yourself you can ask for help by creating an
issue on GitHub. If the problem if specific to your model you may be ask to continue the conversation
through email. 

*UWGeodynamics* is an open source free software and we cannot guarantee that it
is free of bugs. Feel free to signal any strange behaviour by raising an issue (see below section
on how to contribute.)


Contributing
------------

If you want to contribute to the UWGeodynamics projects and make it better, your help is very welcome.

So how can you contribute?

- Found a bug? Submit an issue using the issue tracker here on GitHub
- Have some suggestions? You can create an issue. Just add [Feature Request] in the title.

If you have developed some code and you think that it should be included in UWGeodynamics, you
can create a Pull Request and We will be happy to review it.

How to create a Pull Request (PR)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Create a personal fork of the project on Github.

   You will need a Github Account to do that. Just click on
   the Fork button at the top right corner of this repository.

#. Clone the fork on your local machine. Your remote repo on Github is called origin.

   :code:`git clone https://github.com/your-github-name/UWGeodynamics`

   replacing "your-github-name" with your actual github name...

#. Add the original repository as a remote called upstream.

   :code:`git add remote upstream https://github.com/underworldcode/UWGeodynamics`

#. If you created your fork a while ago be sure to pull upstream changes into your local repository.

   :code:`git pull upstream`

#. Create a new branch to work on! Branch from development!

   :code:`git checkout upstream/development`
   :code:`git checkout -b newFeature`

#. Implement/fix your feature, comment your code.

#. Follow the code style of the project, including indentation.

#. Include some tests or usage cases

#. Add or change the documentation as needed.
   The UWGeodynamics documentation is located in the `docs` directory.

#. Push your branch to your fork on Github, the remote origin.
   :code:`git push origin newFeature`

#. From your fork open a pull request in the correct branch. Target the project's `development`.

Always write your commit messages in the present tense.
Your commit message should describe what the commit, when applied, does to the code – not what you did to the code.


There is no small contribution!


Community driven
----------------

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/lgpl-3.0.en.html>.


.. _binder: https://mybinder.org/v2/gh/rbeucher/UWGeodynamics-binder/master
.. _Underworld: https://github.com/underworldcode/underworld2
.. _Jupyter: http://jupyter.org/
.. _Docker: https://www.docker.com
.. _Docker Hub: https://hub.docker.com/r/underworldcode/uwgeodynamics
.. _Kitematic: https://kitematic.com/
.. _github: https://github.com/underworldcode/UWGeodynamics.git
.. _Pint: https://pint.readthedocs.io/en/latest
