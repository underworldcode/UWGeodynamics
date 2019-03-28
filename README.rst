Underworld Geodynamics project
==============================

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

Getting started
---------------

The full documentation is available on `ReadTheDocs <http://uwgeodynamics.readthedocs.org/>`_

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

1. Create a personal fork of the project on Github.

   You will need a Github Account to do that. Just click on
   the Fork button at the top right corner of this repository.

1. Clone the fork on your local machine. Your remote repo on Github is called origin.

   `git clone https://github.com/your-github-name/UWGeodynamics`

   replacing "your-github-name" with your actual github name...

1. Add the original repository as a remote called upstream.

   `git add remote upstream https://github.com/underworldcode/UWGeodynamics`

1. If you created your fork a while ago be sure to pull upstream changes into your local repository.

   `git pull upstream`

1. Create a new branch to work on! Branch from development!

   `git checkout upstream/development`
   `git checkout -b newFeature`

1. Implement/fix your feature, comment your code.

1. Follow the code style of the project, including indentation.

1. Include some tests or usage cases

1. Add or change the documentation as needed.
   The UWGeodynamics documentation is located in the `docs` directory.

1. Push your branch to your fork on Github, the remote origin.
   `git push origin newFeature`

1. From your fork open a pull request in the correct branch. Target the project's `development`.

Always write your commit messages in the present tense.
Your commit message should describe what the commit, when applied, does to the code – not what you did to the code.


There is no small contribution!


Community driven
----------------

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/lgpl-3.0.en.html>.

Versioning
----------

Current releases (**DOI** citable): 

`DOI <https://zenodo.org/badge/114189389.svg)](https://zenodo.org/badge/latestdoi/114189389>`_

