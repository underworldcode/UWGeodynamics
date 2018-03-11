# Underworld Geodynamics project

[![DOI](https://zenodo.org/badge/114189389.svg)](https://zenodo.org/badge/latestdoi/114189389)

The UWGeodynamics module intents to facilitate rapid prototyping of geodynamics models using Underworld. 
It can be seen as a set of high-level functions within the Underworld ecosystem. 
It is a means to quickly get the user into Underworld modelling and assumes very
little knowledge in coding. The module make some assumptions based on how the user
defines the boundary conditions and the properties of the materials (rocks, phases).
Its simplicity comes with a relatively more rigid workflow (compared to the classic Underworld functions).
However, the user can easily break the high level objects and get back to core
Underworld function at any step of model design.

As we think the low-level interface is more flexible, and in so allows for more complex models,
we strongly encourage users to explore and break the High Level functions.

We hope that the user will naturally move to the low-level functionalities as he
or her gets more confident, and by doing so will access the wide range of 
possibilities offered by Underworld.

![https://github.com/rbeucher/UWGeodynamics/wiki/img/cover.png](https://github.com/rbeucher/UWGeodynamics/wiki/img/cover.png)


## Getting started

For installation information and documentation visit our github [**wiki page**](https://github.com/rbeucher/UWGeodynamics/wiki) which provides several useful notes on how to start using the tool.

The easiest way to get started is with the [Docker container](https://hub.docker.com/r/rbeucher/underworld2_geodynamics/) using [Kitematic](https://docs.docker.com/kitematic/userguide/). Once **Kitematic** is installed on your computer, open it and look for **underworld2_geodynamics** via the *search* menu.

The latest UWGeodynamics version is the one that’s in our Github [repository](https://github.com/rbeucher/UWGeodynamics). Get it using this shell command, which requires Git: 
* `git clone https://github.com/rbeucher/UWGeodynamics.git`


## Installation

### Requirements

An up-to-date working version of Underworld is required.

UWGeodynamics is available via pip:

```
   pip install UWGeodynamics
```

## Community driven

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/lgpl-3.0.en.html>.

## Versioning

Current releases (**DOI** citable): 

[![DOI](https://zenodo.org/badge/114189389.svg)](https://zenodo.org/badge/latestdoi/114189389)

