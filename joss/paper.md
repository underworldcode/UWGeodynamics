---
title: 'UWGeodynamics: A teaching an research tool for numerical geodynamic modelling'
tags:
  - Python
  - geodynamics
  - geophysics
  - geology
  - earth science
  - tectonics
  - rift
  - subduction
  - basin evolution
  - mountains
  - sandbox experiment
  - stokes
authors:
 - name: Romain Beucher
   orcid:0000-0003-3891-5444 
   affiliation: "1"
 - name: Louis Moresi
   orcid: 0000-0003-3685-174X
   affiliation: "1"
 - name: Juian Giordani
   orcid:
   affiliation: "1"
 - name: John Mansour
   orcid:
   affiliation: "2"
 - name: Rebecca Farrington
   orcid: 0000-0002-2594-6965 
   affiliatioin: "1"
 - name: Luke Mondy
   orcid: 0000-0001-7779-509X
   affiliation:  "3"
 - name: Claire Mallard
   orcid: 0000-0003-2595-2414
   affiliation: "3"
 - name: Patrice Rey
   orcid: 0000-0002-1767-8593 
   affiliation: "3"
 - name: Guillaume Duclaux
   orcid: 0000-0002-9512-7252
   affiliation: "4"
 - name: Owen Kaluza
   orcid: 0000-0001-6303-5671
   affiliation: "2"
 - name: Arijit Laik
   orcid: 0000-0002-3484-7985
   affiliation: "5"
 - name: Sara Morón
   orcid: 0000-0002-1270-437
   affiliation: "1"  

affiliations:
 - name: School of Earth Science, The University of Melbourne, Melbourne, Australia
   index: 1
 - name: Monash eResearch Centre, Monash University, Clayton, Australia 
   index: 2
 - name: School of Geosciences, Earthbyte Research Group, The University of Sydney, Australia
   index: 3
 - name: Laboratoire Géoazur, Université Nice Sophia Antipolis, Nice, France
   index: 4
 - name: Department of Earth Science, Faculty of Science, Vrije Universiteit, Amsterdam 
   index: 5

date: 05 October 2018
bibliography: paper.bib
---

# Summary

UWGeo is a high-level python interface to the Underworld API. 
It provides a framework for rapid prototyping of mechanical and 
thermo-mechanical numerical models both in 2D and 3D.

UWGeo functionalities includes:

* Handles units through pint
* semi-automatic scaling through definition of characteristic variables.
* Function to define 2D and 3D shapes from simple boxes to more complex composite shapes.
* Lithostatic pressure calculation
* Isostatic module.
* Partial melt calculation and associated change in viscosity / heat production.
* Viscous, Visco-plastic, Visco-elasto-plastic rheologies implementation.
* Database of common rheologies used in Geodynamics.
* Possibility to import personal rheology databases.
* von Mises, Coulomb and Drucker-prager yield criterion.
* Phase changes
* Frictional boundary condition
* Indentor definition
* Simple Adaptive mesh

UWGeo comes with a series of examples and tutorials such as:

* Simple


# References
