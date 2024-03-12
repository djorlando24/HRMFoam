# HRMFoam

**The Homogeneous Relaxation Model for ESI OpenFOAM with Sigma-Y droplet model.**

Copyright 2006-2010 David Schmidt, Shiva Gopalakrishnan, Kshitij Neroorkar
University of Massachusetts Amherst

Repository maintainer: Dr Daniel Duke
Department of Mechanical & Aerospace Engineering
Monash University, Australia

## Notice to users (including students)

Access to this repository must be approved by the maintainers.
Do NOT share the code in this repository without prior permission.

## License
This file is part of HRMFoam, which includes software from OpenFOAM.

OpenFOAM is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with OpenFOAM; if not, write to the Free Software Foundation,
Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

## Contents of this repository

### HRMFoam versions

#### HRMFoam_OFv2212
HRMFoam source code for ESI OpenFOAM. Currently supports:
- Sigma Y with URANS (k-Epsilon, k-Omega)
- Noncondensible gas transport
- LES, but not with Sigma Y.

Currently NOT supported
- Sigma Y with Large Eddy Simulation (Omega Eqn needs to be updated)

As currently used on LCRC Improv cluster.

#### HRMFoam_fe32

HRMFoam for foam-extend v3.2.
As used on LCRC Blues and Bebop clusters.
This version came from the SVN repo on Cod.

#### HRMFoam_fe32_LESStatic

HRMFoam for foam-extend v3.2, for use with static meshes.
Contains a rough LES implementation of the Omega eqn.
Writes psi2phase field out for calculating isenthalpic Mach number.
As used by DD on LCRC Bebop.

### Property table generators

#### GE_Flash_moleFrac

GEFlash isenthalphic flash evaporation property table generator for ethanol blends.
Writes out the ethanol mole fractions instead of thermal conductivity data.

#### generateTable

REFPROP simple (single component or ideal mixture) property table generator.

### Sample cases

#### Wedge_Throttle_kEps
Sample case for an axisymmetric wedge of a throttled flow of R-134a into air.
Uses k-Epsilon URANS but Sigma Y is disabled.

#### Wedge_Throttle_kEPS+SY
Sample case for an axisymmetric wedge of a throttled flow of R-134a into air.
Uses k-Epsilon URANS with Sigma Y enabled.

#### Wedge_Throttle_LES
Sample case for an axisymmetric wedge of a throttled flow of R-134a into air.
Uses WALE LES, which is technically inappropriate for a 2D geometry but this serves as a technical demonstration that WALE LES works with HRMFoam.
Sigma Y is disabled, and is currently not supported with LES.
