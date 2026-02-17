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

#### HRMFoam_OFv25
HRMFoam source code for ESI OpenFOAM v25xx. 
Modified from previous version to account for API changes in fvPatchField.
Should be backward compatible with OFv24

#### HRMFoam_OFv24
HRMFoam source code for ESI OpenFOAM v22xx up to v24xx.
Edits to the latest version are maintaned concurrently with this one, so both should work.

#### Feature support
Both OFv24 and OFv25 currently supports:
- Sigma Y with URANS (k-Epsilon, k-Omega)
- Noncondensible gas transport
- LES, but not with Sigma Y.

Currently NOT supported
- Sigma Y with Large Eddy Simulation (Omega Eqn needs to be updated)

#### Old Versions
Older versions of HRMFoam are contained in the old_versions subdirectory.

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

Several test cases have been provided, in 2d and 3d geometries.
The 3d cases will take quite some time to run and will need multiple processors.

Git LFS is used to store the grids; to download these you will need to have git-lfs
installed, i.e. `apt-get install git-lfs` on Ubuntu or `brew install git-lfs` on MacOS.
Then run `git lfs install` once.

More information: https://git-lfs.com

#### 2d_r134a/Wedge_Throttle_laminar
Sample case for an axisymmetric wedge of a throttled flow of R-134a.
No turbulence model - minimal example.

#### 2d_r134a/Wedge_Throttle_kEps
Sample case for an axisymmetric wedge of a throttled flow of R-134a into air.
Uses k-Epsilon URANS but Sigma Y is disabled.

#### 2d_r134a/Wedge_Throttle_kEPS+SY
Sample case for an axisymmetric wedge of a throttled flow of R-134a into air.
Uses k-Epsilon URANS with Sigma Y enabled.

#### 2d_r134a/Wedge_Throttle_LES
Sample case for an axisymmetric wedge of a throttled flow of R-134a into air.
Uses WALE LES, which is technically inappropriate for a 2D geometry but this serves as a technical demonstration that WALE LES works with HRMFoam.
Sigma Y is disabled, and is currently not supported with LES.

#### 3d_water/3d_flashingHotWaterNozzle 
A 3d snappyHexMesh generated simulation for flashing hot water into 1 bar ambient pressure through an 0.4x1.0mm nozzle.
Based on Ankit Rawat's PhD experiment.

#### 3d_water/3d_mobyDick
A 3d snappyHexMesh generated simulation for flash evaporation/cavitation of water from ~1.38 bar at ~80C into a partial vacuum in a narrow pipe,
simulating the 'Moby Dick' experiments of Reocreux et al. Based on James Puli's PhD experiments.
