# Mesh generation using snappyHexMesh

    author Daniel Duke <daniel.duke@monash.edu>
    copyright (c) 2026 Daniel Duke
    license GPL-3.0+
    version 1.0.0
    date 20/02/2026

    Multiphase Flow Laboratory
    Monash University, Australia

# User guide

- Source geometry in constant/triSurface/*.stl
Note that the geometry is in mm, which snappyHexMesh prefers.
We will scale it down to m after mesh generation.

- Base mesh resolution is governed by constant/polyMesh/blockMeshDict
Also in mm to match the stl

- Refinements and final resolution is governed by system/snappyHexMeshDict

- run_local.sh will attempt parallel mesh generation. Number of CPUs must be set in system/decomposeParDict.? files (all must match)

- cleanup.sh will wipe away the mesh and start over, leaving the STL files intact.

- Once the mesh is generated, copy it to your simulation case:

        cp -r 3e-08/polyMesh $MYCASE/constant/

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
