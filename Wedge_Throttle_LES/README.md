# HRMFoam sample problem
## Axisymmetric wedge throttle geometry

**WALE (Wall-adapting Local Eddy Simulation)**
note: LES in an axisymmetric 2d grid is technically not physically feasible.
This is simply a technical demonstration that the code will execute properly.

Sigma-Y is OFF. Pls note that as of early 2024, LES + SIGMAY is not yet supported

Noncondensible gas (air) in outlet.
Fluid is saturated R-134a.
setFields has already been run on t=0 to establish initial condition of liquid in the nozzle and air outside. (0/y.gz written).
