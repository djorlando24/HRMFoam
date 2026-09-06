# Sample axisymmetric case -strongly flashing hot water nozzle

d=0.6mm nozzle using grid created by Ankit Rawat (Monash University). 
Grid was made in gmsh, then extruded in OpenFOAM.

Using fluid.dat for water, set h = 0.110000E+07 and p = 60 bar inlet stagn pressure to obtain 526K (253C) compressed liquid at 0.796673E+03 kg/m3.

Back pressure ramps down over 10ms to 1.01 bar.  This can be UNSTABLE for strongly flashing cases and you should probably first try the 'sealConstants' approach for a large outlet plenum filled with air.

To run this case;
1. use ./clean and ./clean0 for fresh start
2. run setFields to create nonuniform initial condition
3. run decomposePar
4. runParallel HRMFoam

