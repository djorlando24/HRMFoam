# This setup script does the following:
# It removes 6 layers of the needle boundary cells. It then decomposes the case
# and adds the layers back with noHydro. This prevents processor boundaries from
# cutting into the 6 uniform layers and it allows for HRMFoam to start with 6 uniform
# motion layers. It also sets up the internal cellZones and faceZones needed for mesh
# motion.

# define the number or processors here
np=24
#define the number of uniform layers here
nUniformLayers=6


# First, input the number of uniform layers into the system/controlDict_addRemoveLayers
nEndTime="$(($nUniformLayers+1))"e-09
echo "$:1,\$s/endTime  .*/endTime  $nEndTime;/" > edit.vim
echo ":wq" >> edit.vim
vi -s edit.vim system/controlDict_addRemoveLayers
rm edit.vim
# and input it into the nUniformLayers.dat file. This is read by the pacMan3 library in the
# noHydro solver, and used when adding the very first layer, so that it knows where to put it.
echo "$nUniformLayers" > constant/nUniformLayers.dat

# next, see if there is a 0orig directory. If not, make one.
if [ ! -d 0orig ]; then
	  cp -r 0 0orig
fi

# next, setup the cellZones and faceZones which are used in pacManTopoFvMesh3Coeffs as frontFaces and movingCells. Here they are called "movingFaces" and "movingCells" 
rm -r constant/polyMesh/sets constant/polyMesh/*Zones* 
setSet -batch setBatch
rm -r constant/polyMesh/sets/*_old*
setsToZones

# next, lower the mesh with noHydro, removing layers
cp constant/pistonStroke.dat_down constant/pistonStroke.dat
cp system/controlDict_addRemoveLayers system/controlDict
cp constant/dynamicMeshDict_addRemoveLayers constant/dynamicMeshDict
rm constant/polyMesh/meshModifiers*
noHydro &> noHydroLog &
wait

# copy that mesh into the constant/polyMesh directory for proper decomposition
mv constant/polyMesh constant/polyMeshOrig
mv $nEndTime/polyMesh constant/
sh tmpRemove
rm -r 0
cp -r 0orig 0

# now decomposePar and clean parallel face zones
decomposePar
mpirun -np $np cleanParallelFaceZones -parallel -overwrite

# then add the layers back using noHydro
# note that the addLayers routine in noHydro has been modified
# to have proper layer thickness for 6 initial layers
rm constant/polyMesh/meshModifiers* processor*/0/polyMesh/meshModifiers*
cp constant/meshModifiers_save constant/polyMesh/meshModifiers
cp constant/pistonStroke.dat_orig constant/pistonStroke.dat
mpirun -np $np noHydro -parallel >> noHydroLog &
wait

# move the new mesh to time 0 and clean case with tmpRemove
rm -r processor*/0
for i in `seq 0 $(($np-1))`;
 do
     echo "Moving processor file" $i
     mv ./processor$i/$nEndTime processor$i/0
 done 

# clean case and setup for HRMFoam
sh tmpRemove
rm constant/polyMesh/meshModifiers* processor*/0/polyMesh/meshModifiers* constant/nUniformLayers.dat
cp constant/meshModifiers_run constant/polyMesh/meshModifiers
cp system/controlDict_run system/controlDict
cp constant/dynamicMeshDict_run constant/dynamicMeshDict

#set the initial conditions
mpirun -np $np funkySetFields -parallel -latestTime -field h -keepPatches -expression '-325457-58913*tanh(4500*pos().x+.0027*4500)'
mpirun -np $np funkySetFields -parallel -latestTime -field p -keepPatches -expression '53000000+47000000*tanh(-4500*pos().x-.0027*4500)'
mpirun -np $np funkySetFields -parallel -latestTime -field y -keepPatches -expression 'faceAverage(fpos().x >= surf(-.00001) ? surf(1.0) : surf(0.))'
