#!/bin/bash
#
#    OPENFOAM SNAPPYHEXMESH SCRIPT
#
#    author Daniel Duke <daniel.duke@monash.edu>
#    copyright (c) 2020 Daniel Duke
#    license GPL-3.0+
#    version 1.0.0
#    date 01/07/2020
#        __   ____________    ___    ______
#       / /  /_  ____ __  \  /   |  / ____/
#      / /    / /   / /_/ / / /| | / /
#     / /___ / /   / _, _/ / ___ |/ /_________
#    /_____//_/   /_/ |__\/_/  |_|\__________/
#
#    Laboratory for Turbulence Research in Aerospace & Combustion (LTRAC)
#    Monash University, Australia

echo "snappyHexMesh local script"
date

# Load OpenFOAM
source ~/OpenFOAM/OpenFOAM-v2412/etc/bashrc
NCPU=8

. $WM_PROJECT_DIR/bin/tools/RunFunctions

# determine surface feature edges
surfaceFeatureExtract | tee log.surfaceFeatureExtract

# base mesh
blockMesh | tee log.blockMesh

# decompose
decomposePar | tee log.decomposePar

# run snappy
mpirun -np $NCPU snappyHexMesh -parallel | tee log.snappyHexMesh

if [ ! -d processor0/2e-08 ] ; then
    echo "I think something went wrong. Stopping"
    exit 1
fi

# Reconstruct mesh
reconstructParMesh -latestTime | tee log.reconstructParMesh

if [ ! -d 2e-08 ] ; then
    echo "I think something went wrong. Stopping"
    exit 1
fi

# transform co-ordinates (move origin & scale)
#transformPoints -translate '(-44 -4.5 -4.5)' >> log.transformPoints.translate
transformPoints -scale '(1e-3 1e-3 1e-3)' >> log.transformPoints.scale # convert base unit to mm

# Check mesh integrity
checkMesh -latestTime | tee log.checkMesh

# Extract subsets of bad geometry for inspection
for setfile in `ls -1 --color=no 2e-08/polyMesh/sets/*Face*` ; do
    faceSet=$(basename $setfile .gz)
    foamToVTK -latestTime -faceSet $faceSet | tee 'log.foamToVTK_'$faceSet
done

echo "Done."
date

