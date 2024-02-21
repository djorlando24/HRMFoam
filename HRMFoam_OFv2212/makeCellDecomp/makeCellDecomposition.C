/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    makeCellDecomposition 

Description
    Under certain conditions one need precise control of how cells are 
    distributed over processors. This tool reads in cell sets named 
    processor0, processor1 ... and creates a cellDecomposition file that
    can be used with the manual decomposition method in decomposePar.

Usage
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "Time.H"
#include "cellSet.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//  Main program:
int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// TO DO: Make the app write a labelList file in the constant directory instead of this:

    labelList procIds(mesh.cells().size(), 0);
    
    volScalarField procID
    (
        IOobject
        (
            "cellDecomposition",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    int ID = 0;
    while(true)
    {
        if (!exists(mesh.time().path()/topoSet::localPath(mesh, "processor"+Foam::name(ID)))) 
        {
            Info << "No more processor cellSets found. Stopping." << endl;
            break;
        }

        Info << "Adding processor" << ID << endl;
        cellSet processorID(mesh, "processor"+Foam::name(ID));
          for
          (
            cellSet::const_iterator iter = processorID.begin();
            iter != processorID.end();
            ++iter
          )
          {
             procID[iter.key()] = ID;
             procIds[iter.key()] = ID;
          }

        ID++;
    }

    procID.write();
    Info << nl << "Wrote volScalarField "
         << procID.name() << " for visualization."
         << endl;


    labelIOList cellDecomposition
    (
        IOobject
        (
            "cellDecomposition",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
         ),
         procIds
     );

    cellDecomposition.write();
    Info << nl << "Wrote decomposition to "
         << cellDecomposition.objectPath()
         << " for use in manual decomposition.\n" << endl;

    
    return 0;
}


// ************************************************************************* //
