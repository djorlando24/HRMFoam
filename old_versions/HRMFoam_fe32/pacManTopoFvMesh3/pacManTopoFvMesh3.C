/*---------------------------------------------------------------------------*\
    Modifications copyright ICON 2013
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "pacManTopoFvMesh3.H"
#include "mapPolyMesh.H"
#include "layerAdditionRemoval.H"
#include "volMesh.H"
#include "transformField.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pacManTopoFvMesh3, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        pacManTopoFvMesh3,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::pacManTopoFvMesh3::calcMotionMask() const
{
    Info<< "Updating vertex markup" << endl;

     tmp<scalarField> tvertexMarkup(new scalarField(allPoints().size(), 0));
    scalarField& vertexMarkup = tvertexMarkup();

    cellZoneID movingCellsID(movingCellsName_, cellZones());

    // In order to do a correct update on a mask on processor boundaries,
    // Detection of moving cells should use patchNeighbourField for
    // processor (not coupled!) boundaries.  This is done by expanding
    // a moving cell set into a field and making sure that processor patch
    // points move in sync.  Not done at the moment, probably best to do
    // using parallel update of pointFields.  HJ, 19/Feb/2011

    // If moving cells are found, perform mark-up
    if (movingCellsID.active())
    {
        // Get cell-point addressing
        const labelListList& cp = cellPoints();

        // Get labels of all moving cells
        const labelList& movingCells = cellZones()[movingCellsID.index()];

        forAll (movingCells, cellI)
        {
            const labelList& curCp = cp[movingCells[cellI]];

            forAll (curCp, pointI)
            {
                vertexMarkup[curCp[pointI]] = 1;
            }
        }
    }

    faceZoneID frontFacesID(frontFacesName_, faceZones());

    if (frontFacesID.active())
    {
        const faceZone& frontFaces = faceZones()[frontFacesID.index()];

        const labelList& mp = frontFaces().meshPoints();

        forAll (mp, mpI)
        {
            vertexMarkup[mp[mpI]] = 1;
        }
    }


    return tvertexMarkup;
}


void Foam::pacManTopoFvMesh3::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if (topoChanger_.size() > 0)
    {
        Info<< "void pacManTopoFvMesh3::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    // Add layer addition/removal interfaces
    topoChanger_.setSize(1);
    label nMods = 0;


    faceZoneID frontFacesID(frontFacesName_, faceZones());

    if (frontFacesID.active())
    {
        const faceZone& frontFaces = faceZones()[frontFacesID.index()];

        if (!frontFaces.empty())
        {
            topoChanger_.set
            (
                nMods,
                new layerAdditionRemoval
                (
                    frontFacesName_ + "Layer",
                    nMods,
                    topoChanger_,
                    frontFacesName_,
                    readScalar
                    (
                        dict_.subDict("front").lookup("minThickness")
                    ),
                    readScalar
                    (
                        dict_.subDict("front").lookup("maxThickness")
                    )
                )
            );

            nMods++;
        }
    }


    topoChanger_.setSize(nMods);

    reduce(nMods, sumOp<label>());

    Info << "Adding " << nMods << " mesh modifiers" << endl;

    // Write mesh and modifiers
    topoChanger_.write();

    // No need to write the mesh - only modifiers are added.
    // HJ, 18/Feb/2011
//     write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pacManTopoFvMesh3::pacManTopoFvMesh3(const IOobject& io)
:
    topoChangerFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    movingCellsName_(dict_.lookup("movingCells")),
    frontFacesName_(dict_.lookup("frontFaces")),
    SBMFPtr_(solidBodyMotionFunction::New(dict_, time())),
    motionMask_(),
    nLayers_(0)
{
    addZonesAndModifiers();
    motionMask_ = calcMotionMask();

    IFstream nLayersFileIn
    (
        time().path()/time().timeName()/"uniform"
        /"nLayers_"+frontFacesName_+"_"+movingCellsName_+".raw"
    );

    if ( nLayersFileIn.good() )
    {
        nLayersFileIn >> nLayers_;
    }
    else
    {
        nLayers_ = 0;

        //check that uniform directory exists, create if not
        fileName uniformDir(time().path()/time().timeName()/"uniform");
        Foam::mkDir(uniformDir);
    
        OFstream nLayersFileOut
        (
            uniformDir/"nLayers_"+frontFacesName_+"_"+movingCellsName_+".raw"
        );

        nLayersFileOut << nLayers_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pacManTopoFvMesh3::~pacManTopoFvMesh3()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pacManTopoFvMesh3::update()
{
    // Store points to recreate mesh motion
    pointField oldPointsNew = allPoints();
    pointField newPoints = allPoints();
    label nOldPoints = newPoints.size();

    // add/remove layers
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();


    bool localMeshChanged = topoChangeMap->morphing();
    bool globalMeshChanged = localMeshChanged;
    reduce(globalMeshChanged, orOp<bool>());


    if (globalMeshChanged)
    {
        //update list of moving points 
        Pout<< "Topology change. Calculating motion point mask" << endl;
        motionMask_ = calcMotionMask();
    }

    if (localMeshChanged)
    {
    //         // Map old points onto the new mesh
    //         pointField mappedOldPointsNew(allPoints().size());
    //         mappedOldPointsNew.map(oldPointsNew, topoChangeMap->pointMap());

    //         movePoints(mappedOldPointsNew);
    //         resetMotion();
    //         setV0();

        // Get new points from preMotion
        newPoints = topoChangeMap().preMotionPoints();
    }
    //     else
    //     {
    //         // No change, use old points
    //         movePoints(oldPointsNew);
    //         resetMotion();
    //         setV0();
    //     }



    //find out if layer was added or removed and update number of layers
    if (globalMeshChanged)
    {
        if (nOldPoints < newPoints.size()) //layer has been added
        {
            nLayers_++;
        }       
        else if (nOldPoints > newPoints.size() && nLayers_ > 0) //layer has been removed
        {                                                      //nLayers can't go negative
            nLayers_--;
        }
    }

    Info << "Number of layers to be made uniform by this topoChanger: " << nLayers_ << endl;

    //for restart-capability write number of layers to time dump directory
    if (time().outputTime())
    {
        //check that uniform directory exists, create if not
        fileName uniformDir(time().path()/time().timeName()/"uniform");
        Foam::mkDir(uniformDir);

        OFstream nLayersFileOut
        (
            uniformDir/"nLayers_"+frontFacesName_+"_"+movingCellsName_+".raw"
        );

        nLayersFileOut << nLayers_ << endl;
    }

    //find mesh points of face zone and for each layer
    labelListList layerPoints(nLayers_+1);

    const faceZoneID frontFacesID(frontFacesName_, faceZones());
    labelList curFaces = faceZones()[frontFacesID.index()];
    labelList curCells = faceZones()[frontFacesID.index()].masterCells();

    //add face zone points to labelListList[0] 
    layerPoints[0] = faceZones()[frontFacesID.index()]().meshPoints();

    //loop over layers
    for (label layerI = 1; layerI <= nLayers_; layerI++)  //? what if nLayers = 0 or 1
    {
        faceList fl (curFaces.size());

        //loop over layer cells
        forAll (curCells, i)
        {
            // find opposing face of cell
            curFaces[i] = cells()[curCells[i]].opposingFaceLabel(curFaces[i], faces());

            // find other cell of face
            if (faceOwner()[curFaces[i]] == curCells[i])
            {
                curCells[i] = faceNeighbour()[curFaces[i]];
            }
            else if (faceNeighbour()[curFaces[i]] == curCells[i])
            {
                curCells[i] = faceOwner()[curFaces[i]];
            }
            else
            {
                FatalErrorIn ("bool Foam::pacManTopoFvMesh3::update()")
                    << "Detected face in layers without owner or neighbour in layer."
                    << abort(FatalError);
            }

            fl[i] = faces()[curFaces[i]];
        }
        
        //add current layer's point indices to layerPoints list
        primitiveFacePatch pp ( fl, points() );
        layerPoints[layerI] =  pp.meshPoints();
    }


//?? restartable

    // Calculate new points using a velocity transformation
    newPoints += motionMask_*
        transform(SBMFPtr_().velocity(), newPoints)*time().deltaT().value();

    //linearly distribute points between layers
    forAll (layerPoints[0], i)
    {
        vector start = newPoints[layerPoints[0][i]];
        vector dir = newPoints[layerPoints[layerPoints.size()-1][i]] - start;

        label layerI = 1;
        while ( layerI < layerPoints.size())  //first and last not needed
        {
            newPoints[layerPoints[layerI][i]] = start + dir * ( scalar(layerI) / nLayers_ ); 
            layerI++;
        }
    }


    //update mesh.
    Info << "Executing mesh motion" << endl;
    movePoints(newPoints);

    return localMeshChanged;
}


// ************************************************************************* //
