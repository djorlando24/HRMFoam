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

#include "injNeedleTopoFvMesh.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "pointField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(injNeedleTopoFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        injNeedleTopoFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::injNeedleTopoFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
    )
    {
        Info<< "void injNeedleTopoFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        if (topoChanger_.size() == 0)
        {
            FatalErrorIn
            (
                "void injNeedleTopoFvMesh::addZonesAndModifiers()"
            )   << "Mesh modifiers not read properly"
                << abort(FatalError);
        }

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    // Add zones
    List<pointZone*> pz(1);
    List<faceZone*> fz(6);
    List<cellZone*> cz(0);

         label nPointZones = 0;
       label nFaceZones = 0;


    // Add an empty zone for cut points

    pz[nPointZones] = new pointZone
    (
        "cutPointZone",
        labelList(0),
        nPointZones,
        pointZones()
    );
    
    nPointZones++;
    // Do face zones for slider

    // Inner slider
    const word innerSliderName(motionDict_.subDict("slider").lookup("inside"));
    const polyPatch& innerSlider =
        boundaryMesh()[boundaryMesh().findPatchID(innerSliderName)];

    labelList isf(innerSlider.size());

    //include face ids of all faces in the patch
    forAll (isf, i)
    {
        isf[i] = innerSlider.start() + i;
    }

    fz[nFaceZones] = new faceZone
    (
        "insideSliderZone",
        isf,
        boolList(innerSlider.size(), false),
        nFaceZones,
        faceZones()
    );

  nFaceZones++;

    // Outer slider
    const word outerSliderName
    (
        motionDict_.subDict("slider").lookup("outside")
    );

    const polyPatch& outerSlider =
        boundaryMesh()[boundaryMesh().findPatchID(outerSliderName)];

    labelList osf(outerSlider.size());

    //include face ids of all faces in the patch
    forAll (osf, i)
    {
        osf[i] = outerSlider.start() + i;
    }

    fz[nFaceZones] = new faceZone
    (
        "outsideSliderZone",
        osf,
        boolList(outerSlider.size(), false),
        nFaceZones,
        faceZones()
    );

   nFaceZones++;

    // Add empty zone for cut faces
    fz[nFaceZones] = new faceZone
    (
        "cutFaceZone",
        labelList(0),
        boolList(0, false),
        nFaceZones,
        faceZones()
    );

   nFaceZones++;


    // Add face zone for layer addition
    const word layerPatchName
    (
        motionDict_.subDict("layer").lookup("needle")
    );

    const polyPatch& layerPatch =
        boundaryMesh()[boundaryMesh().findPatchID(layerPatchName)];

    labelList lpf(layerPatch.size());


    forAll (lpf, i)
    {
        lpf[i] = layerPatch.start() + i;
    }

    fz[nFaceZones] = new faceZone
    (
        "needleLayerZone",
        lpf,
        boolList(layerPatch.size(), true),
        nFaceZones,
        faceZones()
    );

   nFaceZones++;




   //additng ggi zones 

   Info << "Adding zone 1 for ggi" << endl;
   const word Ggi1Name
   (
       motionDict_.subDict("nonConformal").lookup("patch1")
   );

   const polyPatch& GgiPatch1 =
            boundaryMesh()[boundaryMesh().findPatchID(Ggi1Name)];

   labelList GgiPatch1Labels(GgiPatch1.size(), GgiPatch1.start());


   forAll (GgiPatch1Labels, i)
   {
          GgiPatch1Labels[i] += i;
   }


   fz[nFaceZones] =
        new faceZone
           (
                GgiPatch1.name()+"Zone",
                GgiPatch1Labels,
                boolList(GgiPatch1Labels.size(), false),
                nFaceZones,
                faceZones()
           );
   nFaceZones++;


   Info << "Adding zone 2 for ggi" << endl;
 
   const word Ggi2Name
   (
        motionDict_.subDict("nonConformal").lookup("patch2")
   );
   
   const polyPatch& GgiPatch2 =
         boundaryMesh()[boundaryMesh().findPatchID(Ggi2Name)];

   labelList GgiPatch2Labels(GgiPatch2.size(), GgiPatch2.start());


   forAll (GgiPatch2Labels, i)
        {
              GgiPatch2Labels[i] += i;
        }


   fz[nFaceZones] =
      new faceZone
         (
             GgiPatch2.name()+"Zone",
             GgiPatch2Labels,
             boolList(GgiPatch2Labels.size(), false),
             nFaceZones,
            faceZones()
         );
   nFaceZones++;
 


   Info << "Adding point and face zones" << endl;
   fz.setSize(nFaceZones);
   pz.setSize(nPointZones);

   addZones(pz, fz, cz);

   // Add topology modifiers
   topoChanger_.setSize(2);
   topoChanger_.set
   (
    0,
    new slidingInterface
       (
            "needleSlider",
            0,
            topoChanger_,
            "outsideSliderZone",
            "insideSliderZone",
            "cutPointZone",
            "cutFaceZone",
            outerSliderName,
            innerSliderName,
            slidingInterface::INTEGRAL,   // Edge matching algorithm
            true,                         // Attach-detach action
            intersection::VISIBLE         // Projection algorithm
       )
   );

   

   topoChanger_.set
    (
        1,
        new layerAdditionRemoval
        (
            "needleLayer",
            1,
            topoChanger_,
            "needleLayerZone",
            readScalar
            (
                motionDict_.subDict("layer").lookup("minThickness")
            ),
            readScalar
            (
                motionDict_.subDict("layer").lookup("maxThickness")
            )
        )
    );

    // Write mesh and modifiers
   topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
   topoChanger_.write();
   write();
}


void Foam::injNeedleTopoFvMesh::makeLayersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable layering
    forAll (topoChanges, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].disable();
        }
        else
        {
            FatalErrorIn("void injNeedleTopoFvMesh::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::injNeedleTopoFvMesh::makeSlidersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable sliding interface
    forAll (topoChanges, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanges[modI]))
        {
            topoChanges[modI].disable();
        }
        else if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else
        {
            FatalErrorIn("void injNeedleTopoFvMesh::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


bool Foam::injNeedleTopoFvMesh::attached() const
{
    const polyTopoChanger& topoChanges = topoChanger_;

    bool result = false;

    forAll (topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            result =
                result
             || refCast<const slidingInterface>(topoChanges[modI]).attached();
        }
    }

    // Check thal all sliders are in sync (debug only)
    forAll (topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            if
            (
                result
             != refCast<const slidingInterface>(topoChanges[modI]).attached()
            )
            {
                FatalErrorIn("bool injNeedleTopoFvMesh::attached() const")
                    << "Slider " << modI << " named "
                    << topoChanges[modI].name()
                    << " out of sync: Should be" << result
                    << abort(FatalError);
            }
        }
    }

    return result;
}


Foam::tmp<Foam::pointField> Foam::injNeedleTopoFvMesh::newPoints() const
{
    tmp<pointField> tnewPoints
    (
        new pointField(points())
    );

    pointField& np = tnewPoints();

    const word layerPatchName
    (
        motionDict_.subDict("layer").lookup("needle")
    );

    const polyPatch& layerPatch =
        boundaryMesh()[boundaryMesh().findPatchID(layerPatchName)];

    const labelList& patchPoints = layerPatch.meshPoints();


    const Time& t =time();
    
    scalar newVal = table_(t.value());
    scalar oldVal = table_(t.value() - t.deltaT().value());
    

    forAll (patchPoints, ppI)
    {

          np[patchPoints[ppI]] += ( axis_ * (newVal - oldVal) );

    }

    return tnewPoints;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::injNeedleTopoFvMesh::injNeedleTopoFvMesh(const IOobject& io)
:
   topoChangerFvMesh(io),
   motionDict_
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
   axis_(motionDict_.subDict("needleMotion").lookup("axis")),
   profileFile_(motionDict_.subDict("needleMotion").lookup("profileFile")),
   table_(profileFile_)
 {
    
    addZonesAndModifiers();
    // Normalize the axis
    axis_ /= mag(axis_) + VSMALL;
     // Set outOfBounds handling to clamp
    table_.outOfBounds(interpolationTable<scalar>::CLAMP);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::injNeedleTopoFvMesh::~injNeedleTopoFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::injNeedleTopoFvMesh::update()
{
    // Detaching the interface
    if (attached())
    {
        Info << "Decoupling sliding interfaces" << endl;
        makeSlidersLive();

        // Changing topology by hand
        topoChanger_.changeMesh();
    }
    else
    {
        Info << "Sliding interfaces decoupled" << endl;
    }

    // Perform layer action and mesh motion
    makeLayersLive();

    // Changing topology by hand
    {
        autoPtr<mapPolyMesh> topoChangeMap2 = topoChanger_.changeMesh();

        if (topoChangeMap2->morphing())
        {
            if (topoChangeMap2->hasMotionPoints())
            {
                Info << "Topology change; executing pre-motion" << endl;
                movePoints(topoChangeMap2->preMotionPoints());
            }
        }
    }

    // Move points
    movePoints(newPoints());

    // Attach the interface
    Info << "Coupling sliding interfaces" << endl;
    makeSlidersLive();

    // Changing topology by hand
    {
        // Grab old points to correct the motion
        pointField oldPointsNew = oldAllPoints();

        autoPtr<mapPolyMesh> topoChangeMap3 = topoChanger_.changeMesh();

        if (topoChangeMap3->morphing())
        {
            if (debug)
            {
                Info << "Moving points post slider attach" << endl;
            }

            pointField newPoints = allPoints();
            pointField mappedOldPointsNew(newPoints.size());

            mappedOldPointsNew.map(oldPointsNew, topoChangeMap3->pointMap());

            // Solve the correct mesh motion to make sure motion fluxes
            // are solved for and not mapped
            movePoints(mappedOldPointsNew);
            resetMotion();
            setV0();
            movePoints(newPoints);
        }
    }

    return true;
}


// ************************************************************************* //

