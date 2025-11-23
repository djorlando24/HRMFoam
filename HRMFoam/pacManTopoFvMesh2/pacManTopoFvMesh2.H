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

Class
    Foam::pacManTopoFvMesh2

Description
    Sample topoChangerFvMesh that moves a cell zone in a given direction,
    introducing/removing layers around it.

    Cells in the movingCells zone shall be moved given the prescribed velocity
    and will be bounded in "front" and "back" by other cell zones.
    Layer addition/removal interfaces are inserted at boundaries between the
    moving zone and front and back, pointing outside of the moving cell zone

Author and rewrite
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    pacManTopoFvMesh2.C

\*---------------------------------------------------------------------------*/

#ifndef pacManTopoFvMesh2_H
#define pacManTopoFvMesh2_H

#include "topoChangerFvMesh.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class pacManTopoFvMesh2 Declaration
\*---------------------------------------------------------------------------*/

class pacManTopoFvMesh2
:
    public topoChangerFvMesh
{
    // Private data

        //- Dynamic mesh dictionary
        dictionary dict_;

        //- Moving cell zone name
        word movingCellsName_;

        //- Front face zone name
        word frontFacesName_;


        //- Motion control function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        //- Vertex motion mask
        scalarField motionMask_;

        label nLayers_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pacManTopoFvMesh2(const pacManTopoFvMesh2&);

        //- Disallow default bitwise assignment
        void operator=(const pacManTopoFvMesh2&);


        //- Add mixer zones and modifiers
        void addZonesAndModifiers();

        //- Mark motion vertices
        tmp<scalarField> calcMotionMask() const;


public:

    //- Runtime type information
    TypeName("pacManTopoFvMesh2");


    // Constructors

        //- Construct from database
        explicit pacManTopoFvMesh2(const IOobject& io);


    // Destructor

        virtual ~pacManTopoFvMesh2();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
