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
    Foam::solidBodyMotionFunctions::needleMotion2

Description
    Linear oscillation motion function

SourceFiles
    needleMotion2.C

\*---------------------------------------------------------------------------*/

#ifndef needleMotion2_H
#define needleMotion2_H

#include "solidBodyMotionFunction.H"
#include "interpolationTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                      Class needleMotion2 Declaration
\*---------------------------------------------------------------------------*/

class needleMotion2
:
    public solidBodyMotionFunction
{
    // Private data

        //- Motion amplitude
        vector amplitude_;

        //- Motion period
        scalar period_;


    // Private Member Functions

        //- Disallow copy construct
        needleMotion2(const needleMotion2&);

        //- Disallow default bitwise assignment
        void operator=(const needleMotion2&);


        //- Calculate position
        vector calcPosition(const scalar t) const;


      //- Axis of displacement
        vector axis_;

      //- Profile file containing displacement data

        fileName profileFile_;

      //- Interpolation table
        interpolationTable<scalar> table_;

      //- Profile file containing displacement data for side motion

        fileName sideProfileFile_;

      //- Interpolation table for side motion
        interpolationTable<vector> sideTable_;


public:

    //- Runtime type information
    TypeName("needleMotion2");


    // Constructors

        //- Construct from components
        needleMotion2
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    // Destructor

        virtual ~needleMotion2();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Return the solid-body motion velocity
        virtual septernion velocity() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
