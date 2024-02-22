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
    Foam::timeVaryingTransonicTotalPressureFvPatchField

Description
    Foam::timeVaryingTransonicTotalPressureFvPatchField smoothly transitions from a totalPressure BC
    to zero gradient based on Mach number.

SourceFiles
    timeVaryingTransonicTotalPressureFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingTransonicTotalPressureFvPatchScalarField_H
#define timeVaryingTransonicTotalPressureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include  "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class timeVaryingTransonicTotalPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingTransonicTotalPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{

  // total pressure value
  scalar p0_ ;
  interpolationTable<scalar> timeSeries_;

protected:

    // Protected data

        word phiName_;


public:

    //- Runtime type information
    TypeName("timeVaryingTransonicTotalPressure");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingTransonicTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingTransonicTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingTransonicTotalPressureFvPatchField onto a new patch
        timeVaryingTransonicTotalPressureFvPatchScalarField
        (
            const timeVaryingTransonicTotalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingTransonicTotalPressureFvPatchScalarField
        (
            const timeVaryingTransonicTotalPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField > clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new timeVaryingTransonicTotalPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingTransonicTotalPressureFvPatchScalarField
        (
            const timeVaryingTransonicTotalPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new timeVaryingTransonicTotalPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
//        virtual void updateCoeffs();

        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "timeVaryingTransonicTotalPressureFvPatchField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
