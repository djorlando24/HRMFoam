/*---------------------------------------------------------------------------*\
 
License
    Written by David P. Schmidt, UMass Amherst

    This is free software; you can redistribute it and/or modify it
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
    Foam::fixedHeatFluxFvPatchScalarField

Description
    Foam::fixedHeatFluxFvPatchScalarField

SourceFiles
    fixedHeatFluxFvPatchScalarField.C

This boundary condition takes heatFlux as a parameter (positive values means
that the heat flow is directed into the domain) and imposes the correct 
boundary condition on enthalpy.  The value of heat flux is in [W/sq. m].

\*---------------------------------------------------------------------------*/

#ifndef fixedHeatFluxFvPatchScalarFields_H
#define fixedHeatFluxFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fixedHeatFluxFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedHeatFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
       scalarField heatFlux_;


public:

    //- Runtime type information
    TypeName("fixedHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        fixedHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedHeatFluxFvPatchScalarField onto
        //  a new patch
        fixedHeatFluxFvPatchScalarField
        (
            const fixedHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedHeatFluxFvPatchScalarField
        (
            const fixedHeatFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedHeatFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedHeatFluxFvPatchScalarField
        (
            const fixedHeatFluxFvPatchScalarField&,
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
                new fixedHeatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
