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
#include "transonicTotalPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "./propReader/refprop/refprop.H"
#include "./propReader/utrc/utrc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

transonicTotalPressureFvPatchScalarField::transonicTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    p0_(0.0),
    phiName_("phi")
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


transonicTotalPressureFvPatchScalarField::transonicTotalPressureFvPatchScalarField
(
    const transonicTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    p0_(ptf.p0_),
    phiName_(ptf.phiName_)
{}


transonicTotalPressureFvPatchScalarField::transonicTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    p0_(readScalar(dict.lookup("p0"))),
    phiName_("phi")
{


    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(this->refValue());
    }

    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;

    if (dict.found("phi"))
    {
        dict.lookup("phi") >> phiName_;
    }
}


transonicTotalPressureFvPatchScalarField::transonicTotalPressureFvPatchScalarField
(
    const transonicTotalPressureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    p0_(ptf.p0_),
    phiName_(ptf.phiName_)
{}


transonicTotalPressureFvPatchScalarField::transonicTotalPressureFvPatchScalarField
(
    const transonicTotalPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    p0_(ptf.p0_),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void transonicTotalPressureFvPatchScalarField::updateCoeffs()
{
   if (updated())
   {
       return;
   }
 
   fvsPatchField<scalar> phip = this->patch().lookupPatchField
   (
        phiName_,
        reinterpret_cast<const surfaceScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
   );

   const fvPatchField<scalar>& psi = this->patch().lookupPatchField
   (
        "psi",
        reinterpret_cast<const volScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
   );
   

         const fvPatchField<scalar>& rho = this->patch().lookupPatchField
   (
        "rho",
        reinterpret_cast<const volScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
   );
   
 
         const fvPatchField<vector>& Up = this->patch().lookupPatchField
   (
        "U",
        reinterpret_cast<const volVectorField*>(NULL),
        reinterpret_cast<const vector*>(NULL)
   );
     
   // BC = valueFraction*refValue 
   // + (1-valueFraction)*(internalValue+refGrad/deltaCoeff)

   this->refValue() = p0_ - 0.5*rho*(1.0 - pos(phip))*magSqr(Up);
   this->valueFraction() = Foam::max(1 - psi*magSqr(Up),0.0)*pos(phip)+(1-pos(phip));
   this->refGrad() = pTraits<scalar>::zero;
     
   mixedFvPatchScalarField::updateCoeffs();
}


void transonicTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi")
            << phiName_ << token::END_STATEMENT << nl;
    }
    this->refValue().writeEntry("inletValue", os);

    os.writeKeyword("p0")<<p0_<<token::END_STATEMENT <<nl;

    this->writeEntry("value", os);
}

makePatchTypeField(fvPatchScalarField, transonicTotalPressureFvPatchScalarField);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
