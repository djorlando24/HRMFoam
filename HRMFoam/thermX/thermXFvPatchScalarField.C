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
#include "thermXFvPatchScalarField.H"
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

thermXFvPatchScalarField::thermXFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_("phi")
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


thermXFvPatchScalarField::thermXFvPatchScalarField
(
    const thermXFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{}


thermXFvPatchScalarField::thermXFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phi")
{


    if (dict.found("inletValue"))
         this->refValue() = Field<scalar>("inletValue", dict, p.size());
    else
         this->refValue() = pTraits<scalar>::zero;
   

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


thermXFvPatchScalarField::thermXFvPatchScalarField
(
    const thermXFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    phiName_(ptf.phiName_)
{}


thermXFvPatchScalarField::thermXFvPatchScalarField
(
    const thermXFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void thermXFvPatchScalarField::updateCoeffs()
{
   if (updated())
   {
       return;
   }
  
   const thermoBase& tPoint =
        this->db().objectRegistry::lookupObject<thermoBase>("tPoint");

   fvsPatchField<scalar> phip = this->patch().lookupPatchField<surfaceScalarField>(phiName_);

   const fvPatchField<scalar>& P = this->patch().lookupPatchField<volScalarField>("p");

   const fvPatchField<scalar>& H = this->patch().lookupPatchField<volScalarField>("h");

   const fvPatchField<scalar>& X = this->patch().lookupPatchField<volScalarField>("x");
   
   
   Field<scalar>& patchField = *this;
   forAll(patchField,facei)
   {	
        // use total, not partial pressure for property lookup

	tPoint.getProperties( P[facei],H[facei], X[facei]);
         
         //calculate equilibrium x from the fluid.dat table
 
	 patchField[facei] = tPoint.xbar();

   }
              
   this->refValue() = patchField;
   this->valueFraction() = 1- pos(phip);
   this->refGrad() = pTraits<scalar>::zero;
     
   mixedFvPatchScalarField::updateCoeffs();
}


void thermXFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi")
            << phiName_ << token::END_STATEMENT << nl;
    }
    this->refValue().writeEntry("inletValue", os);
    this->writeEntry("value", os);
}

makePatchTypeField(fvPatchScalarField, thermXFvPatchScalarField);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
