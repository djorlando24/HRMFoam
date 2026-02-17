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
#include "isoThermHFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "./propReader/refprop/refprop.H"
#include "./propReader/utrc/utrc.H"
#include "./modelCal/modelCalc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

isoThermHFvPatchScalarField::isoThermHFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    tVal_(0.0),
    phiName_("phi")
{
    this->refValue() = pTraits<scalar>::zero;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = 0.0;
}


isoThermHFvPatchScalarField::isoThermHFvPatchScalarField
(
    const isoThermHFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    tVal_(ptf.tVal_),
    phiName_(ptf.phiName_)
{}


isoThermHFvPatchScalarField::isoThermHFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    tVal_(readScalar(dict.lookup("T"))),
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


isoThermHFvPatchScalarField::isoThermHFvPatchScalarField
(
    const isoThermHFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    tVal_(ptf.tVal_),
    phiName_(ptf.phiName_)
{}


isoThermHFvPatchScalarField::isoThermHFvPatchScalarField
(
    const isoThermHFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    tVal_(ptf.tVal_),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void isoThermHFvPatchScalarField::updateCoeffs()
{
   if (updated())
   {
       return;
   }
   const IOdictionary& thermophysicalProperties = this->db().objectRegistry::lookupObject<IOdictionary>("thermophysicalProperties");
   const refprop& tPoint =
        this->db().objectRegistry::lookupObject<refprop>("tPoint");
   fvsPatchField<scalar> phip = this->patch().lookupPatchField
   (
        phiName_,
        reinterpret_cast<const surfaceScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
   );

   const fvPatchField<scalar>& P = this->patch().lookupPatchField
   (
        "p",
        reinterpret_cast<const volScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
   );

   const fvPatchField<scalar>& X = this->patch().lookupPatchField
   (
        "x",
        reinterpret_cast<const volScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
   );
   const fvPatchField<scalar>& Y = this->patch().lookupPatchField
   (
        "y",
        reinterpret_cast<const volScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
   );
   const fvPatchField<scalar>& T = this->patch().lookupPatchField
   (
        "T",
        reinterpret_cast<const volScalarField*>(NULL),
        reinterpret_cast<const scalar*>(NULL)
   );


   dimensionedScalar hDatum("hDatum", thermophysicalProperties);
   dimensionedScalar cpGas("cpGas", thermophysicalProperties);
   Field<scalar>& patchField = *this;

   forAll(patchField,facei)
   {
	//Bisection method to find required enthalpy.  GLJ 9/2018
        double up=tPoint.hEnd();
	double low=tPoint.hStart();
	double enth=0.0;
	double temp=0.0;
	double tol=1.0;
	int count=0;
	while(tol>0.0001 && count < 100){
		enth=(up+low)*0.5;
		tPoint.getProperties(P[facei],enth,X[facei]);
		//temp=tPoint.temp(); //Pure fuel
		temp=tPoint.temp()*(1-Y[facei])+(enth-hDatum.value())/cpGas.value()*Y[facei]; //All phases
		tol=temp-tVal_;
		(tol > 0) ? up=enth : low=enth;
		tol=tol*tol; //Make sure it's positive
		count ++;
	}
	//Info << "Your wastefulness is " << count << endl;
	patchField[facei] = enth;
   }
   if(min(T) < tVal_ - (tVal_*0.01) || max(T) > tVal_ + (tVal_*0.01)){
		Info << "WARNING: Isothermal BC not converged." << endl;
                Info << "         Ensure fluid.dat captures desired T." << nl << endl;
	}
   this->refValue() = patchField;
   this->valueFraction() = 1; //Doing this as a mixed BC was pretty dumb
   this->refGrad() = pTraits<scalar>::zero;

   mixedFvPatchScalarField::updateCoeffs();
}


void isoThermHFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi")
            << phiName_ << token::END_STATEMENT << nl;
    }
    os.writeKeyword("T")<<tVal_<<token::END_STATEMENT <<nl;
    this->refValue().writeEntry("inletValue", os);
    this->writeEntry("value", os);
}

makePatchTypeField(fvPatchScalarField, isoThermHFvPatchScalarField);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
