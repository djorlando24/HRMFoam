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
#include "thermDensityFvPatchScalarField.H"
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

  thermDensityFvPatchScalarField::thermDensityFvPatchScalarField
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


  thermDensityFvPatchScalarField::thermDensityFvPatchScalarField
  (
   const thermDensityFvPatchScalarField& ptf,
   const fvPatch& p,
   const DimensionedField<scalar, volMesh>& iF,
   const fvPatchFieldMapper& mapper
   )
    :
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
  {}


  thermDensityFvPatchScalarField::thermDensityFvPatchScalarField
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


  thermDensityFvPatchScalarField::thermDensityFvPatchScalarField
  (
   const thermDensityFvPatchScalarField& ptf
   )
    :
    mixedFvPatchScalarField(ptf),
    phiName_(ptf.phiName_)
  {}


  thermDensityFvPatchScalarField::thermDensityFvPatchScalarField
  (
   const thermDensityFvPatchScalarField& ptf,
   const DimensionedField<scalar, volMesh>& iF
   )
    :
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_)
  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


  void thermDensityFvPatchScalarField::updateCoeffs()
  {
    if (updated())
      {
	return;
      }

    const thermoBase& tPoint =
      this->db().objectRegistry::lookupObject<thermoBase>("tPoint");

    const IOdictionary& thermProps =
      this->db().objectRegistry::lookupObject<IOdictionary>("thermophysicalProperties");


    dimensionedScalar   rhoMin("rhoMin", thermProps);


    Switch comp(thermProps.lookup("compressible"));

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

    const fvPatchField<scalar>& H = this->patch().lookupPatchField
      (
       "h",
       reinterpret_cast<const volScalarField*>(NULL),
       reinterpret_cast<const scalar*>(NULL)
       );


    const fvPatchField<scalar>& x = this->patch().lookupPatchField
      (
       "x",
       reinterpret_cast<const volScalarField*>(NULL),
       reinterpret_cast<const scalar*>(NULL)
       );


    const fvPatchField<scalar>& y = this->patch().lookupPatchField
      (
       "y",
       reinterpret_cast<const volScalarField*>(NULL),
       reinterpret_cast<const scalar*>(NULL)
       );

     const fvPatchField<scalar>& pSat = this->patch().lookupPatchField
      (
       "pSat",
       reinterpret_cast<const volScalarField*>(NULL),
       reinterpret_cast<const scalar*>(NULL)
       );

    Field<scalar>& patchField = *this;

    if (comp)
      {
	dimensionedScalar Rgas("Rgas", thermProps);

	dimensionedScalar cpGas("cpGas", thermProps);

	dimensionedScalar hDatum("hDatum", thermProps);

	dimensionedScalar aL("aL", thermProps);

	dimensionedScalar aV("aV", thermProps);

	dimensionedScalar psiv = 1/sqr(aV);
	dimensionedScalar psil = 1/sqr(aL);

	forAll(patchField,facei)
	  {
	    // use total, not partial pressure
	    tPoint.getProperties(P[facei],H[facei], x[facei]);
	    scalar rhoNCgas = P[facei]/((H[facei]-hDatum.value())/cpGas.value())/Rgas.value();
	    scalar rhoVCalc = tPoint.rhov(); //Uncorrected saturation densities
	    scalar rhoLCalc = tPoint.rhol();

	      if(pSat[facei]>P[facei]) //Extrapolate with specific volume to avoid negative rhoV at low pressures
	      {
		scalar vSat= 1/rhoVCalc;
		rhoVCalc= 1/(vSat - sqr(vSat)*psiv.value()*(P[facei]-pSat[facei]));
	      }
	      else rhoVCalc += psiv.value()*(P[facei]-pSat[facei]); //Otherwise extrapolate using density

	      rhoLCalc += psil.value()*(P[facei]-pSat[facei]);
	      rhoVCalc = max(rhoVCalc,rhoMin.value());

	    patchField[facei] = 1.0/(
				     (  ( tPoint.xbar()/max(rhoVCalc,SMALL) + (1-tPoint.xbar())/max(rhoLCalc,SMALL) )
					*  (1-y[facei])
					+ y[facei]/rhoNCgas
					)
				     );

            // DPS--I think this is OK for when we are above the fluid.dat
            // since the pSat value will top out at the same place as the
            // rho value.


	  }
      }
    else
      {

	forAll(patchField,facei)
	  {
	    // use total, not partial pressure
	    tPoint.getProperties(P[facei],H[facei], x[facei]);

	    dimensionedScalar rhoNCgas("rhoNCgas", thermProps);
	    //for incompressible cases, the subcooled density equals the saturated liquid density
	    patchField[facei] = 1.0/(
				     (  ( tPoint.xbar()/max(tPoint.rhov(),SMALL)
					  + (1-tPoint.xbar())/max(tPoint.rhol(),SMALL) )
					*  (1-y[facei])
					+ y[facei]/rhoNCgas.value()
					)
				     );

	  }
      }

    // limit density
    forAll(patchField,facei)
      {
	patchField[facei] = max ( patchField[facei] , rhoMin.value() );
      }

    this->refValue() = patchField;
    this->valueFraction() = 1- pos(phip);
    this->refGrad() = pTraits<scalar>::zero;

    mixedFvPatchScalarField::updateCoeffs();
  }


  void thermDensityFvPatchScalarField::write(Ostream& os) const
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

  makePatchTypeField(fvPatchScalarField, thermDensityFvPatchScalarField);
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
