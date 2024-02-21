/*---------------------------------------------------------------------------*\

License
    Written by David P. Schmidt, UMass Amherst 2014

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

\*---------------------------------------------------------------------------*/

#include "fixedHeatFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedHeatFluxFvPatchScalarField::fixedHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    heatFlux_(p.size(),0.0)
{}


Foam::fixedHeatFluxFvPatchScalarField::fixedHeatFluxFvPatchScalarField
(
    const fixedHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    heatFlux_(ptf.heatFlux_)
{}


Foam::fixedHeatFluxFvPatchScalarField::fixedHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    heatFlux_("heatFlux",dict,p.size())
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::fixedHeatFluxFvPatchScalarField::fixedHeatFluxFvPatchScalarField
(
    const fixedHeatFluxFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    heatFlux_(wbppsf.heatFlux_)
{}


Foam::fixedHeatFluxFvPatchScalarField::fixedHeatFluxFvPatchScalarField
(
    const fixedHeatFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    heatFlux_(wbppsf.heatFlux_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

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

  const fvPatchField<scalar>& K = this->patch().lookupPatchField
      (
       "K",
       reinterpret_cast<const volScalarField*>(NULL),
       reinterpret_cast<const scalar*>(NULL)
       );


    const IOdictionary& thermProps =
      this->db().objectRegistry::lookupObject<IOdictionary>("thermophysicalProperties");


    const dimensionedScalar cpGas("cpGas", thermProps);
    const dimensionedScalar cpV("cpV", thermProps);
    const dimensionedScalar cpL("cpL", thermProps);

    const scalarField cpMix
    (
        (1-x)*(1-y)*cpL.value()
      + x*(1-y)*cpV.value() + y*cpGas.value()
    );

    scalarField alphaMix(K/cpMix);


    if (x.db().foundObject<volScalarField>("alphat"))
    {
      const fvPatchField<scalar>& alphat = this->patch().lookupPatchField
      (
       "alphat",
       reinterpret_cast<const volScalarField*>(NULL),
       reinterpret_cast<const scalar*>(NULL)
       );

      alphaMix += alphat;
    }

    gradient() = heatFlux_ / alphaMix;  // assumes non-zero alphaMix

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::fixedHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    this->heatFlux_.writeEntry("heatFlux",os);

    gradient().writeEntry("gradient", os);

    this->writeEntry("value", os);  // this is just here to make paraFoam happy
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedHeatFluxFvPatchScalarField
    );
}

// ************************************************************************* //
