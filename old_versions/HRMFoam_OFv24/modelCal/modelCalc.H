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
    Foam::modelCalc

Description
    Evaluates the Homogenous Relaxation Model and updates properties from the
    property readers

SourceFiles
    modelCalc.C

\*---------------------------------------------------------------------------*/

#ifndef modelCalc_H
#define modelCalc_H

#include "volFields.H"

#include "./propReader/thermoBase.H"
#include "./propReader/refprop/refprop.H"

#include "./propReader/utrc/utrc.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                Class modelCalc Declaration
\*---------------------------------------------------------------------------*/

class modelCalc
:
    public psiThermo

{
      // Private Data

      const  fvMesh& mesh_;
      dimensionedScalar theta_0_;
      dimensionedScalar pNonDimExp_;
      bool compressible_;
      dimensionedScalar theta_floor_;
      dimensionedScalar alphaFracExp_;
      const Time& runTime_;


      dimensionedScalar Rgas_;
      dimensionedScalar cpGas_;
      dimensionedScalar hDatum_;
      dimensionedScalar psiv_;
      dimensionedScalar psil_;
      dimensionedScalar Pcrit;
      dimensionedScalar rhoMin_;
      const dimensionedScalar rhoL_cap;
      dimensionedScalar  muNCG_;  // viscosity of non-condensable gas

      volScalarField  rhog_;       // a field (works for compressible or incomp.)
      volScalarField  rhoV_;
      volScalarField  rhoL_;
      volScalarField  xbar_;
      volScalarField  pSat_;
      volScalarField  K_;
      volScalarField  alphaFrac_;
      volScalarField  omega_HRM_;
      volScalarField  theta_;
      volScalarField  DxDt_;
      volScalarField  dMdp_;
      volScalarField  M_;
      volScalarField  M_old_;
      volScalarField  rhoBar_;  // density at current p,h but at x=xbar
      volScalarField  rhoTemp_; // Dummy density required for the bThermo model

      volScalarField  Liq_Mol_a_;
      volScalarField  Vap_Mol_a_;
  
      autoPtr<thermoBase> thermPoint1;

    // Private Member Functions

    //- Calculates properties and returns model
    tmp<volScalarField> getModel
    (
        const volScalarField& p,
        const  volScalarField& H,
				const volScalarField& x,
        const volScalarField& y,
			  const volScalarField& rho
    );

    void calcRhoBar
    (
        const volScalarField& x,
        const volScalarField& y,
        const volScalarField& h
    );


public:


    // Constructors

    modelCalc
    (
      dimensionedScalar  theta_0,
      dimensionedScalar  pNonDimExp,
      bool  pcritical,
      bool  compressible,
      dimensionedScalar  theta_floor,
      dimensionedScalar  alphaFracExp,
      const volScalarField&  h,
      volScalarField&  rho,
      volScalarField&  x,
      const volScalarField&  y,
      const volVectorField&  U,
      int conditionalWrite
    );


    // Destructor
    virtual ~modelCalc() = default;

    // Member Functions

       //update all properties and dMdp
  void update(const volScalarField& H, const volScalarField& rho, const volScalarField& x,
	      const volScalarField& y, const volScalarField& dp )
       {

	 M_old_=getModel((p_+dp),H,x,y,rho);
	 M_ = getModel(p_,H,x,y,rho);
	 //Model derivative with respect to pressure
	 dMdp_= (M_old_ - M_)/dp;
       }

       //update model term and properties
  void update_props(const volScalarField& H, const volScalarField& x, const volScalarField& y,
		    const volScalarField& rho)
       {
	 M_ = getModel(p_,H,x,y,rho);
       }

  void update_psi(const volScalarField& rhoEOS, const volScalarField& y, const volScalarField& x, const volScalarField& h)
  {
	 // calculate minimum pressure field based upon rhoMin_ and the ideal gas law
	 volScalarField pMin_(rhoMin_*(h-hDatum_)*Rgas_/cpGas_);
         // update compressibility
         psi_ = rhoEOS*(y/Foam::max(p_,pMin_) + (1-y)*(x*psiv_/(Foam::max(rhoV_,rhoMin_)) + (1-x)*psil_/rhoL_));
         // volScalarField gasPsi = 1/Rgas_/((h-hDatum_)/cpGas_);
         // I suspect this expression, below, is more accurate.  But it is less stable.
         // psi_ = sqr(rho)*( ( x * psiv_/sqr(rhoV_) + (1-x) * psil_/sqr(rhoL_))*(1-y) + y*gasPsi/sqr(rhoNCgas));
	 // try this
	 //psi_ = sqr(rho)*( ( x * psiv_/sqr(rhoV_) + (1-x) * psil_/sqr(rhoL_))*(1-y) + (y/Foam::max(p_,pMin_))/rhog_);
  }

       const volScalarField& rhoL() const
       {
            return rhoL_;
       }

       const volScalarField& rhoV() const
       {
            return rhoV_;
       }



       const volScalarField& xbar() const
       {
            return xbar_;
       }


       const volScalarField& pSat() const
       {
           return pSat_;
       }


       const volScalarField& MTerm() const
       {
           return M_;
       }

       const volScalarField& dMdp() const
       {
           return dMdp_;
       }

      //thermal conductivity
       const volScalarField& K() const
       {
           return K_;
       }

      //Void Fraction
       const volScalarField& alphaFrac() const
       {
           return alphaFrac_;
       }


       const volScalarField& omega_HRM() const
       {
           return omega_HRM_;
       }

      const volScalarField& theta() const
       {
           return theta_;
       }

      const volScalarField& DxDt() const
       {
           return DxDt_;
       }

      const volScalarField& rhog() const
      {
	return rhog_;
      }

  const dimensionedScalar psil() const
  {
    return psil_;
  }

  const dimensionedScalar psiv() const
  {
    return psiv_;
  }

  const dimensionedScalar Rgas() const
  {
    return Rgas_;
  }

  const dimensionedScalar hDatum() const
  {
    return hDatum_;
  }

  const dimensionedScalar cpGas() const
  {
    return cpGas_;
  }

  const volScalarField& rhoBar() const
  {
    return rhoBar_;
  }


  const volScalarField& Liq_Mol_a() const
  {
    return Liq_Mol_a_;  
  }

  const volScalarField& Vap_Mol_a() const
  {
    return Vap_Mol_a_;  
  }


//This returns a dummy variable with density dimensions, since no turbulence model
//calls the rho() function from the thermoPhysical model
      virtual tmp<volScalarField> rho() const;

      void correct();

//Setter function alpha_. This field needs to be updated in order for dynamic
//mesh motion to work, otherwise topology change chokes on MapGeometricFields.
  void set_alpha(const volScalarField& alphaMix) { alpha_ = alphaMix; }


    // Compatibility

    virtual word thermoName() const
    {
        NotImplemented;
        return word::null;
    }

    virtual bool incompressible() const
    {
        return false;
    }

    virtual bool isochoric() const
    {
        return false;
    }

    virtual volScalarField& he()
    {
        NotImplemented;
        return const_cast<volScalarField&>(volScalarField::null());
    }

    virtual const volScalarField& he() const
    {
        NotImplemented;
        return volScalarField::null();
    }

    virtual tmp<volScalarField> he
    (
        const volScalarField& p,
        const volScalarField& T
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> he
    (
        const scalarField& p,
        const scalarField& T,
        const labelList& cells
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> he
    (
        const scalarField& p,
        const scalarField& T,
        label patchi
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> hc() const
    {
        NotImplemented;
        return  nullptr;
    }

    virtual tmp<scalarField> THE
    (
        const scalarField& h,
        const scalarField& p,
        const scalarField& T0,
        const labelList& cells
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> THE
    (
        const scalarField& h,
        const scalarField& p,
        const scalarField& T0,
        const label patchi
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> Cp() const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> Cp
    (
        const scalarField& p,
        const scalarField& T,
        const label patchi
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    //- Heat capacity using pressure and temperature [J/kg/K]
    virtual tmp<scalarField> Cp
    (
        const scalarField& p,
        const scalarField& T,
        const labelList& cells
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> Cv() const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> Cv
    (
        const scalarField& p,
        const scalarField& T,
        const label patchi
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    //- Density from pressure and temperature from EoS
    virtual tmp<scalarField> rhoEoS
    (
        const scalarField& p,
        const scalarField& T,
        const labelList& cells
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> gamma() const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> gamma
    (
        const scalarField& p,
        const scalarField& T,
        const label patchi
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> Cpv() const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> Cpv
    (
        const scalarField& p,
        const scalarField& T,
        const label patchi
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> CpByCv() const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> CpByCv
    (
        const scalarField& p,
        const scalarField& T,
        const label patchi
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> CpByCpv() const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> CpByCpv
    (
        const scalarField& p,
        const scalarField& T,
        const label patchi
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> W() const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> kappa() const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> kappa(const label patchi) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> alphahe() const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> alphahe(const label patchi) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> kappaEff(const volScalarField& alphat) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> kappaEff
    (
        const scalarField& alphat,
        const label patchi
    ) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<volScalarField> alphaEff(const volScalarField& alphat) const
    {
        NotImplemented;
        return nullptr;
    }

    virtual tmp<scalarField> alphaEff
    (
        const scalarField& alphat,
        const label patchi
    ) const
    {
        NotImplemented;
        return nullptr;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
