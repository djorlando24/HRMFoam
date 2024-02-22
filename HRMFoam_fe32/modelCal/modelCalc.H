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
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                Class modelCalc Declaration
\*---------------------------------------------------------------------------*/

class modelCalc
:     
    public basicThermo

{
      // Private Data

      const  fvMesh& mesh_;          
      dimensionedScalar theta_0_;
      dimensionedScalar pNonDimExp_;
      bool pcritical_;
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
      autoPtr<thermoBase> thermPoint1;

    // Private Member Functions

    //- Calculates properties and returns model
  const volScalarField getModel( const  volScalarField& p ,const  volScalarField& H, 
				 const volScalarField &x, const volScalarField &y,
				 const volScalarField& rho ); 

  void calcRhoBar(const volScalarField& x, const volScalarField& y, const volScalarField& h);

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
    ~modelCalc();

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
	 volScalarField pMin_ = rhoMin_*(h-hDatum_)*Rgas_/cpGas_;
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


//This returns a dummy variable with density dimensions, since no turbulence model
//calls the rho() function from the thermoPhysical model       
      virtual tmp<volScalarField> rho() const;

      void correct();

//Setter function alpha_. This field needs to be updated in order for dynamic
//mesh motion to work, otherwise topology change chokes on MapGeometricFields.
  void set_alpha(const volScalarField& alphaMix) { alpha_ = alphaMix; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
