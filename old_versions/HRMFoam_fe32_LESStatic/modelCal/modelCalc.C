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

\*---------------------------------------------------------------------*/

#include "modelCalc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

modelCalc::modelCalc
(
    dimensionedScalar  theta_0,
    dimensionedScalar  pNonDimExp,
    bool  pcritical,
    bool  compressible,
    dimensionedScalar  theta_floor,
    dimensionedScalar  alphaFracExp, 
    const volScalarField& h, 
    volScalarField&  rho,
    volScalarField& x,
    const  volScalarField& y,
    const  volVectorField& U,
    int conditionalWrite
)
:
    basicThermo(U.mesh(),U.mesh()),
    mesh_(U.mesh()),
    theta_0_(theta_0),
    pNonDimExp_(pNonDimExp),
    pcritical_(pcritical),
    compressible_(compressible),
    theta_floor_(theta_floor),
    alphaFracExp_(alphaFracExp),
    runTime_(U.time()),
    Rgas_    ("Rgas_",    dimensionSet(0,  2, -2, -1, 0, 0, 0) ,0.0),
    cpGas_   ("cpGas_",   dimensionSet(0,  2, -2, -1, 0, 0, 0) ,0.0),
    hDatum_  ("hDatum_",  dimensionSet(0,  2, -2,  0, 0, 0, 0) ,0.0),
    psiv_    ("psiv_",    dimensionSet(0, -2,  2,  0, 0, 0, 0) ,0.0),
    psil_    ("psil_",    dimensionSet(0, -2,  2,  0, 0, 0, 0) ,0.0),
    Pcrit    ("Pcrit",    dimensionSet(1, -1, -2,  0, 0, 0, 0) ,0.0),
    rhoMin_  ("rhoMin_",  dimensionSet(1, -3,  0,  0, 0, 0, 0) ,0.0),
    rhoL_cap ("rhoL_cap", dimensionSet(1, -3,  0,  0, 0, 0, 0) ,SMALL),
    muNCG_   ("muNCG",    dimensionSet(1, -1, -1,  0, 0, 0, 0) ,0.0),
    rhog_
    (
        IOobject
        (
            "rhog",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0, 0, 0) ,1.0)
    ),

    rhoV_
    (
        IOobject
        (
            "rhoV",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0, 0, 0) ,0.0)
    ),

    rhoL_
    (
        IOobject
        (
            "rhoL",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0, 0, 0) ,0.0)
    ),


    xbar_
    (
        IOobject
        (
            "xbar",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
        ),
       mesh_,
       dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0.0)
    ),

    pSat_
    (
        IOobject
        (
            "pSat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(1, -1, -2,0,0,0,0),0.0)
    ),



    K_
    (
        IOobject
        (
            "K",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::writeOption(conditionalWrite)
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(1, 1, -3, -1, 0, 0, 0) ,0.0)
    ),

    alphaFrac_
    (
        IOobject
        (
            "alphaFrac",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0.0)
     ), 

    omega_HRM_
    (
        IOobject
        (
            "omega_HRM",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0, 0, 0) ,0.0)
    ),  



    theta_
    (
        IOobject
        (
            "theta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0, 0, 0) ,0.0)
    ),

    DxDt_
    (
        IOobject
        (
            "DxDt",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0) ,0.0)
    ),


    dMdp_
    (
        IOobject
        (
            "dMdp",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(-1, 1, 1, 0, 0, 0, 0) ,0.0)
    ),

    M_
    (
      IOobject
      (
            "M",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
      ),
      mesh_,
      dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0) ,0.0)
    ),

    M_old_
    (
      IOobject
      (
            "M_old",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
      ),
      mesh_,
      dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0) ,0.0)
    ),

    rhoBar_
    (
      IOobject
      (
            "rhoBar",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
      ),
      mesh_,
      dimensionedScalar("zero", dimDensity ,0.0)
    ),

//this is a dummy density required for the bThermo class
    rhoTemp_
   (
       IOobject
        (
            "rhoBThermo",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
       dimensionedScalar("one", dimensionSet(1, -3, 0, 0, 0), 1.0) // will be over-written
    ),

    thermPoint1
    (
       thermoBase::New(
                        U ,
                        IOobject
                        (
                            "tPoint",
                             U.time().constant(),
                             U.db(),
                             IOobject::NO_READ,
                             IOobject::NO_WRITE,
                             true
                         )
                      )
     )
{

   if(pcritical)
        Pcrit.value() = thermPoint1->Pcrit();
   

   // read requisite thermo properties for compressible or incompressible gas

   const IOdictionary& dict =  db().lookupObject<IOdictionary>("thermophysicalProperties");

   rhoMin_ = dimensionedScalar (dict.lookup("rhoMin"));
   cpGas_ =  dimensionedScalar (dict.lookup("cpGas"));
   hDatum_ = dimensionedScalar (dict.lookup("hDatum"));

   if (dict.found("muNCG"))
   {
	     muNCG_ = dimensionedScalar(dict.lookup("muNCG"));
   }
   else
   {
	     muNCG_ = dimensionedScalar("muNCG",dimensionSet(1,-1,-1,0,0,0,0),1.813e-5);
	     Info<<endl<<"WARNING: The non-condensable gas viscosity is not specified in thermophysicalProperties."<<endl;
	     Info<<"Any non-condensable gas will be assumed to have the viscosity of dry air at NTP."<<endl;
  	     Info<<"To specify non-condensable gas viscosity insert the following in thermophysicalProperties:"<<endl;
  	     Info<<"muNCG   muNCG  [1 -1 -1 0 0 0 0] 1.813e-5;"<<endl<<endl;
   }
   if (compressible_)
     {

       Rgas_  = dimensionedScalar (dict.lookup("Rgas"));

       dimensionedScalar aL_  (dict.lookup("aL"));
       dimensionedScalar aV_  (dict.lookup("aV"));

       psiv_ = 1/sqr(aV_);
       psil_ = 1/sqr(aL_);

       // calculate  density
       rhog_ = Foam::max(p_/((h-hDatum_)/cpGas_)/Rgas_,rhoMin_);

     }
   else
     {
       rhog_ = dimensionedScalar(dict.lookup("rhoNCgas"));
     }

   calcRhoBar(x,y,h);

   //these feature overwrites the inital density and quality fields to their equilibrium values
   //It only works starting at time 0. This cannot be done at restarts.
   //If features are not desires, set rhoOverwrite = false and/or xOverwrite in thermophysicalProperties
   Switch rhoOverwrite("yes");
   rhoOverwrite.readIfPresent("rhoOverwrite", dict);
   if (U.time().timeName() == "0" && rhoOverwrite == 1)
   {       
      Info<< 
	"WARNING: Over-riding initial condition and changing the density to the equilibrium value. \n\n";
      rho = rhoBar_;
      rho.correctBoundaryConditions();
   }

   getModel(p_,h,x,y,rho);

   Switch xOverwrite("yes");
   xOverwrite.readIfPresent("xOverwrite", dict);
   if (U.time().timeName() == "0" && xOverwrite == 1)
   {       
      Info<< 
	"WARNING: Over-riding initial condition and changing the quality to the equilibrium value. \n\n";
      x = xbar_;
      x.correctBoundaryConditions();
   }
}
   

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

modelCalc::~modelCalc()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
  const volScalarField modelCalc::getModel(const volScalarField& pres,  const volScalarField& h, 
					   const volScalarField& x, const volScalarField &y,
					   const volScalarField& rho) 
{

  forAll(mesh_.C(), celli)
  {
    
       

    //finally call getProperties with the total, not partial pressure
    thermPoint1->getProperties(pres[celli], h[celli],x[celli]);
      
   
     //Now update the properties 

     xbar_[celli]=thermPoint1->xbar();
     pSat_[celli]=thermPoint1->pSat();
     T_[celli]= thermPoint1->temp()*(1-y[celli]) 
               + (h[celli]-hDatum_.value())/cpGas_.value()*y[celli];
     K_[celli] = thermPoint1->K();
  
     rhoV_[celli]=thermPoint1->rhov();
     rhoL_[celli]=thermPoint1->rhol();
     mu_[celli]=(max(min(y[celli]*rho[celli]/rhog_[celli],1),0))*muNCG_.value()+(1-(max(min(y[celli]*rho[celli]/rhog_[celli],1),0)))*thermPoint1->mu();

     if (compressible_)
       {
	 if(pSat_[celli]>pres[celli]) //Extrapolate with specific volume to avoid negative rhoV at low pressures
	 { 
		scalar vSat= 1/rhoV_[celli];
		rhoV_[celli]= 1/(vSat - sqr(vSat)*psiv_.value()*(pres[celli]-pSat_[celli]));
	 }
	 else rhoV_[celli] =rhoV_[celli] + psiv_.value()*(pres[celli]-pSat_[celli]); //Otherwise extrapolate using density

	 rhoL_[celli] = rhoL_[celli] + psil_.value()*(pres[celli]-pSat_[celli]);
	 rhoV_[celli] = max(rhoV_[celli],rhoMin_.value());
	
       }

     alphaFrac_[celli]=max
                   (
                     min
                     (
		      x[celli]*rhoL_[celli]/(rhoV_[celli]+x[celli]*(rhoL_[celli]-rhoV_[celli])),
                       1.0
                     ),
                     0.0
                   ); 


     if (rhoV_[celli] == 0)
     {
         omega_HRM_[celli]=0;
     }
     else
       {  
	 omega_HRM_[celli]=(y[celli]-1)*
	   (1/rhoV_[celli]-1/rhoL_[celli])*rho[celli]*rho[celli];
     }
   }

   //access boundary elements
   //const fvPatchList& patches = mesh.boundary(); 
   forAll(mesh_.boundary(), patchI)
   {
      const fvPatch& cPatch = mesh_.boundary()[patchI];
  
      // cache constant fields
      const fvPatchScalarField& patchP = pres.boundaryField()[patchI];
      const fvPatchScalarField& patchH = h.boundaryField()[patchI];
      const fvPatchScalarField& patchRho = rho.boundaryField()[patchI];
      const fvPatchScalarField& patchx = x.boundaryField()[patchI];
      const fvPatchScalarField& patchy = y.boundaryField()[patchI];
      const fvPatchScalarField& patchrhog = rhog_.boundaryField()[patchI];

  
      //cache fields to update
      fvPatchScalarField& patchrhoV = rhoV_.boundaryField()[patchI];
      fvPatchScalarField& patchrhoL = rhoL_.boundaryField()[patchI];
      //      fvPatchScalarField& patchpsi = psi_.boundaryField()[patchI];
      fvPatchScalarField& patchmu = mu_.boundaryField()[patchI];
      fvPatchScalarField& patchxbar = xbar_.boundaryField()[patchI];
      fvPatchScalarField& patchpSat = pSat_.boundaryField()[patchI];
      fvPatchScalarField& patchT = T_.boundaryField()[patchI];
      fvPatchScalarField& patchK = K_.boundaryField()[patchI];
      fvPatchScalarField& patchalphaFrac = alphaFrac_.boundaryField()[patchI];
      fvPatchScalarField& patchomega_HRM = omega_HRM_.boundaryField()[patchI];
  
        
  
      forAll(cPatch, faceI)
      {

        patchmu[faceI]=(max(min(patchy[faceI]*patchRho[faceI]/patchrhog[faceI],1),0))*muNCG_.value()
          +(1-(max(min(patchy[faceI]*patchRho[faceI]/patchrhog[faceI],1),0)))*thermPoint1->mu();
            
	// get properties with the total, not partial pressure
	thermPoint1->getProperties(patchP[faceI], patchH[faceI],patchx[faceI]);

	patchxbar[faceI]=thermPoint1->xbar();

        patchpSat[faceI]=thermPoint1->pSat();

        patchT[faceI]=thermPoint1->temp()*(1-patchy[faceI])
	     + (patchH[faceI]-hDatum_.value())/cpGas_.value()*patchy[faceI];

        patchK[faceI]=thermPoint1->K();
 
        patchrhoV[faceI]=thermPoint1->rhov();

        patchrhoL[faceI]=thermPoint1->rhol();

        if (compressible_)
        {
          if(patchpSat[faceI]>patchP[faceI]) //Extrapolate with specific volume to avoid negative rhoV at low pressures
	  { 
		scalar vSat= 1/patchrhoV[faceI];
		patchrhoV[faceI]= 1/(vSat - sqr(vSat)*psiv_.value()*(patchP[faceI]-patchpSat[faceI]));
	  }
	  else patchrhoV[faceI] = patchrhoV[faceI] + psiv_.value()*(patchP[faceI]-patchpSat[faceI]); //Otherwise extrapolate using density

          patchrhoL[faceI] = patchrhoL[faceI] + psil_.value()*(patchP[faceI]-patchpSat[faceI]);
	  patchrhoV[faceI] = max( patchrhoV[faceI], rhoMin_.value() );
        }

        patchalphaFrac[faceI]=max(min((patchx[faceI]*patchrhoL[faceI])
          /(patchrhoV[faceI] + patchx[faceI]*(patchrhoL[faceI]-patchrhoV[faceI]) ),1.0),0.0);  


        if (patchrhoV[faceI] == 0)  //Never triggered?  GLJ
        {
          patchomega_HRM[faceI]=0;
        }
        else
  	{
          patchomega_HRM[faceI]= (patchy[faceI]-1)*	            
            (1/patchrhoV[faceI]-1/patchrhoL[faceI]) *  patchRho[faceI]*patchRho[faceI];                                  
        }
      }
   }

   if (compressible_)  
       rhog_ = Foam::max(pres/((h-hDatum_)/cpGas_)/Rgas_,rhoMin_);
   // no else--incompressible density never needs updating
     
   //non dimensional pressure
   volScalarField  pNonDim_= max
                              (
                               mag((pSat_-pres)/(Pcrit-pSat_)),
                               SMALL
                              );  

   const dimensionedScalar rhoV_floor
                            (
                              "rhoV_floor",
                              dimensionSet(1, -3 ,0, 0,0,0,0),
                              SMALL
                            ); 

   calcRhoBar(x,y,h);

   //D-Z correlation
   theta_= theta_0_
         *pow(max(alphaFrac_,SMALL),alphaFracExp_)
         *pow(pNonDim_,pNonDimExp_);
   //for future considerion,  a correction for large alphaFrac:    *pow(max(1.0-alphaFrac,SMALL),alphaFracExp);

   DxDt_=(xbar_-x)/max(theta_,theta_floor_);  

   return (omega_HRM_*DxDt_/rho) ;
   //    return ( (1-y) * (rho- rhoBar_)/max(theta_,theta_floor_)/rho );
}



//This is a dummy variable with density dimensions. this is needded because this
//function is a pure virtual in basicThermo class
//It is a dummy because none of the turbulence models call this function 
Foam::tmp<Foam::volScalarField> Foam::modelCalc::rho() const
{

   return rhoTemp_;
}


void Foam::modelCalc::correct()
{
//no need to update properties here since they will be update in 
//the main program
}

// calculate equilibrium density, assuming rhog_ is already updated
void Foam::modelCalc::calcRhoBar(const volScalarField& x, const volScalarField& y, const volScalarField& h)
{
  

  forAll(mesh_.C(), celli) 
    {

      //finally call get properties with the total, not partial pressure
      thermPoint1->getProperties(p_[celli],h[celli], x[celli]);
      scalar rhoVCalc = thermPoint1->rhov(); //Uncorrected saturation densities
      scalar rhoLCalc = thermPoint1->rhol();      

      if(compressible_) //Correct saturation densities for compressibility
	{   
	  
	  
	 if(pSat_[celli]>p_[celli]) //Extrapolate with specific volume to avoid negative rhoV at low pressures
	 { 
		scalar vSat = 1/rhoVCalc;
		rhoVCalc = 1/(vSat - sqr(vSat)*psiv_.value()*(p_[celli]-pSat_[celli]));
	 }
	 else rhoVCalc += psiv_.value()*(p_[celli]-pSat_[celli]); //Otherwise extrapolate using density

	 rhoLCalc += psil_.value()*(p_[celli]-pSat_[celli]);
	 rhoVCalc = max(rhoVCalc,rhoMin_.value());
	}
	  rhoBar_[celli] = 1.0/(
	       (  
		( thermPoint1->xbar()/max(rhoVCalc,SMALL) + (1-thermPoint1->xbar())/max(rhoLCalc,SMALL) )
					*  (1-y[celli])
					+ y[celli]/rhog_[celli]
	       )
				);
            // DPS--I think this is OK for when we are above the fluid.dat
            // since the pSat value will top out at the same place as the 
            // rho value.

	  rhoBar_[celli] = Foam::max(rhoBar_[celli],rhoMin_.value());

	

//Patches	   
    }
    forAll(mesh_.boundary(), patchI)
  	 {
    	  const fvPatch& cPatch = mesh_.boundary()[patchI];
	  const fvPatchScalarField& patchP = p_.boundaryField()[patchI];
	  const fvPatchScalarField& patchPSat = pSat_.boundaryField()[patchI];
          const fvPatchScalarField& patchY = y.boundaryField()[patchI];
	  const fvPatchScalarField& patchH = h.boundaryField()[patchI];
	  const fvPatchScalarField& patchX = x.boundaryField()[patchI];
          const fvPatchScalarField& patchRhoG = rhog_.boundaryField()[patchI];

          fvPatchScalarField& patchrhoBar = rhoBar_.boundaryField()[patchI];

          forAll(cPatch, faceI)	
          {
	    thermPoint1->getProperties(patchP[faceI],patchH[faceI],patchX[faceI]);
	    scalar rhoVCalc = thermPoint1->rhov(); //Uncorrected saturation densities
	    scalar rhoLCalc = thermPoint1->rhol(); 

	    if(compressible_) //Correct saturation densities for compressibility
	    {   
	  
	      if(patchPSat[faceI]>patchP[faceI]) //Extrapolate with specific volume to avoid negative rhoV at low pressures
	      { 
		scalar vSat = 1/rhoVCalc;
		rhoVCalc = 1/(vSat - sqr(vSat)*psiv_.value()*(patchP[faceI]-patchPSat[faceI]));
	      }
	      else rhoVCalc += psiv_.value()*(patchP[faceI]-patchPSat[faceI]); //Otherwise extrapolate using density

	      rhoLCalc += psil_.value()*(patchP[faceI]-patchPSat[faceI]);
	      rhoVCalc = max(rhoVCalc,rhoMin_.value());
	    }

            patchrhoBar[faceI]= 1.0/(
				     (  ( thermPoint1->xbar()/max(rhoVCalc,SMALL) 
				      + (1-thermPoint1->xbar())/max(rhoLCalc,SMALL) )
					*  (1-patchY[faceI])
					+ (patchY[faceI]/patchRhoG[faceI])
					)
				     );
	    patchrhoBar[faceI] = max( patchrhoBar[faceI], rhoMin_.value() );
 	  }
       }  

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam



// ************************************************************************* //
