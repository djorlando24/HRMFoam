#include<iostream>
#include<fstream>
#include<cmath>
#include"utrc.H"
#include"thermTran.H"
#include"addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(utrc,0);
addToRunTimeSelectionTable(thermoBase, utrc, dictionary);

utrc::utrc(const volVectorField& U, const IOobject& thermPoint )
:
  thermoBase(U, thermPoint )
{
  getArrayInfo();
  getWorkingArray();
}




//void utrc::getPropertiesUTRC(const double &P, const double &H, const double &X)
void utrc::getProperties(const double& P, const double& H, const double& X ) const
{


       opflag_ = 1; // get calculation
              

       Pinout_ = P/100000.00;
//if pressure goes below the lower bound of the ptessure table equate it to the lower bound
       if(Pinout_ < pow(10.00,P0_) )
          Pinout_ = pow(10.00,P0_);

       Hin_    = H*wtm_/1000;
       xin_    = X;

       interpCppWrapper(
       &opflag_,   //INPUT: integer flag -1:get array size and index strides; 0:initialize; 1:get thermo properties;
       &Pinout_,   //INOUTPUT: double pressure in [bar] (absolute) from 0.1 to 20 as input, , when opflag.eq.2, it is output
       &Hin_,      //INPUT:    double enthalpy in  [J/mol] from -1d5 to -5d4
       &xin_,      //INPUT:    double quality (mass fraction of vapour) [non-D] from 0-0.1
       &wtm_,      //OUTPUT:   double overall molecular weight [g/mol] set when opflag.eq.0
       &rhoV_,     //OUTPUT:   double density of vap [kg/m^3]
       &rhoL_,     //OUTPUT:   double density of liq [kg/m^3]
       &cpV_,      //OUTPUT:   double Cp of vap [J/kg]
       &cpL_,      //OUTPUT:   double Cp of liq [J/kg]
       &muV_,      //OUTPUT:   double viscosity of vap [Pa-S]
       &muL_,      //OUTPUT:   double viscosity of liq [Pa-S]
       &lambdaV_,  //OUTPUT:   double thermo conductivity of vap [W/m-K]
       &lambdaL_,  //OUTPUT:   double thermo conductivity of liq [W/m-K]
       &aV_,       //OUTPUT:   double acoustic speed of vap [m/S]
       &aL_,       //OUTPUT:   double acoustic speed of liq [m/S]
       &TV_,       //OUTPUT: double temp of vap [K]
       &TL_,       //OUTPUT: double temp of liq [K]
       &Teq_,      //OUTPUT: double temp of system in equilibrium (P can differ) [K]
       &SurfTen_,  // OUTPUT: double Surface Tension [N/m]
       &superheat_,   //OUTPUT:   integer 1: superheated, 0: otherwise
       &supercool_,   //OUTPUT:   integer 1: supercooled, 0: otherwise
       &xbar_,     //OUTPUT:   double [non-D]
       &Pbub_,
       &sizep_,       //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       &sizeh_,       //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       &sizex_,       //INOUTPUT: array size in the x dimension, output when opflag=0, else input
       &sizebub_,     //INOUTPUT: array size of bubble point curve data
       &sizexbarp_,   //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       &sizexbarh_,   //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       &dP_,       //INOUTPUT: delta P in thermo-transport data array, output when opflag=0, else input
       &dH_,       //INOUTPUT: delta H in thermo-transport data array, output when opflag=0, else input
       &dx_,       //INOUTPUT: delta X in thermo-transport data array, output when opflag=0, else input
       &dbubH_,    //INOUTPUT: delta H in bubble point data array,     output when opflag=0, else input
       &dxbarP_,   //INOUTPUT: delta P in xbar data array,             output when opflag=0, else input
       &dxbarH_,   //INOUTPUT: delta h in xbar data array,             output when opflag=0, else input
       rhoVD_,    //INOUTPUT: double array (sizep,sizeh,sizex) e.g. (100*60*20) array size passed down
       rhoLD_,    //INOUTPUT: double array (sizep,sizeh,sizex)
       cpVD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       cpLD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       muVD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       muLD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       lambdaVD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       lambdaLD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       aVD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       aLD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       TVD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       TLD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       TeqD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       SurfTenD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       superheatD_,  //INOUTPUT: double array (sizep,sizeh,sizex)
       PbubD_,    //INOUTPUT: double array (sizebub    ) e.g. (100)
       HbubD_,    //INOUTPUT: double array (sizebub    ) e.g. (100)
       xbarD_,    //INOUTPUT: double array (sizexbarp,sizexbarh) e.g. (100*60) array size passed down
       &H0_,       //INOUTPUT: double lower limit of H 
       &H1_,       //INOUTPUT: double upper limit of H
       &P0_,       //INOUTPUT: double lower limit of P
       &P1_,       //INOUTPUT: double upper limit of P
       &x0_,       //INOUTPUT: double lower limit of x
       &x1_,       //INOUTPUT: double upper limit of x
       &ierr_);    //OUTPUT: error flag: 0 normal; -1 opflag out of range(OOR) ;-2 P OOR;-3: H OOR;-4 x OOR;-5: nonuniform interpolation pts
       
       //  conductivity = ((1 - xin)*lambdaL) + (xin*lambdaV);
/*       prop[0]=wtm;      //OUTPUT:   double overall molecular weight [g/mol] set when opflag.eq.0
       prop[1]=rhoV;     //OUTPUT:   double density of vap [kg/m^3]
       prop[2]=rhoL;     //OUTPUT:   double density of liq [kg/m^3]
       prop[3]=cpV;      //OUTPUT:   double Cp of vap [J/kg]
       prop[4]=cpL;      //OUTPUT:   double Cp of liq [J/kg]
       prop[5]=muV;      //OUTPUT:   double viscosity of vap [Pa-S]
       prop[6]=muL;      //OUTPUT:   double viscosity of liq [Pa-S]
       prop[7]=lambdaV;  //OUTPUT:   double thermo conductivity of vap [W/m-K]
       prop[8]=lambdaL;  //OUTPUT:   double thermo conductivity of liq [W/m-K]
       prop[9]=aV;       //OUTPUT:   double acoustic speed of vap [m/S]
       prop[10]=aL;       //OUTPUT:   double acoustic speed of liq [m/S]
       prop[11]=xbar;     //OUTPUT:   double [non-D]
       prop[12]=Pbub*1e5; //OUTPUT:  double bubble point pressure [Pa]
       prop[13]=TL;       //OUTPUT:  double Temperature [K]
*/
        mu_ = (1-X) * muL_ + X * muV_;
        lambda_ =  (1-X) * lambdaL_ + X* lambdaV_;
        psi_ = (1.0-X)/max((aL_*aL_),SMALL) +X/max((aV_*aV_),SMALL) ; 
       return;



}

void utrc::getArrayInfo()  const

{
     
       opflag_ = -1; // get one example calculation
      
       interpCppWrapper(
       &opflag_,   //INPUT: integer flag -1:get array size and index strides; 0:initialize; 1:get thermo properties;
       &Pinout_,   //INOUTPUT: double pressure in [bar] (absolute) from 0.1 to 20 as input, , when opflag.eq.2, it is output
       &Hin_,      //INPUT:    double enthalpy in  [J/mol] from -1d5 to -5d4
       &xin_,      //INPUT:    double quality (mass fraction of vapour) [non-D] from 0-0.1
       &wtm_,      //OUTPUT:   double overall molecular weight [g/mol] set when opflag.eq.0
       &rhoV_,     //OUTPUT:   double density of vap [kg/m^3]
       &rhoL_,     //OUTPUT:   double density of liq [kg/m^3]
       &cpV_,      //OUTPUT:   double Cp of vap [J/kg]
       &cpL_,      //OUTPUT:   double Cp of liq [J/kg]
       &muV_,      //OUTPUT:   double viscosity of vap [Pa-S]
       &muL_,      //OUTPUT:   double viscosity of liq [Pa-S]
       &lambdaV_,  //OUTPUT:   double thermo conductivity of vap [W/m-K]
       &lambdaL_,  //OUTPUT:   double thermo conductivity of liq [W/m-K]
       &aV_,       //OUTPUT:   double acoustic speed of vap [m/S]
       &aL_,       //OUTPUT:   double acoustic speed of liq [m/S]
       &TV_,       //OUTPUT: double temp of vap [K]
       &TL_,       //OUTPUT: double temp of liq [K]
       &Teq_,      //OUTPUT: double temp of system in equilibrium (P can differ) [K]
       &SurfTen_,  // OUTPUT: double Surface Tension [N/m]
       &superheat_,   //OUTPUT:   integer 1: superheated, 0: otherwise
       &supercool_,   //OUTPUT:   integer 1: supercooled, 0: otherwise
       &xbar_,     //OUTPUT:   double [non-D]
       &Pbub_,
       &sizep_,       //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       &sizeh_,       //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       &sizex_,       //INOUTPUT: array size in the x dimension, output when opflag=0, else input
       &sizebub_,     //INOUTPUT: array size of bubble point curve data
       &sizexbarp_,   //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       &sizexbarh_,   //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       &dP_,       //INOUTPUT: delta P in thermo-transport data array, output when opflag=0, else input
       &dH_,       //INOUTPUT: delta H in thermo-transport data array, output when opflag=0, else input
       &dx_,       //INOUTPUT: delta X in thermo-transport data array, output when opflag=0, else input
       &dbubH_,    //INOUTPUT: delta H in bubble point data array,     output when opflag=0, else input
       &dxbarP_,   //INOUTPUT: delta P in xbar data array,             output when opflag=0, else input
       &dxbarH_,   //INOUTPUT: delta h in xbar data array,             output when opflag=0, else input
       rhoVD_,    //INOUTPUT: double array (sizep,sizeh,sizex) e.g. (100*60*20) array size passed down
       rhoLD_,    //INOUTPUT: double array (sizep,sizeh,sizex)
       cpVD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       cpLD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       muVD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       muLD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       lambdaVD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       lambdaLD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       aVD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       aLD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       TVD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       TLD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       TeqD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       SurfTenD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       superheatD_,  //INOUTPUT: double array (sizep,sizeh,sizex)
       PbubD_,    //INOUTPUT: double array (sizebub    ) e.g. (100)
       HbubD_,    //INOUTPUT: double array (sizebub    ) e.g. (100)
       xbarD_,    //INOUTPUT: double array (sizexbarp,sizexbarh) e.g. (100*60) array size passed down
       &H0_,       //INOUTPUT: double lower limit of H 
       &H1_,       //INOUTPUT: double upper limit of H
       &P0_,       //INOUTPUT: double lower limit of P
       &P1_,       //INOUTPUT: double upper limit of P
       &x0_,       //INOUTPUT: double lower limit of x
       &x1_,       //INOUTPUT: double upper limit of x
       &ierr_);    //OUTPUT: error flag: 0 normal; -1 opflag out of range(OOR) ;-2 P OOR;-3: H OOR;-4 x OOR;-5: nonuniform interpolation pts

       std:: cout << "Array Sizes " << sizep_ <<" "<< sizeh_ <<" "<< sizex_ <<" "<< sizebub_ <<" "<< sizexbarp_ <<" "<< sizexbarh_;
       std:: cout << "Delta Sizes " << dP_ <<" "<< dH_ <<" "<< dx_ <<" "<< dbubH_ <<" "<< dxbarP_ <<" "<< dxbarH_;
      
     

       return;


}

void utrc::getWorkingArray()  const

{


       NP_ = sizep_;
       NH_ = sizeh_;
       NX_ = sizex_;
       NBUB_ = sizebub_;

       std::cout <<"Sizes,p,h,x,bub,xp,xh:"<<sizep_<<" "<<sizeh_<<" "<<sizex_<<" "<<sizebub_<<" "<<sizexbarp_<<" "<<sizexbarh_<<std::endl;
       rhoVD_       = new double [NP_*NH_*NX_];
       rhoLD_       = new double [NP_*NH_*NX_];
       cpVD_        = new double [NP_*NH_*NX_];
       cpLD_        = new double [NP_*NH_*NX_];
       muVD_        = new double [NP_*NH_*NX_];
       muLD_        = new double [NP_*NH_*NX_];
       lambdaVD_    = new double [NP_*NH_*NX_];
       lambdaLD_    = new double [NP_*NH_*NX_];
       aVD_         = new double [NP_*NH_*NX_];
       aLD_         = new double [NP_*NH_*NX_];
       TVD_         = new double [NP_*NH_*NX_];
       TLD_         = new double [NP_*NH_*NX_];
       TeqD_        = new double [NP_*NH_*NX_];
       SurfTenD_    = new double [NP_*NH_*NX_];
       superheatD_  = new int    [NP_*NH_*NX_];
       PbubD_       = new double [NBUB_    ];
       HbubD_       = new double [NBUB_    ];
       xbarD_       = new double [sizexbarp_*sizexbarh_   ];


       opflag_ = 0; // get working arrays
      
       
       interpCppWrapper(
       &opflag_,   //INPUT: integer flag -1:get array size and index strides; 0:initialize; 1:get thermo properties;
       &Pinout_,   //INOUTPUT: double pressure in [bar] (absolute) from 0.1 to 20 as input, , when opflag.eq.2, it is output
       &Hin_,      //INPUT:    double enthalpy in  [J/mol] from -1d5 to -5d4
       &xin_,      //INPUT:    double quality (mass fraction of vapour) [non-D] from 0-0.1
       &wtm_,      //OUTPUT:   double overall molecular weight [g/mol] set when opflag.eq.0
       &rhoV_,     //OUTPUT:   double density of vap [kg/m^3]
       &rhoL_,     //OUTPUT:   double density of liq [kg/m^3]
       &cpV_,      //OUTPUT:   double Cp of vap [J/kg]
       &cpL_,      //OUTPUT:   double Cp of liq [J/kg]
       &muV_,      //OUTPUT:   double viscosity of vap [Pa-S]
       &muL_,      //OUTPUT:   double viscosity of liq [Pa-S]
       &lambdaV_,  //OUTPUT:   double thermo conductivity of vap [W/m-K]
       &lambdaL_,  //OUTPUT:   double thermo conductivity of liq [W/m-K]
       &aV_,       //OUTPUT:   double acoustic speed of vap [m/S]
       &aL_,       //OUTPUT:   double acoustic speed of liq [m/S]
       &TV_,       //OUTPUT: double temp of vap [K]
       &TL_,       //OUTPUT: double temp of liq [K]
       &Teq_,      //OUTPUT: double temp of system in equilibrium (P can differ) [K]
       &SurfTen_,  // OUTPUT: double Surface Tension [N/m]
       &superheat_,   //OUTPUT:   integer 1: superheated, 0: otherwise
       &supercool_,   //OUTPUT:   integer 1: supercooled, 0: otherwise
       &xbar_,     //OUTPUT:   double [non-D]
       &Pbub_,
       &sizep_,       //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       &sizeh_,       //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       &sizex_,       //INOUTPUT: array size in the x dimension, output when opflag=0, else input
       &sizebub_,     //INOUTPUT: array size of bubble point curve data
       &sizexbarp_,   //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       &sizexbarh_,   //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       &dP_,       //INOUTPUT: delta P in thermo-transport data array, output when opflag=0, else input
       &dH_,       //INOUTPUT: delta H in thermo-transport data array, output when opflag=0, else input
       &dx_,       //INOUTPUT: delta X in thermo-transport data array, output when opflag=0, else input
       &dbubH_,    //INOUTPUT: delta H in bubble point data array,     output when opflag=0, else input
       &dxbarP_,   //INOUTPUT: delta P in xbar data array,             output when opflag=0, else input
       &dxbarH_,   //INOUTPUT: delta h in xbar data array,             output when opflag=0, else input
       rhoVD_,    //INOUTPUT: double array (sizep,sizeh,sizex) e.g. (100*60*20) array size passed down
       rhoLD_,    //INOUTPUT: double array (sizep,sizeh,sizex)
       cpVD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       cpLD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       muVD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       muLD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       lambdaVD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       lambdaLD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       aVD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       aLD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       TVD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       TLD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       TeqD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       SurfTenD_,      //INOUTPUT: double array (sizep,sizeh,sizex)  
       superheatD_,  //INOUTPUT: double array (sizep,sizeh,sizex)
       PbubD_,    //INOUTPUT: double array (sizebub    ) e.g. (100)
       HbubD_,    //INOUTPUT: double array (sizebub    ) e.g. (100)
       xbarD_,    //INOUTPUT: double array (sizexbarp,sizexbarh) e.g. (100*60) array size passed down
       &H0_,       //INOUTPUT: double lower limit of H 
       &H1_,       //INOUTPUT: double upper limit of H
       &P0_,       //INOUTPUT: double lower limit of P
       &P1_,       //INOUTPUT: double upper limit of P
       &x0_,       //INOUTPUT: double lower limit of x
       &x1_,       //INOUTPUT: double upper limit of x
       &ierr_);    //OUTPUT: error flag: 0 normal; -1 opflag out of range(OOR) ;-2 P OOR;-3: H OOR;-4 x OOR;-5: nonuniform interpolation pts

       return;


}


void utrc::interpCppWrapper(
       int    *opflag_,   //INPUT: integer flag -1:get array size and index strides; 0:initialize; 1:get thermo properties;
                         //                     2:get Pbub for given H; 3:get xbar for P,H
       double* Pinout_,   //INOUTPUT: double pressure in [bar] (absolute) from 0.1 to 20 as input, , when opflag.eq.2, it is output
       double* Hin_,      //INPUT:    double enthalpy in  [J/mol] from -1d5 to -5d4
       double* xin_,      //INPUT:    double quality (mass fraction of vapour) [non-D] from 0-0.1
       double* wtm_,      //OUTPUT:   double overall molecular weight [g/mol] set when opflag.eq.0
       double* rhoV_,     //OUTPUT:   double density of vap [kg/m^3]
       double* rhoL_,     //OUTPUT:   double density of liq [kg/m^3]
       double* cpV_,      //OUTPUT:   double Cp of vap [J/kg]
       double* cpL_,      //OUTPUT:   double Cp of liq [J/kg]
       double* muV_,      //OUTPUT:   double viscosity of vap [Pa-S]
       double* muL_,      //OUTPUT:   double viscosity of liq [Pa-S]
       double* lambdaV_,  //OUTPUT:   double thermo conductivity of vap [W/m-K]
       double* lambdaL_,  //OUTPUT:   double thermo conductivity of liq [W/m-K]
       double* aV_,       //OUTPUT:   double acoustic speed of vap [m/S]
       double* aL_,       //OUTPUT:   double acoustic speed of liq [m/S]
       double *TV_,       //OUTPUT: double temp of vap [K]
       double *TL_,       //OUTPUT: double temp of liq [K]
       double *Teq_,      //OUTPUT: double temp of system in equilibrium (P can differ) [K]
       double *SurfTen_,  // OUTPUT: double Surface Tension [N/m]
       int* superheat_,   //OUTPUT:   integer 1: superheated, 0: otherwise
       int* supercool_,   //OUTPUT:   integer 1: supercooled, 0: otherwise
       double* xbar_,     //OUTPUT:   double [non-D]
       double* Pbub_,
       int* sizep_,       //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       int* sizeh_,       //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       int* sizex_,       //INOUTPUT: array size in the x dimension, output when opflag=0, else input
       int* sizebub_,     //INOUTPUT: array size of bubble point curve data
       int* sizexbarp_,   //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       int* sizexbarh_,   //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       double* dP_,       //INOUTPUT: delta P in thermo-transport data array, output when opflag=0, else input
       double* dH_,       //INOUTPUT: delta H in thermo-transport data array, output when opflag=0, else input
       double* dx_,       //INOUTPUT: delta X in thermo-transport data array, output when opflag=0, else input
       double* dbubH_,    //INOUTPUT: delta H in bubble point data array,     output when opflag=0, else input
       double* dxbarP_,   //INOUTPUT: delta P in xbar data array,             output when opflag=0, else input
       double* dxbarH_,   //INOUTPUT: delta h in xbar data array,             output when opflag=0, else input
       double *rhoVD_,    //INOUTPUT: double array (sizep,sizeh,sizex) e.g. (100*60*20) array size passed down
       double *rhoLD_,    //INOUTPUT: double array (sizep,sizeh,sizex)
       double *cpVD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       double *cpLD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       double *muVD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       double *muLD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       double *lambdaVD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       double *lambdaLD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       double *aVD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       double *aLD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
        double *TVD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       double *TLD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       double *TeqD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       double *SurfTenD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       int *superheatD_,  //INOUTPUT: double array (sizep,sizeh,sizex)
       double *PbubD_,    //INOUTPUT: double array (sizebub    ) e.g. (100)
       double *HbubD_,    //INOUTPUT: double array (sizebub    ) e.g. (100)
       double *xbarD_,    //INOUTPUT: double array (sizexbarp,sizexbarh) e.g. (100*60) array size passed down
       double* H0_,       //INOUTPUT: double lower limit of H 
       double* H1_,       //INOUTPUT: double upper limit of H
       double* P0_,       //INOUTPUT: double lower limit of P
       double* P1_,       //INOUTPUT: double upper limit of P
       double* x0_,       //INOUTPUT: double lower limit of x
       double* x1_,       //INOUTPUT: double upper limit of x
       int* ierr_)   const  //OUTPUT: error flag: 0 normal; -1 opflag out of range(OOR) ;-2 P OOR;-3: H OOR;-4 x OOR;-5: nonuniform interpolation pts

{

       FORTNAME(interp)(
       opflag_,      //INPUT: integer flag -1:get array size and index strides; 0:initialize; 1:get thermo properties;
                         //                     2:get Pbub for given H; 3:get xbar for P,H
       Pinout_,   //INOUTPUT: double pressure in [bar] (absolute) from 0.1 to 20 as input, , when opflag.eq.2, it is output
       Hin_,      //INPUT:    double enthalpy in  [J/mol] from -1d5 to -5d4
       xin_,      //INPUT:    double quality (mass fraction of vapour) [non-D] from 0-0.1
       wtm_,      //OUTPUT:   double overall molecular weight [g/mol] set when opflag.eq.0
       rhoV_,     //OUTPUT:   double density of vap [kg/m^3]
       rhoL_,     //OUTPUT:   double density of liq [kg/m^3]
       cpV_,      //OUTPUT:   double Cp of vap [J/kg]
       cpL_,      //OUTPUT:   double Cp of liq [J/kg]
       muV_,      //OUTPUT:   double viscosity of vap [Pa-S]
       muL_,      //OUTPUT:   double viscosity of liq [Pa-S]
       lambdaV_,  //OUTPUT:   double thermo conductivity of vap [W/m-K]
       lambdaL_,  //OUTPUT:   double thermo conductivity of liq [W/m-K]
       aV_,       //OUTPUT:   double acoustic speed of vap [m/S]
       aL_,       //OUTPUT:   double acoustic speed of liq [m/S]
       TV_,       //OUTPUT: double temp of vap [K]
       TL_,       //OUTPUT: double temp of liq [K]
       Teq_,      //OUTPUT: double temp of system in equilibrium (P can differ) [K]
       SurfTen_,  // OUTPUT: double Surface Tension [N/m]
       superheat_,   //OUTPUT:   integer 1: superheated, 0: otherwise
       supercool_,   //OUTPUT:   integer 1: supercooled, 0: otherwise
       xbar_,     //OUTPUT:   double [non-D]
       Pbub_,
       sizep_,       //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       sizeh_,       //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       sizex_,       //INOUTPUT: array size in the x dimension, output when opflag=0, else input
       sizebub_,     //INOUTPUT: array size of bubble point curve data
       sizexbarp_,   //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       sizexbarh_,   //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       dP_,       //INOUTPUT: delta P in thermo-transport data array, output when opflag=0, else input
       dH_,       //INOUTPUT: delta H in thermo-transport data array, output when opflag=0, else input
       dx_,       //INOUTPUT: delta X in thermo-transport data array, output when opflag=0, else input
       dbubH_,    //INOUTPUT: delta H in bubble point data array,     output when opflag=0, else input
       dxbarP_,   //INOUTPUT: delta P in xbar data array,             output when opflag=0, else input
       dxbarH_,   //INOUTPUT: delta h in xbar data array,             output when opflag=0, else input
       rhoVD_,    //INOUTPUT: double array (sizep,sizeh,sizex) e.g. (100*60*20) array size passed down
       rhoLD_,    //INOUTPUT: double array (sizep,sizeh,sizex)
       cpVD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       cpLD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       muVD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       muLD_,     //INOUTPUT: double array (sizep,sizeh,sizex)
       lambdaVD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       lambdaLD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       aVD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       aLD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
        TVD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       TLD_, //INOUTPUT: double array (sizep,sizeh,sizex)
       TeqD_,      //INOUTPUT: double array (sizep,sizeh,sizex)
       SurfTenD_,      //INOUTPUT: double array (sizep,sizeh,sizex) 
       superheatD_,  //INOUTPUT: double array (sizep,sizeh,sizex)
       PbubD_,    //INOUTPUT: double array (sizebub    ) e.g. (100)
       HbubD_,    //INOUTPUT: double array (sizebub    ) e.g. (100)
       xbarD_,    //INOUTPUT: double array (sizexbarp,sizexbarh) e.g. (100*60) array size passed down
       H0_,       //INOUTPUT: double lower limit of H 
       H1_,       //INOUTPUT: double upper limit of H
       P0_,       //INOUTPUT: double lower limit of P
       P1_,       //INOUTPUT: double upper limit of P
       x0_,       //INOUTPUT: double lower limit of x
       x1_,       //INOUTPUT: double upper limit of x
       ierr_);       //OUTPUT: error flag: 0 normal; -1 opflag out of range(OOR) ;-2 P OOR;-3: H OOR;-4 x OOR;-5: nonuniform interpolation pts


       return;

}

} //end namespace foam

