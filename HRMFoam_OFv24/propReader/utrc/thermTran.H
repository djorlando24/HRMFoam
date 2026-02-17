/* defines the basic CSPComp that
last modified on : 21st may 2007
By: Jerry Lee
Location: UTRC UTC (Copy righted)
actions: included DensePort stuffs from Sophia

*/

//#include <string>
//#include <iostream>
//using namespace std;
#include "fortran_mangle.H"


extern "C"
{
  void FORTNAME(interp)(
       int *opflag_,      //INPUT: integer flag -1:get array size and index strides; 0:initialize; 1:get thermo properties;
                         //                     2:get Pbub for given H; 3:get xbar for P,H
       double *Pinout_,   //INOUTPUT: double pressure in [bar] (absolute) from 0.1 to 20 as input, , when opflag.eq.2, it is output
       double *Hin_,      //INPUT:    double enthalpy in  [J/mol] from -1d5 to -5d4
       double *xin_,      //INPUT:    double quality (mass fraction of vapour) [non-D] from 0-0.1
       double *wtm_,      //OUTPUT:   double overall molecular weight [g/mol] set when opflag.eq.0
       double *rhoV_,     //OUTPUT:   double density of vap [kg/m^3]
       double *rhoL_,     //OUTPUT:   double density of liq [kg/m^3]
       double *cpV_,      //OUTPUT:   double Cp of vap [J/kg]
       double *cpL_,      //OUTPUT:   double Cp of liq [J/kg]
       double *muV_,      //OUTPUT:   double viscosity of vap [Pa-S]
       double *muL_,      //OUTPUT:   double viscosity of liq [Pa-S]
       double *lambdaV_,  //OUTPUT:   double thermo conductivity of vap [W/m-K]
       double *lambdaL_,  //OUTPUT:   double thermo conductivity of liq [W/m-K]
       double *aV_,       //OUTPUT:   double acoustic speed of vap [m/S]
       double *aL_,       //OUTPUT:   double acoustic speed of liq [m/S]
       double *TV_,       //OUTPUT: double temp of vap [K]
       double *TL_,       //OUTPUT: double temp of liq [K]
       double *Teq_,      //OUTPUT: double temp of system in equilibrium (P can differ) [K]
       double *SurfTen_,  // OUTPUT: double Surface Tension [N/m]
       int *superheat_,   //OUTPUT:   integer 1: superheated, 0: otherwise
       int *supercool_,   //OUTPUT:   integer 1: supercooled, 0: otherwise
       double *xbar_,     //OUTPUT:   double [non-D]
       double *Pbub_,
       int *sizep_,       //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       int *sizeh_,       //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       int *sizex_,       //INOUTPUT: array size in the x dimension, output when opflag=0, else input
       int *sizebub_,     //INOUTPUT: array size of bubble point curve data
       int *sizexbarp_,   //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       int *sizexbarh_,   //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       double *dP_,       //INOUTPUT: delta P in thermo-transport data array, output when opflag=0, else input
       double *dH_,       //INOUTPUT: delta H in thermo-transport data array, output when opflag=0, else input
       double *dx_,       //INOUTPUT: delta X in thermo-transport data array, output when opflag=0, else input
       double *dbubH_,    //INOUTPUT: delta H in bubble point data array,     output when opflag=0, else input
       double *dxbarP_,   //INOUTPUT: delta P in xbar data array,             output when opflag=0, else input
       double *dxbarH_,   //INOUTPUT: delta h in xbar data array,             output when opflag=0, else input
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
       double *H0_,       //INOUTPUT: double lower limit of H 
       double *H1_,       //INOUTPUT: double upper limit of H
       double *P0_,       //INOUTPUT: double lower limit of P
       double *P1_,       //INOUTPUT: double upper limit of P
       double *x0_,       //INOUTPUT: double lower limit of x
       double *x1_,       //INOUTPUT: double upper limit of x
       int *ierr_);       //OUTPUT: error flag: 0 normal; -1 opflag out of range(OOR) ;-2 P OOR;-3: H OOR;-4 x OOR;-5: nonuniform interpolation pts

};
