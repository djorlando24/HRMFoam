#ifndef utstat
#define utstat
#include"./propReader/thermoBase.H"
#include"IOdictionary.H"
#include"fortran_mangle.H"

//using namespace std;
//this is a derived class which provides properties for
//JP8
namespace Foam
{

class utrc
: public thermoBase
{
    
private:

//UTRC property table variables

     mutable  int     opflag_;
     mutable  double  Pinout_,  Hin_,  xin_;
      
    mutable double  wtm_, rhoV_,  rhoL_,cpV_,cpL_, muV_,  muL_, lambdaV_, lambdaL_, aV_;
    mutable double  aL_, TV_,  TL_, Teq_, SurfTen_;
    mutable int  superheat_, supercool_, *superheatD_;
    mutable double  xbar_, Pbub_, mu_, lambda_, psi_ ;

           
    mutable   int  sizep_,sizeh_, sizex_, sizebub_, sizexbarp_, sizexbarh_, ierr_;
    mutable   double  dP_, dH_, dx_, dbubH_, dxbarP_, dxbarH_, *rhoVD_, *rhoLD_, *cpVD_, *cpLD_, *muVD_, *muLD_, *lambdaVD_, *lambdaLD_;
    mutable   double *aVD_, *aLD_, *TVD_, *TLD_, *TeqD_, *SurfTenD_, *PbubD_, *HbubD_, *xbarD_;
    mutable   double  H0_, H1_,  P0_,  P1_,  x0_,  x1_;

       long error_;
// local variables
     mutable  int ierrout_;
     mutable  int NP_,NH_,NX_,NBUB_;

     //- Disallow default bitwise copy construct
             utrc(const utrc&);
     //
     //- Disallow default bitwise assignment
             void operator=(const utrc&);
                    


//END UTRC variables


  void   interpCppWrapper(    // returns the value of ierr
       int  *opflag_,      //INPUT: integer flag -1:get array size and index strides; 0:initialize; 1:get thermo properties;
                         //                     2:get Pbub for given H; 3:get xbar for P,H
       double  *Pinout_,   //INOUTPUT: double pressure in [bar] (absolute) from 0.1 to 20 as input, , when opflag.eq.2, it is output
       double  *Hin_,      //INPUT:    double enthalpy in  [J/mol] from -1d5 to -5d4
       double  *xin_,      //INPUT:    double quality (mass fraction of vapour) [non-D] from 0-0.1
       double  *wtm_,      //OUTPUT:   double overall molecular weight [g/mol] set when opflag.eq.0
       double  *rhoV_,     //OUTPUT:   double density of vap [kg/m^3]
       double  *rhoL_,     //OUTPUT:   double density of liq [kg/m^3]
       double  *cpV_,      //OUTPUT:   double Cp of vap [J/kg]
       double  *cpL_,      //OUTPUT:   double Cp of liq [J/kg]
       double  *muV_,      //OUTPUT:   double viscosity of vap [Pa-S]
       double  *muL_,      //OUTPUT:   double viscosity of liq [Pa-S]
       double  *lambdaV_,  //OUTPUT:   double thermo conductivity of vap [W/m-K]
       double  *lambdaL_,  //OUTPUT:   double thermo conductivity of liq [W/m-K]
       double  *aV_,       //OUTPUT:   double acoustic speed of vap [m/S]
       double  *aL_,       //OUTPUT:   double acoustic speed of liq [m/S]
       double *TV_,       //OUTPUT: double temp of vap [K]
       double *TL_,       //OUTPUT: double temp of liq [K]
       double *Teq_,      //OUTPUT: double temp of system in equilibrium (P can differ) [K]
       double *SurfTen_,  // OUTPUT: double Surface Tension [N/m]
       int  *superheat_,   //OUTPUT:   integer 1: superheated, 0: otherwise
       int  *supercool_,   //OUTPUT:   integer 1: supercooled, 0: otherwise
       double  *xbar_,     //OUTPUT:   double [non-D]
       double *Pbub_,
       int  *sizep_,       //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       int  *sizeh_,       //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       int  *sizex_,       //INOUTPUT: array size in the x dimension, output when opflag=0, else input
       int  *sizebub_,     //INOUTPUT: array size of bubble point curve data
       int  *sizexbarp_,   //INOUTPUT: array size in the p dimension, output when opflag=0, else input
       int  *sizexbarh_,   //INOUTPUT: array size in the h dimension, output when opflag=0, else input
       double  *dP_,       //INOUTPUT: delta P in thermo-transport data array, output when opflag=0, else input
       double  *dH_,       //INOUTPUT: delta H in thermo-transport data array, output when opflag=0, else input
       double  *dx_,       //INOUTPUT: delta X in thermo-transport data array, output when opflag=0, else input
       double  *dbubH_,    //INOUTPUT: delta H in bubble point data array,     output when opflag=0, else input
       double  *dxbarP_,   //INOUTPUT: delta P in xbar data array,             output when opflag=0, else input
       double  *dxbarH_,   //INOUTPUT: delta h in xbar data array,             output when opflag=0, else input
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
       double  *H0_,       //INOUTPUT: double lower limit of H 
       double  *H1_,       //INOUTPUT: double upper limit of H
       double  *P0_,       //INOUTPUT: double lower limit of P
       double  *P1_,       //INOUTPUT: double upper limit of P
       double  *x0_,       //INOUTPUT: double lower limit of x
       double  *x1_,       //INOUTPUT: double upper limit of x
       int  *ierr_)   const;       //OUTPUT: error flag: 0 normal; -1 opflag out of range(OOR) ;-2 P OOR;-3: H OOR;-4 x OOR;-5: nonuniform interpolation pts
// END interpCppWrapper declaration

public:

    TypeName("utrc");

 //Constructor
    utrc(const volVectorField& U, const IOobject& thermPoint);
    

 //Member functions
    void getArrayInfo() const ;
    
    void getWorkingArray()  const;

 //requires P,H and X to interpolate
    void getProperties(const double& P, const double& H, const double& X)  const;


    inline double xbar() const
    { 
       return xbar_;
    }

    inline double rho() const
    { 
       return  ( 1.0/( xbar_*(1.0/max(rhoV_,SMALL) -1.0/max(rhoL_,SMALL) ) + 1/max(rhoL_,SMALL)))    ;
    }

    inline double rhov() const
    {  
       return rhoV_;
    }

    inline double rhol() const
    {  
       return rhoL_;
    }

    inline double temp() const
    { 
       return TL_;
    }

    inline double visc_l() const
    { 
       return   muL_;
    }

    inline double visc_v() const
    {  
       return   muV_;
    }
   
    

    inline double cond_l() const
    { 
        return lambdaL_;
    }

    inline double cond_v() const
    {  
       return lambdaV_;
    }

    inline double mu() const
    {
       return mu_;
    }     

    inline double K() const
    {
          return lambda_;
    }
 
    inline double cp_V() const
    {  
       return cpV_;
    }


    inline double cp_L() const
    { 
       return cpL_;
    }

    inline double a_L() const
    {
       return aL_;
    }

    inline double a_V() const
    { 
       return aV_;
    }

     inline double psi() const
    { 
       return  psi_;
    }

 
    inline double Pcrit() const
    {
       //should be read in from the data file thermtran.dat 
       return 151.1E05;
    }

   inline double pSat() const
   {  
       return  (Pbub_*1e5) ;
   }

   inline int NCall() const
   {
       return 2;
   }

};

}//end namespace foam

#endif

