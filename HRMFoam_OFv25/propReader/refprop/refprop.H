#ifndef refstat
#define refstat
#include "./propReader/thermoBase.H"
#include "IOdictionary.H"

// Length of char buffer for reading fluid.dat
#define FLUID_DAT_BUFFER_LEN 255

//using namespace std;
//this is the derived class which reads the property table for refprop
//and performs interpolations
namespace Foam
{


class refprop
: public thermoBase
{

private:

    double *pressure_, *h_, *xbar_, *rho_, *rhov_, *rhol_, *temp_, *visc_, *cond_, *Liq_Mol_a_, *Vap_Mol_a_, *sos_, *pSat_ ;
    char junk[FLUID_DAT_BUFFER_LEN+1];
    long points,mark,hres,pres;

    mutable double xbar_cal, rho_cal, rhov_cal, rhol_cal, temp_cal, visc_cal, cond_cal, Liq_Mol_a_cal, Vap_Mol_a_cal, sos_cal, pSat_cal;
    double  Pc;

    double p_st,p_end,h_st,h_end,step_h,step_p;

    mutable double fac_x,fac_y;
    mutable long indp, indh;
    mutable long NW, NE, N, S, SW,SE;


//- Disallow default bitwise copy construct
    refprop(const refprop&);

//- Disallow default bitwise assignment
    void operator=(const refprop&);

    void errorCheck(const long error) const;

public:
//public variables

  TypeName("refprop");

  refprop(const volVectorField& U, const IOobject& thermPoint);

  double interp(double*)  const;
  //member functions
    
  void getProperties(const double& P, const double& H, const double& dummy) const ;


  inline double xbar() const
  { 
    return xbar_cal;
  }
  
      
  inline double rho() const
  { 
    return rho_cal;
  }
   
  inline double rhov() const
  {
    return rhov_cal;
  }

  inline double rhol() const
  {  
    return rhol_cal;
  }
   
  inline double temp() const
  { 
    return temp_cal;
  }

  inline double mu() const
  {  
    return visc_cal;
  }

  inline double K() const
  {  
    return cond_cal;
  }

  inline double sos() const
  {  
    return sos_cal;
  }

    inline double Liq_Mol_a() const
  {
    return Liq_Mol_a_cal;
  }

  inline double Vap_Mol_a() const
  {
    return Vap_Mol_a_cal;
  }


  inline double psi() const  // warning:  refprop does not currently produce usable speed of sound information
  {  
    return 1.0/(sos()*sos());
  }

  inline double pSat() const
  { 
    return pSat_cal;
  }

  inline double Pcrit() const
  {
    return Pc;
  }

  inline int NCall() const
  {
    return 1;
  }

  void printTableSize();

  inline double pStart() const
  {
    return p_st;
  }

  inline double pEnd() const
  {
    return p_end;
  }

  inline double hStart() const
  {
    return h_st;
  }

  inline double hEnd() const
  {
    return h_end;
  }

};

}//end namespace foam


#endif



