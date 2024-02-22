#include<iostream>
#include<fstream>
#include"math.h"
#include"refprop.H"
#include "addToRunTimeSelectionTable.H"
//using namespace std;

namespace Foam
{

defineTypeNameAndDebug(refprop,0);
addToRunTimeSelectionTable(thermoBase, refprop, dictionary);
                        

refprop::refprop(const volVectorField& U, const IOobject& thermPoint )
: thermoBase(U, thermPoint)
{
  std::ifstream infile("fluid.dat");

  if(!infile.good())
  {
   errorCheck(1);
  }
     

  infile.getline(junk,FLUID_DAT_BUFFER_LEN);
  infile >> junk >> p_st >> junk >> p_end >> junk >> step_p;
  infile >> junk >> h_st >> junk >> h_end >> junk >> step_h;
  infile >> junk >> points >> junk >> pres >> junk >> hres;
  infile >> junk >> Pc;
  infile.getline(junk,FLUID_DAT_BUFFER_LEN);
  infile.getline(junk,FLUID_DAT_BUFFER_LEN);

  if(!infile.good())
  {  
   errorCheck(3);
  }

  pressure_=new double[points];
  h_ = new double[points];
  xbar_=new double[points];
  rho_ = new double[points];
  rhov_=new double[points];
  rhol_=new double[points];
  temp_=new double[points];
  visc_=new double[points];
  cond_=new double[points];
  cv_=new double[points];
  cp_=new double[points];
  sos_=new double[points];
  pSat_= new double[points]; 

  
  for(long i=0;i<points;i++)
  {
      infile >> pressure_[i] >> h_[i] >> xbar_[i] >> rho_[i] >> rhov_[i] >> rhol_[i] 
             >> temp_[i] >> visc_[i] >> cond_[i] >> cv_[i] >> cp_[i] >> sos_[i] >>pSat_[i];
      infile.getline(junk,FLUID_DAT_BUFFER_LEN);    
  }
        

   if(!infile.good()) 
   { 
      errorCheck(2);
   }
  //p_st and h_st are the lowest pressure and enthalpy in the fluid.dat
  //p_end and h_end are the highest
  //step_p and step_h are the step sizes in the table
  //pres and hres are the number of p points and h points
  //points is pres * hres
 
}


//this function is called in the file enthfuncs and does the interpolations to get the
//properties from the refprop table
void refprop::getProperties(const double& P,const double& H , const double& dummy ) const
{
  
 
  indp =std::min(long(pres-2),std::max(long(floor((P - p_st) / step_p)),long (0)));
  indh = std::min(long(hres-2),std::max(long(floor((H - h_st) / step_h)),long (0)));
  fac_y = std::min(double(1),std::max(double(0),(P - pressure_[indp*hres+indh]) / step_p));
  fac_x = std::min(double(1),std::max(double(0),(H - h_[indp*hres+indh]) / step_h));
   
  NW = indp*hres+indh ;
  NE = indp*hres+indh+1;
  SW = (indp+1)*hres+indh;
  SE = (indp+1)*hres+indh+1;

 
  xbar_cal=  interp(xbar_);
   
  temp_cal=  interp(temp_); 

  rho_cal=  interp(rho_);

  rhol_cal=  interp(rhol_);  

  rhov_cal=  interp(rhov_);

  visc_cal=  interp(visc_);    
 
  cond_cal=  interp(cond_);

  cv_cal=  interp(cv_);

  cp_cal=  interp(cp_);   

  sos_cal= interp(sos_); 
 
 
  pSat_cal= interp(pSat_); 
  

 if (P < p_st)
   {
     xbar_cal=1;
   }

  //unit change
  visc_cal=visc_cal/1e6;

}

inline double refprop::interp(double* var) const
{ 

  return  (var[NW] + fac_y*(var[SW] - var[NW]) + fac_x*(var[NE]-var[NW] + fac_y*(var[SE] - var[SW] - var[NE] + var[NW])));
 
 } 


void refprop::printTableSize()
  {
   std::cout << "fluid.dat table Size"<<std::endl;
   std::cout << p_st<<" "<<p_end<<" "<<step_p<<std::endl;
   std::cout << h_st<<" "<<h_end<<" "<<step_h<<std::endl;
   std::cout <<points<<" "<<pres<<" "<<hres<<std::endl;
  }

void refprop::errorCheck(const long error)  const
{

//Checking for errors in the fluid.dat file
  if (error == 1)
   {
       FatalIOErrorIn
       (
           "refprop::refprop()","COULD NOT FIND PROPERTY FILE :: fluid.dat"
        )   << exit(FatalIOError);
   }
  else if(error == 2)
    {
       FatalErrorIn
       (
           "refprop::refprop"
	   ) 
	   << "ERROR IN PROPERTY FILE::fluid.dat \n"
           << "MUST HAVE ALL OF  THE FOLLOWING VARIABLES \n"
           << "Pressure \n"<<"Enthalpy\n"<<"Quality\n"<<"Density\n"<<"DensityVapor\n"<<"DensityLiquid\n"<<"Temperature\n"
           << "Viscosity\n"<<"Conductivity\n"<<"Cv(KJ/KgK)\n"<<"Cp(KJ/KgK)\n"<<"SoS(m/s)\n"<<"PSat(Pa)"
           <<  exit(FatalError);
   }
   else if (error == 3)
   {
       FatalIOErrorIn
       (
           "refprop::refprop","CRITICAL PRESSURE MISSING IN PROPERTY FILE ::fluid.dat"
        )   << exit(FatalIOError);
   } 


}

} //end namespace foam
