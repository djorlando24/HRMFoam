#include <fstream>
#include <iostream>
#include <cmath>
#include<iomanip>
using namespace std;

#include"Parameters.H"
#include"Input.H"
#include "FuelVap.H"

//Constructor- gets the temperature and pure vapor pressure tables and system temperature
//and enthalpy
Fuel_Props::Fuel_Props(Input Fuel):
gam_a(0.0),gam_g(0.0),PVap(0.0),PDew(0.0),
hfg(0.0), Temp(0.0), rho_l(0.0), rho_v(0.0),
P_alcohol(0.0), P_gasoline(0.0),
VapMolFract(0.0), Tref(0.0), Href(0.0), x_a(Fuel.x_a),
x_g(Fuel.x_g), MolWt(Fuel.MixMol), CP(Fuel.CP), Mu(Fuel.Mu), 
Cond(Fuel.Cond)
 {
  
  int i;
  points = Fuel.points;
  P_a = new double[points];  
  P_g = new double[points];
  T   = new double[points];


  for (i=0;i<points;i++)
   {   T[i] = Fuel.T[i];
      P_a[i] =Fuel.P_a[i];
      P_g[i] =Fuel.P_g[i];
  } 
   
//  x_a = Fuel.x_a;
//  x_g = Fuel.x_g;
//  MolWt = Fuel.MixMol;
//  CP = Fuel.CP;
//  Mu = Fuel.Mu;
//  Cond = Fuel.Cond;
 }

void Fuel_Props::InterpCall(double Temperature)
{ 
 
   double step_T = (T[points-1] - T[0])/(points-1);


  double temp_1 = Temperature - T[points-1];
  double temp_2 = Temperature - T[0];
   long upper,lower;


  if(( temp_1<=0 ) && (temp_2>=0))
       lower = floor((Temperature - T[0])/step_T);
  else if (temp_1 >0)
       lower = points-2;
  else if (temp_2 < 0)
       lower = 0; 
 
   upper = lower+1;
 
   double temp_P_a = Interp(T[lower], T[upper], P_a[lower], P_a[upper], Temperature);  
   double temp_P_g = Interp(T[lower], T[upper], P_g[lower], P_g[upper], Temperature);  


//calculate gamma a and gamma g at the given temperature
   ActCoeff(Temperature);
//calculate vapor pressure and dewpoint pressure
   PVap = PVapCalculator(temp_P_a , temp_P_g);
   PDew = DewPCalculator(temp_P_a , temp_P_g);  
//calculate enthalpy of vaporization 
   hfg =  HfgCall();
//pure component vapor pressures at the given temperature
   P_alcohol = temp_P_a;
   P_gasoline = temp_P_g;
   
}


double Fuel_Props::Interp(double A1, double A2, double B1, double B2, double a)
{
   return (B1-B2)/(A1-A2)*(a-A1)+B1;
   
}

//Vapor pressure calculation using the Wilson equation parameters from Pumphrey et al.
//Introduction to chemical engineering thermodynamics by smith
double Fuel_Props::PVapCalculator(double P_a, double P_g)  
{
    double VaporPressure;
   
//    ActCoeff(P_a, P_g, gam_a, gam_g);
   //vapor pressure accordin to Pumphrey at al.
    VaporPressure = x_a * gam_a * P_a + x_g * gam_g * P_g;
   
   return VaporPressure;

}

//Dew Point calculation from the Wilson equation parameters from Pumphrey et al.
double Fuel_Props::DewPCalculator(double P_a, double P_g)  
{
    double DPPressure;
  
   // ActCoeff(P_a, P_g, gam_a, gam_g);
   //vapor pressure accordin to Pumphrey at al.
    DPPressure = 1.00/(x_a/(gam_a*P_a) + x_g/(gam_g*P_g)); 
      
   return DPPressure;
}



//acitivity coefficient calculator
void Fuel_Props::ActCoeff(double Temperature)
{
//setting parameters of the Wilsons equation - from Pumphrey et al
    double    G_12 = 0.1665;
    double    G_21 = 0.3527;

//For Methanol-Toluene mixture 
//  double    G_12 =   106.86/40.714 *  exp(- 1702.3230/(1.98721*Temperature));
//  double    G_21 =   40.714/106.86 *   exp(- 186.7324/(1.98721*Temperature));


   double temp_a, temp_g;

  //from differentiating Wilsons equation w.r.t. mole fration to give activity coeffs
  //we get the following equations
    temp_a = -log(x_a + G_12* x_g) + x_g * (G_12/(x_a+G_12*x_g) - G_21/(G_21*x_a+x_g));
  
    temp_g = -log(x_g + G_21* x_a) - x_a * (G_12/(x_a+G_12*x_g) - G_21/(G_21*x_a+x_g));
  
    gam_a = exp(temp_a);
    gam_g = exp(temp_g);
    
}


 
//Calls the vapor pressure calculator for the first and last point
//in the table (for slope calculation)
//Calls the subroutine which calculates the enthalpy of vaporization
double Fuel_Props::HfgCall()
{
   long i;
   double PVap1, TVap1, PVap2, TVap2;


//first point in the table
    i=0;
    PVap1 = PVapCalculator(P_a[i], P_g[i]);     
    TVap1 = T[i];   

//last point in the table 
    i = points-1;

    PVap2 = PVapCalculator(P_a[i], P_g[i]);
    TVap2 = T[i];    
    
//subroutine which calculates the enthalpy of vaporization 
    return  HfgCalculator(PVap1,PVap2,TVap1,TVap2);
        
}


//Calculates the enthalpy of vaporization  from clausius clapeyrn equation
double  Fuel_Props::HfgCalculator(double P1, double P2, double T1, double T2)
{
  double slope;

//Clausius clapeyron equation according to Kenneth Kar et al.
  slope = (log(P2) - log(P1))/(1.0/T2-1.0/T1);

//hfg in j/mol
  hfg =  -slope *R ;
  
return hfg;
}

     
