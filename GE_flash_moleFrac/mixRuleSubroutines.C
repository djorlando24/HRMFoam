#include <fstream>
#include <iostream>
#include <cmath>
#include<iomanip>
using namespace std;


#include"Parameters.H"
#include "mixRule.H"

//CONSTRUCTOR################################################


MixRule::MixRule(double P, double T, double Mol_a, double Mol_gas, double Eth_a, double Eth_b, double Gas_a, double Gas_b, double u, double w):
Pvap(P), Tvap(T) , u_m(u), w_m(w)
{

  b_m = Mol_a * Eth_b + Mol_gas * Gas_b;
 
//ignoring binary interation parameter
  a_m =pow((Mol_a* sqrt(Eth_a)+ Mol_gas * sqrt(Gas_a)),2);

//  u_m =1.0;
//  w_m =0.0;


  CubicEqParam();
  cubic();  
  SatDensity();
} 

//DONE CONSTRUCTOR##########################################



//analytical solution of cubic equation to find three roots
void MixRule::cubic()
{ 

  double f,g,h,i,j,k,l,m,n,p, R, S, T,U;
  f=( (3.0*c/a) - ( b*b/(a*a)))/3.0;

  g=((2.0*b*b*b/(a*a*a)) - (9.0*b*c/(a*a)) + (27.0*d/a))/27.0;
  
  h= (g*g/4.0) + (f*f*f/27.0);

  if (h<= 0.0)
    {
     i = pow(((g*g/4.0)-h),0.5);
     j = pow(i, 1/3.0);

     k= acos(-g/(2.0*i));
     l= j*-1;
     m= cos(k/3.0);
     n= pow(3,0.5) * sin(k/3.0);
     p= (b/(3.0*a))*-1;
    
     x1= 2.0*j*cos(k/3.0)-(b/(3.0*a));
     x2= l*(m+n)+p;
     x3= l*(m-n)+p;
    }

else if(h>0)
    {
     cout<<" ONLY ONE REAL ROOT "<<endl;      
     
     R = -(g/2.0)+pow(h,0.5);
     if(R<0) 
        S =- pow(-R,1.0/3.0);
     else
        S = pow(R, 1.0/3.0)  ;      

     T = -(g/2.0) - pow(h,0.5);

     if(T<0) 
          U = -pow(-T,1.0/3.0);
     else
          U =pow(T,1.0/3.0);
  
     x1 = (S+U)-(b/(3.0*a));

     x2=1e10;
     x3=1e10;
    }

minmaxFinder();
}


//find the maximum and minimum roots which correspond to the compressibility 
void MixRule::minmaxFinder()
{
  maxRoot = max(x1,(max(x2,x3)));
  
  if(x1<0)  minRoot = min(x2,x3);
  else if (x2<0) minRoot = min(x1,x3);
  else minRoot = min(x1,x2); 

}

//Calculate the molar volumes from the compressibility
void MixRule::SatDensity()
{
     

  SatVol_vap = maxRoot*R*Tvap/Pvap;
  
  SatVol_liq = minRoot*R*Tvap/Pvap;    


// SatVol_liq is in m^3/mol
  SatDen_liq = 1/(SatVol_liq);  

  SatDen_vap = 1/(SatVol_vap);  
 
}



//Calculates the 
void MixRule::CubicEqParam()
{
//cubic equation coefficients 
//az^3+bz^2+cz+d=0 
  
  double tempA, tempB;

  tempA = a_m*Pvap/(R*R*Tvap*Tvap);
  tempB = b_m*Pvap/(R*Tvap);

  a= 1.0;
  b= -( 1+ tempB - u_m*tempB) ;
  c= (tempA +w_m*tempB*tempB-u_m*tempB- u_m* tempB*tempB) ; 
  d= -tempA*tempB - w_m*tempB*tempB - w_m * tempB*tempB*tempB ;


}




      
