#include <fstream>
#include <iostream>
#include <cmath>
#include<iomanip>

#include "NewtonRaphson.H"
using namespace std;


double Newton(double guess, double z_a, double z_g, double K_a, double K_g)
{
 double x_new, x,error;

 x_new =guess;
 long num=0;
 do
  {
     x=x_new;
     x_new = x - (funct(x,z_a,z_g,K_a,K_g)/deriv(x,z_a,z_g,K_a,K_g)) ;
 
     error = abs((x-x_new)/x);   

//     cout<<"Newton  "<< x <<" "<<x_new<<" "<<error<<" "<<  funct(x,z_a,z_g,K_a,K_g)<<" "<< deriv(x,z_a,z_g,K_a,K_g) <<endl; 
    
     num = num+1;
  }
  while ((error > 1e-6))  ;
//  while ((error > 1e-3) && num <10)  ;

  return x_new;
}

double funct(double x, double z_a, double z_g, double K_a, double K_g)
{

  return z_a*(K_a-1.0)/(1.0+x*(K_a-1.0)) + z_g*(K_g-1.0)/(1.0+x*(K_g-1.0)) ;
}

//derivative expression from introduction to chemical thermodynamics by Smith et al.
double deriv(double x, double z_a, double z_g, double K_a, double K_g)
{
 return   -z_a*pow((K_a-1.0),2.0)/ pow((1.0+x*(K_a-1.0)),2.0) - z_g*pow((K_g-1.0),2.0)/pow((1.0+x*(K_g-1.0)),2.0)  ;
}


