#include <fstream>
#include <iostream>
#include <cmath>
#include<iomanip>
#include "Parameters.H"
#include "NewtonRaphson.H"

using namespace std;



double RachfordRice(double x_a, double K_a, double K_g)
{
  
  double x_g = 1-x_a;
//  double K_a = P_a*gam_a /pressure ;
//  double K_g = P_g* gam_g /pressure ;
 double guess= 0.001;


  double MolFract;


//  if (K_a <1.0 && K_g<1.0) 
 //     MolFract = 0.0;
 // else if (K_a >1.0 && K_g>1.0)
 //     MolFract =1.0;
 // else 

 //      cout<< x_a<<"  "<<x_g<<" "<<K_a<<"  "<<K_g<<endl;  
      
      MolFract = Newton(guess,x_a, x_g, K_a, K_g); 
 //      cout<<MolFract<<" "<< K_a<<" "<<K_g<<endl;

     
// double MolFract = (1-x_a*K_a-x_g*K_g)/(1.00-K_a-K_g+K_a*K_g);
  //return MolFract;  
  

    return MolFract ;     

}


double EnthalpyCalc(double T, double MolWt, double Cp, double Tref, double Href)
{
//In J/K/g from the paper by Nan
//    double Cp = 2.22;
//Converting to J/K/mol
    Cp = Cp * MolWt;
 
    return Cp*(T-Tref) + Href;     
     
}


//Rackett equation to calculate the density 
double Rackett(double x_a,double x_g, double Tc_a, double Tc_g,double Pc_a,double Pc_g,
               double Vc_a, double Vc_g, double ZRA_a,double ZRA_g, double  T_sys)
{
      double ZRA_M,  temp,temp2, temp1, temp3, phi_a, phi_g, Tcm, Tr, Vm, k12;

      ZRA_M = x_a*ZRA_a + x_g*ZRA_g;
      temp =Vc_a*x_a+ Vc_g*x_g;
      phi_a = x_a*Vc_a/temp;
      phi_g = x_g*Vc_g/temp;

       
      k12 = 1 - 8.0 * pow((Vc_a*Vc_g),0.5)/pow((pow(Vc_a,1.0/3.0)+ pow(Vc_g, 1.0/3.0)),3.0);  

       
      Tcm = phi_a*phi_a * Tc_a + phi_g*phi_g * Tc_g + 2.0* (1.0-k12)* phi_a*phi_g*pow((Tc_a*Tc_g),0.5);
      Tr = T_sys/Tcm;

      temp1 = 1.0+ pow((1-Tr), 2.0/7.0);
         
      temp2 =   pow(ZRA_M, temp1);
   
      temp3 = x_a*Tc_a/Pc_a+x_g*Tc_g/Pc_g;
      
      Vm = R* temp2* temp3;
       

      return 1/Vm; 
}     
