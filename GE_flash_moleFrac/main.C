#include <fstream>
#include <iostream>
#include <cmath>
#include<iomanip>
#include"stdlib.h"
#include"EqOfSt.H"
#include"mixRule.H"
#include"Input.H"
#include "FuelVap.H"

double RachfordRice( double, double, double);
double EnthalpyCalc( double, double, double, double, double);
double Fvalue(double, double, double , Fuel_Props, double& );
double Rackett(double, double,double, double ,double, double, double, double, double, double,double);     
using namespace std;


//MAIN PROGRAM#####################################################################
int main()
{

//  double x_alcohol;
  char junk[100]; 

  double P_sys, H_start, H_end, H_sys, Tref, Href, P_start, P_end; 
  double  stepH, stepP;
  long pres, hres, resolution, H_iter, P_iter;

//Reader for pressure and Enthalpy############################################################################
  ifstream infile("PH.inp");

  infile.getline(junk,199);
  infile>>Tref>>Href;
  infile.getline(junk,199);
  infile>>junk>>P_start>>junk>>P_end>>junk>>stepP;
  infile.getline(junk,199);
  infile>>junk>>H_start>>junk>>H_end>>junk>>stepH;

//Done reading pressure and enthalpy#######################################################################
   
   pres = (P_end - P_start)/stepP +1; 
   hres = (H_end - H_start)/stepH +1;
   resolution = pres*hres; 

//Write header for  VaporPressure.dat##########################################################################
  ofstream outfile("fluid.dat");
  outfile.precision(10);

//Reads in all the input files
//PureVap.inp
//ThermoProps.inp

   Input Fuel;
//mole fraction of alcohol

  double Mol_alc = Fuel.x_a;
  double Mol_gas = Fuel.x_g ;    
   
  cout<<" MolWt "<<Fuel.MixMol<<endl;


//Original liquid density of the mixture######################################################################
  Fuel_Props Gasahol_orig(Fuel);

  Gasahol_orig.InterpCall(Tref);
 
//Calculating  Volume Fraction of Ethanol in the mixture######################################################
//In the initial mixture at reference temperature 
#include "VolFractEth.H"
//done calculating volume fraction of Ethanol #############################################################################
//done calculatiing the original liquid density of the mixture ############################################################

  outfile<<"GEFlash Ethanol Fuel = E"<<VolFract*100<<" with "<<Fuel.x_a<<"  moles ethanol" <<endl;
  outfile<<"P_start "<<P_start<<" P_end "<<P_end<<" stepP "<<stepP<<endl;
  outfile<<"H_start "<<H_start<<" H_end "<<H_end<<" stepH "<<stepH<<endl;
  outfile<<"Points "<<resolution<<" pres "<<pres<<" Hres "<<hres<<endl;
  outfile<<"CriticalPressure  2.56E6    "<<endl;

  outfile<<"# Pressure(Pa) Enthalpy(j/kg) VapMassFract  Density(kg/m^3)  RhoV(kg/m^3)   RhoL(kg/m^3)  Temp(K) Viscosity(microPa.s)  "
            <<"  Conductivity(W/m.K)    Liq. Alc. Frac (-)  Gas Alc. Frac (-)  SOS(m/s)   PBUB(Pa)  "<<endl;


    cout<<" Final Enthalpy (J/kg) " <<Href + (Fuel.CP* (Fuel.T[Fuel.points -1] - Tref ) *1000.00  ) <<endl;    
//converting the enthalpy values to j/mol to be used in the code
H_start = H_start*Fuel.MixMol/1000.00;
Href = Href*Fuel.MixMol/1000.00;
H_end = H_end*Fuel.MixMol/1000.00;
stepH = stepH*Fuel.MixMol/1000.00;


//Start pressure loop
P_iter =0;
P_sys = P_start;
do
 {      

//Start Enthalpy loop
      H_sys = H_start; 
      H_iter = 0;
      do
    {
    
         //Gasahol will calculate the Pdew, PVap   
         //Gasahol will have all the final properties
           Fuel_Props Gasahol(Fuel);
         
         //(A)Initial guess of temperature  ###########################################
            double h_l,h_v; 
            double h_flash;
            double funct_value;
            double deriv_value;
            double VapMolFract=0.0;
            double error;
         
            Gasahol.Tref = Tref;
            Gasahol.Href = Href;
         
            double T_a = Gasahol.T[0] ;
            double T_b = Gasahol.T[Gasahol.points-1] ;
            
            double T_sys;

//for the end temperature, calculate the  vapor pressure########################################
//and dew point pressure because the vaporpressure needs to be for the subcooled temperature
            double Temp_end = (H_sys - Gasahol.Href)/(Gasahol.CP * Fuel.MixMol)+Tref;   
             
            if((Temp_end)>T_b)
                cout<<"End temperature specified is higher than limit of "<<T_b<<endl;

//the vapor presssure required in HRMFoam is for the subcooled temperature because it is an isenthalpic process
            Gasahol.InterpCall(Temp_end);
            double PVap_Tend = Gasahol.PVap;        
 
            double F_a =  Fvalue(T_a, P_sys, H_sys,Gasahol,VapMolFract);
         
            double F_b =  Fvalue(T_b, P_sys,H_sys,Gasahol,VapMolFract);
         
            double F_c;    

            
            if (F_a*F_b >0) 
             { cout<<"out of Range, Exiting"<<endl;
               cout<<"Function value at end points  "<<F_a<<" "<<F_b<<endl;
	       cout<<"P="<<P_sys<<", H="<<H_sys/Fuel.MixMol*1000.00<<endl;
               exit(1);
             }
                 
         //start phflsh loop
	   int nphflsh=0;
           do 
             {
                    
                  T_sys = (T_a+T_b)/2;  
              
                  VapMolFract =0;
                  F_c = Fvalue(T_sys, P_sys, H_sys,Gasahol, VapMolFract);    
              
                  if ((F_a *F_c) > 0.0)
                    {    F_a = F_c;
                         T_a = T_sys;
                    }
                 else if((F_b*F_c)> 0.0)
                    { 
                        F_b = F_c;
                        T_b = T_sys;
                    }
              
 
                   error = abs(T_a-T_b)/T_a ; 
		   nphflsh++;
		   if (nphflsh>1000) {
			   cout<<"Maximum number of iterations exceeded."<<endl;
			   cout<<"P="<<P_sys<<", H="<<H_sys/Fuel.MixMol*1000.00<<endl;
			   exit(1);
		   }
	     }
             while(error >1e-9);
            
           Gasahol.InterpCall(T_sys);
         
         
         //(C) Calculate the EOS parameters for the two pure components###################################
            EqOfSt Ethanol(Fuel.Tcrit_a, Fuel.Pcrit_a, T_sys, Fuel.ActFt_a, Fuel.Flag_eos(), Fuel.K1_a);
         
            EqOfSt Gasoline(Fuel.Tcrit_g, Fuel.Pcrit_g, T_sys, Fuel.ActFt_g, Fuel.Flag_eos(), Fuel.K1_g);
         
         
         //(D) Applying the mixing rule and calculating the mixture saturation density#####################
         
         //Calculate the Saturation density of the fuel
         

          MixRule Blend(Gasahol.PVap, T_sys, Mol_alc, Mol_gas, Ethanol.a_eos(), Ethanol.b_eos(), Gasoline.a_eos(), Gasoline.b_eos(),
                          Ethanol.u_eos(), Ethanol.w_eos());
         
         //densities are converted to kg/m^3 
            Gasahol.rho_l = Blend.SatDen_liq * Fuel.MixMol * 1E-3 ;
            Gasahol.rho_v = Blend.SatDen_vap * Fuel.MixMol * 1E-3 ;
         //Done calculating Saturated densities############################################################
         
           if(Fuel.Flag_eos() ==6) 
         //modified Rackett equation
            { 
         //Density is converted to kg/m^3
                  Gasahol.rho_l = Rackett(Mol_alc, Mol_gas, Fuel.Tcrit_a, Fuel.Tcrit_g, Fuel.Pcrit_a, Fuel.Pcrit_g, Fuel.Vcrit_a, Fuel.Vcrit_g,
                                         Fuel.ZRA_a, Fuel.ZRA_g, T_sys) * Fuel.MixMol * 1E-3;
            }
         
         //Calculating mass fraction of Vapor
         
              double K_a = Gasahol.P_alcohol * Gasahol.gam_a/P_sys;
              double K_g = Gasahol.P_gasoline * Gasahol.gam_g/P_sys;
              double  Vap_Mol_a, Liq_Mol_a, Vap_Mol_g, Liq_Mol_g; 
         
         
             if (P_sys >Gasahol.PVap) 
         //subcooled, the liquid mole fraction is equal to the total mole fraction and
         //vapor mole fractions are zero
                {   Vap_Mol_a = 0.0;
                    Vap_Mol_g = 0.0;    
                    Liq_Mol_a = Mol_alc;
                    Liq_Mol_g = Mol_gas;
		    cout<<"    Subcooled";
                 }
            else if (P_sys < Gasahol.PDew)       
         //superheated, the vapor mole fraction is equal to the total mole fractionand
         //liquid mole fractions are zero
                {  Vap_Mol_a = Mol_alc;
                   Vap_Mol_g = Mol_gas;
                   Liq_Mol_a = 0.0;       
                  Liq_Mol_g = 0.0;
		  cout<<"    Superheated";
                }
            else
               { //Vapor moleFract of ethanol
                 Vap_Mol_a = Gasahol.x_a*K_a/(1+VapMolFract*(K_a -1));
         //Liquid moleFract of ethonal
                 Liq_Mol_a = Vap_Mol_a/K_a;
                 Vap_Mol_g = Gasahol.x_g*K_g/(1+VapMolFract*(K_g -1));
                 Liq_Mol_g = Vap_Mol_g/K_g;
		 cout<<"      Mixture";
              }  
         
         
         //Molecular Weight of Liquid
               double MolWt_liq = Liq_Mol_a * Fuel.MolWt_a + Liq_Mol_g * Fuel.MolWt_g;
               double MolWt_vap = Vap_Mol_a * Fuel.MolWt_a + Vap_Mol_g * Fuel.MolWt_g;
              
               double VapMassFract = VapMolFract*MolWt_vap/( VapMolFract * MolWt_vap + (1-VapMolFract)* MolWt_liq); 
         
               double density = Gasahol.rho_l*Gasahol.rho_v/(Gasahol.rho_v*(1-VapMassFract) + VapMassFract * Gasahol.rho_l);
         
         //P in Pa, H in j/kg, xbar in kg/kg, rhoV and rhoL in kg/m^3,  
	       outfile<<fixed<<" "<<P_sys<<"     "<<H_sys/Fuel.MixMol*1000.00<<"  "<<VapMassFract<<"  "<<density<<"   " <<Gasahol.rho_v<<"  "
		  <<Gasahol.rho_l<<"  "<<T_sys<< " "<<Fuel.Mu*1e6<<"  "<<Fuel.Cond<<"   "<<Liq_Mol_a<<"   "<<Vap_Mol_a << "      1e-15   "<< PVap_Tend<<endl;

 
         //step in H
           H_sys = H_sys +stepH;
           H_iter++;    
    	   cout<< ", H=" << H_sys/Fuel.MixMol*1000.0 << " J/kg" <<endl;
    
	    }
             while(H_iter < hres);
 
           P_sys = P_sys +stepP;
           P_iter++;
	   cout<< "P=" << P_sys<<" Pa"<<endl;  
 }
 while(P_iter < pres);
  
  return 0;
 }


//END MAIN PROGRAM######################################################################


//
double Fvalue(double T_sys, double P_sys, double H_sys,Fuel_Props Gasahol, double& VapMolFract)
{
      double funct_value, K_a, K_g, h_l, h_v;
  

         Gasahol.InterpCall(T_sys);

//calculate h_l and h_v for that temperature
//  
     
         h_l = EnthalpyCalc(T_sys, Gasahol.MolWt, Gasahol.CP, Gasahol.Tref, Gasahol.Href);

         h_v = h_l+ Gasahol.hfg;
  
      
//(E) Calculate the molFraction from the Rachford Rice Equation##################################################
 
         if (P_sys >Gasahol.PVap) 
            VapMolFract =0.0;             
         else if (P_sys < Gasahol.PDew)       
            VapMolFract = 1.0;
         else
           { 
       
            K_a = Gasahol.P_alcohol * Gasahol.gam_a/P_sys;
            K_g = Gasahol.P_gasoline * Gasahol.gam_g/P_sys;
            
            //cout<< P_sys<<" "<< K_a<<" "<<K_g<<endl;
            VapMolFract =  RachfordRice(Gasahol.x_a, K_a, K_g);
           }
          if (VapMolFract > 1) 
             VapMolFract =1.0;
          else if (VapMolFract<0) 
             VapMolFract =0.0;

//Done calculating the molfraction of the gasahol###########################################################
//calculate the enthalpy from the flash calculation


     funct_value = VapMolFract*h_v+(1-VapMolFract)*h_l-H_sys;
   
   return funct_value;
}


