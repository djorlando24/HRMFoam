#include <fstream>
#include <iostream>
#include <cmath>
#include<iomanip>
using namespace std;

#include "mixRule.H"
#include "Input.H" 


//CONSTRUCTOR###############################################################
//constructor reads in the experimental data table
Input::Input()
{
   char junk[200];
//file with pure fluid properties
  ifstream infile("PureVap.inp");
  infile>>points;
  infile.getline(junk,199);
  infile.getline(junk,199);
   
  P_a = new double[points];  
  P_g = new double[points];
  T   = new double[points];


//setting the vapor pressure of the ethanol and gasoline from experiments 
//Kpa
//

 for(long i=0;i<points;i++)
   {   infile>>T[i]>>P_a[i]>>P_g[i];
      infile.getline(junk,199);
   }


  infile.close(); 
  ifstream thermo("ThermoProps.inp");
  thermo.getline(junk,199);
  thermo>>Flag;


  thermo.getline(junk,199);
  thermo.getline(junk,199);
  thermo.getline(junk,199);

  thermo>>Tcrit_a>>Pcrit_a>>Vcrit_a>>MolWt_a>>x_a; 

  thermo.getline(junk,199);
  thermo.getline(junk,199);
  thermo.getline(junk,199);
  thermo>>ActFt_a>>K1_a>>ZRA_a>>CP_a>>Mu_a>>Cond_a;

  thermo.getline(junk,199);
  thermo.getline(junk,199);
  thermo.getline(junk,199);
  thermo.getline(junk,199);

  thermo>>Tcrit_g>>Pcrit_g>>Vcrit_g>>MolWt_g>>x_g;
  thermo.getline(junk,199);
  thermo.getline(junk,199);
  thermo.getline(junk,199);
  thermo>>ActFt_g>>K1_g>>ZRA_g>>CP_g>>Mu_g>>Cond_g;

  thermo.close();

//calculate the vapor pressure and dewpoint pressur of all the experimental data points
//  PVapCall();
   MixMol = MolWt_a*x_a + MolWt_g * x_g;


 //specific heat, viscosity, conductivity of the mixture is a molar average
   CP = CP_a*x_a + CP_g*x_g;
   Mu = Mu_a*x_a + Mu_g*x_g;
   Cond =Cond_a*x_a + Cond_g*x_g; 

}
//END CONSTRUCTOR##################################################################

//MEMBER FUNCTIONS#################################################################



