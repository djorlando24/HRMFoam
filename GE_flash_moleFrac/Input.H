class Input
{
private:

  char junk[100];
  long Flag;

public:
//mole fractions of alcohol and gasoline
  double x_a, x_g;  
  
//vapor pressures of pure gasoline and pure ethanol and of the mixture
  double *P_a, *P_g, *T ;

//number of temperatures with vapor pressure data
  long  points;

// thermo properties for alcohol and gasoline
  double Tcrit_a,Tcrit_g, Pcrit_a, Pcrit_g, MolWt_a, MolWt_g, Vcrit_a, Vcrit_g;
  double MixMol, ActFt_a, ActFt_g, K1_a, K1_g, ZRA_a, ZRA_g;
  double CP_a, CP_g, CP, Mu_a, Mu_g, Mu, Cond_a, Cond_g, Cond;

  double Flag_eos()
   { return Flag;
   }


  Input();
};


