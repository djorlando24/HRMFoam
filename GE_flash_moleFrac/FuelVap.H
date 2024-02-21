class Fuel_Props
{

private:
 
void InterpCall();

//wilson equation coefficients  
  double G_12, G_21;

double PVapCalculator( double , double);
double DewPCalculator( double , double);
void ActCoeff(double);
double Interp(double, double, double, double, double);
double converter();
double HfgCall();
double HfgCalculator(double , double , double , double);

public:
 Fuel_Props(Input);
 void  InterpCall(double); 
 double x_a, x_g, MolWt;
 double gam_a, gam_g;
 double PVap , CP, Mu, Cond;
 double PDew;
 double hfg, Temp;
 double rho_l, rho_v;
double *T, *P_a, *P_g;
 double P_alcohol, P_gasoline;
 double VapMolFract, Tref, Href;
 long points;
};

