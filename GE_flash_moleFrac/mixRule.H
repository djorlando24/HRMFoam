
//calculates the mixture coefficients using the mixing rule
//and the saturated liquid and vapour densities
class MixRule
{
private:
double a_m, b_m, u_m,w_m;
double  maxRoot, minRoot;
double x1,x2,x3;
double a, b, c, d;
double SatVol_liq, SatVol_vap;

void minmaxFinder();
void cubic();
void CubicEqParam();
void SatDensity();


public:

double Pvap, Tvap;

MixRule(double, double, double, double,double,
        double, double, double, double, double );

double SatDen_liq,SatDen_vap;
};

