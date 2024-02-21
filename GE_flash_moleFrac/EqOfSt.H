

//calculates the cubic equation parameters for each component
class EqOfSt
{
private:
double a, b, u, w;

public:
double Tvap, Pcrit, Tcrit, M;
double ActFt, K1;
long Flag;
EqOfSt(double,  double, double, double, long, double);

void EOS();
double a_eos()
   {
    return a;
    }
double b_eos()
   {
    return b;
    }

double u_eos()
   {
    return u;
    }
double w_eos()
   {
    return w;
    }



//EqOfSt();
 
};

