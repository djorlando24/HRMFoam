
#include <fstream>
#include <iostream>
#include <cmath>
#include<iomanip>

using namespace std;

int main()
{
  double i;
  ofstream outfile("pistonStroke.dat");
  outfile<<"("<<endl;
  i=0.0;

  do 
{
   outfile<<"( "<< i <<"  "<<(0.18E-3 * sin(2.00 * 3.141592654 * 1/(0.0045) * i))<<" ) "<<endl; 
   i= i+0.00005; 
  
} while (i<= 0.00225);

   outfile<<"(100 0 )"<<endl;
   outfile<<")"<<endl;


return 0;
}
