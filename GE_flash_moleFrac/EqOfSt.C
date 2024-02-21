#include <fstream>
#include <iostream>
#include <cmath>
#include<iomanip>
using namespace std;

#include"Parameters.H"
#include "EqOfSt.H"

EqOfSt::EqOfSt(double TC, double PC, double T, double AF, long Fl, double K_1):
Tcrit(TC), Pcrit(PC),  Tvap(T), ActFt(AF), Flag(Fl), K1(K_1)  
{


double fw, Tr, K0, alpha;
switch(Flag)
{
//Redlich-Kwong equation of state from Reid 
 case 1:
//    cout<<"RK EOS" <<endl;
    b = 0.08664*R*Tcrit/Pcrit;
    a = 0.42748*R*R*pow(Tcrit,2.5)/(Pcrit* pow(Tvap,0.5));  
//For RK and SRK
    u =1.0;
    w =0.0;
    break;  

// Souave Redlich-Kwong equation of state from Reid 
 case 2:
//    cout<<"SRK EOS"<<endl;
    b = 0.08664*R*Tcrit/Pcrit;
    fw =  0.48+1.574*ActFt -0.176*ActFt*ActFt;
    a = 0.42748*R*R*Tcrit*Tcrit/(Pcrit)* pow((1+fw*(1-pow(Tvap/Tcrit,0.5))),2.00);  
//For RK and SRK
    u =1.0;
    w =0.0;
    break;

//Peng Robinson 
  case 3:
//     cout<<"Peng Robinson "<<endl;
    b = 0.0778*R*Tcrit/Pcrit;
    fw =  0.37464+1.54226*ActFt - 0.26992*ActFt*ActFt;
    a = 0.45724*R*R*Tcrit*Tcrit/(Pcrit)* pow((1+fw*(1-pow(Tvap/Tcrit,0.5))),2.00);  
//For Peng Robinson
    u =2.0;
    w =-1.0;   
     break;

//Peng Robinson Gasem
  case 4:
 //   cout<<"Peng-Robinson-Gasem "<<endl;
    b = 0.0778*R*Tcrit/Pcrit;
    fw =  0.134+0.508*ActFt-0.0467*ActFt*ActFt;
    Tr = Tvap/Tcrit;
    a = 0.45724*R*R*Tcrit*Tcrit/(Pcrit)* exp( (2.0+0.836*Tr)*(1-pow(Tr,fw)));  
//For Peng Robinson
   u =2.0;
   w =-1.0;   
   break;

  case 5:
//    cout<<"Peng Robinson - Stryjek and Vera "<<endl;
    b = 0.077796*R*Tcrit/Pcrit;
    K0 = 0.378893+ 1.4897153* ActFt- 0.17131848*ActFt*ActFt+ 0.019655*ActFt*ActFt*ActFt;
    Tr = Tvap/Tcrit;
    fw = K0 + K1*(1.0+ pow(Tr,0.5))*(0.7-Tr) ;
    alpha =  pow((1.0+ fw*(1.0 - pow(Tr,0.5))),2);
    a = 0.457235 *(R*R*Tcrit*Tcrit)*alpha/Pcrit;     
   break;

  case 6:
 //   cout<<"Peng Robinson - Stryjek and Vera for Vapor Density "<<endl;
    b = 0.077796*R*Tcrit/Pcrit;
    K0 = 0.378893+ 1.4897153* ActFt- 0.17131848*ActFt*ActFt+ 0.019655*ActFt*ActFt*ActFt;
    Tr = Tvap/Tcrit;
    fw = K0 + K1*(1.0+ pow(Tr,0.5))*(0.7-Tr) ;
    alpha =  pow((1.0+ fw*(1.0 - pow(Tr,0.5))),2);
    a = 0.457235 *(R*R*Tcrit*Tcrit)*alpha/Pcrit;     
   break;


}

}

