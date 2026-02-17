#ifndef thermoB
#define thermoB
#include "volFields.H"
#include "IOdictionary.H"
#include "runTimeSelectionTables.H"

//using namespace std;
//this is the abstract base class for the property table 
//readers 
namespace Foam
{


class thermoBase 
: 
  public IOdictionary

{

//private:

   //Disallow default bitwise copy construct
   thermoBase(const thermoBase&);

   //Disallow default bitwise assignment
   void operator=(const thermoBase&);


public:
//public variables

   TypeName("thermoBase");

   
   thermoBase(const volVectorField& U, const IOobject& thermPoint);  

   //selector
   static autoPtr<thermoBase> New
   (
           const volVectorField& U,
           const IOobject& thermPoint
   );


   //declare runtime selection table
   declareRunTimeSelectionTable
   (
        autoPtr,
        thermoBase,
        dictionary,
        (
             const volVectorField& U,
             const IOobject& thermPoint
        ),
        (U, thermPoint)
   );
   

   //Member functions
   virtual void getProperties(const double&, const double&, const double&)const=0;

   virtual double rho() const = 0 ;

   virtual double xbar() const = 0 ;

   virtual double rhov() const = 0 ;
 
   virtual double rhol() const = 0 ;

   virtual double Pcrit() const=0;

   virtual double psi() const=0;

   virtual double mu() const=0;

   virtual double pSat() const=0;

   virtual double K() const=0;

   virtual double temp() const=0;

   virtual int NCall() const=0;

   virtual double Liq_Mol_a() const=0  ;
  
   virtual double Vap_Mol_a() const=0  ;

};

}//end namespace foam


#endif



