#include<iostream>
#include<fstream>
#include"math.h"
#include"thermoBase.H"
//using namespace std;

namespace Foam
{

defineTypeNameAndDebug(thermoBase,0);
defineRunTimeSelectionTable(thermoBase, dictionary);

                        

thermoBase::thermoBase(const volVectorField& U, const IOobject& thermPoint)
:
IOdictionary(thermPoint)
{}


} //end namespace foam


