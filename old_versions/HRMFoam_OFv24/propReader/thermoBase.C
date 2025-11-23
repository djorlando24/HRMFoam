#include"thermoBase.H"

namespace Foam
{

defineTypeNameAndDebug(thermoBase,0);
defineRunTimeSelectionTable(thermoBase, dictionary);



thermoBase::thermoBase(const volVectorField& U, const IOobject& thermPoint)
:
IOdictionary(thermPoint)
{}


} //end namespace foam


