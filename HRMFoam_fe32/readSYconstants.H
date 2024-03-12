Info<< "Reading SigmaY model constants\n" << endl;

IOdictionary sigmayConstants
(
 IOobject
 (
  "sigmayConstants",
  runTime.constant(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::NO_WRITE
  )
);

dimensionedScalar Sc("Sc",dimless,VGREAT);
if (sigmayConstants.found("Sc")) 
{
  Sc = dimensionedScalar(sigmayConstants.lookup("Sc"));
}

Info<< "Turbulent Schmidt number = "<< Sc.value() << endl;

Switch solveSigma("no");
solveSigma.readIfPresent("solveSigma", sigmayConstants);

dimensionedScalar Cmu("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0.09);
if (sigmayConstants.found("Cmu")) 
  {
    Cmu = dimensionedScalar(sigmayConstants.lookup("Cmu"));
  }

dimensionedScalar sigmaST("zero", dimensionSet(1, 0, -2, 0, 0, 0, 0) ,2.89e-2);
if (sigmayConstants.found("surfaceTension")) 
  {
    sigmaST = dimensionedScalar(sigmayConstants.lookup("surfaceTension"));
  }

    dimensionedScalar C4("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,5.5);
if (sigmayConstants.found("C4")) 
  {
    C4 = dimensionedScalar(sigmayConstants.lookup("C4"));
  }

    dimensionedScalar C5("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,1.0);
if (sigmayConstants.found("C5")) 
  {
    C5 = dimensionedScalar(sigmayConstants.lookup("C5"));
  }

    dimensionedScalar C6("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,1.0);
if (sigmayConstants.found("C6")) 
  {
    C6 = dimensionedScalar(sigmayConstants.lookup("C6"));
  }

// To make HRMFoam backwards compatable and avoid redeclaration of 
// Sigma constants in every iteration, Sigma constants will be declared here.
// The read / write options will be conditional so that non-sigma runs don't see
// Sigma fields.

int conditionalWrite = 1;

if (solveSigma)
{
  conditionalWrite = 0;
}

/*volScalarField Sigma
(
  IOobject
  (
    "Sigma",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::writeOption(conditionalWrite)
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0.0)
);*/


volScalarField Omega
(
  IOobject
  (
    "Omega",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::writeOption(conditionalWrite)
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(1, -2, 0, 0, 0, 0, 0) ,0.0)
);

volScalarField Ybar
(
  IOobject
  (
    "Ybar",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::writeOption(conditionalWrite)
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0.0)
);

volScalarField lt
(
  IOobject
  (
    "lt",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 1, 0, 0, 0, 0, 0) ,0.0)
);

volScalarField d32
(
  IOobject
  (
    "d32",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::writeOption(conditionalWrite)
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 1, 0, 0, 0, 0, 0) ,0.0),
  zeroGradientFvPatchScalarField::typeName
);

volScalarField A
(
  IOobject
  (
    "A",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0) ,0.0)
);

volScalarField a
(
  IOobject
  (
    "a",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0) ,0.0)
);

    /*    volScalarField rEq
    (
        IOobject
        (
            "rEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
            ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 1, 0, 0, 0, 0, 0) ,0.0)
    );*/

volScalarField sEq
(
  IOobject
  (
    "sEq",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0, 0, 0) ,0.0)
);

volScalarField sMin
(
  IOobject
  (
    "sMin",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0, 0, 0) ,0.0)
);

    volScalarField sInit
    (
        IOobject
        (
            "sInit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -1, -1, 0, 0, 0, 0), 0.0)
	);

volScalarField Sigma
(
 IOobject
 (
  "Sigma",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::writeOption(conditionalWrite)
 ),
 mesh,
 dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0)
 );
  
