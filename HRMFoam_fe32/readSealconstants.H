Info<< "Reading Seal constants\n" << endl;



IOdictionary sealConstants//Reading Sealing constants, if the file is present---CM
(
 IOobject
 (
  "sealConstants",
  runTime.constant(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::NO_WRITE
  )
);

//Set default values for sealing constants----CM
dimensionedScalar drag("drag",dimensionSet(0,0,1,0,0,0,0),VGREAT);//drag constant
dimensionedScalar slRelaxation("slRelaxation",dimensionSet(0,0,1,0,0,0,0),VGREAT);//time relaxation factor
point lowerCorner(VSMALL,VSMALL,VSMALL);//bounding boxlowecorner
point upperCorner(SMALL,SMALL,SMALL);//bounding box upper
dimensionedScalar triggerDisp("triggerDisp",dimensionSet(0,0,0,0,0,0,0),VGREAT); //lift to toggle sealing 

//If the seal constants file is found, set up constants
bool haveSeal;
if (sealConstants.headerOk()) {
  haveSeal = true;
  Info << "Seal constants dictionary found.  Continuing sealing setup." << endl;

  slRelaxation = dimensionedScalar(sealConstants.lookup("slRelaxation"));
  Info<< "Sealing time relaxation factor = "<< slRelaxation.value() << endl;

  drag = dimensionedScalar(sealConstants.lookup("drag"));
  Info<< "Drag constant for sealing = "<< drag.value() <<endl;

  lowerCorner = point(sealConstants.lookup("lowerCorner"));
  Info<< "Lower corner coordinate for sealing bounding box: "<<lowerCorner <<endl;

  upperCorner = point(sealConstants.lookup("upperCorner"));
  Info<< "Upper corner coordinate for sealing bounding box: "<<upperCorner <<endl;

  triggerDisp = dimensionedScalar(sealConstants.lookup("triggerDisp"));
  Info<< "Trigger Mesh Displacement = "<<triggerDisp.value() <<endl;

}
else {
  haveSeal = false;
  Info << "No seal constants dictionary found.  Sealing deactivated." << endl;//deactivating sealing
}
