#include "C12MagneticField.h"
#include "TVector3.h"
#include "CLHEP/Units/SystemOfUnits.h"

C12MagneticField::C12MagneticField(bool use_toroid, bool use_solenoid  ) : 
   fUseToroid(use_toroid), fUseSolenoid(use_solenoid), 
   fScaleToroid(1.0), fScaleSolenoid(1.0)
{
   if(fUseSolenoid) fSolenoidField.ReadMap();
   if(fUseToroid)   fToroidField.ReadMap();
}
//______________________________________________________________________________

C12MagneticField::~C12MagneticField()
{ }
//______________________________________________________________________________

G4ThreeVector C12MagneticField::GetFieldValue( const G4double Point[4]) const
{
   using namespace CLHEP;
   TVector3 pos(Point[0]/cm, Point[1]/cm, Point[2]/cm);
   TVector3 B = {0.0,0.0,0.0};
   if(fUseSolenoid) B += fSolenoidField.GetField(pos);
   if(fUseToroid) B += fToroidField.GetField(pos);
   return G4ThreeVector(B.x()*tesla, B.y()*tesla, B.z()*tesla);
}
//______________________________________________________________________________

void  C12MagneticField::GetFieldValue( const G4double Point[4], G4double *Bfield ) const
{
   using namespace CLHEP;
   G4ThreeVector B = GetFieldValue(Point);
   Bfield[0] = B.x()*tesla;
   Bfield[1] = B.y()*tesla;
   Bfield[2] = B.z()*tesla;
}
//______________________________________________________________________________
