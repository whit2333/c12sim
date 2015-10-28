#include "C12MagneticField.h"
#include "TVector3.h"

C12MagneticField::C12MagneticField()
{ }
//______________________________________________________________________________

C12MagneticField::~C12MagneticField()
{ }
//______________________________________________________________________________

void  C12MagneticField::GetFieldValue( const G4double Point[4], G4double *Bfield ) const
{
   TVector3 pos(Point[0],Point[1],Point[2]);
   TVector3 B = fToroidField.GetField(pos) + fSolenoidField.GetField(pos);
   Bfield[0] = B.x();
   Bfield[1] = B.y();
   Bfield[2] = B.z();
}
//______________________________________________________________________________

