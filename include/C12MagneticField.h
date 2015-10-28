#ifndef C12MagneticField_HH
#define C12MagneticField_HH 1

#include "G4MagneticField.hh"
#include "ToroidField.h"
#include "SolenoidField.h"

class C12MagneticField : public G4MagneticField {

   protected:
      clas12::mag::ToroidField   fToroidField;
      clas12::mag::SolenoidField fSolenoidField;

      C12MagneticField();
      virtual ~C12MagneticField();
      virtual void  GetFieldValue( const G4double Point[4], G4double *Bfield ) const override;

};

#endif
