#ifndef C12MagneticField_HH
#define C12MagneticField_HH 1

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "ToroidField.h"
#include "SolenoidField.h"

class C12MagneticField : public G4MagneticField {
   private:
      bool   fUseToroid;
      bool   fUseSolenoid;
      double fScaleToroid;
      double fScaleSolenoid;

   protected:
      clas12::mag::ToroidField   fToroidField;
      clas12::mag::SolenoidField fSolenoidField;

   public:
      C12MagneticField(bool use_toroid = false, bool use_solenoid = true, double tor = 1.0, double sol = 1.0 );
      virtual ~C12MagneticField();
      G4ThreeVector GetFieldValue( const G4double Point[4]) const;
      virtual void  GetFieldValue( const G4double Point[4], G4double *Bfield ) const override;

};

#endif
