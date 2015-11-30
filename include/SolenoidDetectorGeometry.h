#ifndef SolenoidDetectorGeometry_HH
#define SolenoidDetectorGeometry_HH 1

#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"

class SolenoidDetectorGeometry {

   private:

      G4VSolid          * fSolenoid_solid;
      G4Material        * fSolenoid_mat;
      G4LogicalVolume   * fSolenoid_log;
      G4VPhysicalVolume * fSolenoid_phys;
      G4ThreeVector       fSolenoid_pos;

   public:
      SolenoidDetectorGeometry();
      ~SolenoidDetectorGeometry();

      void BuildLogicalVolumes();
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother);
};


#endif
