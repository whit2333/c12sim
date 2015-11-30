#ifndef TorusDetectorGeometry_HH
#define TorusDetectorGeometry_HH 1

#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"

class TorusDetectorGeometry {

   private:

      G4VSolid          * fTorus_solid;
      G4Material        * fTorus_mat;
      G4LogicalVolume   * fTorus_log;
      G4VPhysicalVolume * fTorus_phys;
      G4ThreeVector       fTorus_pos;

   public:
      TorusDetectorGeometry();
      ~TorusDetectorGeometry();

      void BuildLogicalVolumes();
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother);
};


#endif
