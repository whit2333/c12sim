#ifndef DriftChamberDetectorGeometry_HH
#define DriftChamberDetectorGeometry_HH

#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Region.hh"
#include "DriftChamberSensitiveDetector.h"
#include "SensitiveRegionDetector.h"

#include "DCWire.h"

class DriftChamberDetectorGeometry {

   private:

      G4VSolid * fHexWireVolume_solid;

      std::array<G4VSolid*,3>        fRegions_solid;
      std::array<G4LogicalVolume*,3> fRegions_log;

      std::array<G4VSolid*,3>        fRegionContainers_solid;
      std::array<G4LogicalVolume*,3> fRegionContainers_log;

      //std::array<G4LogicalVolume*,3> fEmptyRegions_log;
      //std::array<G4VSolid*,3>        fClippingRegions_solid;

      std::array<G4VSolid*,3>        fEndPlates_solid;
      std::array<G4LogicalVolume*,3> fEndPlates_log;
      std::array<G4ThreeVector,3>    fLeftEndPlates_pos;
      std::array<G4ThreeVector,3>    fRightEndPlates_pos;

      std::array< G4RotationMatrix,6>   fLeftEndPlates_rot;
      std::array< G4RotationMatrix,6>   fRightEndPlates_rot;

      G4Region * fRegion_g4Region = nullptr;

   public:

      DriftChamberSensitiveDetector * fSensitiveDetector;
      SensitiveRegionDetector       * fSensitiveDetector2;

      G4Material * fGasMaterial;

      bool   fBuildWireUnitCells = true;

   public:
      DriftChamberDetectorGeometry();
      ~DriftChamberDetectorGeometry();

      void BuildLogicalVolumes();
      
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother, int sec, int region);
      G4VPhysicalVolume * PlaceParallelPhysicalVolume(G4LogicalVolume * mother, int sec, int region);

};

#endif

