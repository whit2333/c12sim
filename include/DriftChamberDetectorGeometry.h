#ifndef DriftChamberDetectorGeometry_HH
#define DriftChamberDetectorGeometry_HH

#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Material.hh"
#include "DriftChamberSensitiveDetector.h"
#include "SensitiveRegionDetector.h"

#include "DCWire.h"

class DriftChamberDetectorGeometry {

   private:

      G4VSolid * fHexWireVolume_solid;

      std::array<G4LogicalVolume*,3> fRegions_log;
      std::array<G4LogicalVolume*,3> fEmptyRegions_log;
      std::array<G4VSolid*,3> fRegions_solid;

   public:

      DriftChamberSensitiveDetector * fSensitiveDetector;
      SensitiveRegionDetector       * fSensitiveDetector2;

      G4VSolid * fRegion1_solid;
      G4VSolid * fRegion2_solid;
      G4VSolid * fRegion3_solid;

      G4LogicalVolume * fRegion1_log;
      G4LogicalVolume * fRegion2_log;
      G4LogicalVolume * fRegion3_log;

      G4LogicalVolume * fEmptyRegion1_log;
      G4LogicalVolume * fEmptyRegion2_log;
      G4LogicalVolume * fEmptyRegion3_log;

      G4Material * fGasMaterial;

   public:
      DriftChamberDetectorGeometry();
      ~DriftChamberDetectorGeometry();

      void BuildLogicalVolumes();
      
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother, int sec, int region);
      G4VPhysicalVolume * PlaceParallelPhysicalVolume(G4LogicalVolume * mother, int sec, int region);

};

#endif

