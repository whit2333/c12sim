#ifndef MicromegasVertexTrackerDetectorGeometry_HH
#define MicromegasVertexTrackerDetectorGeometry_HH

#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "SensitiveRegionDetector.h"


class MicromegasVertexTrackerDetectorGeometry {

   private:

      //G4VSolid * fHexWireVolume_solid;

      //std::array<G4LogicalVolume*,3> fRegions_log;
      //std::array<G4LogicalVolume*,3> fEmptyRegions_log;
      //std::array<G4VSolid*,3> fRegions_solid;
      //std::array<G4VSolid*,3> fClippingRegions_solid;

   public:

      //DriftChamberSensitiveDetector * fSensitiveDetector;
      //SensitiveRegionDetector       * fSensitiveDetector2;

      G4Material * fMMGas_mat;
      G4Material * fPeek_mat;
      G4Material * fMMMylar_mat;

      G4VSolid          * fFMT_solid;
      G4LogicalVolume   * fFMT_log;
      G4VPhysicalVolume * fFMT_phys;
      G4Material        * fFMT_mat;
      G4ThreeVector       fFMT_pos;

      G4VSolid          * fFMT_Peek_solid;
      G4LogicalVolume   * fFMT_Peek_log;
      G4VPhysicalVolume * fFMT_Peek_phys;
      G4Material        * fFMT_Peek_mat;
      G4ThreeVector       fFMT_Peek_pos;

      G4VSolid          * fFMT_Drift_V_solid;
      G4LogicalVolume   * fFMT_Drift_V_log;
      G4VPhysicalVolume * fFMT_Drift_V_phys;
      G4Material        * fFMT_Drift_V_mat;
      G4ThreeVector       fFMT_Drift_V_pos;

      G4VSolid          * fFMT_Drift_W_solid;
      G4LogicalVolume   * fFMT_Drift_W_log;
      G4VPhysicalVolume * fFMT_Drift_W_phys;
      G4Material        * fFMT_Drift_W_mat;
      G4ThreeVector       fFMT_Drift_W_pos;

      G4VSolid          * fFMT_Gas2_V_solid;
      G4LogicalVolume   * fFMT_Gas2_V_log;
      G4VPhysicalVolume * fFMT_Gas2_V_phys;
      G4Material        * fFMT_Gas2_V_mat;
      G4ThreeVector       fFMT_Gas2_V_pos;

   public:
      MicromegasVertexTrackerDetectorGeometry();
      ~MicromegasVertexTrackerDetectorGeometry();

      void BuildLogicalVolumes();
      
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother);
      G4VPhysicalVolume * PlaceParallelPhysicalVolume(G4LogicalVolume * mother);

};


#endif

