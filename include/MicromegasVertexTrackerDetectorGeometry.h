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
      std::array<G4ThreeVector,3> fOffsets = {{{0,0,0},{0,0,0},{0,0,0}}};

   public:

      //DriftChamberSensitiveDetector * fSensitiveDetector;
      //SensitiveRegionDetector       * fSensitiveDetector2;

      G4Material * fFMT_mat     = nullptr;
      G4Material * fMMGas_mat   = nullptr;
      G4Material * fPeek_mat    = nullptr;
      G4Material * fMMMylar_mat = nullptr;
      G4Material * fStrip_mat   = nullptr;
      G4Material * fMesh_mat    = nullptr;
      G4Material * fPCBoard_mat = nullptr;
      G4Material * fEpoxy_mat   = nullptr;

      G4VSolid          * fFMT_solid;
      G4LogicalVolume   * fFMT_log;
      G4VPhysicalVolume * fFMT_phys;
      G4ThreeVector       fFMT_pos;

      G4VSolid          * fFMT_Peek_solid;
      G4LogicalVolume   * fFMT_Peek_log;
      G4ThreeVector       fFMT_Peek_pos;

      G4VSolid          * fFMT_Drift_V_solid;
      G4LogicalVolume   * fFMT_Drift_V_log;
      G4ThreeVector       fFMT_Drift_V_pos;

      G4VSolid          * fFMT_Drift_W_solid;
      G4LogicalVolume   * fFMT_Drift_W_log;
      G4ThreeVector       fFMT_Drift_W_pos;

      G4VSolid          * fFMT_Strips_V_solid;
      G4LogicalVolume   * fFMT_Strips_V_log;
      G4ThreeVector       fFMT_Strips_V_pos;

      G4VSolid          * fFMT_Strips_W_solid;
      G4LogicalVolume   * fFMT_Strips_W_log;
      G4ThreeVector       fFMT_Strips_W_pos;

      G4VSolid          * fFMT_Epoxy_solid;
      G4LogicalVolume   * fFMT_Epoxy_log;
      G4ThreeVector       fFMT_Epoxy_pos;

      G4VSolid          * fFMT_PCBoard_V_solid;
      G4LogicalVolume   * fFMT_PCBoard_V_log;
      G4ThreeVector       fFMT_PCBoard_V_pos;

      G4VSolid          * fFMT_PCBoard_W_solid;
      G4LogicalVolume   * fFMT_PCBoard_W_log;
      G4ThreeVector       fFMT_PCBoard_W_pos;

      G4VSolid          * fFMT_Gas2_V_solid;
      G4LogicalVolume   * fFMT_Gas2_V_log;
      G4ThreeVector       fFMT_Gas2_V_pos;

      G4VSolid          * fFMT_Gas1_V_solid;
      G4LogicalVolume   * fFMT_Gas1_V_log;
      G4ThreeVector       fFMT_Gas1_V_pos;

      G4VSolid          * fFMT_Gas2_W_solid;
      G4LogicalVolume   * fFMT_Gas2_W_log;
      G4ThreeVector       fFMT_Gas2_W_pos;

      G4VSolid          * fFMT_Gas1_W_solid;
      G4LogicalVolume   * fFMT_Gas1_W_log;
      G4ThreeVector       fFMT_Gas1_W_pos;

      G4VSolid          * fFMT_Mesh_V_solid;
      G4LogicalVolume   * fFMT_Mesh_V_log;
      G4ThreeVector       fFMT_Mesh_V_pos;

      G4VSolid          * fFMT_Mesh_W_solid;
      G4LogicalVolume   * fFMT_Mesh_W_log;
      G4ThreeVector       fFMT_Mesh_W_pos;


   public:
      MicromegasVertexTrackerDetectorGeometry();
      ~MicromegasVertexTrackerDetectorGeometry();

      void BuildLogicalVolumes();
      
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother);
      G4VPhysicalVolume * PlaceParallelPhysicalVolume(G4LogicalVolume * mother);

};


#endif

