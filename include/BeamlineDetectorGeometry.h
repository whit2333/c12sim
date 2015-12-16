#ifndef BeamlineDetectorGeometry_HH
#define BeamlineDetectorGeometry_HH 1

#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

class BeamlineDetectorGeometry {

   protected:

      G4VSolid          * fMollerShieldTube_solid;
      G4Material        * fMollerShieldTube_mat;
      G4LogicalVolume   * fMollerShieldTube_log;
      G4VPhysicalVolume * fMollerShieldTube_phys;
      G4ThreeVector       fMollerShieldTube_pos;
      G4VisAttributes   * fMollerShieldTube_vis;

      G4VSolid          * fMollerShieldTubeW_solid;
      G4Material        * fMollerShieldTubeW_mat;
      G4LogicalVolume   * fMollerShieldTubeW_log;
      G4VPhysicalVolume * fMollerShieldTubeW_phys;
      G4ThreeVector       fMollerShieldTubeW_pos;
      G4VisAttributes   * fMollerShieldTubeW_vis;

      G4VSolid          * fMollerShieldTubeA_solid;
      G4Material        * fMollerShieldTubeA_mat;
      G4LogicalVolume   * fMollerShieldTubeA_log;
      G4VPhysicalVolume * fMollerShieldTubeA_phys;
      G4ThreeVector       fMollerShieldTubeA_pos;
      G4VisAttributes   * fMollerShieldTubeA_vis;

      G4VSolid          * fMollerShieldTubeV_solid;
      G4Material        * fMollerShieldTubeV_mat;
      G4LogicalVolume   * fMollerShieldTubeV_log;
      G4VPhysicalVolume * fMollerShieldTubeV_phys;
      G4ThreeVector       fMollerShieldTubeV_pos;
      G4VisAttributes   * fMollerShieldTubeV_vis;

      G4VSolid          * fMollerShieldCone_solid;
      G4Material        * fMollerShieldCone_mat;
      G4LogicalVolume   * fMollerShieldCone_log;
      G4VPhysicalVolume * fMollerShieldCone_phys;
      G4ThreeVector       fMollerShieldCone_pos;
      G4VisAttributes   * fMollerShieldCone_vis;

      G4VSolid          * fMollerShieldConeA_solid;
      G4Material        * fMollerShieldConeA_mat;
      G4LogicalVolume   * fMollerShieldConeA_log;
      G4VPhysicalVolume * fMollerShieldConeA_phys;
      G4ThreeVector       fMollerShieldConeA_pos;
      G4VisAttributes   * fMollerShieldConeA_vis;

      G4VSolid          * fMollerShieldConeV_solid;
      G4Material        * fMollerShieldConeV_mat;
      G4LogicalVolume   * fMollerShieldConeV_log;
      G4VPhysicalVolume * fMollerShieldConeV_phys;
      G4ThreeVector       fMollerShieldConeV_pos;
      G4VisAttributes   * fMollerShieldConeV_vis;

      G4VSolid          * fBeamline_solid;
      G4Material        * fBeamline_mat;
      G4LogicalVolume   * fBeamline_log;
      G4VPhysicalVolume * fBeamline_phys;
      G4ThreeVector       fBeamline_pos;
      G4VisAttributes   * fBeamline_vis;

      G4VSolid          * fBeamlinePipe_solid;
      G4Material        * fBeamlinePipe_mat;
      G4LogicalVolume   * fBeamlinePipe_log;
      G4VPhysicalVolume * fBeamlinePipe_phys;
      G4ThreeVector       fBeamlinePipe_pos;
      G4VisAttributes   * fBeamlinePipe_vis;

      G4VSolid          * fBeamlinePipe2_solid;
      G4Material        * fBeamlinePipe2_mat;
      G4LogicalVolume   * fBeamlinePipe2_log;
      G4VPhysicalVolume * fBeamlinePipe2_phys;
      G4ThreeVector       fBeamlinePipe2_pos;
      G4VisAttributes   * fBeamlinePipe2_vis;

      G4VSolid          * fBeamlinePipe3_solid;
      G4Material        * fBeamlinePipe3_mat;
      G4LogicalVolume   * fBeamlinePipe3_log;
      G4VPhysicalVolume * fBeamlinePipe3_phys;
      G4ThreeVector       fBeamlinePipe3_pos;
      G4VisAttributes   * fBeamlinePipe3_vis;

   public:
      BeamlineDetectorGeometry();
      ~BeamlineDetectorGeometry();

      void BuildLogicalVolumes();
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother);
};


#endif
