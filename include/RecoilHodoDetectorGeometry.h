#ifndef RecoilHodoDetectorGeometry_HH
#define RecoilHodoDetectorGeometry_HH

#include "G4SystemOfUnits.hh"
#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"

class SiPMSD;
class RecoilScintSensitiveDetector;

/** Recoil Chamber.
 * 
 *  Note: Each "layer" consists of multiple radial wires. Typically
 *  there is only 3 wires radailly with the middle one being a sense wire.
 *  
 */
class RecoilHodoDetectorGeometry {

   private:

      std::array<double, 24>   fScintLightOutput;
      std::array<double, 24>   fWavelengths;
      std::array<double, 24>   fPhotonEnergies;
      std::array<double, 24>   fRefIndex;

      RecoilScintSensitiveDetector * fScint1_det = nullptr;
      RecoilScintSensitiveDetector * fScint2_det = nullptr;
      
      bool    fUseOpticalPhotons = true;
      double  fRefScint          = 1000.0/CLHEP::MeV;

      G4ThreeVector       fScint1_pos;
      G4Material        * fScint1_mat   = nullptr;
      G4VSolid          * fScint1_solid = nullptr;
      G4LogicalVolume   * fScint1_log   = nullptr;

      G4ThreeVector       fScint2_pos;
      G4Material        * fScint2_mat   = nullptr;
      G4VSolid          * fScint2_solid = nullptr;
      G4LogicalVolume   * fScint2_log   = nullptr;

      G4ThreeVector       fPhotonDet1_pos ;
      G4Material        * fPhotonDet1_mat   = nullptr;
      G4VSolid          * fPhotonDet1_solid = nullptr;
      G4LogicalVolume   * fPhotonDet1_log   = nullptr;
      G4VPhysicalVolume * fPhotonDet1_phys  = nullptr;
      SiPMSD            * fPhotonDet1_det   = nullptr;

      G4ThreeVector       fPhotonDet2_pos;
      G4Material        * fPhotonDet2_mat   = nullptr;
      G4VSolid          * fPhotonDet2_solid = nullptr;
      G4LogicalVolume   * fPhotonDet2_log   = nullptr;
      G4VPhysicalVolume * fPhotonDet2_phys  = nullptr;
      SiPMSD            * fPhotonDet2_det   = nullptr;

      G4ThreeVector       fPhotonDet3_pos ;
      G4Material        * fPhotonDet3_mat   = nullptr;
      G4VSolid          * fPhotonDet3_solid = nullptr;
      G4LogicalVolume   * fPhotonDet3_log   = nullptr;
      G4VPhysicalVolume * fPhotonDet3_phys  = nullptr;
      SiPMSD            * fPhotonDet3_det   = nullptr;

      G4ThreeVector       fPhotonDet4_pos;
      G4Material        * fPhotonDet4_mat   = nullptr;
      G4VSolid          * fPhotonDet4_solid = nullptr;
      G4LogicalVolume   * fPhotonDet4_log   = nullptr;
      G4VPhysicalVolume * fPhotonDet4_phys  = nullptr;
      SiPMSD            * fPhotonDet4_det   = nullptr;

      const int fNScint                 = 60;
      std::array<G4ThreeVector          , 60> fScint1_positions;
      std::array<G4VPhysicalVolume*     , 60> fScint1_physicals;
      std::array<G4LogicalBorderSurface*, 60> fScint1_borders;

      std::array<G4ThreeVector          , 60> fScint2_positions;
      std::array<G4VPhysicalVolume*     , 60> fScint2_physicals;
      std::array<G4LogicalBorderSurface*, 60> fScint2_borders;

   public:

      double fScintWrapThickness  =  2.54*0.004*CLHEP::cm; 
      double fScintGap            =  2.54*0.002*CLHEP::cm; 

      double fScintLength         = 30.0*CLHEP::cm;

      double fScint1Thickness      =  2.0*CLHEP::mm;
      double fScint2Thickness      = 50.0*CLHEP::mm;

      double fPhotonDetThickness =  25.4*0.00002*CLHEP::mm;

      double fInnerRadius        = 8.0*CLHEP::cm;
      double fScintDeltaTheta    = 6.0*CLHEP::degree; // angle subtended by bar 


   public:
      RecoilHodoDetectorGeometry();
      ~RecoilHodoDetectorGeometry();

      void BuildMaterials();
      void BuildLogicalVolumes();
      G4VPhysicalVolume * PlacePhysicalVolume(
            G4LogicalVolume   * mother, 
            G4VPhysicalVolume * mother_phys,
            G4VPhysicalVolume * adjacent_phys = 0 );


};

#endif

