#ifndef RecoilHodoDetectorGeometry_HH
#define RecoilHodoDetectorGeometry_HH

#include "G4SystemOfUnits.hh"
#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"

class BeamTestSD;

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

      
      G4ThreeVector       fRadiator_pos;
      G4Material        * fRadiator_mat   = nullptr;
      G4VSolid          * fRadiator_solid = nullptr;
      G4LogicalVolume   * fRadiator_log   = nullptr;
      //G4VPhysicalVolume * fRadiator_phys  = nullptr;

      G4ThreeVector       fPhotonDet1_pos ;
      G4Material        * fPhotonDet1_mat   = nullptr;
      G4VSolid          * fPhotonDet1_solid = nullptr;
      G4LogicalVolume   * fPhotonDet1_log   = nullptr;
      G4VPhysicalVolume * fPhotonDet1_phys  = nullptr;
      BeamTestSD        * fPhotonDet1_det   = nullptr;

      G4ThreeVector       fPhotonDet2_pos;
      G4Material        * fPhotonDet2_mat   = nullptr;
      G4VSolid          * fPhotonDet2_solid = nullptr;
      G4LogicalVolume   * fPhotonDet2_log   = nullptr;
      G4VPhysicalVolume * fPhotonDet2_phys  = nullptr;
      BeamTestSD        * fPhotonDet2_det   = nullptr;

   public:

      double fScintWrapThickness =  2.54*0.004*CLHEP::cm; 
      double fScintGap           =  2.54*0.002*CLHEP::cm; 
      double fScintWidth         =  1.0*CLHEP::cm;
      double fScintLength        = 40.0*CLHEP::cm;
      double fScintThickness     =  3.2*CLHEP::cm;
      double ffPhotonDetThickness =  25.4*0.00002*CLHEP::mm;

      double fInnerRadius        = 20*CLHEP::cm;
      double fScintDeltaTheta    = 6.0*CLHEP::degree; // angle subtended by bar 

      const int fNScint     = 60;
      std::array<G4ThreeVector     , 60>      fScint_positions;
      std::array<G4VPhysicalVolume*, 60>      fScint_physicals;
      std::array<G4LogicalBorderSurface*, 60> fScint_borders;

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

