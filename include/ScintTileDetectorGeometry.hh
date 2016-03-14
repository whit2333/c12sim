#ifndef ScintTileDetectorGeometry_HH
#define ScintTileDetectorGeometry_HH

#include "G4SystemOfUnits.hh"
#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4RotationMatrix.hh"

class BeamTestSD;
class SiPMSD;
class ScintTileSensitiveDetector;

/** Recoil Chamber.
 * 
 *  Note: Each "layer" consists of multiple radial wires. Typically
 *  there is only 3 wires radailly with the middle one being a sense wire.
 *  
 */
class ScintTileDetectorGeometry {

   private:

      bool    fUseOpticalPhotons = true;
      double  fRefScint          = 1000.0/CLHEP::MeV;

      std::array<double, 33>   fScintLightOutput;
      std::array<double, 33>   fScintLightOutput_slow;
      std::array<double, 33>   fScintWavelengths;
      std::array<double, 33>   fScintPhotonEnergies;
      std::array<double, 33>   fScintRefIndex;

      const int fNScintPerStrip         = 10;
      const int fNScint                 = 60;
      std::array<G4ThreeVector     ,      600> fScint_positions;
      std::array<G4RotationMatrix  ,      600> fScint_rotations;
      std::array<G4VPhysicalVolume*,      600> fScint_physicals;
      std::array<G4LogicalBorderSurface*, 600> fScint_borders;


      G4ThreeVector       fSiPM_pos;
      G4RotationMatrix    fSiPM_rot ;
      G4Material        * fSiPM_mat   = nullptr;
      G4VSolid          * fSiPM_solid = nullptr;
      G4LogicalVolume   * fSiPM_log   = nullptr;
      G4VPhysicalVolume * fSiPM_phys  = nullptr;
      SiPMSD            * fSiPM_det   = nullptr;

      G4ThreeVector       fScint_pos;
      G4Material        * fScint_mat   = nullptr;
      G4VSolid          * fScint_solid = nullptr;
      G4LogicalVolume   * fScint_log   = nullptr;
      G4VPhysicalVolume * fScint_phys  = nullptr;
      ScintTileSensitiveDetector * fScint_det = nullptr;

      G4ThreeVector       fWLSFiber_pos ;
      G4RotationMatrix    fWLSFiber_rot ;
      G4Material        * fWLSFiber_mat   = nullptr;
      G4VSolid          * fWLSFiber_solid = nullptr;
      G4LogicalVolume   * fWLSFiber_log   = nullptr;
      std::array<G4VPhysicalVolume*, 4>      fWLSFibers_phys;
      G4VPhysicalVolume * fWLSFiber_phys  = nullptr;

      G4ThreeVector       fCurvedWLSF_pos ;
      G4RotationMatrix    fCurvedWLSF_rot ;
      G4Material        * fCurvedWLSF_mat   = nullptr;
      G4VSolid          * fCurvedWLSF_solid = nullptr;
      G4LogicalVolume   * fCurvedWLSF_log   = nullptr;
      std::array<G4VPhysicalVolume*, 4>      fCurvedWLSFs_phys;
      G4VPhysicalVolume * fCurvedWLSF_phys  = nullptr;

      G4ThreeVector       fCurvedFiber_pos ;
      G4RotationMatrix    fCurvedFiber_rot ;
      G4Material        * fCurvedFiber_mat   = nullptr;
      G4VSolid          * fCurvedFiber_solid = nullptr;
      G4LogicalVolume   * fCurvedFiber_log   = nullptr;
      std::array<G4VPhysicalVolume*, 4>      fCurvedFibers_phys;
      G4VPhysicalVolume * fCurvedFiber_phys  = nullptr;

      //G4ThreeVector       fPhotonDet1_pos ;
      //G4Material        * fPhotonDet1_mat   = nullptr;
      //G4VSolid          * fPhotonDet1_solid = nullptr;
      //G4LogicalVolume   * fPhotonDet1_log   = nullptr;
      //G4VPhysicalVolume * fPhotonDet1_phys  = nullptr;
      //BeamTestSD        * fPhotonDet1_det   = nullptr;

      //G4ThreeVector       fPhotonDet2_pos;
      //G4Material        * fPhotonDet2_mat   = nullptr;
      //G4VSolid          * fPhotonDet2_solid = nullptr;
      //G4LogicalVolume   * fPhotonDet2_log   = nullptr;
      //G4VPhysicalVolume * fPhotonDet2_phys  = nullptr;
      //BeamTestSD        * fPhotonDet2_det   = nullptr;

   public:

      double fRadius             = 9.2*CLHEP::cm;
      double fScintGap           = 2.54*0.002*CLHEP::cm; 
      double fScintWidth         = 10.0*CLHEP::mm;
      double fScintLength        = 30.0*CLHEP::mm;
      double fScintThickness     = 10.0*CLHEP::mm;

      double fPhotonDetThickness = 25.4*0.00002*CLHEP::mm;

      double fFiberCurveRadius   = 15.0*CLHEP::cm;

      double fWLSFiberDiameter   = 1.0*CLHEP::mm;
      double fWLSFiberLength     = 4.0*CLHEP::cm;

      double fWLSFiberEmbedAngle = 0.0;

      double fSiPMThickness      = 1.0*CLHEP::mm;

   public:
      ScintTileDetectorGeometry();
      ~ScintTileDetectorGeometry();

      void BuildMaterials();
      void BuildLogicalVolumes();
      G4VPhysicalVolume * PlacePhysicalVolume(
            G4LogicalVolume   * mother, 
            G4VPhysicalVolume * mother_phys,
            G4VPhysicalVolume * adjacent_phys = 0 );


};

#endif

