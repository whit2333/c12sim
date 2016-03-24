#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "DriftChamberDetectorGeometry.h"
#include "DriftChamberSensitiveDetector.h"
#include "ECDetectorGeometry.h"
#include "RecoilChamberDetectorGeometry.h"
#include "RecoilHodoDetectorGeometry.h"
#include "RecoilHodoDetectorGeometry2.h"
#include "RecoilHodoDetectorGeometry3.h"
#include "HTCCDetectorGeometry.h"
#include "SolenoidDetectorGeometry.h"
#include "TorusDetectorGeometry.h"
#include "BeamlineDetectorGeometry.h"
#include "MicromegasVertexTrackerDetectorGeometry.h"
#include "SiliconVertexTrackerDetectorGeometry.h"
using MVTDetectorGeometry = MicromegasVertexTrackerDetectorGeometry;
using SVTDetectorGeometry = SiliconVertexTrackerDetectorGeometry;

#include "G4SystemOfUnits.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class B1DetectorMessenger;
class G4Material;
class G4VSolid;
class FakeSD;
class G4UserLimits;
#include "G4ThreeVector.hh"
#include "G4String.hh"

/// Detector construction class to define materials and geometry.

class B1DetectorConstruction : public G4VUserDetectorConstruction
{
   private:

      G4double world_x                       ;
      G4double world_y                       ;
      G4double world_z                       ;

      //G4double radiator_thickness            ;
      //G4double collimator_target_center_gap  ;
      //G4double outer_collimator_ID           ;
      //G4double outer_collimator_OD           ;
      //G4double collimator_ID                 ;
      //G4double collimator_OD                 ;
      //G4double collimator_diameter           ;
      //G4double collimator_z_end              ;
      //G4double radiator_collimator_gap ;
      //G4double collimator_length   ;
      //G4double beampipe_length   ;
      //G4double beampipe_diameter ;
      //G4double radiator_diameter ;
      //G4double scoring_diameter  ;
      //G4double scoring_length    ;
      //G4double window_diameter   ;
      //G4double window_thickness  ;
      //G4double scoring2_diameter ;
      //G4double scoring2_length   ;


      //--Parameters of the magnetic field
      G4double   fieldValue = 5*tesla;
      G4double   acc_inter = 0.1*mm;
      //--


      //--Parameters of the target
      G4double   innerRadiusOfTheTarget = 0.0*mm;
      G4double   outerRadiusOfTheTarget = 6.*mm;
      G4double   fTargetLength = 20.0*cm;//5.0*cm;
      G4double   startAngleOfTheTarget = 0.*deg;
      G4double   spanningAngleOfTheTarget = 360.*deg;
      G4double   target_posX = 0.*mm;
      G4double   target_posY = 0.*mm;
      G4double   target_posZ = 0.*mm;

      //--Kapton foil around the target
      G4double      innerRadiusOfTheKapton = 6.0*mm;
      G4double      outerRadiusOfTheKapton = 6.028*mm;
      G4double      hightOfTheKapton = 200.*mm;
      G4double      startAngleOfTheKapton = 0.*deg;
      G4double      spanningAngleOfTheKapton = 360.*deg;

      //--Parameters of space around the target
      G4double  innerRadiusOfAround = 6.028*mm;
      G4double  outerRadiusOfAround = 29.990*mm;
      G4double  hightOfAround = 200.*mm;
      G4double  startAngleOfAround = 0.*deg;
      G4double spanningAngleOfAround = 360.*deg;

      //--Mylar foil around the clear space
      G4double   innerRadiusOfTheOclKapton = 29.996*mm;
      G4double   outerRadiusOfTheOclKapton = 30.000*mm;
      G4double   hightOfTheOclKapton = 200.*mm;
      G4double   startAngleOfTheOclKapton = 0.*deg;
      G4double   spanningAngleOfTheOclKapton = 360.*deg;
      //--


      //--Mylar foil around the gas detector
      G4double     innerRadiusOfTheOKapton = 79.995*mm;
      G4double     outerRadiusOfTheOKapton = 80.0*mm;
      G4double     hightOfTheOKapton = 200.*mm;
      G4double     startAngleOfTheOKapton = 0.*deg;
      G4double     spanningAngleOfTheOKapton = 360.*deg;
      //--

      //--Parameters of the scintillators
      G4double  innerRadiusOfTheSiDetector = 80.*mm;
      G4double  outerRadiusOfTheSiDetector = 140.*mm;
      G4double  hightOfTheSiDetector = 200.*mm;
      G4double  startAngleOfTheSiDetector = 0.*deg;
      G4double  spanningAngleOfTheSiDetector = 360.*deg;
      G4double  SiDetector_posX = 0.*mm;
      G4double  SiDetector_posY = 0.*mm;
      G4double  SiDetector_posZ = 0.*mm; 

      G4Material* Air;
      G4Material* Ar;
      G4Material* Silicon;
      G4Material* Scinti;
      G4Material* Lead;
      G4Material* Kapton;
      G4Material* He_Target;
      G4Material* He_ClearS;
      G4Material* Xe_varPT;
      G4Material* isobutane;
      G4Material* Vacuum;
      G4Material* Deuterium;
      G4Material* CO2;
      G4Material* He10CO2;
      G4Material* HeiC4H10;
      G4Material* Tungsten; 
      G4Material* Mylar;
      G4Material* LH2;

      G4Element* H;
      G4Element* C;
      G4Element* N;
      G4Element* O;
      G4Element* elHe;
      G4Element* elXe;
      G4Element* ele_D;
      G4Element* elW;
      G4UserLimits* fStepLimit; // pointer to user step limits

   protected:

      G4LogicalVolume     * fScoringVolume;
      B1DetectorMessenger * fMessenger;
      G4String    fCollimatorMatName;

      bool fHasBeenBuilt;


   private:

      G4Material        * world_mat   ;
      G4VSolid          * world_solid ;
      G4LogicalVolume   * world_log   ;
      G4VPhysicalVolume * world_phys  ;

      HTCCDetectorGeometry           * fHTCC          = nullptr;
      RecoilChamberDetectorGeometry  * fRecoilChamber = nullptr;
      RecoilHodoDetectorGeometry     * fRecoilHodo    = nullptr;
      RecoilHodoDetectorGeometry2    * fRecoilHodo2   = nullptr;
      RecoilHodoDetectorGeometry3    * fRecoilHodo3   = nullptr;
      DriftChamberDetectorGeometry   * fDriftChamber  = nullptr;
      SolenoidDetectorGeometry       * fSolenoid      = nullptr;
      TorusDetectorGeometry          * fTorus         = nullptr;
      BeamlineDetectorGeometry       * fBeamline      = nullptr;
      MVTDetectorGeometry            * fMVT           = nullptr;
      SVTDetectorGeometry            * fSVT           = nullptr;
      ECDetectorGeometry             * fEC            = nullptr;

   protected:
      void DefineMaterials();

   public:
      B1DetectorConstruction();
      virtual ~B1DetectorConstruction();

      virtual G4VPhysicalVolume* Construct();

      void ToggleCherenkov(bool l);

      void PrintConfigInfo() const;

      void CalculatePositions();

      void Rebuild();

      //G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

   protected:
};
//______________________________________________________________________________

#endif

