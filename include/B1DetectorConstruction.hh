#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "DriftChamberDetectorGeometry.h"
#include "DriftChamberSensitiveDetector.h"
#include "RecoilChamberDetectorGeometry.h"
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
      G4double radiator_thickness            ;
      G4double collimator_target_center_gap  ;
      G4double outer_collimator_ID           ;
      G4double outer_collimator_OD           ;
      G4double collimator_ID                 ;
      G4double collimator_OD                 ;
      G4double collimator_diameter           ;
      G4double collimator_z_end              ;
      G4double radiator_collimator_gap ;
      G4double collimator_length   ;
      G4double beampipe_length   ;
      G4double beampipe_diameter ;
      G4double radiator_diameter ;
      G4double scoring_diameter  ;
      G4double scoring_length    ;
      G4double window_diameter   ;
      G4double window_thickness  ;
      G4double scoring2_diameter ;
      G4double scoring2_length   ;


      //--Parameters of the magnetic field
      G4double   fieldValue = 5*tesla;
      G4double   acc_inter = 0.1*mm;
      //--


      //--Parameters of the target
      G4double   innerRadiusOfTheTarget = 0.0*mm;
      G4double   outerRadiusOfTheTarget = 6.*mm;
      G4double   hightOfTheTarget = 200.*mm;
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
      //--
      //--


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

      G4LogicalVolume   * drift_chamber_log   ;
      G4VPhysicalVolume * drift_chamber_phys  ;
      G4LogicalVolume   * htcc_log   ;
      G4VPhysicalVolume * htcc_phys  ;

      FakeSD * scoring_det;
      FakeSD * scoring2_det;

      G4Material        * world_mat   ;
      G4VSolid          * world_solid ;
      G4LogicalVolume   * world_log   ;
      G4VPhysicalVolume * world_phys  ;

      G4ThreeVector       beampipe_pos;
      G4Material        * beampipe_mat   ;
      G4VSolid          * beampipe_solid ;
      G4LogicalVolume   * beampipe_log   ;
      G4VPhysicalVolume * beampipe_phys  ;
      G4ThreeVector       radiator_pos;
      G4Material        * radiator_mat  ;
      G4VSolid          * radiator_solid;
      G4LogicalVolume   * radiator_log  ;
      G4VPhysicalVolume * radiator_phys ;
      G4ThreeVector       collimator_pos;
      G4Material        * collimator_mat ;  
      G4VSolid          * collimator_solid ;
      G4LogicalVolume   * collimator_log ;  
      G4VPhysicalVolume * collimator_phys;  
      G4ThreeVector       collimator2_pos;
      G4Material        * collimator2_mat ;  
      G4VSolid          * collimator2_solid ;
      G4LogicalVolume   * collimator2_log ;  
      G4VPhysicalVolume * collimator2_phys;  
      G4ThreeVector       outer_collimator_pos;
      G4Material        * outer_collimator_mat ;  
      G4VSolid          * outer_collimator_solid ;
      G4LogicalVolume   * outer_collimator_log ;  
      G4VPhysicalVolume * outer_collimator_phys;  
      G4ThreeVector       scoring_pos;
      G4Material        * scoring_mat   ;
      G4VSolid          * scoring_solid ;
      G4LogicalVolume   * scoring_log   ;
      G4VPhysicalVolume * scoring_phys  ;
      G4ThreeVector       window_pos;
      G4Material        * window_mat  ;
      G4VSolid          * window_solid;
      G4LogicalVolume   * window_log  ;
      G4VPhysicalVolume * window_phys ;
      G4ThreeVector       scoring2_pos;
      G4Material        * scoring2_mat  ;
      G4VSolid          * scoring2_solid;
      G4LogicalVolume   * scoring2_log  ;
      G4VPhysicalVolume * scoring2_phys ;

   protected:

      //int nCherenkovPMT;
      //PMTArraySensitiveDetector * fCherenkov_PMTArray;

      //int nDriftChamberVols;
      //DriftChamberSensitiveDetector * fDriftChamber_SD;

      RecoilChamberDetectorGeometry  * fRecoilChamber;
      DriftChamberDetectorGeometry   * fDriftChamber;
      G4LogicalVolume   * wire_hex_log   ;

   protected:
      void DefineMaterials();

   public:
      B1DetectorConstruction();
      virtual ~B1DetectorConstruction();

      virtual G4VPhysicalVolume* Construct();

      void SetRadiatorMaterial(G4String);
      void SetCollimatorMaterial(G4String);

      G4double GetRadiatorCollimatorGap() const {return radiator_collimator_gap;}
      G4double GetCollimatorLength() const {return collimator_length;}

      void     SetRadiatorCollimatorGap(G4double l) ;
      void     SetCollimatorLength(G4double l) ;
      void     SetInnerCollimatorOD(G4double l) ;

      void     PrintConfigInfo() const;

      void CalculatePositions();

      void Rebuild();

      G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

   protected:
};
//______________________________________________________________________________

#endif

