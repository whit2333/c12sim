#include "B1ParallelWorldConstruction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"    
#include "G4SDManager.hh"    
#include "G4Colour.hh"    
#include "G4VisAttributes.hh"    
#include <string>
#include "FakeSD.hh"    
#include "SimulationManager.h"
#include "RecoilChamberDetectorGeometry.h"
#include "DriftChamberDetectorGeometry.h"
#include "HTCCDetectorGeometry.h"

B1ParallelWorldConstruction ::B1ParallelWorldConstruction(G4String& parallelWorldName) :
   G4VUserParallelWorld(parallelWorldName),
   fConstructed(false),
   fNeedsRebuilt(true),
   fSpanDistance(40.0*cm),
   fDirection(0,0,1.0),
   fStartingPoint(0,0,0),
   fDet_size(50.0*cm)
{
   for(int i = 0; i < fNplanes; i++) {
      fDet_solid[i]  = nullptr;
      fDet_log[i]    = nullptr;
      fDet_phys[i]   = nullptr;
      fDet_det[i]    = nullptr;
      fDet_vis[i]    = nullptr;
   }
   std::cout << "Parallel world ctor" <<std::endl;
}
//______________________________________________________________________________

B1ParallelWorldConstruction::~B1ParallelWorldConstruction()
{ }
//______________________________________________________________________________

void B1ParallelWorldConstruction::Construct()
{
   if( !fNeedsRebuilt && fConstructed) return;
   fConstructed = true;

   std::cout << "constriing Parallel world " <<std::endl;

   bool    checkOverlaps    = false;
   int     natoms           = 0;
   int     ncomponents      = 0;
   double  A                = 0.0;
   double  Z                = 0.0;
   double  thickness        = 0.0;
   double  surface_z        = 0.0;
   double  red              = 0.0;
   double  green            = 0.0;
   double  blue             = 0.0;
   double  alpha            = 0.0;
   double  density          = 0.0;
   double  pressure         = 0.0;
   double  temperature      = 0.0;

   G4VSolid          * scoring_solid  = nullptr;
   G4LogicalVolume   * scoring_log    = nullptr;
   G4VPhysicalVolume * scoring_phys   = nullptr;
   FakeSD            * scoring_det    = nullptr;
   G4VisAttributes   * scoring_vis    = nullptr;


   // --------------------------------------------------------------
   // World
   //
   G4VPhysicalVolume * ghostWorld   = GetWorld();
   G4LogicalVolume   * worldLogical = ghostWorld->GetLogicalVolume();

   // --------------------------------------------------------------
   // 
   double        scoring_length = 1.0*um;
   double        step_size   = fSpanDistance/double(fNplanes-1);
   G4ThreeVector scoring_pos = fStartingPoint;
   G4ThreeVector step        = fDirection;
   step.setMag(step_size);

   //   red       = 177.0/256.0;
   //   green     = 104.0/256.0;
   //   blue      = 177.0/256.0;
   //   alpha     = 0.4;

   if(!scoring_det) scoring_det = new FakeSD("/FakeSD0");
   G4SDManager::GetSDMpointer()->AddNewDetector(scoring_det);

   for (int i = 0; i<fNplanes; i++) {
      double radius = 5.0*cm + double(i)*10.0*cm;

      G4String solid_name = G4String("scoring_tube_solid_") + std::to_string(i);
      fDet_solid[i] = new G4Tubs(solid_name, radius, radius + 1.0*mm, 5.0*m, 0,  360.0*CLHEP::degree);
      solid_name = G4String("scoring_tube_log_") + std::to_string(i);
      fDet_log[i]   = new G4LogicalVolume(fDet_solid[i], 0, solid_name,0, 0,0);
      solid_name = G4String("scoring_tube_phys_") + std::to_string(i);
      fDet_phys[i]  = new G4PVPlacement(0,G4ThreeVector(0.0,0,4.0*m), fDet_log[i], solid_name,worldLogical,false,i,false);
      fDet_log[i]->SetSensitiveDetector(scoring_det);
   }

   // ------------------------------------------------------------------------
   // Drift Chamber
   // ------------------------------------------------------------------------
   //DriftChamberDetectorGeometry * fDriftChamber = SimulationManager::GetInstance()->GetDriftDetectorGeometry();
   //// Sectors
   //for(int i = 1; i<=6; i++ ) {

   //   // Region I
   //   fDriftChamber->PlaceParallelPhysicalVolume( worldLogical, i, 1);

   //   // Region II
   //   fDriftChamber->PlaceParallelPhysicalVolume( worldLogical, i, 2);

   //   // Region III
   //   fDriftChamber->PlaceParallelPhysicalVolume( worldLogical, i, 3);
   //}

   // ------------------------------------------------------------------------
   // Recoil Chamber
   // ------------------------------------------------------------------------
   //RecoilChamberDetectorGeometry * fRecoilChamber = SimulationManager::GetInstance()->GetRecoilDetectorGeometry();
   //fRecoilChamber->He10CO2   = He10CO2;
   //fRecoilChamber->HeiC4H10  = HeiC4H10;
   //fRecoilChamber->Tungsten  = Tungsten; 
   //fRecoilChamber->Mylar     = Mylar;
   //fRecoilChamber->PlaceParallelPhysicalVolume( worldLogical);

   // ------------------------------------------------------------------------
   // HTCC
   // ------------------------------------------------------------------------
   //HTCCDetectorGeometry * fHTCC = SimulationManager::GetInstance()->GetHTCCDetectorGeometry();
   //fHTCC->He10CO2   = He10CO2;
   //fHTCC->HeiC4H10  = HeiC4H10;
   //fHTCC->Tungsten  = Tungsten; 
   //fHTCC->Mylar     = Mylar;
   //fHTCC->PlaceParallelPhysicalVolume( worldLogical);

   fNeedsRebuilt = false;

}
//______________________________________________________________________________

