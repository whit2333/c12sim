#include "MicromegasVertexTrackerDetectorGeometry.h"

#include "G4VisAttributes.hh"
#include "G4GenericTrap.hh"
#include "G4Box.hh"
#include "G4IntersectionSolid.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4UserLimits.hh"
#include "G4Polyhedra.hh"
#include "G4Tubs.hh"
#include <algorithm>

MicromegasVertexTrackerDetectorGeometry::MicromegasVertexTrackerDetectorGeometry()
{
   using namespace clas12::geo;
   using namespace CLHEP;

   //G4SDManager* SDMan = G4SDManager::GetSDMpointer();
   //fSensitiveDetector = new DriftChamberSensitiveDetector("DriftChamber",6*6*6*112);
   //SDMan->AddNewDetector(fSensitiveDetector);

   // Huge box for geometric intersection with trapezoid  so that the placement
   // of the region is easier . 
   G4VSolid * temp_box    = new G4Box("DC_placement_box_solid",5*cm,5*cm,5*cm);
   G4VSolid * temp_region = 0;

// FMT | root |Forward Micromegas Vertex Tracker | 0*mm 0*mm 326.250*mm | 0*deg 0*deg 0*deg | aaaaff | Tube | 20.000*mm 250.000*mm 29*mm 0*deg 360*deg | Air | no | 1 | 1 | 1 | 0 | 0 | no | no | no 
   fFMT_solid = new G4Tubs("fFMT_solid",20.000*mm, 250.000*mm, 29*mm, 0*deg, 360*deg);
   fFMT_pos   = {0*mm, 0*mm, 326.250*mm};

// FMT_Peek1 | FMT | PEEK, Layer 1, | 0*mm 0*mm -21*mm | 0*deg 0*deg 0*deg | 666666 | Tube | 25.000*mm 32.500*mm 5*mm 0*deg 360*deg | peek | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
   fFMT_Peek_solid = new G4Tubs("fFMT_Peek_solid",25.000*mm, 32.500*mm, 5*mm, 0*deg, 360*deg);
   fFMT_Peek_pos   = {0*mm, 0*mm, -21*mm };

// FMT_Drift_V_L1 | FMT | Drift V, Layer 1, | 0*mm 0*mm -26.823*mm | 0*deg 0*deg 0*deg | fff600 | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | MMMylar | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
   fFMT_Drift_V_solid = new G4Tubs("fFMT_Drift_V_solid",32.500*mm, 215.000*mm, 0.05*mm, 0*deg, 360*deg);
   fFMT_Drift_V_pos   = {0*mm, 0*mm, -26.823*mm };

// FMT_Gas2_V_L1 | FMT | Gas2 V, Layer 1, | 0*mm 0*mm -24.273*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 2.5*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | fmt | fmt |superlayer manual 1 type manual 1 segment manual 1 strip manual 1 
   fFMT_Gas2_V_solid = new G4Tubs("fFMT_Gas2_V_solid",32.500*mm, 215.000*mm, 2.5*mm, 0*deg, 360*deg);
   fFMT_Gas2_V_pos   = {0*mm, 0*mm, -24.273*mm };

// FMT_Mesh_V_L1 | FMT | Mesh V, Layer 1, | 0*mm 0*mm -21.758*mm | 0*deg 0*deg 0*deg | 252020 | Tube | 32.500*mm 215.000*mm 0.015*mm 0*deg 360*deg | MMMesh | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas1_V_L1 | FMT | Gas1 V, Layer 1, | 0*mm 0*mm -21.679*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 0.064*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Strips_V_L1 | FMT | Strips V, Layer 1, | 0*mm 0*mm -21.6075*mm | 0*deg 0*deg 0*deg | 353540 | Tube | 32.500*mm 215.000*mm 0.0075*mm 0*deg 360*deg | MMStrips | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_PCBoard_V_L1 | FMT | PC Board V, Layer 1, | 0*mm 0*mm -21.55*mm | 0*deg 0*deg 0*deg | 0000ff | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | pcBoardMaterial | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Epoxy1 | FMT | Epoxy, Layer 1, | 0*mm 0*mm -21*mm | 0*deg 0*deg 0*deg | e200e1 | Tube | 32.500*mm 215.000*mm 0.5*mm 0*deg 360*deg | epoxy | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
// FMT_PCBoard_W_L1 | FMT | PC Board W, Layer 1, | 0*mm 0*mm -20.45*mm | 0*deg 0*deg 0*deg | 0000ff | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | pcBoardMaterial | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Strips_W_L1 | FMT | Strips W, Layer 1, | 0*mm 0*mm -20.3925*mm | 0*deg 0*deg 0*deg | 353540 | Tube | 32.500*mm 215.000*mm 0.0075*mm 0*deg 360*deg | MMStrips | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas1_W_L1 | FMT | Gas1 W, Layer 1, | 0*mm 0*mm -20.321*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 0.064*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Mesh_W_L1 | FMT | Mesh W, Layer 1, | 0*mm 0*mm -20.242*mm | 0*deg 0*deg 0*deg | 252020 | Tube | 32.500*mm 215.000*mm 0.015*mm 0*deg 360*deg | MMMesh | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas2_W_L1 | FMT | Gas2 W, Layer 1, | 0*mm 0*mm -17.727*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 2.5*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | fmt | fmt |superlayer manual 1 type manual 2 segment manual 1 strip manual 1 
// FMT_Drift_W_L1 | FMT | Drift W, Layer 1, | 0*mm 0*mm -15.177*mm | 0*deg 0*deg 0*deg | fff600 | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | MMMylar | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Peek2 | FMT | PEEK, Layer 2, | 0*mm 0*mm 0*mm | 0*deg 0*deg 0*deg | 666666 | Tube | 25.000*mm 32.500*mm 5*mm 0*deg 360*deg | peek | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
// FMT_Drift_V_L2 | FMT | Drift V, Layer 2, | 0*mm 0*mm -5.823*mm | 0*deg 0*deg 0*deg | fff600 | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | MMMylar | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas2_V_L2 | FMT | Gas2 V, Layer 2, | 0*mm 0*mm -3.273*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 2.5*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | fmt | fmt |superlayer manual 2 type manual 1 segment manual 1 strip manual 1 
// FMT_Mesh_V_L2 | FMT | Mesh V, Layer 2, | 0*mm 0*mm -0.758*mm | 0*deg 0*deg 0*deg | 252020 | Tube | 32.500*mm 215.000*mm 0.015*mm 0*deg 360*deg | MMMesh | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas1_V_L2 | FMT | Gas1 V, Layer 2, | 0*mm 0*mm -0.679*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 0.064*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Strips_V_L2 | FMT | Strips V, Layer 2, | 0*mm 0*mm -0.6075*mm | 0*deg 0*deg 0*deg | 353540 | Tube | 32.500*mm 215.000*mm 0.0075*mm 0*deg 360*deg | MMStrips | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_PCBoard_V_L2 | FMT | PC Board V, Layer 2, | 0*mm 0*mm -0.55*mm | 0*deg 0*deg 0*deg | 0000ff | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | pcBoardMaterial | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Epoxy2 | FMT | Epoxy, Layer 2, | 0*mm 0*mm 0*mm | 0*deg 0*deg 0*deg | e200e1 | Tube | 32.500*mm 215.000*mm 0.5*mm 0*deg 360*deg | epoxy | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
// FMT_PCBoard_W_L2 | FMT | PC Board W, Layer 2, | 0*mm 0*mm 0.55*mm | 0*deg 0*deg 0*deg | 0000ff | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | pcBoardMaterial | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Strips_W_L2 | FMT | Strips W, Layer 2, | 0*mm 0*mm 0.6075*mm | 0*deg 0*deg 0*deg | 353540 | Tube | 32.500*mm 215.000*mm 0.0075*mm 0*deg 360*deg | MMStrips | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas1_W_L2 | FMT | Gas1 W, Layer 2, | 0*mm 0*mm 0.679*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 0.064*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Mesh_W_L2 | FMT | Mesh W, Layer 2, | 0*mm 0*mm 0.758*mm | 0*deg 0*deg 0*deg | 252020 | Tube | 32.500*mm 215.000*mm 0.015*mm 0*deg 360*deg | MMMesh | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas2_W_L2 | FMT | Gas2 W, Layer 2, | 0*mm 0*mm 3.273*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 2.5*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | fmt | fmt |superlayer manual 2 type manual 2 segment manual 1 strip manual 1 
// FMT_Drift_W_L2 | FMT | Drift W, Layer 2, | 0*mm 0*mm 5.823*mm | 0*deg 0*deg 0*deg | fff600 | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | MMMylar | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Peek3 | FMT | PEEK, Layer 3, | 0*mm 0*mm 21*mm | 0*deg 0*deg 0*deg | 666666 | Tube | 25.000*mm 32.500*mm 5*mm 0*deg 360*deg | peek | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
// FMT_Drift_V_L3 | FMT | Drift V, Layer 3, | 0*mm 0*mm 15.177*mm | 0*deg 0*deg 0*deg | fff600 | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | MMMylar | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas2_V_L3 | FMT | Gas2 V, Layer 3, | 0*mm 0*mm 17.727*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 2.5*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | fmt | fmt |superlayer manual 3 type manual 1 segment manual 1 strip manual 1 
// FMT_Mesh_V_L3 | FMT | Mesh V, Layer 3, | 0*mm 0*mm 20.242*mm | 0*deg 0*deg 0*deg | 252020 | Tube | 32.500*mm 215.000*mm 0.015*mm 0*deg 360*deg | MMMesh | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas1_V_L3 | FMT | Gas1 V, Layer 3, | 0*mm 0*mm 20.321*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 0.064*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Strips_V_L3 | FMT | Strips V, Layer 3, | 0*mm 0*mm 20.3925*mm | 0*deg 0*deg 0*deg | 353540 | Tube | 32.500*mm 215.000*mm 0.0075*mm 0*deg 360*deg | MMStrips | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_PCBoard_V_L3 | FMT | PC Board V, Layer 3, | 0*mm 0*mm 20.45*mm | 0*deg 0*deg 0*deg | 0000ff | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | pcBoardMaterial | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Epoxy3 | FMT | Epoxy, Layer 3, | 0*mm 0*mm 21*mm | 0*deg 0*deg 0*deg | e200e1 | Tube | 32.500*mm 215.000*mm 0.5*mm 0*deg 360*deg | epoxy | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
// FMT_PCBoard_W_L3 | FMT | PC Board W, Layer 3, | 0*mm 0*mm 21.55*mm | 0*deg 0*deg 0*deg | 0000ff | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | pcBoardMaterial | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Strips_W_L3 | FMT | Strips W, Layer 3, | 0*mm 0*mm 21.6075*mm | 0*deg 0*deg 0*deg | 353540 | Tube | 32.500*mm 215.000*mm 0.0075*mm 0*deg 360*deg | MMStrips | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas1_W_L3 | FMT | Gas1 W, Layer 3, | 0*mm 0*mm 21.679*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 0.064*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Mesh_W_L3 | FMT | Mesh W, Layer 3, | 0*mm 0*mm 21.758*mm | 0*deg 0*deg 0*deg | 252020 | Tube | 32.500*mm 215.000*mm 0.015*mm 0*deg 360*deg | MMMesh | no | 1 | 1 | 1 | 1 | 0 | no | no | no 
// FMT_Gas2_W_L3 | FMT | Gas2 W, Layer 3, | 0*mm 0*mm 24.273*mm | 0*deg 0*deg 0*deg | e10000 | Tube | 32.500*mm 215.000*mm 2.5*mm 0*deg 360*deg | MMGas | no | 1 | 1 | 1 | 1 | 0 | fmt | fmt |superlayer manual 3 type manual 2 segment manual 1 strip manual 1 
// FMT_Drift_W_L3 | FMT | Drift W, Layer 3, | 0*mm 0*mm 26.823*mm | 0*deg 0*deg 0*deg | fff600 | Tube | 32.500*mm 215.000*mm 0.05*mm 0*deg 360*deg | MMMylar | no | 1 | 1 | 1 | 1 | 0 | no | no | no 

}
//______________________________________________________________________________

MicromegasVertexTrackerDetectorGeometry::~MicromegasVertexTrackerDetectorGeometry()
{
}
//______________________________________________________________________________

void MicromegasVertexTrackerDetectorGeometry::BuildLogicalVolumes()
{

   bool check_overlaps = false;

   using namespace CLHEP;
   using namespace clas12::geo;

   G4NistManager * matman = G4NistManager::Instance();


   //G4Element     * Ar     = new G4Element("Argon", "Ar", /*z    = */18, /*a            = */ 39.95*g/mole);
   fMMGas_mat   = new G4Material("fMMGas_mat", /* density = */ 1.70335*mg/cm3, /*nel = */ 3);
   fMMGas_mat->AddMaterial(matman->FindOrBuildMaterial("G4_Ar"),     95.0*perCent);
   fMMGas_mat->AddMaterial(matman->FindOrBuildMaterial("G4_H"),  0.86707*perCent);
   fMMGas_mat->AddMaterial(matman->FindOrBuildMaterial("G4_C"),  4.13293*perCent);

   //mmstrips  |          ft micromegas strips  |      6.72  |         1  |                                           G4_Cu 1  | none  | none  | none  | none  | none  | none  | none  
   //   |   -1  |   -1  |   -1  |   -1  |   -1  
   //   mmmesh  |            ft micromegas mesh  |      8.02  |         5  |G4_Mn 0.02 G4_Si 0.01 G4_Cr 0.19 G4_Ni 0.10 G4_Fe 0.68  | none  | none  | none  | none  | none  | none  | no
   //   ne  |   -1  |   -1  |   -1  |   -1  |   -1  
   //   mmgas  |             ft micromegas gas  |0.00170335  |         3  |          G4_Ar 0.95 G4_H 0.0086707 G4_C 0.0413293  | none  | none  | none  | none  | none  | none  | none  
   //   |   -1  |   -1  |   -1  |   -1  |   -1  
   //   mmmylar  | ft micromegas mylar 1.40g/cm3  |       1.4  |         3  |         G4_H 0.041958 G4_C 0.625017 G4_O 0.333025  | none  | none  | none  | none  | none  | none  | none  
   //   |   -1  |   -1  |   -1  |   -1  |   -1  
   fMMMylar_mat   = new G4Material("fMMMylar_mat", /* density = */ 1.4*mg/cm3, /*nel = */ 3);
   fMMMylar_mat->AddMaterial(matman->FindOrBuildMaterial("G4_H"),   4.1958*perCent);
   fMMMylar_mat->AddMaterial(matman->FindOrBuildMaterial("G4_C"),  62.5017*perCent);
   fMMMylar_mat->AddMaterial(matman->FindOrBuildMaterial("G4_O"),  33.3025*perCent);
   //   epoxy  |              micromegas epoxy  |      1.16  |         4  |                                 C 15 H 32 N 2 O 4  | none  | none  | none  | none  | none  | none  | none  |   -1  |   -1  |   -1  |   -1  |   -1  
   //      MMStrips  |  micromegas strips are copper  |      6.72  |         1  |                                           G4_Cu 1  | none  | none  | none  | none  | none  | none  | none  |   -1  |   -1  |   -1  |   -1  |   -1  
   //      pcBoardMaterial  |         bst pc board material  |     1.860  |         3  |                      G4_Fe 0.3 G4_C 0.4 G4_Si 0.3  | none  | none  | none  | none  | none  | none  | none  |   -1  |   -1  |   -1  |   -1  |   -1  
   //      peek  |                      bst peek  |      1.31  |         3  |                                     C 19 H 12 O 3  | none  | none  | none  | none  | none  | none  | none  |   -1  |   -1  |   -1  |   -1  |   -1  

   fPeek_mat   = new G4Material("fPeek_mat", /* density = */ 1.31*g/cm3, /*nel = */ 3);
   fPeek_mat->AddMaterial(matman->FindOrBuildMaterial("G4_C"),     19.0*perCent);
   fPeek_mat->AddMaterial(matman->FindOrBuildMaterial("G4_H"),     12.0*perCent);
   fPeek_mat->AddMaterial(matman->FindOrBuildMaterial("G4_O"),      3.0*perCent);

   G4VisAttributes * vs = new G4VisAttributes(G4Colour(0.5,0.0,0.5));
   //vs->SetDaughtersInvisible(true);
   vs->SetForceWireframe(true);
   //vs->SetForceSolid(true);

   G4VisAttributes * vs2 = new G4VisAttributes(G4Colour(0.3,0.5,0.2));//G4VisAttributes::GetInvisible());//
   //vs2->SetDaughtersInvisible(true);
   vs2->SetForceWireframe(true);

   G4VisAttributes * vs3 = new G4VisAttributes(G4Colour(0.8,0.1,0.1));//G4VisAttributes::GetInvisible());//;
   //vs3->SetDaughtersInvisible(true);
   vs3->SetForceSolid(true);

   G4VisAttributes * vs_odd = new G4VisAttributes(G4Colour(0.8,0.2,0.0));//G4VisAttributes::GetInvisible());//
   vs_odd->SetForceWireframe(true);
   //vs_odd->SetForceSolid(true);

   G4VisAttributes * vs_even = new G4VisAttributes(G4Colour(0.0,0.2,0.8));//G4VisAttributes::GetInvisible());//
   vs_even->SetForceWireframe(true);
   //vs_even->SetForceSolid(true);

   fFMT_mat = matman->FindOrBuildMaterial("G4_AIR");

   fFMT_log = new G4LogicalVolume(fFMT_solid, fFMT_mat, "fFMT_log");
   fFMT_log->SetVisAttributes(vs);

   fFMT_Peek_log = new G4LogicalVolume(fFMT_Peek_solid, fPeek_mat, "fFMT_Peek_log");
   fFMT_Peek_log->SetVisAttributes(vs2);
   fFMT_Peek_phys = new G4PVPlacement(
         0, fFMT_Peek_pos ,
         fFMT_Peek_log,          // its logical volume
         "fFMT_Peek_phys", // its name
         fFMT_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         false);                        // check for overlaps

   fFMT_Drift_V_log = new G4LogicalVolume(fFMT_Drift_V_solid, fMMMylar_mat, "fFMT_Drift_V_log");
   fFMT_Drift_V_log->SetVisAttributes(vs3);
   fFMT_Drift_V_phys = new G4PVPlacement(
         0, fFMT_Drift_V_pos ,
         fFMT_Drift_V_log,          // its logical volume
         "fFMT_Drift_V_phys", // its name
         fFMT_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         false);                        // check for overlaps

   fFMT_Gas2_V_log = new G4LogicalVolume(fFMT_Gas2_V_solid, fMMGas_mat, "fFMT_Gas2_V_log");
   fFMT_Gas2_V_log->SetVisAttributes(vs_odd);
   fFMT_Gas2_V_phys = new G4PVPlacement(
         0, fFMT_Gas2_V_pos ,
         fFMT_Gas2_V_log,          // its logical volume
         "fFMT_Gas2_V_phys", // its name
         fFMT_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         false);                        // check for overlaps

}
//______________________________________________________________________________

G4VPhysicalVolume * MicromegasVertexTrackerDetectorGeometry::PlacePhysicalVolume(G4LogicalVolume * mother)
{
   using namespace clas12::geo;
   
   fFMT_phys = new G4PVPlacement(
         0, fFMT_pos ,
         fFMT_log,          // its logical volume
         "fFMT_phys", // its name
         mother,                       // its mother (logical) volume
         false,                        // no boolean operations
         1,                     // its copy number
         false);                        // check for overlaps

   return fFMT_phys;
}
//______________________________________________________________________________

G4VPhysicalVolume * MicromegasVertexTrackerDetectorGeometry::PlaceParallelPhysicalVolume(G4LogicalVolume * mother)
{
   using namespace clas12::geo;
   return nullptr;
}
//______________________________________________________________________________

