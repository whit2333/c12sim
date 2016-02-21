#include "SiliconVertexTrackerDetectorGeometry.h"

#include "G4VisAttributes.hh"
#include "G4GenericTrap.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4IntersectionSolid.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4UserLimits.hh"
#include "G4Polyhedra.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <algorithm>

SiliconVertexTrackerDetectorGeometry::SiliconVertexTrackerDetectorGeometry()
{
   using namespace clas12::geo;
   using namespace CLHEP;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();
   //fSensitiveDetector = new DriftChamberSensitiveDetector("DriftChamber",6*6*6*112);
   //SDMan->AddNewDetector(fSensitiveDetector);
   //fSensitiveDetector2 = new SensitiveRegionDetector("DriftChamberRegion",6*3);
   //SDMan->AddNewDetector(fSensitiveDetector2);

   double extra = 0.001*mm;

   fBST_solid[0] = new G4Box( G4String("fBST_solid_0"),
         fRohacell71Thickness/2.0 + fSiLayerThickness ,fBSTModuleWidth/2.0 + extra, fBST_length[0]/2.0 + fBSTReadoutLength/2.0 + extra);
   fBST_solid[1] = new G4Box( G4String("fBST_solid_1"),
         fRohacell71Thickness/2.0 + fSiLayerThickness ,fBSTModuleWidth/2.0 + extra, fBST_length[1]/2.0 + fBSTReadoutLength/2.0 + extra);
   fBST_solid[2] = new G4Box( G4String("fBST_solid_2"),
         fRohacell71Thickness/2.0 + fSiLayerThickness ,fBSTModuleWidth/2.0 + extra, fBST_length[2]/2.0 + fBSTReadoutLength/2.0 + extra);
   fBST_solid[3] = new G4Box( G4String("fBST_solid_3"),
         fRohacell71Thickness/2.0 + fSiLayerThickness ,fBSTModuleWidth/2.0 + extra, fBST_length[3]/2.0 + fBSTReadoutLength/2.0 + extra);

   fBSTSilicon_solid[0] = new G4Box("fBSTSilicon_solid_0", fSiLayerThickness/2.0, fBSTModuleWidth/2.0, fBST_length[0]/2.0 );
   fBSTSilicon_solid[1] = new G4Box("fBSTSilicon_solid_1", fSiLayerThickness/2.0, fBSTModuleWidth/2.0, fBST_length[1]/2.0 );
   fBSTSilicon_solid[2] = new G4Box("fBSTSilicon_solid_2", fSiLayerThickness/2.0, fBSTModuleWidth/2.0, fBST_length[2]/2.0 );
   fBSTSilicon_solid[3] = new G4Box("fBSTSilicon_solid_3", fSiLayerThickness/2.0, fBSTModuleWidth/2.0, fBST_length[3]/2.0 );

   fBSTRohacell71_solid[0] = new G4Box("fBSTRohacell71_solid_0", fRohacell71Thickness/2.0, fBSTModuleWidth/2.0, fBST_length[0]/2.0 );
   fBSTRohacell71_solid[1] = new G4Box("fBSTRohacell71_solid_1", fRohacell71Thickness/2.0, fBSTModuleWidth/2.0, fBST_length[1]/2.0 );
   fBSTRohacell71_solid[2] = new G4Box("fBSTRohacell71_solid_2", fRohacell71Thickness/2.0, fBSTModuleWidth/2.0, fBST_length[2]/2.0 );
   fBSTRohacell71_solid[3] = new G4Box("fBSTRohacell71_solid_3", fRohacell71Thickness/2.0, fBSTModuleWidth/2.0, fBST_length[3]/2.0 );

   fFST_solid[0] = new G4Tubs("fFST_solid_0", 
         fFSTInnerRadus - extra, 
         fFSTSiliconRadus + fBSTReadoutLength/2.0 + extra, 
         fSiLayerThickness + fRohacell71Thickness/2.0 + extra,
         0.0, 360*degree);

   fFST_solid[1] = new G4Tubs("fFST_solid_1", 
         fFSTInnerRadus - extra, 
         fFSTSiliconRadus + fBSTReadoutLength/2.0 + extra, 
         fSiLayerThickness + fRohacell71Thickness/2.0 + extra,
         0.0, 360*degree);
   fFST_solid[2] = new G4Tubs("fFST_solid_2", 
         fFSTInnerRadus - extra, 
         fFSTSiliconRadus + fBSTReadoutLength/2.0 + extra, 
         fSiLayerThickness + fRohacell71Thickness/2.0 + extra,
         0.0, 360*degree);

   fFSTSilicon_solid[0] = new G4Tubs("fFSTSilicon_solid_0", 
         fFSTInnerRadus - extra, 
         fFSTSiliconRadus , 
         fSiLayerThickness/2.0 ,
         0.0, 360*degree);
   fFSTSilicon_solid[1] = new G4Tubs("fFSTSilicon_solid_1", 
         fFSTInnerRadus - extra, 
         fFSTSiliconRadus , 
         fSiLayerThickness/2.0 ,
         0.0, 360*degree);
   fFSTSilicon_solid[2] = new G4Tubs("fFSTSilicon_solid_2", 
         fFSTInnerRadus - extra, 
         fFSTSiliconRadus , 
         fSiLayerThickness/2.0 ,
         0.0, 360*degree);


   fFSTRohacell71_solid[0] = new G4Tubs("fFSTRohacell71_solid_0", 
         fFSTInnerRadus - extra, 
         fFSTSiliconRadus , 
         fRohacell71Thickness/2.0 ,
         0.0, 360*degree);
   fFSTRohacell71_solid[1] = new G4Tubs("fFSTRohacell71_solid_1", 
         fFSTInnerRadus - extra, 
         fFSTSiliconRadus , 
         fRohacell71Thickness/2.0 ,
         0.0, 360*degree);
   fFSTRohacell71_solid[2] = new G4Tubs("fFSTRohacell71_solid_2", 
         fFSTInnerRadus - extra, 
         fFSTSiliconRadus , 
         fRohacell71Thickness/2.0 ,
         0.0, 360*degree);


   G4NistManager * matman = G4NistManager::Instance();
   //G4Material    * Si    = matman->FindOrBuildMaterial("G4_Si");
   G4Element     * Si     = new G4Element("Silicon", "Si", /*z = */14, /*a = */ 28.084*g/mole);
   fSiliconMaterial       = new G4Material("SiliconMaterial", /* density = */ fSiDensity, /*nel = */ 1);
   fSiliconMaterial->AddElement(Si, 100*perCent);

   fRohacell71Material = new G4Material("SiliconMaterial", /* density = */ fRohacell71Density, /*nel = */ 4);
   fRohacell71Material->AddMaterial( matman->FindOrBuildMaterial("G4_C"), 0.6465);
   fRohacell71Material->AddMaterial( matman->FindOrBuildMaterial("G4_H"), 0.0784);
   fRohacell71Material->AddMaterial( matman->FindOrBuildMaterial("G4_N"), 0.0839);
   fRohacell71Material->AddMaterial( matman->FindOrBuildMaterial("G4_O"), 0.1912);
   //G4_C 0.6465 G4_H 0.0784 G4_N 0.0839 G4_O 0.1912
}
//______________________________________________________________________________

SiliconVertexTrackerDetectorGeometry::~SiliconVertexTrackerDetectorGeometry()
{
}
//______________________________________________________________________________

void SiliconVertexTrackerDetectorGeometry::BuildLogicalVolumes()
{

   bool check_overlaps = false;

   using namespace CLHEP;
   using namespace clas12::geo;
   using namespace clas12::geo::DC;

   G4VisAttributes * vs_odd = new G4VisAttributes(G4Colour(0.8,0.2,0.0));//G4VisAttributes::GetInvisible());//
   vs_odd->SetForceWireframe(true);
   //vs_odd->SetForceSolid(true);

   G4VisAttributes * vs_even = new G4VisAttributes(G4Colour(0.0,0.2,0.8));//G4VisAttributes::GetInvisible());//
   vs_even->SetForceWireframe(true);
   //vs_even->SetForceSolid(true);

   G4NistManager * matman = G4NistManager::Instance();
   G4Material    * Air    = matman->FindOrBuildMaterial("G4_AIR");

   // -----------------------------
   // BST Module containers
   fBST_log[0] = new G4LogicalVolume(fBST_solid[0], Air, "fBST_log_0");
   fBST_log[1] = new G4LogicalVolume(fBST_solid[1], Air, "fBST_log_1");
   fBST_log[2] = new G4LogicalVolume(fBST_solid[2], Air, "fBST_log_2");
   fBST_log[3] = new G4LogicalVolume(fBST_solid[3], Air, "fBST_log_3");

   G4VisAttributes * vs = new G4VisAttributes(G4Colour(0.5,0.0,0.5));
   //vs->SetDaughtersInvisible(true);
   vs->SetForceWireframe(true);
   //vs->SetForceSolid(true);
   fBST_log[0]->SetVisAttributes(vs);
   fBST_log[1]->SetVisAttributes(vs);
   fBST_log[2]->SetVisAttributes(vs);
   fBST_log[3]->SetVisAttributes(vs);

   // -----------------------------
   // BST Silicon
   fBSTSilicon_log[0] = new G4LogicalVolume(fBSTSilicon_solid[0], fSiliconMaterial, "fBSTSilicon_log_0");
   fBSTSilicon_log[1] = new G4LogicalVolume(fBSTSilicon_solid[1], fSiliconMaterial, "fBSTSilicon_log_1");
   fBSTSilicon_log[2] = new G4LogicalVolume(fBSTSilicon_solid[2], fSiliconMaterial, "fBSTSilicon_log_2");
   fBSTSilicon_log[3] = new G4LogicalVolume(fBSTSilicon_solid[3], fSiliconMaterial, "fBSTSilicon_log_3");

   G4VisAttributes * vs2 = new G4VisAttributes(G4Colour(0.3,0.5,0.2,0.4));//G4VisAttributes::GetInvisible());//
   //vs2->SetDaughtersInvisible(true);
   //vs2->SetForceWireframe(true);
   vs2->SetForceSolid(true);
   fBSTSilicon_log[0]->SetVisAttributes(vs2);
   fBSTSilicon_log[1]->SetVisAttributes(vs2);
   fBSTSilicon_log[2]->SetVisAttributes(vs2);
   fBSTSilicon_log[3]->SetVisAttributes(vs2);

   // -----------------------------
   // BST Rohacell71
   fBSTRohacell71_log[0] = new G4LogicalVolume(fBSTRohacell71_solid[0], fRohacell71Material, "fBSTRohacell71_log_0");
   fBSTRohacell71_log[1] = new G4LogicalVolume(fBSTRohacell71_solid[1], fRohacell71Material, "fBSTRohacell71_log_1");
   fBSTRohacell71_log[2] = new G4LogicalVolume(fBSTRohacell71_solid[2], fRohacell71Material, "fBSTRohacell71_log_2");
   fBSTRohacell71_log[3] = new G4LogicalVolume(fBSTRohacell71_solid[3], fRohacell71Material, "fBSTRohacell71_log_3");

   G4VisAttributes * vs3 = new G4VisAttributes(G4Colour(0.8,0.1,0.1));//G4VisAttributes::GetInvisible());//;
   //vs3->SetDaughtersInvisible(true);
   vs3->SetForceWireframe(true);
   //vs3->SetForceSolid(true);
   fBSTRohacell71_log[0]->SetVisAttributes(vs3);
   fBSTRohacell71_log[1]->SetVisAttributes(vs3);
   fBSTRohacell71_log[2]->SetVisAttributes(vs3);
   fBSTRohacell71_log[3]->SetVisAttributes(vs3);

   for(int iLayer = 0; iLayer<4; iLayer++){

      // Rohacell
      new G4PVPlacement(
            0,
            G4ThreeVector(0,0,fBSTReadoutLength/2.0),
            fBSTRohacell71_log[iLayer],          // its logical volume
            Form("Rohacell71_phys_%d",iLayer), // its name
            fBST_log[iLayer],                       // its mother (logical) volume
            false,                        // no boolean operations
            iLayer,                         // its copy number
            false);                       // check for overlaps

      // Inner silicon
      new G4PVPlacement(
            0,
            G4ThreeVector(-fSiLayerThickness/2.0-fRohacell71Thickness/2.0,0,fBSTReadoutLength/2.0),
            fBSTSilicon_log[iLayer],          // its logical volume
            Form("Silicon_inner_phys_%d",iLayer), // its name
            fBST_log[iLayer],                       // its mother (logical) volume
            false,                        // no boolean operations
            2*iLayer,                         // its copy number
            false);                       // check for overlaps

      // Outer silicon
      new G4PVPlacement(
            0,
            G4ThreeVector(fSiLayerThickness/2.0 + fRohacell71Thickness/2.0,0,fBSTReadoutLength/2.0),
            fBSTSilicon_log[iLayer],          // its logical volume
            Form("Silicon_outer_phys_%d",iLayer), // its name
            fBST_log[iLayer],                       // its mother (logical) volume
            false,                        // no boolean operations
            2*iLayer+1,                         // its copy number
            false);                       // check for overlaps
   }

   // -----------------------------
   // FST Module containers
   fFST_log[0] = new G4LogicalVolume(fFST_solid[0], Air, "fFST_log_0");
   fFST_log[1] = new G4LogicalVolume(fFST_solid[1], Air, "fFST_log_1");
   fFST_log[2] = new G4LogicalVolume(fFST_solid[2], Air, "fFST_log_2");

   fFST_log[0]->SetVisAttributes(vs);
   fFST_log[1]->SetVisAttributes(vs);
   fFST_log[2]->SetVisAttributes(vs);

   // -----------------------------
   // FST Silicon
   fFSTSilicon_log[0] = new G4LogicalVolume(fFSTSilicon_solid[0], fSiliconMaterial, "fFSTSilicon_log_0");
   fFSTSilicon_log[1] = new G4LogicalVolume(fFSTSilicon_solid[1], fSiliconMaterial, "fFSTSilicon_log_1");
   fFSTSilicon_log[2] = new G4LogicalVolume(fFSTSilicon_solid[2], fSiliconMaterial, "fFSTSilicon_log_2");

   fFSTSilicon_log[0]->SetVisAttributes(vs2);
   fFSTSilicon_log[1]->SetVisAttributes(vs2);
   fFSTSilicon_log[2]->SetVisAttributes(vs2);

   // -----------------------------
   // FST Rohacell71
   fFSTRohacell71_log[0] = new G4LogicalVolume(fFSTRohacell71_solid[0], fRohacell71Material, "fFSTRohacell71_log_0");
   fFSTRohacell71_log[1] = new G4LogicalVolume(fFSTRohacell71_solid[1], fRohacell71Material, "fFSTRohacell71_log_1");
   fFSTRohacell71_log[2] = new G4LogicalVolume(fFSTRohacell71_solid[2], fRohacell71Material, "fFSTRohacell71_log_2");

   fFSTRohacell71_log[0]->SetVisAttributes(vs3);
   fFSTRohacell71_log[1]->SetVisAttributes(vs3);
   fFSTRohacell71_log[2]->SetVisAttributes(vs3);

   for(int iLayer = 0; iLayer<3; iLayer++){

      // Rohacell
      new G4PVPlacement(
            0,
            G4ThreeVector(0,0,0.0),
            fFSTRohacell71_log[iLayer],          // its logical volume
            Form("FST_Rohacell71_phys_%d",iLayer), // its name
            fFST_log[iLayer],                       // its mother (logical) volume
            false,                        // no boolean operations
            iLayer,                         // its copy number
            false);                       // check for overlaps

      // Inner silicon
      new G4PVPlacement(
            0,
            G4ThreeVector(-fSiLayerThickness/2.0 -fRohacell71Thickness/2.0,0,0.0),
            fFSTSilicon_log[iLayer],          // its logical volume
            Form("FST_Silicon_inner_phys_%d",iLayer), // its name
            fFST_log[iLayer],                       // its mother (logical) volume
            false,                        // no boolean operations
            2*iLayer,                         // its copy number
            false);                       // check for overlaps

      // Outer silicon
      new G4PVPlacement(
            0,
            G4ThreeVector(fSiLayerThickness/2.0 + fRohacell71Thickness/2.0,0,0.0),
            fFSTSilicon_log[iLayer],          // its logical volume
            Form("FST_Silicon_outer_phys_%d",iLayer), // its name
            fFST_log[iLayer],                       // its mother (logical) volume
            false,                        // no boolean operations
            2*iLayer+1,                         // its copy number
            false);                       // check for overlaps
   }
}
//______________________________________________________________________________

G4VPhysicalVolume * SiliconVertexTrackerDetectorGeometry::PlacePhysicalVolume(G4LogicalVolume * mother )
{
   using namespace clas12::geo;
   using namespace CLHEP;

   for(int iLayer = 0; iLayer < 4; iLayer++) {

      double angle = 360.0*degree/fNModules[iLayer];

      for(int iMod = 0; iMod < fNModules[iLayer]; iMod++){

         G4RotationMatrix rot = G4RotationMatrix::IDENTITY;
         rot.rotateZ( angle*double(iMod) );

         G4ThreeVector    pos = {
            fBST_radius[iLayer] + fRohacell71Thickness/2.0 + fSiLayerThickness,
            0.0,
            fBSTDownstreamEnd_displacment[iLayer] - fBST_length[iLayer]/2 - fBSTReadoutLength/2.0
         };
         pos.rotateZ( angle*double(iMod) );

         G4VPhysicalVolume * phys = new G4PVPlacement(
               G4Transform3D(rot,pos),
               fBST_log[iLayer],          // its logical volume
               Form("fBST_phys_%d_%d",iLayer,iMod), // its name
               mother,                       // its mother (logical) volume
               false,                        // no boolean operations
               iLayer*fNModules[iLayer]+iMod,                          // its copy number
               false);                       // check for overlaps

      }

   }

   for(int iLayer = 0; iLayer <3; iLayer++) {
         new G4PVPlacement(
               0,
               G4ThreeVector(0,0,fFST_displacement[iLayer]),
               fFST_log[iLayer],          // its logical volume
               Form("fFST_phys_%d",iLayer), // its name
               mother,                       // its mother (logical) volume
               false,                        // no boolean operations
               iLayer,                          // its copy number
               false);                       // check for overlaps
   }


   return nullptr;
}
//______________________________________________________________________________

