#include "SiVertexTrackerDetectorGeometry.h"

#include <cmath>
#include "CLHEP/Units/PhysicalConstants.h"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

#include "G4GenericTrap.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVDivision.hh"


SiVertexTrackerDetectorGeometry::SiVertexTrackerDetectorGeometry()
{
   using namespace clas12::geo;
   using namespace CLHEP;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();

   fSensitiveDetector = new SiTrackerSensitiveDetector("SiVertexTracker",6*6*6*112);
   SDMan->AddNewDetector(fSensitiveDetector);

   fSiVertexTrackerGeometry.Print();

   G4NistManager* nist = G4NistManager::Instance();
   fAir_mat            = nist->FindOrBuildMaterial("G4_AIR");
   fSi_mat             = nist->FindOrBuildMaterial("G4_Si");
}
//______________________________________________________________________________

SiVertexTrackerDetectorGeometry::~SiVertexTrackerDetectorGeometry()
{ }
//______________________________________________________________________________

void SiVertexTrackerDetectorGeometry::BuildLogicalVolumes()
{
   using namespace CLHEP;
   using namespace eic::geo;
   bool check_overlaps = false;

   //G4VisAttributes * vs = new G4VisAttributes(G4Colour(0.0,0.9,0.1));

   for(int i = 0; i<fSiVertexTrackerGeometry.NLayers; i++) {

      std::string solid_name  = "SiVertexTracker_container_solid_";
      solid_name             += std::to_string(i);

      fLayerContainer_solid[i] = new G4Box(
            solid_name.c_str(), 
            fSiVertexTrackerGeometry.fLayerTotalThickness.at(i)/2.0,
            fSiVertexTrackerGeometry.fLayerWidth.at(i)/2.0,
            fSiVertexTrackerGeometry.fLayerLength.at(i)/2.0 );

      fLayerContainer_mat[i] = fAir_mat;

      std::string log_name  = "SiVertexTracker_container_log_";
      log_name           += std::to_string(i);

      fLayerContainer_log[i] = new G4LogicalVolume(
            fLayerContainer_solid[i],
            fLayerContainer_mat[i],
            log_name.c_str() );

      G4VisAttributes * vs1 = new G4VisAttributes(G4Colour(0.8,0.2,0.0));//G4VisAttributes::GetInvisible());//
      vs1->SetForceWireframe(true);
      //vs_even->SetForceSolid(true);
      fLayerContainer_log[i]->SetVisAttributes(vs1);


      // -----------------------------------------------

      solid_name  = "SiVertexTracker_Si_solid_";
      solid_name += std::to_string(i);

      fLayerSi_solid[i] = new G4Box(
            solid_name.c_str(), 
            fSiVertexTrackerGeometry.fLayerTotalThickness.at(i)/2.0,
            fSiVertexTrackerGeometry.fLayerWidth.at(i)/2.0,
            fSiVertexTrackerGeometry.fLayerLength.at(i)/2.0
            );

      fLayerSi_mat[i] = fSi_mat;

      log_name     = "SiVertexTracker_Si_log_";
      log_name    += std::to_string(i);

      fLayerSi_log[i] = new G4LogicalVolume(
            fLayerContainer_solid[i],
            fLayerContainer_mat[i],
            log_name.c_str() );

      G4VisAttributes * vs2 = new G4VisAttributes(G4Colour(0.0+0.2*double(i),0.2,0.8-0.2*double(i),0.5));
      //vs_odd->SetForceWireframe(true);
      vs2->SetForceSolid(true);
      fLayerSi_log[i]->SetVisAttributes(vs2);

      // Not this solid's dimensions will be changed by G4PVDivion
      solid_name     = "Si_slice_solid_";
      solid_name    += std::to_string(i);

      G4Box* slice_solid = new G4Box(solid_name.c_str(),
            fSiVertexTrackerGeometry.fLayerSiThickness.at(i)/2.0,
            fSiVertexTrackerGeometry.fTPixelPitch.at(i)/2.0,
            fSiVertexTrackerGeometry.fZPixelPitch.at(i)/2.0 );

      log_name     = "Si_slice_log_";
      log_name    += std::to_string(i);

      G4LogicalVolume * slice_log = new G4LogicalVolume(slice_solid, fSi_mat, log_name.c_str(),0,0,0);

      // Not this solid's dimensions will be changed by G4PVDivion
      solid_name     = "Si_pixel_solid_";
      solid_name    += std::to_string(i);

      G4Box* pixel_solid = new G4Box(solid_name.c_str(),
            fSiVertexTrackerGeometry.fLayerSiThickness.at(i)/2.0,
            fSiVertexTrackerGeometry.fTPixelPitch.at(i)/2.0,
            fSiVertexTrackerGeometry.fZPixelPitch.at(i)/2.0 );

      log_name     = "Si_pixel_log_";
      log_name    += std::to_string(i);

      G4LogicalVolume * pixel_log = new G4LogicalVolume(pixel_solid, fSi_mat, log_name.c_str(),0,0,0);

      // --------------------------------
      // divide into slices along Z
      std::string  div_name;
      div_name  = "Si_slice_division_";
      div_name += std::to_string(i);

      G4PVDivision * Si_Z_div = new G4PVDivision(
            div_name.c_str(),
            slice_log,
            fLayerSi_log[i],
            kZAxis,
            fSiVertexTrackerGeometry.fZPixelPitch.at(i),
            0 );

      div_name  = "Si_pixel_division_";
      div_name += std::to_string(i);

      G4PVDivision * Si_T_div = new G4PVDivision(
            div_name.c_str(),
            pixel_log,
            slice_log,
            kYAxis,
            fSiVertexTrackerGeometry.fTPixelPitch.at(i),
            0 );

      pixel_log->SetVisAttributes(G4VisAttributes::GetInvisible());
      slice_log->SetVisAttributes(G4VisAttributes::GetInvisible());

      // --------------------------------
      // phys name
      std::string  phys_name;
      phys_name  = "Si_phys_";
      phys_name += std::to_string(i);
      G4VPhysicalVolume * phys = new G4PVPlacement(
            0, G4ThreeVector(0,0,0),
            fLayerSi_log[i], 
            phys_name.c_str(),
            fLayerContainer_log[i],
            false, 0, check_overlaps);
   }

}
//______________________________________________________________________________

G4VPhysicalVolume * SiVertexTrackerDetectorGeometry::PlaceParallelPhysicalVolume(G4LogicalVolume * mother )
{
   using namespace eic::geo;
   bool check_overlaps = false;

   return nullptr;
}
//______________________________________________________________________________

G4VPhysicalVolume * SiVertexTrackerDetectorGeometry::PlacePhysicalVolume(G4LogicalVolume * mother )
{
   // Real world geometry
   BuildLogicalVolumes();

   using namespace eic::geo;
   bool check_overlaps      = false;
   std::string  phys_name   = "";
   G4VPhysicalVolume * phys =  0;

   for(int i = 0; i<fSiVertexTrackerGeometry.NLayers; i++) {
      for(int j = 0; j<fSiVertexTrackerGeometry.fLayerNLadders.at(i); j++) {

         phys_name  = "SiLadder_phys_";
         phys_name += std::to_string(i);
         phys_name += "_";
         phys_name += std::to_string(j);

         new G4PVPlacement(
               G4Transform3D(
                  fSiVertexTrackerGeometry.GetLayerRotation(i,j),
                  fSiVertexTrackerGeometry.GetLayerPosition(i,j)
                  ),
               fLayerContainer_log[i], 
               phys_name.c_str(),
               mother,
               false, 
               fSiVertexTrackerGeometry.GetLadderNumber(i,j),
               check_overlaps);
      }
   }

   return phys;
}
//______________________________________________________________________________


