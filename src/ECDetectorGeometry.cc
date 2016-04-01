#include "ECDetectorGeometry.h"

#include "G4VisAttributes.hh"
#include "G4GenericTrap.hh"
#include "G4Box.hh"
#include "G4TwoVector.hh"
#include "G4IntersectionSolid.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4UserLimits.hh"
#include "G4Polyhedra.hh"
#include <algorithm>

ECDetectorGeometry::ECDetectorGeometry()
{
   using namespace clas12::geo;
   using namespace CLHEP;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();

   fSensitiveDetector = new ECSensitiveDetector("ECRegion",6);
   SDMan->AddNewDetector(fSensitiveDetector);

   //std::vector<G4TwoVector> points = RegionTrapPoints(2);
   //std::cout << clas12::geo::RegionTrapWidth(2) << std::endl;
   //std::for_each(points.begin(), points.end(), [](G4TwoVector p){std::cout << p << std::endl;});

   // Huge box for geometric intersection with trapezoid  so that the placement
   // of the region is easier . 
   G4VSolid * temp_box    = new G4Box("DC_placement_box_solid",5*m,5*m,5*m);
   G4VSolid * temp_region = 0;

      //  4 _________ 3
      //    \ front /
      //     \     /
      //      \   /
      //       \_/
      //      1   2
   std::cout <<  (y_center(2)) << std::endl;
   std::cout << -dy_center(1) << std::endl;
   std::cout <<  fa_params.at(2) + fa_params.at(4)*double(1-1) << std::endl; 
   std::cout << double(1-1)*w_U(1) << std::endl;
   std::cout << w_U(1) << std::endl;

   fContainer_points = {
      {-fd2_param/2.0 + x_corner_U1(39,1),   1.0*y_corner_U1(39, 1) },
      { fd2_param/2.0 + x_corner_U2(39,1),   1.0*y_corner_U1(39, 1) },
      { fd2_param/2.0 + x_corner_U2(39,36),  1.0*y_corner_U2(39,36) },
      {-fd2_param/2.0 + x_corner_U1(39,36),  1.0*y_corner_U2(39,36) },

      {-fd2_param/2.0 + x_corner_U1(39,1),  1.0*y_corner_U1(39, 1) },
      { fd2_param/2.0 + x_corner_U2(39,1),  1.0*y_corner_U1(39, 1) },
      { fd2_param/2.0 + x_corner_U2(39,36), 1.0*y_corner_U2(39,36) },
      {-fd2_param/2.0 + x_corner_U1(39,36), 1.0*y_corner_U2(39,36) }
   };
   for(auto p: fContainer_points) {
      std::cout << " x = " << p.x()/cm << ", y = " << p.y()/cm << " cm \n";
   }

   // ---------------------------------------------
   // Region 1
   fContainer_solid                = new G4GenericTrap("R1Containertrap_solid",  39.0*fLayerHalfThickness, fContainer_points);
   //fRegionContainers_solid[0] = new G4IntersectionSolid("R1Container_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::ContainerTrapWidth(1)) );

   //temp_region       = new G4GenericTrap("R1trap_solid",  clas12::geo::RegionTrapWidth(1), clas12::geo::RegionTrapPoints(1));
   //fRegions_solid[0]       = new G4IntersectionSolid("R1_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::RegionTrapOffset(1)) );

   ////temp_box    = new G4Box("clipping_DC_placement_box_solid",1*m,1*m,1*m);
   ////temp_region    = new G4GenericTrap("clippingR1trap_solid",  clas12::geo::RegionTrapWidth(1), clas12::geo::RegionTrapPoints(1));
   ////fClippingRegion1_solid = new G4IntersectionSolid("ClippingRI_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::RegionTrapOffset(1)) );

   //// ---------------------------------------------
   //// Region 2
   //temp_region                = new G4GenericTrap("R2Containertrap_solid",  clas12::geo::ContainerTrapWidth(2), clas12::geo::ContainerTrapPoints(2));
   //fRegionContainers_solid[1] = new G4IntersectionSolid("R2Container_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::ContainerTrapWidth(2)) );

   //temp_region    = new G4GenericTrap("R2trap_solid",  clas12::geo::RegionTrapWidth(2), clas12::geo::RegionTrapPoints(2));
   //fRegions_solid[1] = new G4IntersectionSolid("RI_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::RegionTrapOffset(2)) );

   //// ---------------------------------------------
   //// Region 3
   //temp_region                = new G4GenericTrap("R3Containertrap_solid",  clas12::geo::ContainerTrapWidth(3), clas12::geo::ContainerTrapPoints(3));
   //fRegionContainers_solid[2] = new G4IntersectionSolid("R3Container_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::ContainerTrapWidth(3)) );

   //temp_region    = new G4GenericTrap("R3trap_solid",  clas12::geo::RegionTrapWidth(3), clas12::geo::RegionTrapPoints(3));
   //fRegions_solid[2] = new G4IntersectionSolid("RI_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::RegionTrapOffset(3)) );

   //// ---------------------------------------------
   //// Endplates

   //G4VSolid * endplate_test  = 0;//new G4Box("endplatetest", 1.0*cm/2.0, 30.0*cm/2.0, 5.0*cm/2.0);
   //fEndPlates_solid[0] = new G4GenericTrap("R1Endplate_solid",  clas12::geo::DC::EndPlateThickness.at(0)/2.0, clas12::geo::EndplateTrapPoints(1));
   //fEndPlates_solid[1] = new G4GenericTrap("R2Endplate_solid",  clas12::geo::DC::EndPlateThickness.at(1)/2.0, clas12::geo::EndplateTrapPoints(2));
   //fEndPlates_solid[2] = new G4GenericTrap("R3Endplate_solid",  clas12::geo::DC::EndPlateThickness.at(2)/2.0, clas12::geo::EndplateTrapPoints(3));

   //// ---------------------------------------------
   //// Materials
   G4NistManager * matman = G4NistManager::Instance();
   G4Element     * Ar     = new G4Element("Argon", "Ar", /*z    = */18, /*a            = */ 39.95*g/mole);
   fGasMaterial           = new G4Material("DC_gas", /* density = */ 1.8*mg/cm3, /*nel = */ 3);
   fGasMaterial->AddElement(Ar, 90*perCent);
   fGasMaterial->AddMaterial(matman->FindOrBuildMaterial("G4_O"),  6.6*perCent);
   fGasMaterial->AddMaterial(matman->FindOrBuildMaterial("G4_C"),  3.4*perCent);
}
//______________________________________________________________________________

ECDetectorGeometry::~ECDetectorGeometry()
{
}
//______________________________________________________________________________

void ECDetectorGeometry::BuildLogicalVolumes()
{

   bool check_overlaps = false;

   using namespace CLHEP;
   using namespace clas12::geo;
   using namespace clas12::geo::DC;

   G4NistManager * matman = G4NistManager::Instance();
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

   // ------------------------------
   // Region Containers
   fContainer_log = new G4LogicalVolume(fContainer_solid, fGasMaterial, "EC_Container_log_0");
   fContainer_log->SetSensitiveDetector(fSensitiveDetector);
   fContainer_log->SetVisAttributes(vs);

}
//______________________________________________________________________________

G4VPhysicalVolume * ECDetectorGeometry::PlacePhysicalVolume(G4LogicalVolume * mother, int sec, int region )
{
   using namespace clas12::geo;
   int index    = region-1;
   int grouping = (sec-1)*3 + (region-1);

   //fEndPlates_log[0] = new G4LogicalVolume(fEndPlates_solid[0], fGasMaterial, "EndPlates_log");

   G4VPhysicalVolume * phys = nullptr;

   //G4VPhysicalVolume * phys = new G4PVPlacement(
   //      G4Transform3D(
   //         RegionRotation(sec,region),
   //         G4ThreeVector(0,0,20.0*CLHEP::cm)+ RegionTranslation(sec, region)
   //         ),
   //      fEmptyRegions_log[index],          // its logical volume
   //      Form("region%d_phys",region), // its name
   //      mother,                       // its mother (logical) volume
   //      false,                        // no boolean operations
   //      sec,                          // its copy number
   //      false);                       // check for overlaps

   return phys;
}
//______________________________________________________________________________

G4VPhysicalVolume * ECDetectorGeometry::PlaceParallelPhysicalVolume(G4LogicalVolume * mother, int sec, int region )
{
   using namespace clas12::geo;
   int index    = region-1;
   int grouping = (sec-1)*3 + (region-1);

   //std::cout << "sec      : " << sec << std::endl;
   //std::cout << "region   : " << region << std::endl;
   //std::cout << "index    : " << index << std::endl;
   //std::cout << "grouping : " << grouping << std::endl;
   double angle = double(sec-1)*60.0*CLHEP::degree;
   G4RotationMatrix rot2 = G4RotationMatrix::IDENTITY;
   rot2.rotateZ(-90.0*CLHEP::degree + angle);
   G4RotationMatrix rot = G4RotationMatrix::IDENTITY;
   rot.rotateX(-25.0*CLHEP::degree);
   rot.rotateZ(-90.0*CLHEP::degree + angle);
   G4VPhysicalVolume * phys = new G4PVPlacement(
         G4Transform3D(
            rot,
            rot2*(fL1_pos + fS_pos)
         ),
         fContainer_log,          // its logical volume
         Form("region%d_phys",region), // its name
         mother,                       // its mother (logical) volume
         false,                        // no boolean operations
         sec,                     // its copy number
         false);                        // check for overlaps

   return phys;
}
//______________________________________________________________________________

