#include "DriftChamberDetectorGeometry.h"

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
#include <algorithm>

DriftChamberDetectorGeometry::DriftChamberDetectorGeometry()
{
   using namespace clas12::geo;
   using namespace CLHEP;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();
   fSensitiveDetector = new DriftChamberSensitiveDetector("DriftChamber",6*6*6*112);
   SDMan->AddNewDetector(fSensitiveDetector);

   fSensitiveDetector2 = new SensitiveRegionDetector("DriftChamberRegion",6*3);
   SDMan->AddNewDetector(fSensitiveDetector2);

   //std::vector<G4TwoVector> points = RegionTrapPoints(2);
   //std::cout << clas12::geo::RegionTrapWidth(2) << std::endl;
   //std::for_each(points.begin(), points.end(), [](G4TwoVector p){std::cout << p << std::endl;});

   // Huge box for geometric intersection with trapezoid  so that the placement
   // of the region is easier . 
   G4VSolid * temp_box    = new G4Box("DC_placement_box_solid",5*m,5*m,5*m);
   G4VSolid * temp_region = 0;

   // RI
   //temp_region    = new G4Box("R1trap_solid_BOX", 4.0*m, 4.0*m, 1.0*m); 
   temp_region    = new G4GenericTrap("R1trap_solid",  clas12::geo::RegionTrapWidth(1), clas12::geo::RegionTrapPoints(1));
   fRegion1_solid = new G4IntersectionSolid("RI_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::RegionTrapOffset(1)) );

   //temp_box    = new G4Box("clipping_DC_placement_box_solid",1*m,1*m,1*m);
   temp_region    = new G4GenericTrap("clippingR1trap_solid",  clas12::geo::RegionTrapWidth(1), clas12::geo::RegionTrapPoints(1));
   fClippingRegion1_solid = new G4IntersectionSolid("ClippingRI_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::RegionTrapOffset(1)) );

   // RII
   temp_region    = new G4GenericTrap("R2trap_solid",  clas12::geo::RegionTrapWidth(2), clas12::geo::RegionTrapPoints(2));
   fRegion2_solid = new G4IntersectionSolid("RI_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::RegionTrapOffset(2)) );


   temp_region    = new G4GenericTrap("R3trap_solid",  clas12::geo::RegionTrapWidth(3), clas12::geo::RegionTrapPoints(3));
   fRegion3_solid = new G4IntersectionSolid("RI_solid", temp_box, temp_region, 0,  G4ThreeVector(0.0,0.0,clas12::geo::RegionTrapOffset(3)) );

   fRegions_solid[0] = fRegion1_solid;
   fRegions_solid[1] = fRegion2_solid;
   fRegions_solid[2] = fRegion3_solid;

   fClippingRegions_solid[0] = fClippingRegion1_solid;
   fClippingRegions_solid[1] = fClippingRegion1_solid;
   fClippingRegions_solid[2] = fClippingRegion1_solid;
   //fRegion1_solid = new G4GenericTrap("RItrap",  clas12::geo::RegionTrapWidth(1), clas12::geo::RegionTrapPoints(1));
   //fRegion2_solid = new G4GenericTrap("RIItrap", clas12::geo::RegionTrapWidth(2), clas12::geo::RegionTrapPoints(2));
   //fRegion3_solid = new G4GenericTrap("RIIItrap",clas12::geo::RegionTrapWidth(3), clas12::geo::RegionTrapPoints(3));

   fRegion1_log = 0;
   fRegion2_log = 0;
   fRegion3_log = 0;

   fEmptyRegion1_log = 0;
   fEmptyRegion2_log = 0;
   fEmptyRegion3_log = 0;

   G4NistManager * matman = G4NistManager::Instance();
   G4Element     * Ar     = new G4Element("Argon", "Ar", /*z    = */18, /*a            = */ 39.95*g/mole);
   fGasMaterial           = new G4Material("DC_gas", /* density = */ 1.8*mg/cm3, /*nel = */ 3);
   fGasMaterial->AddElement(Ar, 90*perCent);
   fGasMaterial->AddMaterial(matman->FindOrBuildMaterial("G4_O"),  6.6*perCent);
   fGasMaterial->AddMaterial(matman->FindOrBuildMaterial("G4_C"),  3.4*perCent);

}
//______________________________________________________________________________

DriftChamberDetectorGeometry::~DriftChamberDetectorGeometry()
{
}
//______________________________________________________________________________

void DriftChamberDetectorGeometry::BuildLogicalVolumes()
{

   bool check_overlaps = false;

   using namespace CLHEP;
   using namespace clas12::geo;
   using namespace clas12::geo::DC;
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

   fRegion1_log = new G4LogicalVolume(fRegion1_solid, fGasMaterial, "DC_Region1_log");
   fRegion2_log = new G4LogicalVolume(fRegion2_solid, fGasMaterial, "DC_Region2_log");
   fRegion3_log = new G4LogicalVolume(fRegion3_solid, fGasMaterial, "DC_Region3_log");

   fRegion1_log->SetSensitiveDetector(fSensitiveDetector2);
   fRegion2_log->SetSensitiveDetector(fSensitiveDetector2);
   fRegion3_log->SetSensitiveDetector(fSensitiveDetector2);

   fRegion1_log->SetVisAttributes(vs);
   fRegion2_log->SetVisAttributes(vs);
   fRegion3_log->SetVisAttributes(vs);

   fRegions_log[0] = fRegion1_log;
   fRegions_log[1] = fRegion2_log;
   fRegions_log[2] = fRegion3_log;

   fEmptyRegion1_log = new G4LogicalVolume(fRegion1_solid, fGasMaterial, "DC_EmptyRegion1_log");
   fEmptyRegion2_log = new G4LogicalVolume(fRegion2_solid, fGasMaterial, "DC_EmptyRegion2_log");
   fEmptyRegion3_log = new G4LogicalVolume(fRegion3_solid, fGasMaterial, "DC_EmptyRegion3_log");

   //fEmptyRegion1_log->SetSensitiveDetector(fSensitiveDetector2);
   //fEmptyRegion2_log->SetSensitiveDetector(fSensitiveDetector2);
   //fEmptyRegion3_log->SetSensitiveDetector(fSensitiveDetector2);

   fEmptyRegion1_log->SetVisAttributes(vs_odd);
   fEmptyRegion2_log->SetVisAttributes(vs_even);
   fEmptyRegion3_log->SetVisAttributes(vs3);

   fEmptyRegions_log[0] = fEmptyRegion1_log;
   fEmptyRegions_log[1] = fEmptyRegion2_log;
   fEmptyRegions_log[2] = fEmptyRegion3_log;

   G4RotationMatrix* yRot = new G4RotationMatrix;  // Rotates X and Z axes only
   yRot->rotateX(M_PI/3.*rad);                     // Rotates 60 degrees

   G4RotationMatrix* yRot2 = new G4RotationMatrix;  // Rotates X and Z axes only
   yRot2->rotateX(2.0*M_PI/3.*rad);                     // Rotates 120 degrees

   G4RotationMatrix* yRot3 = new G4RotationMatrix;  // Rotates X and Z axes only
   yRot3->rotateX(M_PI*rad/6.0);                     // Rotates 30 degrees

   G4ThreeVector zTrans(0, 0, 0);

   const int WiresPerLayer   = 112;
   const int WiresPerSL      = 6*112;
   const int WiresPerSector  = 6*6*112;
   const int TotalWires      = 6*6*6*112;

   // -----------------------------------------------------------------------------
   for( int super_layer = 1; super_layer <=6; super_layer++) {

      for( int layer = 1; layer<=6; layer++ ) {

         for(int i = 1; i<=112; i++ ) {

            // Has to be really long for some reason, otherwise there is a seg fault...
            G4double hex_length = WireLength(super_layer,layer,i) - 3.0*cm;
            double factor = 1.0/2.0;//TMath::Sqrt(3.0);//2.0;
            double zPlane[] = {-hex_length/2.0,hex_length/2.0};
            double rInner[] = {0.0,0.0};
            double rOuter[] = {factor*LayerWireSpacing[super_layer-1], factor*LayerWireSpacing[super_layer-1] };

            // This operation puts the hex tube in a box with the proper orientation ( with flat sides against in-row adjacent)
            //G4VSolid* unionMoved3 = new G4IntersectionSolid(Form("BoxCylinderMoved3%d",super_layer), subtraction_box3, unionMoved2,      yRot3, zTrans);
            G4VSolid * hex_polyhedra = new G4Polyhedra("hex_polyhedra",0.0,360.0*degree, 6, 2, zPlane, rInner, rOuter );
            //G4VSolid * hex_polyhedra = new G4Box("hex_polyhedra_BOX",factor*LayerWireSpacing[super_layer-1]/2.0, factor*LayerWireSpacing[super_layer-1]/2.0, hex_length/2.0 - 1.0*cm );
            G4VSolid * wire_hex_solid = hex_polyhedra;//unionMoved3;

            int sl_ind     = super_layer-1;
            int lay_ind    = layer-1;
            int wire_ind   = i-1;

            int channel = /*WiresPerSector*sec_ind +*/ WiresPerSL*sl_ind + WiresPerLayer*lay_ind + wire_ind;
            //int channel = 112*(layer-1) + 6*112*((super_layer-1)%2) + (i-1);

            //std::cout << " sl_ind " << sl_ind << std::endl;
            //std::cout << " lay_ind " << lay_ind << std::endl;
            //std::cout << " wire_ind " << wire_ind << std::endl;
            //std::cout << " sec_chan " << channel << std::endl;
            //std::cout << " hex_length " << hex_length << std::endl;

            G4RotationMatrix * aRot  = new G4RotationMatrix();
            G4RotationMatrix * noRot = new G4RotationMatrix();
            (*aRot) = clas12::geo::LayerWireRotation(super_layer);
            aRot->rotateY(90.0*degree);
            //aRot->rotateZ(-clas12::geo::DC::ThetaStereo[sl_ind/2]);
            G4ThreeVector zero_trans(0.0,0.0,0.0);

            G4ThreeVector wire_trans = clas12::geo::ToWireMidPlane(super_layer,layer,i);
            wire_trans -= clas12::geo::WireStereoShift(super_layer,layer,i);

            //G4ThreeVector wire_trans_subtract = 0.1*wire_trans;

            //if(layer==1 && super_layer==1 && i==88){
            //   std::cout << wire_trans <<std::endl;
            //}

            // Wire placement happens when creating the solid by chopping off the ends
            // of the hex tube so that it fits inside the sector/region 
            // NOTE : this appears to have a bug in it!
            //G4VSolid* wire_solid  = new G4IntersectionSolid(
            //      Form("sl%d_%d_%d_solid",super_layer,layer,i), 
            //      fRegions_solid[SuperLayerRegionIndex[super_layer-1]],  
            //      wire_hex_solid, 
            //      aRot,
            //      wire_trans-wire_trans_subtract
            //      );

            G4LogicalVolume  * wire_log = new G4LogicalVolume(
                  wire_hex_solid,//wire_solid,
                  fGasMaterial, 
                  Form("sl%d_%d_%d_log",super_layer,layer,i)
                  );

            wire_log->SetSensitiveDetector( fSensitiveDetector ) ;

            G4double maxStep = 0.1*mm; // forces many steps on the order of the average length for creating ion pair
            //wire_log->SetUserLimits(new G4UserLimits(maxStep));

            //check_overlaps = true;
            //if( super_layer >=5 ) check_overlaps = true;
            //check_overlaps = false;
            //if(super_layer%2 == 0 ) {
            //   if(i>60) {
            //      check_overlaps = true;
            //   }
            //}

            G4VPhysicalVolume * phys = new G4PVPlacement(
                  aRot, wire_trans,
                  wire_log,          // its logical volume
                  Form("sl%d_%d_%d_phys",super_layer,layer,i),
                  fRegions_log[SuperLayerRegionIndex[super_layer-1]],  
                  false,                        // no boolean operations
                  channel,                           // its copy number
                  check_overlaps ); // surface check

            // Only visualize one sector (otherwise painfully slow)
            //if(sector == 1 ) {
            if(TMath::Abs(channel%112-10) <= 5 ) {
               if(super_layer%2 == 0 ) {
                  wire_log->SetVisAttributes(vs_even);
                  //wire_log->SetVisAttributes(G4VisAttributes::GetInvisible());
               } else {
                  wire_log->SetVisAttributes(vs_odd);
                  //wire_log->SetVisAttributes(G4VisAttributes::GetInvisible());
               }
            }else {
               wire_log->SetVisAttributes(G4VisAttributes::GetInvisible());
            }

               //wire_log->SetVisAttributes(vs_odd);
            //} else {
            //   wire_log->SetVisAttributes(G4VisAttributes::GetInvisible());
            //}

         }
      }
   }




}
//______________________________________________________________________________

G4VPhysicalVolume * DriftChamberDetectorGeometry::PlacePhysicalVolume(G4LogicalVolume * mother, int sec, int region )
{
   using namespace clas12::geo;
   int index    = region-1;
   int grouping = (sec-1)*3 + (region-1);

   G4VPhysicalVolume * phys = new G4PVPlacement(
         G4Transform3D(
            RegionRotation(sec,region),
            G4ThreeVector(0,0,20.0*CLHEP::cm)+ RegionTranslation(sec, region)
            ),
         fEmptyRegions_log[index],          // its logical volume
         Form("region%d_phys",region), // its name
         mother,                       // its mother (logical) volume
         false,                        // no boolean operations
         sec,                          // its copy number
         false);                       // check for overlaps

   return phys;
}
//______________________________________________________________________________

G4VPhysicalVolume * DriftChamberDetectorGeometry::PlaceParallelPhysicalVolume(G4LogicalVolume * mother, int sec, int region )
{
   using namespace clas12::geo;
   int index    = region-1;
   int grouping = (sec-1)*3 + (region-1);

   //std::cout << "sec      : " << sec << std::endl;
   //std::cout << "region   : " << region << std::endl;
   //std::cout << "index    : " << index << std::endl;
   //std::cout << "grouping : " << grouping << std::endl;

   G4VPhysicalVolume * phys = new G4PVPlacement(
         G4Transform3D(
            RegionRotation(sec,region),
            G4ThreeVector(0,0,20.0*CLHEP::cm) + RegionTranslation(sec, region)
         ),
         fRegions_log[index],          // its logical volume
         Form("region%d_phys",region), // its name
         mother,                       // its mother (logical) volume
         false,                        // no boolean operations
         sec,                     // its copy number
         false);                        // check for overlaps

   return phys;
}
//______________________________________________________________________________

