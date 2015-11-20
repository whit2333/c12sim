#include "RecoilChamberDetectorGeometry.h"
#include <cmath>

#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4GenericTrap.hh"
#include "G4Box.hh"
#include "G4IntersectionSolid.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4UserLimits.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTrap.hh"

#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"

#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "FakeSD.hh"
#include "G4ExtrudedSolid.hh"


RecoilChamberDetectorGeometry::RecoilChamberDetectorGeometry()
{
   using namespace clas12::geo;
   using namespace CLHEP;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();

   fSensitiveDetector = new RecoilChamberSensitiveDetector("RecoilChamber",6*6*6*112);
   SDMan->AddNewDetector(fSensitiveDetector);

   for(int i = 0; i< fNLayers;i++){

      // The inner radius for a given layer of wires
      double layer_innner_radius = (innerRadiusOfTheGasDetector + DeltaR + NsLay*DeltaR*i);
      int    NWiresLay           = int(CLHEP::pi*layer_innner_radius / DeltaP);
      NWiresLay                  = 2*NWiresLay; // to ensure an even number of cells
      double PhiWire             = (2.*CLHEP::pi/ NWiresLay);
      double sign = 1.0;
      if(i%2==0) sign=1.;
      else sign=-1.;

      fLayerRadius[i]     = layer_innner_radius;
      fLayerNwires[i]     = NWiresLay;
      fDeltaPhi[i]        = PhiWire;
      fLayerSteroAngle[i] = steAng*sign;
   }
}
//______________________________________________________________________________

RecoilChamberDetectorGeometry::~RecoilChamberDetectorGeometry()
{
}
//______________________________________________________________________________

G4ThreeVector RecoilChamberDetectorGeometry::GetIntersectionPoint(
      const G4ThreeVector x0, const G4ThreeVector x1,
      const G4ThreeVector p0, const G4ThreeVector norm )
{

   G4ThreeVector w = x0 - p0;
   G4ThreeVector u = x1 - x0;
   double si = -1.0*(norm.dot(w)/(norm.dot(u)));

   G4ThreeVector sterm = si*u;
   G4ThreeVector Psi = p0 + w + sterm  ;

   return Psi;
}
//______________________________________________________________________________


G4VSolid * RecoilChamberDetectorGeometry::BuildWireSolid(int layer, int radial_wire_number)
{

   // The first layer is offset from the nominal inner radius by DeltaP
   // union of four pieces where the sensitive wire is in the middle
   //  _______ 
   //  \__|__/
   //   \_|_/


   double sign = 1.0;
   if(layer%2==0) sign=1.;
   else sign=-1.;

   double overlap_size = 0.1*mm;
   double PhiWire = fDeltaPhi.at(layer);

   //std::cout << "Layer  " << layer << std::endl;
   //std::cout << "radial " << radial_wire_number << std::endl;

   G4double Rlay  = fLayerRadius.at(layer) + DeltaR*double(radial_wire_number);
   G4double Rlay2 = fLayerRadius.at(layer) + DeltaR*double(radial_wire_number+1);

   //std::cout << "Rlay " << Rlay << std::endl;
   //std::cout << "Rlay2 " << Rlay2 << std::endl;
   
   G4double xw0_L = (Rlay * cos(-0.5*PhiWire));
   G4double yw0_L = (Rlay * sin(-0.5*PhiWire));
   G4double xw1_L = (Rlay * cos(0.5*PhiWire));
   G4double yw1_L = (Rlay * sin(0.5*PhiWire));

   //std::cout << "xw0_L " << xw0_L << std::endl;
   //std::cout << "yw0_L " << yw0_L << std::endl;
   //std::cout << "xw1_L " << xw1_L << std::endl;
   //std::cout << "yw1_L " << yw1_L << std::endl;

   G4double xv0_L = (Rlay2 * cos(-0.5*PhiWire));
   G4double yv0_L = (Rlay2 * sin(-0.5*PhiWire));
   G4double xv1_L = (Rlay2 * cos(0.5*PhiWire));
   G4double yv1_L = (Rlay2 * sin(0.5*PhiWire));

   //std::cout << "xv0_L " << xv0_L << std::endl;
   //std::cout << "yv0_L " << yv0_L << std::endl;
   //std::cout << "xv1_L " << xv1_L << std::endl;
   //std::cout << "yv1_L " << yv1_L << std::endl;

   G4ThreeVector p0_L { xw0_L, yw0_L, -lengthOfTheWire/2.0 };
   G4ThreeVector p1_L { xw1_L, yw1_L, -lengthOfTheWire/2.0 };
   G4ThreeVector p2_L { xv0_L, yv0_L, -lengthOfTheWire/2.0 };
   G4ThreeVector p3_L { xv1_L, yv1_L, -lengthOfTheWire/2.0 };
   G4ThreeVector p4_L { xw0_L, yw0_L,  lengthOfTheWire/2.0 };
   G4ThreeVector p5_L { xw1_L, yw1_L,  lengthOfTheWire/2.0 };
   G4ThreeVector p6_L { xv0_L, yv0_L,  lengthOfTheWire/2.0 };
   G4ThreeVector p7_L { xv1_L, yv1_L,  lengthOfTheWire/2.0 };

   G4ThreeVector pivot_w0_L { xw0_L, yw0_L, 0.0 };
   G4ThreeVector pivot_w1_L { xw1_L, yw1_L, 0.0 };
   G4ThreeVector pivot_v0_L { xv0_L, yv0_L, 0.0 };
   G4ThreeVector pivot_v1_L { xv1_L, yv1_L, 0.0 };

   G4ThreeVector endplate_0 { 0,0, -lengthOfTheWire/2.0};
   G4ThreeVector endplate_0_norm { 0,0, 1.0};

   G4ThreeVector endplate_1 { 0,0,  lengthOfTheWire/2.0};
   G4ThreeVector endplate_1_norm { 0,0, 1.0};

   G4RotationMatrix rot0;
   rot0.rotate( fLayerSteroAngle.at(layer), G4ThreeVector(xw0_L, yw0_L,0.) );

   G4RotationMatrix rot1;
   rot1.rotate( fLayerSteroAngle.at(layer), G4ThreeVector(xw1_L, yw1_L,0.) );

   G4ThreeVector q0_L = GetIntersectionPoint( pivot_w0_L, rot0*p0_L, endplate_0, endplate_0_norm);
   G4ThreeVector q1_L = GetIntersectionPoint( pivot_w1_L, rot1*p1_L, endplate_0, endplate_0_norm);
   G4ThreeVector q2_L = GetIntersectionPoint( pivot_v0_L, rot0*p2_L, endplate_0, endplate_0_norm);
   G4ThreeVector q3_L = GetIntersectionPoint( pivot_v1_L, rot1*p3_L, endplate_0, endplate_0_norm);
   G4ThreeVector q4_L = GetIntersectionPoint( pivot_w0_L, rot0*p4_L, endplate_1, endplate_1_norm);
   G4ThreeVector q5_L = GetIntersectionPoint( pivot_w1_L, rot1*p5_L, endplate_1, endplate_1_norm);
   G4ThreeVector q6_L = GetIntersectionPoint( pivot_v0_L, rot0*p6_L, endplate_1, endplate_1_norm);
   G4ThreeVector q7_L = GetIntersectionPoint( pivot_v1_L, rot1*p7_L, endplate_1, endplate_1_norm);

   G4TwoVector pp0_L { q0_L.x(), q0_L.y() };
   G4TwoVector pp1_L { q1_L.x(), q1_L.y() };
   G4TwoVector pp2_L { q2_L.x(), q2_L.y() };
   G4TwoVector pp3_L { q3_L.x(), q3_L.y() };
   G4TwoVector pp4_L { q4_L.x(), q4_L.y() };
   G4TwoVector pp5_L { q5_L.x(), q5_L.y() };
   G4TwoVector pp6_L { q6_L.x(), q6_L.y() };
   G4TwoVector pp7_L { q7_L.x(), q7_L.y() };

   //std::cout << "pp0_L " << pp0_L <<std::endl;
   //std::cout << "pp1_L " << pp1_L <<std::endl;
   //std::cout << "pp2_L " << pp2_L <<std::endl;
   //std::cout << "pp3_L " << pp3_L <<std::endl;
   //std::cout << "pp4_L " << pp4_L <<std::endl;
   //std::cout << "pp5_L " << pp5_L <<std::endl;
   //std::cout << "pp6_L " << pp6_L <<std::endl;
   //std::cout << "pp7_L " << pp7_L <<std::endl;

   std::vector<G4TwoVector> trap_points_L = {
      pp0_L, pp1_L, pp3_L,  pp2_L, 
      pp4_L, pp5_L, pp7_L, pp6_L
   };

   G4VSolid * trap_1 = new  G4GenericTrap("trap1", lengthOfTheWire/2.0, trap_points_L);

   return trap_1;
}
//______________________________________________________________________________

void RecoilChamberDetectorGeometry::BuildUnitCells()
{

   for(G4int tlay=0; tlay<fNLayers; tlay++){
      // Determine the radius of this layer of wires
      // The first layer is offset from the nominal inner radius by DeltaP
      G4double Rtlay = fLayerRadius[tlay];//(innerRadiusOfTheGasDetector + DeltaR + NsLay*DeltaR*tlay);
      //G4int      NWiresLay = int(CLHEP::pi*Rtlay / DeltaP);
      //NWiresLay = 2*NWiresLay; // to ensure an even number of wires
      //fLayerNwires[tlay] = NWiresLay;
      G4double PhiWire = fDeltaPhi[tlay];//(2.*CLHEP::pi/ NWiresLay);

      // union of four pieces where the sensitive wire is in the middle
      //  _______ 
      //  \__|__/
      //   \_|_/

      G4Tubs   * container     = new G4Tubs("container", 0,outerRadiusOfTheGasDetector, hightOfTheGasDetector, 0.0, 360.0*deg);
      G4VSolid * combined      = 0;
      G4VSolid * combined_temp = 0;
      G4VSolid * combined_temp1 = 0;
      G4VSolid * combined_temp2 = 0;

      double sign = 1.0;
      if(tlay%2==0) sign=1.;
      else sign=-1.;

      fWireVolume_solid[tlay][0] = BuildWireSolid(tlay,0);
      fWireVolume_solid[tlay][1] = BuildWireSolid(tlay,0);
      fWireVolume_solid[tlay][2] = BuildWireSolid(tlay,1);
      fWireVolume_solid[tlay][3] = BuildWireSolid(tlay,1);


      for(int i = 0; i<4; i++){
         std::string log_name = "wire_volume_log_" + std::to_string(tlay) + "_part" + std::to_string(i);
         fWireVolume_log[tlay][i]   = new G4LogicalVolume(fWireVolume_solid[tlay][i], HeiC4H10, log_name);
      }
   }
}
//______________________________________________________________________________
void  RecoilChamberDetectorGeometry::PlaceCells(G4LogicalVolume * mother, int layer, double z_rotation, int wire_number ){

   bool check_overlaps = false;

   G4RotationMatrix* rotationMatrix1 = 0;
   G4RotationMatrix* rotationMatrix2 = 0;
   G4double PhiWire = fDeltaPhi[layer];

   rotationMatrix1 = new G4RotationMatrix();
   rotationMatrix1->rotateZ( -0.5*PhiWire + z_rotation);
   fWireVolume_rot[layer][0] = rotationMatrix1;
   fWireVolume_rot[layer][2] = rotationMatrix1;

   rotationMatrix2 = new G4RotationMatrix();
   rotationMatrix2->rotateZ( 0.5*PhiWire + z_rotation);
   fWireVolume_rot[layer][1] = rotationMatrix2;
   fWireVolume_rot[layer][3] = rotationMatrix2;

   G4VisAttributes * vs = 0;
   if(layer%2==0) {
      vs = new G4VisAttributes(G4Colour(0.0,0.9,0.1));
   } else {
      vs = new G4VisAttributes(G4Colour(0.9,0.0,0.1));
   }
   vs->SetForceWireframe(true);
   //vs->SetForceSolid(true);

   for(int i = 0; i<4; i++) {
      std::string wire_name = "wire_volume_" + std::to_string(wire_number) + "_part" + std::to_string(i);
      fWireVolume_log[layer][i]->SetSensitiveDetector(fSensitiveDetector);
      fWireVolume_log[layer][i]->SetVisAttributes(vs);//G4VisAttributes::GetInvisible());

      new G4PVPlacement(G4Transform3D(*(fWireVolume_rot[layer][i]),G4ThreeVector(0,0,0.)),
            fWireVolume_log[layer][i], wire_name,
            mother,
            false,
            wire_number,
            check_overlaps);
   }
}
//______________________________________________________________________________


void RecoilChamberDetectorGeometry::BuildLogicalVolumes()
{
   using namespace CLHEP;
   using namespace clas12::geo;
   using namespace clas12::geo::DC;
   //G4VisAttributes * vs = new G4VisAttributes(G4Colour(0.0,0.9,0.1));
}
//______________________________________________________________________________

G4VPhysicalVolume * RecoilChamberDetectorGeometry::PlaceParallelPhysicalVolume(G4LogicalVolume * mother )
{
   using namespace clas12::geo;

   bool check_overlaps = false;

   BuildUnitCells();

   //--Creating the geometry
   G4ThreeVector gasDetector_pos = G4ThreeVector(gasDetector_posX, gasDetector_posY, gasDetector_posZ);

   G4Tubs* gasDetector = new G4Tubs("GasDetectorPara", 
         innerRadiusOfTheGasDetector,outerRadiusOfTheGasDetector, 
         hightOfTheGasDetector, 
         startAngleOfTheGasDetector, spanningAngleOfTheGasDetector);

   G4LogicalVolume* logicGasDetector =                         
      new G4LogicalVolume(gasDetector, HeiC4H10, "GasDetectorPara_log");            

   G4VPhysicalVolume * gasDetector_phys = new G4PVPlacement(0, gasDetector_pos, logicGasDetector, 
         "gasDetectorPara_phys", mother, false, 0, true);            

   G4VisAttributes * GasDetectorVisAtt = new G4VisAttributes(G4Colour(0.3,0.1,0.1));
   GasDetectorVisAtt->SetForceWireframe(true);
   //GasDetectorVisAtt->SetDaughtersInvisible(true);
   logicGasDetector->SetVisAttributes(GasDetectorVisAtt);

   //--Step max
   G4double maxStep = 2.0*mm;
   G4UserLimits * fStepLimit = new G4UserLimits(maxStep);
   fStepLimit->SetMaxAllowedStep(maxStep);
   //   fStepLimit->SetUserMaxTrackLength(maxStep);
   logicGasDetector->SetUserLimits(fStepLimit);

   //--Making it a sensitive detector
   G4SDManager* SDman = G4SDManager::GetSDMpointer();
   G4String gasDetectorSDname = "/mydet/GasDetector";
   //--

   int wire_number = 0;

   for(G4int tlay=0; tlay<fNLayers; tlay++){

      bool placed_one = false;
      // Determine the radius of this layer of wires
      // The first layer is offset from the nominal inner radius by DeltaP
      G4double   Rtlay = fLayerRadius[tlay];//(innerRadiusOfTheGasDetector + DeltaP + NsLay*DeltaP*tlay);
      G4int      NWiresLay = fLayerNwires[tlay];//int(CLHEP::pi*Rtlay / DeltaP);
      //NWiresLay = 2*NWiresLay; // to ensure an even number of wires
      G4double PhiWire = fDeltaPhi[tlay];//(2.*CLHEP::pi/ NWiresLay);

      for(G4int slay=0;slay<NsLay;slay++){

         //G4double Rlay = Rtlay + DeltaR*slay;

         for(int wi=0;wi<NWiresLay/2.0;wi++){  
            

            if( (!placed_one) ) {

               PlaceCells( logicGasDetector, tlay, 2.0*double(wi)*PhiWire, wire_number );

               //placed_one=true;
               //std::string wire_name = "wire_volume_" + std::to_string(wire_number);

               //new G4PVPlacement(G4Transform3D(*rotationMatrix,G4ThreeVector(0,0,0.)),
               //      fWireVolume_log[tlay], wire_name, logicGasDetector,
               //      false,wire_number,check_overlaps);
               wire_number++;
            }
         }   
      } // slay
   } // tlay

   return gasDetector_phys;
}
//______________________________________________________________________________

G4VPhysicalVolume * RecoilChamberDetectorGeometry::PlacePhysicalVolume(G4LogicalVolume * mother )
{
   // Real world geometry

   using namespace clas12::geo;

   G4ThreeVector gasDetector_pos = G4ThreeVector(gasDetector_posX, gasDetector_posY, gasDetector_posZ);

   G4Tubs* gasDetector = new G4Tubs("GasDetector", 
         innerRadiusOfTheGasDetector,outerRadiusOfTheGasDetector, 
         hightOfTheGasDetector, 
         startAngleOfTheGasDetector, spanningAngleOfTheGasDetector);

   G4LogicalVolume* logicGasDetector =                         
      new G4LogicalVolume(gasDetector, HeiC4H10, "GasDetector");            

   G4VPhysicalVolume * gasDetector_phys = new G4PVPlacement(0, gasDetector_pos, logicGasDetector, 
         "GasDetector", mother, false, 0, true);            

   G4VisAttributes * GasDetectorVisAtt = new G4VisAttributes(G4Colour(0.3,0.1,0.1));
   GasDetectorVisAtt->SetForceWireframe(true);
   logicGasDetector->SetVisAttributes(GasDetectorVisAtt);

   //--Step max
   G4double maxStep = 2.0*mm;
   G4UserLimits * fStepLimit = new G4UserLimits(maxStep);
   fStepLimit->SetMaxAllowedStep(maxStep);
   //   fStepLimit->SetUserMaxTrackLength(maxStep);
   logicGasDetector->SetUserLimits(fStepLimit);

   //--Making it a sensitive detector
   G4SDManager* SDman = G4SDManager::GetSDMpointer();

   G4String gasDetectorSDname = "/mydet/GasDetector";
   //GasDetectorSD * gasDetectorSD = new GasDetectorSD(gasDetectorSDname);
   //SDman->AddNewDetector(gasDetectorSD);
   //logicGasDetector->SetSensitiveDetector(gasDetectorSD);
   //--

   // --Placing the wires
   G4Tubs* wiresolid = new G4Tubs("Wire", innerRadiusOfTheWire,outerRadiusOfTheWire, lengthOfTheWire/2.0, startAngleOfTheWire, spanningAngleOfTheWire);
   G4LogicalVolume* wires_even_log = new G4LogicalVolume(wiresolid, Tungsten, "wires_even_log");       
   G4LogicalVolume* wires_odd_log = new G4LogicalVolume(wiresolid, Tungsten, "wires_odd_log");       

			
   //bool checkOverlaps = false;
   //double sign = 1.0;
   //int wire_element;
   //for(G4int tlay=0;tlay<NTLay;tlay++){

   //   G4double Rtlay = (innerRadiusOfTheGasDetector + DeltaP + NsLay*DeltaP*tlay)*mm;
   //   G4int NWiresLay = (G4int) (CLHEP::pi*Rtlay / DeltaP);
   //   NWiresLay = 2*NWiresLay; // to ensure an even number of wires
   //   G4double PhiWire = (2.*CLHEP::pi/ NWiresLay)*rad;

   //   for(G4int slay=0;slay<NsLay;slay++){

   //      G4double Rlay = Rtlay + DeltaR*slay;

   //      for(G4double wi=0;wi<NWiresLay;wi++){  
   //         G4double xw = (Rlay * cos(wi*PhiWire))*mm;
   //         G4double yw = (Rlay * sin(wi*PhiWire))*mm;

   //         if(tlay%2==0) sign=1.;
   //         else sign=-1.;
   //         G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
   //         rotationMatrix->rotate(sign*steAng,G4ThreeVector(xw,yw,0.));

   //         G4LogicalVolume * logicWires = wires_odd_log;
   //         if(tlay%2 == 0 )  logicWires = wires_even_log;

   //         std::string wire_name = "wire_element_" + std::to_string(wire_element);

   //         new G4PVPlacement(rotationMatrix,G4ThreeVector(xw,yw,0.),logicWires,
   //               wire_name,logicGasDetector,
   //               false,(tlay+1)*10000+(slay+1)*1000+wi,checkOverlaps);
   //         wire_element++;
   //      } // wi		     
   //   } // slay
   //} // tlay
   ////--

   //G4VisAttributes * wires_even_vis = new G4VisAttributes(G4Colour(1.,1.0,0.0));
   //wires_even_vis->SetForceSolid(true);
   //wires_even_log->SetVisAttributes(wires_even_vis);

   //G4VisAttributes * wires_odd_vis = new G4VisAttributes(G4Colour(140./256.,222./256.,242./256.));
   //wires_odd_vis->SetForceSolid(true);
   //wires_odd_log->SetVisAttributes(wires_odd_vis);

   return gasDetector_phys;
}
//______________________________________________________________________________


