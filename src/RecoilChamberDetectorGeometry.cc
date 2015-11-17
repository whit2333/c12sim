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

}
//______________________________________________________________________________

RecoilChamberDetectorGeometry::~RecoilChamberDetectorGeometry()
{
}
//______________________________________________________________________________

void RecoilChamberDetectorGeometry::BuildUnitCells()
{

   for(G4int tlay=0; tlay<NTLay; tlay++){
      // Determine the radius of this layer of wires
      // The first layer is offset from the nominal inner radius by DeltaP
      G4double Rtlay = (innerRadiusOfTheGasDetector + DeltaP + NsLay*DeltaP*tlay);
      fLayerRadius[tlay] = Rtlay;
      G4int      NWiresLay = int(CLHEP::pi*Rtlay / DeltaP);
      NWiresLay = 2*NWiresLay; // to ensure an even number of wires
      fLayerNwires[tlay] = NWiresLay;
      G4double PhiWire = (2.*CLHEP::pi/ NWiresLay);
      fDeltaPhi[tlay] = PhiWire;

      // union of four pieces where the sensitive wire is in the middle
      //  _______ 
      //  \__|__/
      //   \_|_/

      G4Tubs   * container     = new G4Tubs("container", 0,outerRadiusOfTheGasDetector, hightOfTheGasDetector, 0.0, 360.0*deg);
      G4VSolid * combined      = 0;
      G4VSolid * combined_temp = 0;
      G4VSolid * combined_temp1 = 0;
      G4VSolid * combined_temp2 = 0;

      //G4UnionSolid * solid = 0;
      for(G4int slay=0;slay<NsLay-1;slay++){

         double sign = 1.0;
         if(tlay%2==0) sign=1.;
         else sign=-1.;

         G4double Rlay  = Rtlay + DeltaR*slay;
         G4double Rlay2 = Rtlay + DeltaR*(slay+1);

         G4double xw0 = (Rlay * cos(-0.5*PhiWire));
         G4double yw0 = (Rlay * sin(-0.5*PhiWire));
         G4double xw1 = (Rlay * cos(0.5*PhiWire));
         G4double yw1 = (Rlay * sin(0.5*PhiWire));

         G4double xv0 = (Rlay2 * cos(-0.5*PhiWire));
         G4double yv0 = (Rlay2 * sin(-0.5*PhiWire));
         G4double xv1 = (Rlay2 * cos(0.5*PhiWire));
         G4double yv1 = (Rlay2 * sin(0.5*PhiWire));

         G4double xc0 = (Rlay + DeltaR/2.0)*cos(-0.5*PhiWire);
         G4double yc0 = (Rlay + DeltaR/2.0)*sin(-0.5*PhiWire);
         G4double xc1 = (Rlay + DeltaR/2.0)*cos( 0.5*PhiWire);
         G4double yc1 = (Rlay + DeltaR/2.0)*sin( 0.5*PhiWire);

         G4ThreeVector p0 { xw0, yw0, -lengthOfTheWire/2.0 };
         G4ThreeVector p1 { xw1, yw1, -lengthOfTheWire/2.0 };
         G4ThreeVector p2 { xv0, yv0, -lengthOfTheWire/2.0 };
         G4ThreeVector p3 { xv1, yv1, -lengthOfTheWire/2.0 };
         G4ThreeVector p4 { xw0, yw0,  lengthOfTheWire/2.0 };
         G4ThreeVector p5 { xw1, yw1,  lengthOfTheWire/2.0 };
         G4ThreeVector p6 { xv0, yv0,  lengthOfTheWire/2.0 };
         G4ThreeVector p7 { xv1, yv1,  lengthOfTheWire/2.0 };

         G4RotationMatrix rot0;
         rot0.rotate( sign*steAng, G4ThreeVector(xw0, yw0,0.) );

         G4RotationMatrix rot1;
         rot1.rotate( sign*steAng, G4ThreeVector(xw1, yw1,0.) );

         G4ThreeVector q0 = rot0*p0;
         G4ThreeVector q1 = rot1*p1;
         G4ThreeVector q2 = rot0*p2;
         G4ThreeVector q3 = rot1*p3;
         G4ThreeVector q4 = rot0*p4;
         G4ThreeVector q5 = rot1*p5;
         G4ThreeVector q6 = rot0*p6;
         G4ThreeVector q7 = rot1*p7;

         G4TwoVector pp0 { q0.x(), q0.y() };
         G4TwoVector pp1 { q1.x(), q1.y() };
         G4TwoVector pp2 { q2.x(), q2.y() };
         G4TwoVector pp3 { q3.x(), q3.y() };
         G4TwoVector pp4 { q4.x(), q4.y() };
         G4TwoVector pp5 { q5.x(), q5.y() };
         G4TwoVector pp6 { q6.x(), q6.y() };
         G4TwoVector pp7 { q7.x(), q7.y() };
         //std::cout << q0 << "\n";
         //std::cout << q1 << "\n";
         //std::cout << q2 << "\n";
         //std::cout << q3 << "\n";

         std::vector<G4TwoVector> trap_points = { pp0, pp1, pp3,  pp2, pp4, pp5, pp7,pp6 };

         G4VSolid * trap_1 = new  G4GenericTrap("trap1", lengthOfTheWire/2.0, trap_points);
         G4VSolid * trap_2 = new  G4GenericTrap("trap2", lengthOfTheWire/2.0, trap_points);

         G4RotationMatrix* rotationMatrix1 = 0;
         G4RotationMatrix* rotationMatrix2 = 0;

         rotationMatrix1 = new G4RotationMatrix();
         rotationMatrix1->rotateZ( -0.5*PhiWire);

         rotationMatrix2 = new G4RotationMatrix();
         rotationMatrix2->rotateZ( 0.5*PhiWire);

         if( combined == 0 ) {

            combined      = new G4IntersectionSolid("combined",container, trap_1, rotationMatrix1, G4ThreeVector(0,0,0.)); 
            combined_temp = new G4UnionSolid("combined2",combined,        trap_2, rotationMatrix2, G4ThreeVector(0,0,0.)); 

         } else {

            combined_temp1 = new G4UnionSolid("combined_temp1",combined_temp, trap_1, rotationMatrix1, G4ThreeVector(0,0,0.)); 
            combined_temp2 = new G4UnionSolid("combined_temp2",combined_temp1, trap_2, rotationMatrix2, G4ThreeVector(0,0,0.)); 

         }

      } // slay

      fWireVolume_solid[tlay] = combined_temp2;
      fWireVolume_log[tlay]   = new G4LogicalVolume(fWireVolume_solid[tlay], HeiC4H10, "Unitcell");
      fWireVolume_log[tlay]->SetVisAttributes(G4VisAttributes::GetInvisible());
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
   GasDetectorVisAtt->SetDaughtersInvisible(true);
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
   double sign = 1.0;
   for(G4int tlay=0; tlay<NTLay; tlay++){

      // Determine the radius of this layer of wires
      // The first layer is offset from the nominal inner radius by DeltaP
      G4double   Rtlay = (innerRadiusOfTheGasDetector + DeltaP + NsLay*DeltaP*tlay);
      G4int      NWiresLay = int(CLHEP::pi*Rtlay / DeltaP);
      NWiresLay = 2*NWiresLay; // to ensure an even number of wires
      G4double PhiWire = (2.*CLHEP::pi/ NWiresLay);

      for(G4int slay=0;slay<NsLay;slay++){

         G4double Rlay = Rtlay + DeltaR*slay;

         for(int wi=0;wi<NWiresLay;wi++){  

            G4double xw = (Rlay * cos(double(wi)*PhiWire));
            G4double yw = (Rlay * sin(double(wi)*PhiWire));

            if(tlay%2==0) sign=1.;
            else sign=-1.;

            G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
            rotationMatrix->rotateZ(-1.0*double(wi)*PhiWire);
            //rotationMatrix->rotate(sign*steAng,G4ThreeVector(xw,yw,0.));

            //G4LogicalVolume * logicWires = wires_odd_log;
            //if(tlay%2 == 0 )  logicWires = wires_even_log;

            if( ((wi+1)%2 == 0) && ( (slay+1)%2==0 ) ) {

               std::string wire_name = "wire_volume_" + std::to_string(wire_number);

               fWireVolume_log[tlay]->SetSensitiveDetector(fSensitiveDetector);
               new G4PVPlacement(rotationMatrix,G4ThreeVector(0,0,0.),
                     fWireVolume_log[tlay], wire_name, logicGasDetector,
                     false,wire_number,false);
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


