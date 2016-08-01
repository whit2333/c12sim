#include "RecoilChamberDetectorGeometry.h"

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


RecoilChamberDetectorGeometry::RecoilChamberDetectorGeometry()
{
   using namespace clas12::geo;
   using namespace CLHEP;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();

   fSensitiveDetector = new RecoilChamberSensitiveDetector("RecoilChamber",6*6*6*112);
   SDMan->AddNewDetector(fSensitiveDetector);

   //fRCGeometry.Print();
}
//______________________________________________________________________________

RecoilChamberDetectorGeometry::~RecoilChamberDetectorGeometry()
{ }
//______________________________________________________________________________

G4VSolid * RecoilChamberDetectorGeometry::BuildWireSolid(int layer, int subcell)
{
   G4RotationMatrix aROT;
   aROT.rotateZ(5.0*CLHEP::degree);
   std::string name = "trap" + std::to_string(layer) + "_sub" + std::to_string(subcell);
   std::vector<G4TwoVector> trap_points = fRCGeometry.GetSubCellTrapPoints(layer,subcell);
   G4VSolid * trap  = new  G4GenericTrap(name.c_str(), fRCGeometry.WireLength/2.0, trap_points);
   G4VSolid * box   = new G4Box("container_box", 4.0*mm, 4.0*mm, fRCGeometry.WireLength/2.0);
   G4VSolid * trap0 = new G4IntersectionSolid("trap0",box,trap,fRCGeometry.GetFirstCellTransform(layer,subcell));
         //G4Transform3D(
         //   aROT ,
         //   G4ThreeVector(0,0,0)
         //)* G4Transform3D(
         //   G4RotationMatrix(fRCGeometry.GetFirstCellRotation(layer,subcell)),
         //   -1.0*fRCGeometry.GetFirstCellPosition(layer,subcell)
         //)
         //);
   return trap0;
}
//______________________________________________________________________________

void RecoilChamberDetectorGeometry::BuildUnitCells()
{
   for(int lay = 0; lay<fRCGeometry.NLayers; lay++){
      //  _______ 
      //  \__|__/
      //   \_|_/

      fWireVolume_solid[lay][0] = BuildWireSolid(lay,0);
      fWireVolume_solid[lay][1] = BuildWireSolid(lay,1);
      fWireVolume_solid[lay][2] = BuildWireSolid(lay,2);
      fWireVolume_solid[lay][3] = BuildWireSolid(lay,3);

      //--Step max
      G4double maxStep = 1.0*mm;
      G4UserLimits * fStepLimit = new G4UserLimits(maxStep);
      fStepLimit->SetMaxAllowedStep(maxStep);
      //fStepLimit->SetUserMaxTrackLength(maxStep);

      for(int i = 0; i<4; i++){

         std::string log_name = "wire_volume_log_" 
            + std::to_string(lay) + "_part" + std::to_string(i);

         fWireVolume_log[lay][i] = new G4LogicalVolume(
               fWireVolume_solid[lay][i], 
               HeiC4H10, 
               log_name);
         fWireVolume_log[lay][i]->SetUserLimits(fStepLimit);
      }
   }
}
//______________________________________________________________________________

void  RecoilChamberDetectorGeometry::PlaceCells(G4LogicalVolume * mother, int layer, double z_rotation, int wire_number ){

   bool check_overlaps = false;

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
      //fWireVolume_log[layer][i]->SetVisAttributes(G4VisAttributes::GetInvisible());

      new G4PVPlacement(
            //   0,G4ThreeVector(0,0,0),
            G4Transform3D(
               fRCGeometry.GetSenseWireRotation(wire_number),
               fRCGeometry.GetSenseWirePosition(wire_number)),
            //G4Transform3D(*(fWireVolume_rot[layer][i]),G4ThreeVector(0,0,0.)),
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
         "gasDetectorPara_phys", mother, false, 0, check_overlaps);            

   G4VisAttributes * GasDetectorVisAtt = new G4VisAttributes(G4Colour(0.3,0.1,0.1));
   GasDetectorVisAtt->SetForceWireframe(true);
   //GasDetectorVisAtt->SetDaughtersInvisible(true);
   logicGasDetector->SetVisAttributes(GasDetectorVisAtt);


   //--Making it a sensitive detector
   G4SDManager* SDman = G4SDManager::GetSDMpointer();
   G4String gasDetectorSDname = "/mydet/GasDetector";

   int wire_number = 0;

   for(G4int lay=0; lay<fRCGeometry.NLayers; lay++){

      bool placed_one = false;

      G4double PhiWire = fRCGeometry.fCellDeltaPhi.at(lay);

      for(int wi=0;wi<fRCGeometry.fNCells.at(lay); wi++){  

         if( (!placed_one) ) {

            PlaceCells( logicGasDetector, lay, double(wi)*PhiWire, wire_number );

            //placed_one=true;
            //std::string wire_name = "wire_volume_" + std::to_string(wire_number);
            //new G4PVPlacement(G4Transform3D(*rotationMatrix,G4ThreeVector(0,0,0.)),
            //      fWireVolume_log[lay], wire_name, logicGasDetector,
            //      false,wire_number,check_overlaps);
            wire_number++;
         }
      }   
   }

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
   G4Tubs* wiresolid = new G4Tubs("Wire", innerRadiusOfTheWire,outerRadiusOfTheWire, fRCGeometry.WireLength/2.0, startAngleOfTheWire, spanningAngleOfTheWire);
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


