#ifndef RecoilChamberDetectorGeometry_HH
#define RecoilChamberDetectorGeometry_HH

#include <array>

#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "RecoilChamberSensitiveDetector.h"
#include "SensitiveRegionDetector.h"

#include "DCWire.h"
#include "RCGeometry.h"

/** Recoil Chamber.
 * 
 *  Note: Each "layer" consists of multiple radial wires. Typically
 *  there is only 3 wires radailly with the middle one being a sense wire.
 *  
 */
class RecoilChamberDetectorGeometry {

   private:

      std::array<std::array<G4VSolid*,         4 >, 8> fWireVolume_solid;
      std::array<std::array<G4LogicalVolume*,  4 >, 8> fWireVolume_log;
      std::array<std::array<G4VPhysicalVolume*,4 >, 8> fWireVolume_phys;
      std::array<std::array<G4RotationMatrix* ,4 >, 8> fWireVolume_rot;

      clas12::geo::RCGeometry          fRCGeometry;

   public:

      //--Parameters of the gas detector
      G4double  innerRadiusOfTheGasDetector   = 30.00*mm;
      G4double  outerRadiusOfTheGasDetector   = 79.995*mm;
      G4double  hightOfTheGasDetector         = 200.*mm;
      G4double  startAngleOfTheGasDetector    = 0.*deg;
      G4double  spanningAngleOfTheGasDetector = 360.*deg;
      G4double   gasDetector_posX             = 0.*mm;
      G4double   gasDetector_posY             = 0.*mm;
      G4double   gasDetector_posZ             = 0.*mm;

      //--Parameters of the wires
      G4double  innerRadiusOfTheWire   = 0.00*mm;
      G4double  outerRadiusOfTheWire   = 0.04*mm;
      G4double  lengthOfTheWire        = 30.*cm;   // not the "hight"
      G4double  startAngleOfTheWire    = 0.*deg;
      G4double  spanningAngleOfTheWire = 360.*deg; 
      G4double  DeltaP                 = 2.0*mm; // *desired* wire separation around the circumference. 
                                                 // This number is only approximate if it does not divide the circumference
      //G4double  DeltaR   = 2.0*mm;
      //G4double  fNLayers = 7.0;//NTLay  = 1.0;//4.; // Number of T layers
      //G4double  NsLay  = 3.0;//5.; // Number of s layers  // wire planes (cylinders) per layer
      //G4double  steAng = 10.0*deg;

      G4Material* He10CO2;
      G4Material* HeiC4H10;
      G4Material* Tungsten; 
      G4Material* Mylar;

      RecoilChamberSensitiveDetector * fSensitiveDetector;
      SensitiveRegionDetector       * fSensitiveDetector2;

   public:
      RecoilChamberDetectorGeometry();
      ~RecoilChamberDetectorGeometry();

      void BuildLogicalVolumes();
      void BuildUnitCells();
      G4VSolid * BuildWireSolid(int layer, int radial_wire_number);
      void  PlaceCells(G4LogicalVolume * mother, int layer, double z_rotation, int wire_number );

      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother);
      G4VPhysicalVolume * PlaceParallelPhysicalVolume(G4LogicalVolume * mother);

};

#endif

