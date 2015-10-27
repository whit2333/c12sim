#ifndef RecoilChamberDetectorGeometry_HH
#define RecoilChamberDetectorGeometry_HH

#include "G4SystemOfUnits.hh"
#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Material.hh"
#include "RecoilChamberSensitiveDetector.h"
#include "SensitiveRegionDetector.h"

#include "DCWire.h"

class RecoilChamberDetectorGeometry {

   private:

      std::array<G4VSolid*,         8> fWireVolume_solid;
      std::array<G4LogicalVolume*,  8> fWireVolume_log;
      std::array<G4VPhysicalVolume*,8> fWireVolume_phys;

      std::array<double,8>             fLayerRadius;
      std::array<double,8>             fLayerNwires;
      std::array<double,8>             fDeltaPhi;  

   public:

      //--Parameters of the gas detector
      G4double  innerRadiusOfTheGasDetector = 30.00*mm;
      G4double  outerRadiusOfTheGasDetector = 79.995*mm;
      G4double  hightOfTheGasDetector = 200.*mm;
      G4double  startAngleOfTheGasDetector = 0.*deg;
      G4double  spanningAngleOfTheGasDetector = 360.*deg;
      G4double   gasDetector_posX = 0.*mm;
      G4double   gasDetector_posY = 0.*mm;
      G4double   gasDetector_posZ = 0.*mm; 
      //--

      //--Parameters of the wires
      G4double  innerRadiusOfTheWire   = 0.00*mm;
      G4double  outerRadiusOfTheWire   = 0.04*mm;
      G4double  lengthOfTheWire        = 30.*cm;   // not the "hight"
      G4double  startAngleOfTheWire    = 0.*deg;
      G4double  spanningAngleOfTheWire = 360.*deg; 
      G4double  DeltaP                 = 2.0*mm; // wire separation around the circumference. 
                                                 // This number is only approximate if it does not divide the circumference
      G4double  DeltaR = 2.0*mm;
      G4double  NTLay  = 6.0;//4.; // Number of T layers
      G4double  NsLay  = 3.0;//5.; // Number of s layers 
      G4double  steAng = 1.0*deg;
      //--

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
      
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother);
      G4VPhysicalVolume * PlaceParallelPhysicalVolume(G4LogicalVolume * mother);

};

#endif

