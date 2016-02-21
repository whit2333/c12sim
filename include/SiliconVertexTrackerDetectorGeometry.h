#ifndef SiliconVertexTrackerDetectorGeometry_HH
#define SiliconVertexTrackerDetectorGeometry_HH

#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "CLHEP/Units/PhysicalConstants.h"

//#include "DriftChamberSensitiveDetector.h"
//#include "SensitiveRegionDetector.h"
#include "DCWire.h"
using namespace CLHEP;

class SiliconVertexTrackerDetectorGeometry {

   private:

      double fRohacell71Density   = 75.0*kg/m3;
      double fRohacell71Thickness = 2.5*mm;
      double fSiDensity           = 2.57*g/cm3;
      double fSiLayerThickness    = 0.320*mm;

      double fBSTModuleWidth      = 42.0*mm;
      double fBSTReadoutLength    = 41.0*mm;

      std::array<double          ,4> fBST_length =  {{ 223.25*mm, 223.25*mm, 334.875*mm, 334.875*mm }};
      std::array<double          ,4> fBST_radius =  {{  65.0*mm,   93.0*mm,   120.0*mm, 161.250*mm }};
      std::array<G4VSolid*       ,4> fBST_solid;
      std::array<G4LogicalVolume*,4> fBST_log;
      std::array<int             ,4> fNModules   = {{ 10, 14, 18, 24 }};
      std::array<double          ,4> fBSTDownstreamEnd_displacment =  {{ 73.13*mm, 111.9*mm, 170.1*mm, 170.1*mm }};

      std::array<G4LogicalVolume*,4> fBSTSilicon_log;
      std::array<G4VSolid*       ,4> fBSTSilicon_solid;

      std::array<G4LogicalVolume*,4> fBSTRohacell71_log;
      std::array<G4VSolid*       ,4> fBSTRohacell71_solid;

      double fFSTSiliconOuterWidth  = 79.615*mm;
      double fFSTSiliconMidWidth    = 26.041*mm;
      double fFSTSiliconInnerWidth  =  7.786*mm;

      double fFSTInnerRadus         = 17.0*mm;
      double fFSTSiliconRadus       = 181.0*mm;

      std::array<double,          3> fFST_length = {{ 164.0*mm, 164.0*mm, 164.0*mm }};
      std::array<double,          3> fFST_displacement = {{ 190.1*mm, 210.1*mm, 230.1*mm }};
      std::array<G4VSolid*       ,3> fFST_solid;
      std::array<G4LogicalVolume*,3> fFST_log;

      std::array<G4LogicalVolume*,4> fFSTSilicon_log;
      std::array<G4VSolid*       ,4> fFSTSilicon_solid;

      std::array<G4LogicalVolume*,4> fFSTRohacell71_log;
      std::array<G4VSolid*       ,4> fFSTRohacell71_solid;

      G4Material * fSiliconMaterial = nullptr;
      G4Material * fRohacell71Material = nullptr;

   public:

      //DriftChamberSensitiveDetector * fSensitiveDetector;
      //SensitiveRegionDetector       * fSensitiveDetector2;


   public:
      SiliconVertexTrackerDetectorGeometry();
      ~SiliconVertexTrackerDetectorGeometry();

      void BuildLogicalVolumes();
      
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother);

};

#endif

