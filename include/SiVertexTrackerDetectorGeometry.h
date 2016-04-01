#ifndef SiVertexTrackerDetectorGeometry_HH
#define SiVertexTrackerDetectorGeometry_HH

#include <array>

#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "SiTrackerSensitiveDetector.h"
#include "SensitiveRegionDetector.h"

#include "SiVertexTrackerGeometry.h"

/** Recoil Chamber.
 * 
 *  Note: Each "layer" consists of multiple radial wires. Typically
 *  there is only 3 wires radailly with the middle one being a sense wire.
 *  
 */
class SiVertexTrackerDetectorGeometry {

   private:

      std::array<G4VSolid*,3>         fLayerContainer_solid;
      std::array<G4LogicalVolume*,3>  fLayerContainer_log;
      std::array<G4Material*,3>       fLayerContainer_mat;

      std::array<G4VSolid*,3>         fLayerSi_solid;
      std::array<G4LogicalVolume*,3>  fLayerSi_log;
      std::array<G4Material*,3>       fLayerSi_mat;

      //std::array<std::array<G4VPhysicalVolume*,4>, 8> fWireVolume_phys;
      //std::array<std::array<G4RotationMatrix* ,4>, 8> fWireVolume_rot;

      eic::geo::SiVertexTrackerGeometry  fSiVertexTrackerGeometry;

      G4Material * fSi_mat  = nullptr;
      G4Material * fAir_mat = nullptr;

   public:


      SiTrackerSensitiveDetector * fSensitiveDetector;
      SensitiveRegionDetector    * fSensitiveDetector2;

   public:
      SiVertexTrackerDetectorGeometry();
      ~SiVertexTrackerDetectorGeometry();

      void BuildLogicalVolumes();
      //G4VSolid * BuildLayer(int layer, int radial_wire_number);
      //void  PlaceCells(G4LogicalVolume * mother, int layer, double z_rotation, int wire_number );

      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother);
      G4VPhysicalVolume * PlaceParallelPhysicalVolume(G4LogicalVolume * mother);

};

#endif

