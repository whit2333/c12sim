#include "TorusDetectorGeometry.h"

#include "G4Polycone.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"

TorusDetectorGeometry::TorusDetectorGeometry()
{
   using namespace CLHEP;
//  Torus |  root |   Torus |    0*cm 0.0*cm -898.093*mm |    0*deg 0*deg 0*deg | 0000ff  |  Polycone |0.0*deg 360.0*deg 14*counts 536.25*mm 390*mm 390*mm 390*mm 390*mm 390*mm 390*mm 390*mm 390*mm 390*mm 390*mm 390*mm 390*mm 784*mm 1020*mm 1020*mm 1038*mm 1038*mm 1020*mm 1020*mm 1038*mm 1038*mm 1020*mm 1020*mm 1038*mm 1038*mm 1020*mm 1020*mm 0*mm 425*mm 440.7*mm 546.4*mm 562.1*mm 814.92*mm 830.62*mm 936.32*mm 952.02*mm 1204.84*mm 1220.54*mm 1326.24*mm 1341.94*mm 1783.44*mm |  G4_Fe |   no | 1  | 1  | 1  |  1  |  1  |   no |   no |     no 
   int zplanes = 14; // number of planes in z directions
   double rInner[] = {
      536.25*mm, 390*mm, 390*mm, 390*mm, 390*mm,
      390*mm, 390*mm, 390*mm, 390*mm, 390*mm,
      390*mm, 390*mm, 390*mm, 784*mm
   };
   double rOuter[] = { 
      1020*mm, 1020*mm, 1038*mm, 1038*mm, 1020*mm,
      1020*mm, 1038*mm, 1038*mm, 1020*mm, 1020*mm,
      1038*mm, 1038*mm, 1020*mm, 1020*mm
   };
   double zPlane[] = { 
      0*mm, 425*mm, 440.7*mm, 546.4*mm, 562.1*mm,
      814.92*mm, 830.62*mm, 936.32*mm, 952.02*mm, 1204.84*mm,
      1220.54*mm, 1326.24*mm, 1341.94*mm, 1783.44*mm
   };

   fTorus_solid = new G4Polycone("fTorus_solid",            ///< name
         0,  ///< Initial Phi starting angle
         360*deg,  ///< Total Phi angle
         zplanes,        ///< Number of z planes
         zPlane,         ///< z coordinate of corners
         rInner,         ///< Tangent distance to inner surface
         rOuter);        ///< Tangent distance to outer surface

   G4NistManager * matman = G4NistManager::Instance();
   fTorus_mat = matman->FindOrBuildMaterial("G4_Fe");

   fTorus_pos = {0*cm,0.0*cm,-898.093*mm};
}
//______________________________________________________________________________

TorusDetectorGeometry::~TorusDetectorGeometry()
{ }
//______________________________________________________________________________

void TorusDetectorGeometry::BuildLogicalVolumes()
{
   fTorus_log = new G4LogicalVolume(fTorus_solid, fTorus_mat, "fTorus_log");
}
//______________________________________________________________________________

G4VPhysicalVolume * TorusDetectorGeometry::PlacePhysicalVolume(G4LogicalVolume * mother)
{
   G4VPhysicalVolume * phys = new G4PVPlacement(
         0, fTorus_pos,
         fTorus_log,          // its logical volume
         "fTorus_phys", // its name
         mother,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         false);                        // check for overlaps
   return phys;
}
//______________________________________________________________________________

