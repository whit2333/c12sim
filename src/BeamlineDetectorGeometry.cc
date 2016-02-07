#include "BeamlineDetectorGeometry.h"

#include "G4Polycone.hh"
#include "G4Tubs.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"

BeamlineDetectorGeometry::BeamlineDetectorGeometry()
{
   using namespace CLHEP;
   G4NistManager * matman = G4NistManager::Instance();

   // -------------------------------
   // noft 

// moller_shield_tube | root |flanged tube from moller cone to torus | 0 0 0 | 0 0 0 | ff8883 | Polycone |0.0*deg 360*deg 6*counts 0.0*mm 0.0*mm 0.0*mm 0.0*mm 0.0*mm 0.0*mm 90*mm 90*mm 80*mm 80*mm 90*mm 90*mm 2267.217*mm 2282.217*mm 2282.217*mm 2647.817*mm 2647.817*mm 2662.817*mm | G4_STAINLESS-STEEL | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
   int zplanes = 6; // number of planes in z directions
   double rInner[] = { 0.0*mm, 0.0*mm, 0.0*mm, 0.0*mm, 0.0*mm, 0.0*mm };
   double rOuter[] = { 90*mm, 90*mm, 80*mm, 80*mm, 90*mm, 90*mm };
   double zPlane[] = { 2267.217*mm, 2282.217*mm, 2282.217*mm, 2647.817*mm, 2647.817*mm, 2662.817*mm };

   fMollerShieldTube_mat = matman->FindOrBuildMaterial("G4_STAINLESS-STEEL");
   fMollerShieldTube_pos = {0*cm,0.0*cm,0.0*mm};
   fMollerShieldTube_solid = new G4Polycone("fMollerShieldTube_solid",            ///< name
         0,  ///< Initial Phi starting angle
         360*degree,  ///< Total Phi angle
         zplanes,        ///< Number of z planes
         zPlane,         ///< z coordinate of corners
         rInner,         ///< Tangent distance to inner surface
         rOuter);        ///< Tangent distance to outer surface

//moller_shield_tube_w | moller_shield_tube |tungsten tube from moller cone to torus | 0.0*cm 0.0*cm 2465.017*mm | 0 0 0 | ffff9b | Tube | 0.0*mm 60*mm 197.8*mm 0.*deg 360.*deg | beamline_W | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
   fMollerShieldTubeW_mat = matman->FindOrBuildMaterial("G4_W");
   fMollerShieldTubeW_pos = {0.0*cm, 0.0*cm, 2465.017*mm};
   fMollerShieldTubeW_solid = new G4Tubs("fMollerShieldTubeW_solid",0.0*mm, 60*mm, 197.8*mm, 0.*degree, 360.*degree);

//moller_shield_tube_a |moller_shield_tube_w |aluminum tube from moller cone to torus | 0.0*cm 0.0*cm 0.0*mm | 0 0 0 | ffffff | Tube | 0.0*mm 41*mm 197.8*mm 0.*deg 360.*deg | G4_Al | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
   fMollerShieldTubeA_mat = matman->FindOrBuildMaterial("G4_Al");
   fMollerShieldTubeA_pos = {0.0*cm, 0.0*cm, 0*mm};
   fMollerShieldTubeA_solid = new G4Tubs("fMollerShieldTubeA_solid",0.0*mm, 41.0*mm, 197.8*mm , 0.*degree, 360.*degree);

//moller_shield_tube_v |moller_shield_tube_a |vacumm tube from moller cone to torus | 0.0*cm 0.0*cm 0.0*mm | 0 0 0 | aaffff | Tube | 0.0*mm 30*mm 197.8*mm 0.*deg 360.*deg | G4_Galactic | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
   fMollerShieldTubeV_mat = matman->FindOrBuildMaterial("G4_Galactic");
   fMollerShieldTubeV_pos = {0.0*cm, 0.0*cm, 0*mm};
   fMollerShieldTubeV_solid = new G4Tubs("fMollerShieldTubeV_solid",0.0*mm, 30.0*mm, 197.8*mm , 0.*degree, 360.*degree);

// moller_shield_cone | root | tungsten moller cone | 0 0 0 | 0 0 0 | ffff9b | Polycone |0.0*deg 360*deg 3*counts 0.0*mm 0.0*mm 0.0*mm 30*mm 90*mm 90*mm 420*mm 1390*mm 2267.117*mm | beamline_W | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
   int zplanes2 = 3; // number of planes in z directions
   double rInner2[] = { 0.0*mm, 0.0*mm, 0.0*mm };
   double rOuter2[] = { 30*mm, 90*mm, 90*mm };
   double zPlane2[] = { 420*mm, 1390*mm, 2267.117*mm };

   fMollerShieldCone_mat = matman->FindOrBuildMaterial("G4_W");
   fMollerShieldCone_pos = {0*cm,0.0*cm,0.0*mm};
   fMollerShieldCone_solid = new G4Polycone("fMollerShieldCone_solid",            ///< name
         0,  ///< Initial Phi starting angle
         360*degree,  ///< Total Phi angle
         zplanes2,        ///< Number of z planes
         zPlane2,         ///< z coordinate of corners
         rInner2,         ///< Tangent distance to inner surface
         rOuter2);        ///< Tangent distance to outer surface

//moller_shield_cone_a | moller_shield_cone | aluminum moller cone | 0 0 0 | 0 0 0 | ffffff | Polycone |0.0*deg 360*deg 3*counts 0.0*mm 0.0*mm 0.0*mm 24*mm 41*mm 41*mm 420*mm 1390*mm 2267.117*mm | G4_Al | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
   int   zplanes3 = 3; // number of planes in z directions
   double rInner3[] = { 0.0*mm, 0.0*mm, 0.0*mm };
   double rOuter3[] = { 24*mm, 41*mm, 41*mm};
   double zPlane3[] = { 420*mm, 1390*mm, 2267.117*mm };

   fMollerShieldConeA_mat = matman->FindOrBuildMaterial("G4_Al");
   fMollerShieldConeA_pos = {0*cm,0.0*cm,0.0*mm};
   fMollerShieldConeA_solid = new G4Polycone("fMollerShieldConeA_solid",            ///< name
         0,  ///< Initial Phi starting angle
         360*degree,  ///< Total Phi angle
         zplanes3,        ///< Number of z planes
         zPlane3,         ///< z coordinate of corners
         rInner3,         ///< Tangent distance to inner surface
         rOuter3);        ///< Tangent distance to outer surface

//moller_shield_cone_v |moller_shield_cone_a | vacuum moller cone | 0 0 0 | 0 0 0 | aaffff | Polycone |0.0*deg 360*deg 3*counts 0.0*mm 0.0*mm 0.0*mm 23*mm 30*mm 30*mm 420*mm 1390*mm 2267.117*mm | G4_Galactic | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
   int   zplanes4 = 3; // number of planes in z directions
   double rInner4[] = { 0.0*mm, 0.0*mm, 0.0*mm };
   double rOuter4[] = { 23*mm, 30*mm, 40*mm};
   double zPlane4[] = { 420*mm, 1390*mm, 2267.117*mm };

   fMollerShieldConeV_mat = matman->FindOrBuildMaterial("G4_Galactic");
   fMollerShieldConeV_pos = {0*cm,0.0*cm,0.0*mm};
   fMollerShieldConeV_solid = new G4Polycone("fMollerShieldConeV_solid",            ///< name
         0,  ///< Initial Phi starting angle
         360*degree,  ///< Total Phi angle
         zplanes4,        ///< Number of z planes
         zPlane4,         ///< z coordinate of corners
         rInner4,         ///< Tangent distance to inner surface
         rOuter4);        ///< Tangent distance to outer surface
   // -------------------------------
   // ft
   //pipe_after_torus_ring | root | Pipe after Torus Ring | 0*mm 0.0*mm 6135.917*mm | 0 0 0 | 993333 | Tube | 135*mm 175*mm 1000*mm 0.0*deg 360*deg | beamline_W | no | 1 | 1 | 1 | 1 | 1 | no | no | no 

   //pipe_after_torus_ring_shield | root | Shielding after Torus Ring | 0*mm 0.0*mm 5051.917*mm | 0 0 0 | 0000ff | Polycone |0.0*deg 360*deg 4*counts 62*mm 62*mm 62*mm 62*mm 170*mm 170*mm 110*mm 110*mm 0*mm 25*mm 25.1*mm 300*mm | beamline_W | no | 1 | 1 | 1 | 1 | 1 | no | no | no 

   //beamline_pipe_shielding | root | Pipe after Torus Ring | 0*mm 0.0*mm 8663.317*mm | 0 0 0 | 999966 | Tube | 0*mm 60*mm 6000*mm 0.0*deg 360*deg | beamline_W | no | 1 | 1 | 1 | 1 | 1 | no | no | no 

   // beamline_pipe |beamline_pipe_shielding |Aluminum Pipe after Torus Ring | 0*mm 0.0*mm 0*mm | 0 0 0 | 87AFC7 | Tube | 0*mm 40*mm 6000*mm 0.0*deg 360*deg | G4_Al | no | 1 | 1 | 1 | 1 | 1 | no | no | no 
   //fBeamline_solid = new G4Tubs(0*mm 40*mm 6000*mm 0.0*deg 360*deg);

   // beamline_vacuum | beamline_pipe | Vacuum Pipe after Torus Ring | 0*mm 0.0*mm 0*mm | 0 0 0 | aaffff | Tube | 0*mm 37*mm 6000*mm 0.0*deg 360*deg | Vacuum | no | 1 | 1 | 1 | 1 | 1 | no | no | no 

}
//______________________________________________________________________________

BeamlineDetectorGeometry::~BeamlineDetectorGeometry()
{ }
//______________________________________________________________________________

void BeamlineDetectorGeometry::BuildLogicalVolumes()
{

   fMollerShieldTube_log = new G4LogicalVolume(fMollerShieldTube_solid, fMollerShieldTube_mat, "fMollerShieldTube_log");
   fMollerShieldTube_vis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.5));
   //fMollerShieldTube_vis->SetForceWireframe(true);
   fMollerShieldTube_log->SetVisAttributes(fMollerShieldTube_vis);

   fMollerShieldTubeW_log = new G4LogicalVolume(fMollerShieldTubeW_solid, fMollerShieldTubeW_mat, "fMollerShieldTubeW_log");
   fMollerShieldTubeW_vis = new G4VisAttributes(G4Colour(0.0, 0.9, 0.1, 0.3));
   fMollerShieldTubeW_log->SetVisAttributes(fMollerShieldTubeW_vis);

   fMollerShieldTubeA_log = new G4LogicalVolume(fMollerShieldTubeA_solid, fMollerShieldTubeA_mat, "fMollerShieldTubeA_log");
   fMollerShieldTubeA_vis = new G4VisAttributes(G4Colour(0.9, 0.0, 0.1, 0.3));
   fMollerShieldTubeA_log->SetVisAttributes(fMollerShieldTubeA_vis);

   fMollerShieldTubeV_log = new G4LogicalVolume(fMollerShieldTubeV_solid, fMollerShieldTubeV_mat, "fMollerShieldTubeV_log");
   fMollerShieldTubeV_vis = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0, 1.0));
   fMollerShieldTubeV_log->SetVisAttributes(fMollerShieldTubeV_vis);

   fMollerShieldCone_log  = new G4LogicalVolume(fMollerShieldCone_solid,  fMollerShieldCone_mat,  "fMollerShieldCone_log");
   fMollerShieldCone_vis = new G4VisAttributes(G4Colour(0.0, 0.9, 0.1, 0.3));
   fMollerShieldCone_log->SetVisAttributes(fMollerShieldCone_vis);

   fMollerShieldConeA_log = new G4LogicalVolume(fMollerShieldConeA_solid, fMollerShieldConeA_mat, "fMollerShieldConeA_log");
   fMollerShieldConeA_vis = new G4VisAttributes(G4Colour(0.9, 0.0, 0.1, 0.3));
   fMollerShieldConeA_log->SetVisAttributes(fMollerShieldConeA_vis);

   fMollerShieldConeV_log = new G4LogicalVolume(fMollerShieldConeV_solid, fMollerShieldConeV_mat, "fMollerShieldConeV_log");
   fMollerShieldConeV_vis = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0, 1.0));
   fMollerShieldConeV_log->SetVisAttributes(fMollerShieldConeV_vis);

}
//______________________________________________________________________________

G4VPhysicalVolume * BeamlineDetectorGeometry::PlacePhysicalVolume(G4LogicalVolume * mother)
{

   //--------------------------------------
   fMollerShieldConeV_phys = new G4PVPlacement(
         0, fMollerShieldConeV_pos,
         fMollerShieldConeV_log,
         "fMollerShieldConeV_phys",
         fMollerShieldConeA_log,
         false,
         0,
         false);

   fMollerShieldConeA_phys = new G4PVPlacement(
         0, fMollerShieldConeA_pos,
         fMollerShieldConeA_log,
         "fMollerShieldConeA_phys",
         fMollerShieldCone_log,
         false,
         0,
         false);

   fMollerShieldCone_phys = new G4PVPlacement(
         0, fMollerShieldCone_pos,
         fMollerShieldCone_log,
         "fMollerShieldCone_phys",
         mother,
         false,
         0,
         false);

   //--------------------------------------

   fMollerShieldTubeV_phys = new G4PVPlacement(
         0, fMollerShieldTubeV_pos,
         fMollerShieldTubeV_log,
         "fMollerShieldTubeV_phys",
         fMollerShieldTubeA_log,
         false,
         0,
         false);

   fMollerShieldTubeA_phys = new G4PVPlacement(
         0, fMollerShieldTubeA_pos,
         fMollerShieldTubeA_log,
         "fMollerShieldTubeA_phys",
         fMollerShieldTubeW_log,
         false,
         0,
         false);

   fMollerShieldTubeW_phys = new G4PVPlacement(
         0, fMollerShieldTubeW_pos,
         fMollerShieldTubeW_log,
         "fMollerShieldTubeW_phys",
         fMollerShieldTube_log,
         false,
         0,
         false);

   fMollerShieldTube_phys = new G4PVPlacement(
         0, fMollerShieldTube_pos,
         fMollerShieldTube_log,          // its logical volume
         "fMollerShieldTube_phys", // its name
         mother,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         false);                        // check for overlaps

   return fMollerShieldTube_phys;
}
//______________________________________________________________________________

