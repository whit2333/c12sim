#include "B1DetectorConstruction.hh"
#include "B1DetectorMessenger.hh"

#include "SimulationManager.h"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "C12MagneticField.h"
#include "G4TransportationManager.hh"
#include "G4CachedMagneticField.hh"

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
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "FakeSD.hh"
#include "DriftChamberDetectorGeometry.h"
#include "HTCCDetectorGeometry.h"
#include "G4ExtrudedSolid.hh"
#include "G4GDMLParser.hh"
#include "G4AutoDelete.hh"
#include "G4ChordFinder.hh"

#ifdef USE_CADMESH
#include "CADMesh.hh"
#endif

B1DetectorConstruction::B1DetectorConstruction() : 
   G4VUserDetectorConstruction(), 
   world_x( 15.0*m ),
   world_y( 15.0*m ),
   world_z( 20.0*m ),
   fHasBeenBuilt(false)
{
   fMessenger = new B1DetectorMessenger(this);
}
//___________________________________________________________________

B1DetectorConstruction::~B1DetectorConstruction()
{
   delete fMessenger;
   delete Scinti;
   delete Silicon;
   delete Ar;
   delete Lead;
   delete Air;
   delete Kapton;
   delete He_Target;
   delete He_ClearS;
   delete Xe_varPT;

   delete O;
   delete N;
   delete C;
   delete H;
   delete elHe;
   delete elXe;
}
//______________________________________________________________________________

void B1DetectorConstruction::Rebuild()
{
   B1DetectorConstruction::Construct();
   G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
//______________________________________________________________________________

void B1DetectorConstruction::CalculatePositions()
{
   //beampipe_pos     = { 0, 0, -beampipe_length/2.0 - radiator_thickness/2.0 };
   //radiator_pos     = { 0, 0, 0.0 };
   //collimator_pos   = { 0, 0, collimator_length/4.0 + radiator_thickness/2.0 + radiator_collimator_gap };
   //collimator2_pos   = { 0, 0, 3.0*collimator_length/4.0 + radiator_thickness/2.0 + radiator_collimator_gap };
   //outer_collimator_pos   = { 0, 0, collimator_length/2.0 + radiator_thickness/2.0 + radiator_collimator_gap };
   //collimator_z_end = collimator_length + radiator_thickness/2.0 + radiator_collimator_gap;
   //scoring_pos      = { 0, 0, radiator_thickness/2.0 + radiator_collimator_gap/2.0 };
   //window_pos       = { 0, 0, collimator_z_end - window_thickness/2.0 };
   //scoring2_pos     = { 0, 0, collimator_z_end + collimator_target_center_gap };
}
//______________________________________________________________________________

void B1DetectorConstruction::DefineMaterials()
{
   G4NistManager* nist = G4NistManager::Instance();

   G4double a       = 0.0;
   G4double z       = 0.0;
   G4double density = 0.0;
   G4int    nel     = 0;

   // -----------------------------------------------
   // Elements 
   H = new G4Element("Hydrogen", "H", z=1., a=  1.01*g/mole);
   C = new G4Element("Carbon",   "C", z=6., a= 12.01*g/mole);
   N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
   O = new G4Element("Oxygen",   "O", z=8., a= 16.00*g/mole);
   elW = new G4Element("elW","W", z=74., a= 183.84*g/mole);
   elHe  = new G4Element("Helium","He" , 2., a = 4.003*g/mole);

   // -----------------------------------------------
   // Elemental Materials
   Lead     = new G4Material("Lead", z    = 82., a = 207.19*g/mole, density = 11.35*g/cm3);
   Ar       = new G4Material("ArgonGas",z = 18., a = 39.95*g/mole, density  = 1.782*mg/cm3);
   Silicon  = new G4Material("Silicon", z = 14., a = 28.09*g/mole, density  = 2.33*g/cm3);
   Tungsten = new G4Material("Tungsten",density = 19.3*g/cm3,1);
   Tungsten->AddElement(elW,1);
   elXe  = new G4Element("Xenon","Xe" , 54., a = 131.293*g/mole);

   // -----------------------------------------------
   // Air
   Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
   Air->AddElement(N, 70.*perCent);
   Air->AddElement(O, 30.*perCent);

   // -----------------------------------------------
   // Mylar
   Mylar = nist->FindOrBuildMaterial("G4_MYLAR");

   // -----------------------------------------------
   // Scintillator plastic
   Scinti = new G4Material("Scintillator", density= 1.032*g/cm3, nel=2);
   Scinti->AddElement(C, 9);
   Scinti->AddElement(H, 10);

   ////////////////////
   //----  ----//
   ////////////////////

   // -----------------------------------------------
   // Helium for the target (5 atm)
   G4int    ncomponents;
   G4double fractionmass;
   G4double temperature, pressure;

   G4double a_noUnit     = 4.003;
   G4double pre_noUnit   = 3.0;
   G4double tempe_noUnit = 298.;
   pressure              = 3.0*atmosphere;
   temperature           = 298.*kelvin;
   density     = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; //0.164*kg/m3 at 1 atm;
   He_Target = new G4Material("He", density, ncomponents=1, kStateGas,temperature,pressure);
   He_Target->AddElement(elHe, fractionmass=1.);

   // -----------------------------------------------
   // Helium for the clear space around the target (1 atm)
   pressure    = 1.0*atmosphere;
   pre_noUnit = 1.0;
   temperature = 298.*kelvin;
   tempe_noUnit = 298.;
   G4double d_He = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; //0.164*kg/m3 at 1 atm;
   He_ClearS = new G4Material("He_ClearS", d_He, ncomponents=1, kStateGas,temperature,pressure);
   He_ClearS->AddElement(elHe, fractionmass=1.);
   //G4cout << "He density is " << d_He/(kg/m3) << "." << G4endl;

   // -----------------------------------------------
   // Xenon
   a_noUnit     = 131.293;
   pressure     = 0.3*atmosphere;
   pre_noUnit   = 0.3;
   temperature  = 298.*kelvin;
   tempe_noUnit = 298.;
   density      = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; //5.37*kg/m3 at 1 atm;
   Xe_varPT = new G4Material("Xe", density, ncomponents=1,
         kStateGas,temperature,pressure);
   Xe_varPT->AddElement(elXe, fractionmass=1.);
   //G4cout << "Xe density is " << density/(kg/m3) << "." << G4endl;

   // -----------------------------------------------
   // iso-Butane
   a = 58.122*g/mole;
   a_noUnit = 58.122;
   pressure    = 1.0*atmosphere;
   pre_noUnit = 1.0;        
   temperature = 298.*kelvin;
   tempe_noUnit = 298.;
   G4double d_iC4H10  = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; 
   isobutane = new G4Material("isoC4H10",d_iC4H10,ncomponents=2,kStateGas,temperature,pressure);
   isobutane->AddElement(C,4);
   isobutane->AddElement(H,10);
   //G4cout << "Isobutane density is " << density/(kg/m3) << "." << G4endl;


   // -----------------------------------------------
   // Deuterium
   a = 2.0141*g/mole;
   a_noUnit = 2.0141;
   pressure = 7.0*atmosphere; 
   pre_noUnit = 7.0; 
   temperature = 298.*kelvin;
   tempe_noUnit = 298.;
   G4Isotope * iso_H2 = new G4Isotope("H2",1,2,a);
   ele_D = new G4Element("Deuterium Atom","D",1);
   ele_D->AddIsotope(iso_H2, 100.*perCent);
   density = (2.*a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3;
   Deuterium = new G4Material("DeuteriumGas", density, ncomponents=1, kStateGas, temperature, pressure); 
   Deuterium->AddElement(ele_D,100.*perCent);

   // -----------------------------------------------
   // LH2
   density   = 0.07085*g/cm3; //g/cm^3
   LH2 = new G4Material("LH2", density, ncomponents=1); 
   LH2->AddElement(H,100.*perCent);

   ////////////////
   //-- Vacuum --//
   ////////////////
   //   Vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole,universe_mean_density,
   //                           kStateGas, 2.73*kelvin, 3.e-18*pascal);


   // -----------------------------------------------
   // Kapton
   Kapton = new G4Material("Kapton",density = 1.42*g/cm3, nel=4);
   Kapton->AddElement(H, fractionmass = 0.0273);
   Kapton->AddElement(C, fractionmass = 0.7213);
   Kapton->AddElement(N, fractionmass = 0.0765);
   Kapton->AddElement(O, fractionmass = 0.1749);

   // -----------------------------------------------
   // CO2
   a_noUnit = 44.01;
   pressure = 1.*atmosphere;
   pre_noUnit = 1.;
   temperature = 298.*kelvin;
   tempe_noUnit = 298.;
   G4double d_CO2 = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3;
   CO2 = new G4Material("CO2", d_CO2, nel=2,
         kStateGas,temperature,pressure);
   CO2->AddElement(C,1);
   CO2->AddElement(O,2);
   //G4cout << "C02 density is " << d_CO2/(kg/m3) << "." << G4endl;

   // -----------------------------------------------
   // Gas mixtures

   // 90% He - 10% CO2
   density = (d_He*90./100. + d_CO2*10./100.);
   He10CO2 = new G4Material("He10CO2"  , density, ncomponents=2);
   He10CO2->AddElement(elHe, fractionmass = 0.9);
   He10CO2->AddMaterial(CO2,    fractionmass = 0.1);
   //G4cout << "He10C02 density is " << density/(kg/m3) << "." << G4endl;

   // 90% He - 10% iC4H10
   density = (d_He*90./100. + d_iC4H10*10./100.);
   HeiC4H10 = new G4Material("HeC4H10"  , density, ncomponents=2);
   HeiC4H10->AddElement(elHe,   fractionmass = 0.9);
   HeiC4H10->AddMaterial(isobutane, fractionmass = 0.1);
   //G4cout << "HeiC4H10 density is " << density/(kg/m3) << "." << G4endl;

}
//______________________________________________________________________________

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
   //std::cout << "B1DetectorConstruction::Construct()" << std::endl;

   DefineMaterials();

   // Get nist material manager

   G4NistManager* nist = G4NistManager::Instance();
   SimulationManager * simMan =  SimulationManager::GetInstance();

   CalculatePositions();

   bool    checkOverlaps    = false;
   int     natoms           = 0;
   int     ncomponents      = 0;
   double  A                = 0.0;
   double  Z                = 0.0;
   double  thickness        = 0.0;
   double  surface_z        = 0.0;
   double  red              = 0.0;
   double  green            = 0.0;
   double  blue             = 0.0;
   double  alpha            = 0.0;
   double  density          = 0.0;
   double  pressure         = 0.0;
   double  temperature      = 0.0;


   // ------------------------------------------------------------------------
   // World
   // ------------------------------------------------------------------------
   density     = universe_mean_density;
   pressure    = 1.e-7*bar;
   temperature = 0.1*kelvin;
   red         = 0.0/256.0;
   green       = 200.0/256.0;
   blue        = 0.0/256.0;
   alpha       = 0.4;

   world_mat   = nist->FindOrBuildMaterial("G4_AIR");
   if(!world_solid) world_solid = new G4Box( "World", 0.5*world_x, 0.5 * world_y, 0.5 * world_z );
   if(!world_log)   world_log   = new G4LogicalVolume( world_solid, world_mat, "world_log" );
   if(!world_phys)  world_phys  = new G4PVPlacement( 0, G4ThreeVector(), world_log, "world_phys", 0, false, 0, checkOverlaps );

   G4Colour            world_color {red, green, blue, alpha }; 
   G4VisAttributes   * world_vis   = new G4VisAttributes(world_color);
   //(*world_vis) = G4VisAttributes::GetInvisible();
   world_vis->SetForceWireframe(true);
   world_log->SetVisAttributes(world_vis);

   //G4UserLimits* userLimits = new G4UserLimits(1.0*cm);
   //world_log->SetUserLimits(userLimits);


   // ------------------------------------------------------------------------
   // Drift chamber
   // ------------------------------------------------------------------------
   //std::cout << " Drift chamber construction \n";
   //fDriftChamber = SimulationManager::GetInstance()->GetDriftDetectorGeometry();
   //fDriftChamber->BuildLogicalVolumes();
   // Sectors
   //for(int i = 1; i<=6; i++ ) {
      // Region I
      //fDriftChamber->PlaceParallelPhysicalVolume( world_log, i, 1);
      //if(i==1) fDriftChamber->PlacePhysicalVolume( world_log, i, 1);
      // Region II
      //fDriftChamber->PlaceParallelPhysicalVolume( world_log, i, 2);
      //fDriftChamber->PlacePhysicalVolume( world_log, i, 2);
      // Region III
      //fDriftChamber->PlaceParallelPhysicalVolume( world_log, i, 3);
      //fDriftChamber->PlacePhysicalVolume( world_log, i, 3);
   //}

   // ------------------------------------------------------------------------
   // HTCC  
   // ------------------------------------------------------------------------
   //std::cout << " HTCC construction \n";
   //fHTCC = SimulationManager::GetInstance()->GetHTCCDetectorGeometry();
   //fHTCC->BuildLogicalVolumes();
   //fHTCC->SetGasIndexOfRefraction(false);
   //fHTCC->PlacePhysicalVolume( world_log, 1, 1);
   //fHTCC->PlacePhysicalVolume( world_log, 2, 2);
   //fHTCC->PlacePhysicalVolume( world_log, 3, 1);
   //fHTCC->PlacePhysicalVolume( world_log, 4, 2);
   //fHTCC->PlacePhysicalVolume( world_log, 5, 1);
   //fHTCC->PlacePhysicalVolume( world_log, 6, 2);

   // ------------------------------------------------------------------------
   // Electromagnetic Calorimeter (EC)
   // ------------------------------------------------------------------------
   //std::cout << " EC construction \n";
   //fEC = SimulationManager::GetInstance()->GetECDetectorGeometry();
   //fEC->BuildLogicalVolumes();
   // Sectors
   //for(int i = 1; i<=6; i++ ) {
   //   fEC->PlaceParallelPhysicalVolume( world_log, i, 1);
   //}

   // ------------------------------------------------------------------------
   // EIC SiVertexTracker
   // ------------------------------------------------------------------------
   //std::cout << " Si Vertex Tracker construction \n";
   //fSiVertexTracker = SimulationManager::GetInstance()->GetSiVertexTrackerDetectorGeometry();
   //fSiVertexTracker->PlacePhysicalVolume( world_log );
   //std::cout << " done \n";

   // ------------------------------------------------------------------------
   // Recoil Chamber
   // ------------------------------------------------------------------------
   std::cout << " Recoil chamber construction \n";
   fRecoilChamber = SimulationManager::GetInstance()->GetRecoilDetectorGeometry();
   fRecoilChamber->He10CO2   = He10CO2;
   fRecoilChamber->HeiC4H10  = HeiC4H10;
   fRecoilChamber->Tungsten  = Tungsten; 
   fRecoilChamber->Mylar     = Mylar;
   fRecoilChamber->PlaceParallelPhysicalVolume( world_log);


   // ------------------------------------------------------------------------
   // Recoil Hodoscope
   // ------------------------------------------------------------------------
   //std::cout << " Recoil Hodo construction ...\n";
   //std::cout << " RecoilHodoscopeDetectorGeometry \n";
   //std::cout << " This version consists of two layers of scintillator boars (one thin, one thick)\n";
   //std::cout << " The photons are detected at both ends\n";
   //fRecoilHodo = SimulationManager::GetInstance()->GetRecoilHodoDetectorGeometry();
   //fRecoilHodo->BuildLogicalVolumes();
   //fRecoilHodo->PlacePhysicalVolume( world_log, world_phys);

   //std::cout << " RecoilHodoscopeDetectorGeometry2 \n";
   //std::cout << " This version consists of one thin scintillator layer and a tiled thick outer layer\n";
   //fRecoilHodo2 = SimulationManager::GetInstance()->GetRecoilHodoDetectorGeometry2();
   //fRecoilHodo2->BuildLogicalVolumes();
   //fRecoilHodo2->PlacePhysicalVolume( world_log, world_phys);

   std::cout << " Recoil Hodo construction ...\n";
   std::cout << " RecoilHodoscopeDetectorGeometry3 \n";
   std::cout << " This version consists of two layers of scintillator bars (one thin, one thick)\n";
   std::cout << " However the outer layer is segmented and read from the outside with one SiPM\n";
   fRecoilHodo3 = SimulationManager::GetInstance()->GetRecoilHodoDetectorGeometry3();
   fRecoilHodo3->BuildLogicalVolumes();
   fRecoilHodo3->PlacePhysicalVolume( world_log, world_phys);

#ifdef USE_CADMESH
   auto alert_frame_mesh = CADMesh::TessellatedMesh::FromPLY(ClasDigi_SHARE_DIR "/alert/geometry_files/alert_frame.ply");
   alert_frame_mesh->SetScale( mm );
   alert_frame_mesh->SetOffset( chamber_port_flange_offset );
   auto alert_frame_solid = alert_frame_mesh->GetSolid();
#endif

   // ------------------------------------------------------------------------
   // Solenoid Geometry
   // ------------------------------------------------------------------------
   fSolenoid = SimulationManager::GetInstance()->GetSolenoidDetectorGeometry();
   fSolenoid->BuildLogicalVolumes();
   fSolenoid->PlacePhysicalVolume( world_log );

   // ------------------------------------------------------------------------
   // Torus Geometry
   // ------------------------------------------------------------------------
   //fTorus = SimulationManager::GetInstance()->GetTorusDetectorGeometry();
   //fTorus->BuildLogicalVolumes();
   //fTorus->PlacePhysicalVolume( world_log );

   // ------------------------------------------------------------------------
   // Beamline
   // ------------------------------------------------------------------------
   //std::cout << " beamline construction \n";
   //fBeamline = SimulationManager::GetInstance()->GetBeamlineDetectorGeometry();
   //fBeamline->BuildLogicalVolumes();
   //fBeamline->PlacePhysicalVolume( world_log );

   // ------------------------------------------------------------------------
   // MVT
   // ------------------------------------------------------------------------
   //fMVT = SimulationManager::GetInstance()->GetMVTDetectorGeometry();
   //fMVT->BuildLogicalVolumes();
   //fMVT->PlacePhysicalVolume( world_log );

   // ------------------------------------------------------------------------
   // MVT
   // ------------------------------------------------------------------------
   //fSVT = SimulationManager::GetInstance()->GetSVTDetectorGeometry();
   //fSVT->BuildLogicalVolumes();
   //fSVT->PlacePhysicalVolume( world_log );

   // ------------------------------------------------------------------------
   // Magnetic field
   // ------------------------------------------------------------------------
   //G4UniformMagField* magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue)); // create a field
   bool use_torus    = true;
   bool use_solenoid = true;
   if( simMan->GetSolenoidFieldScale() == 0.0 ){
      use_solenoid = false;
   }
   if( simMan->GetToroidFieldScale()   == 0.0 ){
      use_torus    = false;
   }
   if( use_torus || use_solenoid ) {
      std::cout << "Loading Field Maps... (this may take a minute) \n";
   }
   //std::cout << " simMan->GetSolenoidFieldScale() = " << simMan->GetSolenoidFieldScale() << std::endl;
   //std::cout << " simMan->GetToroidFieldScale()   = " << simMan->GetToroidFieldScale()   << std::endl;
   C12MagneticField * magField = new C12MagneticField(
         use_torus,                     use_solenoid ,
         simMan->GetToroidFieldScale(), simMan->GetSolenoidFieldScale() );

   G4double                 distanceConst = 2.5 * millimeter; 
   G4MagneticField * pCachedMagField= new G4CachedMagneticField(  magField,  distanceConst); 

   G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
   fieldMgr->SetDetectorField(pCachedMagField);
   fieldMgr->CreateChordFinder(pCachedMagField);
   //G4bool forceToAllDaughters = true;
   //world_log->SetFieldManager(fieldMgr, forceToAllDaughters);

   // change the accuracy of volume intersection

   // Relative accuracy values:
   G4double minEps = 1.0e-6;  //   Minimum & value for smallest steps
   G4double maxEps = 1.0e-5;  //   Maximum & value for largest steps

   fieldMgr->GetChordFinder()->SetDeltaChord(0.001);//1.0e-5);
   fieldMgr->SetMinimumEpsilonStep( minEps );
   fieldMgr->SetMaximumEpsilonStep( maxEps );
   fieldMgr->SetDeltaOneStep( 0.001);//1.0e-6 * mm );  // 0.5 micrometer
   //G4cout << "EpsilonStep: set min= " << minEps << " max= " << maxEps << G4endl;

   // Register the field and its manager for deleting
   //G4AutoDelete::Register(magField);
   //G4AutoDelete::Register(fieldMgr);

   double pos[4] = {0.0,0.0,0.0,0.0};
   G4ThreeVector myfield = magField->GetFieldValue(pos);
   //G4cout << "magnetic field at origin : ( " 
   //   << myfield[0] << " , " 
   //   << myfield[1] << " , " 
   //   << G4BestUnit(myfield[2],"Magnetic flux density") 
   //   << ")" << G4endl;;

   // ------------------------------------------------------------------------
   // Target
   // ------------------------------------------------------------------------
   G4ThreeVector target_pos {target_posX, target_posY, target_posZ};
   G4Tubs* target = new G4Tubs("target_solid", innerRadiusOfTheTarget,outerRadiusOfTheTarget, fTargetLength/2.0, 0.0, 360.0*CLHEP::degree);
   G4Material * target_mat = He_Target; //LH2;//Deuterium;
   G4LogicalVolume* logicTarget = new G4LogicalVolume(target, target_mat,"target_log");
   new G4PVPlacement(0, target_pos,   logicTarget,  "Target_phys", world_log,    
         false,           //no boolean operation
         0,               //copy number
         checkOverlaps);  //overlaps checking
   //Definition of visualisation attributes
   //Instantiation of a set of visualization attributes with cyan colour
   G4VisAttributes * TargetVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.,0.3));
   // Set the forced wireframe style 
   //TargetVisAtt->SetForceWireframe(true);
   // Assignment of the visualization attributes to the logical volume
   logicTarget->SetVisAttributes(TargetVisAtt);


   // ------------------------------------------------------------------------
   // beam vacuum  
   // ------------------------------------------------------------------------
   //density     = universe_mean_density;
   //pressure    = 1.e-7*bar;
   //temperature = 0.1*kelvin;
   //red       = 0.0/256.0;
   //green     = 0.0/256.0;
   //blue      = 192.0/256.0;
   //alpha     = 0.4;

   //if(!beampipe_mat)   beampipe_mat   = new G4Material("beampipe_mat", /*z=*/1.0, /*a=*/1.01*g/mole, density, kStateGas,temperature,pressure);
   //if(!beampipe_solid) beampipe_solid  = new G4Tubs("beampipe_solid", 0.0, beampipe_diameter/2.0, beampipe_length/2.0, 0.0, 360.*deg );
   //if(!beampipe_log  ) beampipe_log   = new G4LogicalVolume(beampipe_solid, beampipe_mat,"beampipe_log");
   //if(!beampipe_phys ) beampipe_phys  = new G4PVPlacement(0,beampipe_pos, beampipe_log, "beampipe_phys",world_log,false,0,checkOverlaps);                                  
   //G4Colour            beampipe_color {red, green, blue, alpha };   // Gray 
   //G4VisAttributes   * beampipe_vis   = new G4VisAttributes(beampipe_color);
   //beampipe_log->SetVisAttributes(beampipe_vis);

   //G4UserLimits * scoring_limits = new G4UserLimits(0.004*um);
   //scoring_log->SetUserLimits(scoring_limits);
   G4ThreeVector kapton_pos = G4ThreeVector(target_posX, target_posY, target_posZ);
   //if(!beampipe_mat)   beampipe_mat   = new G4Material("beampipe_mat", /*z=*/1.0, /*a=*/1.01*g/mole, density, kStateGas,temperature,pressure);
   //if(!beampipe_solid) beampipe_solid  = new G4Tubs("beampipe_solid", 0.0, beampipe_diameter/2.0, beampipe_length/2.0, 0.0, 360.*deg );
   //if(!beampipe_log  ) beampipe_log   = new G4LogicalVolume(beampipe_solid, beampipe_mat,"beampipe_log");
   //if(!beampipe_phys ) beampipe_phys  = new G4PVPlacement(0,beampipe_pos, beampipe_log, "beampipe_phys",world_log,false,0,checkOverlaps);                                  
   //G4Colour            beampipe_color {red, green, blue, alpha };   // Gray 
   //G4VisAttributes   * beampipe_vis   = new G4VisAttributes(beampipe_color);
   //beampipe_log->SetVisAttributes(beampipe_vis);

   //G4UserLimits * scoring_limits = new G4UserLimits(0.004*um);
   //scoring_log->SetUserLimits(scoring_limits);

   //G4NistManager* man = nist;//G4NistManager::Instance();

   // ------------------------------------------------------------------------
   // Kapton foil around the target
   // ------------------------------------------------------------------------
   //G4ThreeVector kapton_pos = G4ThreeVector(target_posX, target_posY, target_posZ);
   G4Tubs* kapton_cyl = new G4Tubs("KaptonCylinder", innerRadiusOfTheKapton, outerRadiusOfTheKapton, hightOfTheKapton, startAngleOfTheKapton, spanningAngleOfTheKapton);

   G4LogicalVolume* logicKapton =                         
      new G4LogicalVolume(kapton_cyl,          //its solid
            Kapton,              //its material
            "KaptonCylinder");   //its name

   new G4PVPlacement(0,                     //no rotation
         kapton_pos,              //at position
         logicKapton,             //its logical volume
         "KaptonCylinder",        //its name
         world_log,    //its mother  volume
         false,                   //no boolean operation
         0,                       //copy number
         checkOverlaps);          //overlaps checking

   // Definition of visualisation attributes
   // Instantiation of a set of visualization attributes with cyan colour
   G4VisAttributes * KaptonVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
   // Set the forced wireframe style 
   KaptonVisAtt->SetForceWireframe(true);
   // Assignment of the visualization attributes to the logical volume
   logicKapton->SetVisAttributes(KaptonVisAtt);


   // ------------------------------------------------------------------------
   // Clear Space - Parameters of space around the target
   // ------------------------------------------------------------------------
   G4ThreeVector around_pos = G4ThreeVector(target_posX, target_posY, target_posZ);
   G4Tubs* around = new G4Tubs("ClearSpace", innerRadiusOfAround, outerRadiusOfAround, hightOfAround, startAngleOfAround, spanningAngleOfAround);
   G4LogicalVolume* logicAround =                         
      new G4LogicalVolume(around,              //its solid
            He_ClearS,           //its material
            "ClearSpace");       //its name
   new G4PVPlacement(0,                     //no rotation
         around_pos,              //at position
         logicAround,             //its logical volume
         "ClearSpace",            //its name
         world_log,    //its mother  volume
         false,                   //no boolean operation
         0,                       //copy number
         checkOverlaps);          //overlaps checking

   // Definition of visualisation attributes
   // Instantiation of a set of visualization attributes with cyan colour
   G4VisAttributes * AroundVisAtt = new G4VisAttributes(G4Colour(1.,0.,1.));
   // Set the forced wireframe style 
   AroundVisAtt->SetForceWireframe(true);
   // Assignment of the visualization attributes to the logical volume
   logicAround->SetVisAttributes(AroundVisAtt);


   // ------------------------------------------------------------------------
   // Mylar layer around clear space
   // ------------------------------------------------------------------------
   G4ThreeVector Oclkapton_pos = G4ThreeVector(target_posX, target_posY, target_posZ);

   G4Tubs* Oclkapton_cyl = new G4Tubs("OclKaptonCylinder",
         innerRadiusOfTheOclKapton, outerRadiusOfTheOclKapton,
         hightOfTheOclKapton, startAngleOfTheOclKapton,
         spanningAngleOfTheOclKapton);

   G4LogicalVolume* logicOclKapton =                         
      new G4LogicalVolume(Oclkapton_cyl,          //its solid
            Mylar,                  //its material
            "OclKaptonCylinder");   //its name

   new G4PVPlacement(0,                     //no rotation
         Oclkapton_pos,              //at position
         logicOclKapton,             //its logical volume
         "OclKaptonCylinder",        //its name
         world_log,    //its mother  volume
         false,                   //no boolean operation
         0,                       //copy number
         checkOverlaps);          //overlaps checking

   // Definition of visualisation attributes
   // Instantiation of a set of visualization attributes with cyan colour
   G4VisAttributes * OclKaptonVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
   // Set the forced wireframe style 
   KaptonVisAtt->SetForceWireframe(true);
   // Assignment of the visualization attributes to the logical volume
   logicOclKapton->SetVisAttributes(OclKaptonVisAtt);
   
   // ------------------------------------------------------------------------
   // Outisde Mylar layer
   // ------------------------------------------------------------------------
   G4ThreeVector okapton_pos = G4ThreeVector(target_posX, target_posY, target_posZ);
   G4Tubs* okapton_cyl = new G4Tubs("OKaptonCylinder", 
         innerRadiusOfTheOKapton, outerRadiusOfTheOKapton, 
         hightOfTheOKapton, startAngleOfTheOKapton, spanningAngleOfTheOKapton);
   G4LogicalVolume* logicOKapton = new G4LogicalVolume(okapton_cyl, Mylar, "OKaptonCylinder");
   new G4PVPlacement(0,                     //no rotation
         okapton_pos,              //at position
         logicOKapton,             //its logical volume
         "OKaptonCylinder",        //its name
         world_log,    //its mother  volume
         false,                   //no boolean operation
         0,                       //copy number
         checkOverlaps);          //overlaps checking

   // Definition of visualisation attributes
   // Instantiation of a set of visualization attributes with cyan colour
   G4VisAttributes * OKaptonVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
   // Set the forced wireframe style 
   OKaptonVisAtt->SetForceWireframe(true);
   // Assignment of the visualization attributes to the logical volume
   logicOKapton->SetVisAttributes(OKaptonVisAtt);

   fHasBeenBuilt = true;

   int run_number = SimulationManager::GetInstance()->GetRunNumber();
   G4GDMLParser parser;
   G4String gdml_name = "data/log/c12geometry_" + std::to_string(run_number) + G4String(".gdml");
   //parser.Write(gdml_name,world_phys);

   return world_phys;
}
//___________________________________________________________________

void     B1DetectorConstruction::ToggleCherenkov(bool l)
{   
   fHTCC = SimulationManager::GetInstance()->GetHTCCDetectorGeometry();
   fHTCC->SetGasIndexOfRefraction(l);
   if(fHasBeenBuilt) Rebuild();
}
//______________________________________________________________________________

void     B1DetectorConstruction::PrintConfigInfo() const
{   

   if(fHasBeenBuilt) {
      std::cout << "------------------------------------------------------------------------\n";
      //std::cout 
      //   << " Geometry parameters :  \n"
      //<< "                   weight : " << collimator_log->GetMass(true)/kg << " kg\n"
      //<< "                 material : " << fCollimatorMatName               << "\n"
      //<< "                   length : " << collimator_length/cm             << " cm\n"
      //<< "           collimator ID  : " << collimator_ID/cm                 << " cm\n"
      //<< "           collimator OD  : " << collimator_OD/cm                 << " cm\n"
      //<< " radiator collimator gap  : " << radiator_collimator_gap/cm       << " cm\n"
      //<< "  collimator target dist  : " << collimator_target_center_gap/cm  << " cm\n"
      //<< "      radiator thickness  : " << radiator_thickness/cm            << " cm\n";
   } else {
      std::cout << " detector not built yet" << std::endl;
   }

}
//______________________________________________________________________________

