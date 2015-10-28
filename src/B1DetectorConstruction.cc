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
#include "G4ExtrudedSolid.hh"

//___________________________________________________________________


B1DetectorConstruction::B1DetectorConstruction() : 
   G4VUserDetectorConstruction(), 
   world_x                      ( 10.0*m                  ),
   world_y                      ( 10.0*m                  ),
   world_z                      ( 16.0*m                  ),
   radiator_thickness           ( 6.0*mm                 ),
   collimator_target_center_gap ( 4.0*cm                 ),
   collimator_ID                ( 1.0*cm                 ),
   collimator_OD                ( 2.54*cm                ),
   outer_collimator_ID          ( 2.54*cm                ),
   outer_collimator_OD          ( 4.0*2.54*cm            ),
   collimator_diameter          ( 8.0*cm            ),
   collimator_z_end             ( 0.0*cm            ),
   radiator_collimator_gap      ( 1.0*mm                 ),
   collimator_length            ( 4.0*cm                 ),
   beampipe_length              ( 1.0*cm ),
   beampipe_diameter            ( 8.0*cm  ),
   radiator_diameter            ( 8.0*cm  ),
   scoring_diameter             ( 20.0*cm ),
   scoring_length               ( 0.01*mm    ),
   window_diameter              ( 1.0*cm  ),
   window_thickness             ( 8.0*mm  ),
   scoring2_diameter            ( 20.0*cm ),
   scoring2_length              ( 0.01*mm    ),
   fScoringVolume               ( 0),
   fHasBeenBuilt(false)
{
   fMessenger = new B1DetectorMessenger(this);
   fCollimatorMatName = "G4_Cu";
   scoring_det      = 0;
   scoring2_det     = 0;

   world_mat        = 0;
   world_solid      = 0;
   world_log        = 0;
   world_phys       = 0;

   beampipe_mat     = 0;
   beampipe_solid   = 0;
   beampipe_log     = 0;
   beampipe_phys    = 0;
   radiator_mat     = 0;
   radiator_solid   = 0;
   radiator_log     = 0;
   radiator_phys    = 0;
   collimator_mat   = 0;
   collimator_solid = 0;
   collimator_log   = 0;
   collimator_phys  = 0;
   collimator2_mat   = 0;
   collimator2_solid = 0;
   collimator2_log   = 0;
   collimator2_phys  = 0;
   outer_collimator_mat   = 0;
   outer_collimator_solid = 0;
   outer_collimator_log   = 0;
   outer_collimator_phys  = 0;
   scoring_mat      = 0;
   scoring_solid    = 0;
   scoring_log      = 0;
   scoring_phys     = 0;
   window_mat       = 0;
   window_solid     = 0;
   window_log       = 0;
   window_phys      = 0;
   scoring2_mat     = 0;
   scoring2_solid   = 0;
   scoring2_log     = 0;
   scoring2_phys    = 0;

   drift_chamber_log     = 0;
   drift_chamber_phys    = 0;
   htcc_log     = 0;
   htcc_phys    = 0;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();

   //fCherenkov_PMTArray = new PMTArraySensitiveDetector("PMTArray",6*8);
   //SDMan->AddNewDetector(fCherenkov_PMTArray);

   //nCherenkovPMT = 0;
   //nDriftChamberVols = 0;
   //fDriftChamber_SD =  new DriftChamberSensitiveDetector("DCIonPairSD",6*8);

   wire_hex_log = 0;
   fDriftChamber = 0;//new DriftChamberDetectorGeometry();
   fRecoilChamber = 0;
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
   beampipe_pos     = { 0, 0, -beampipe_length/2.0 - radiator_thickness/2.0 };
   radiator_pos     = { 0, 0, 0.0 };
   collimator_pos   = { 0, 0, collimator_length/4.0 + radiator_thickness/2.0 + radiator_collimator_gap };
   collimator2_pos   = { 0, 0, 3.0*collimator_length/4.0 + radiator_thickness/2.0 + radiator_collimator_gap };
   outer_collimator_pos   = { 0, 0, collimator_length/2.0 + radiator_thickness/2.0 + radiator_collimator_gap };
   collimator_z_end = collimator_length + radiator_thickness/2.0 + radiator_collimator_gap;
   scoring_pos      = { 0, 0, radiator_thickness/2.0 + radiator_collimator_gap/2.0 };
   window_pos       = { 0, 0, collimator_z_end - window_thickness/2.0 };
   scoring2_pos     = { 0, 0, collimator_z_end + collimator_target_center_gap };
}
//______________________________________________________________________________

void B1DetectorConstruction::DefineMaterials()
{
   ///////////////////////
   //---- Materials ----//
   ///////////////////////

   G4double a, z, density;
   G4int nel;

   H = new G4Element("Hydrogen", "H", z=1., a=  1.01*g/mole);
   C = new G4Element("Carbon",   "C", z=6., a= 12.01*g/mole);
   N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
   O = new G4Element("Oxygen",   "O", z=8., a= 16.00*g/mole);

   Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
   Air->AddElement(N, 70.*perCent);
   Air->AddElement(O, 30.*perCent);

   Lead = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

   Ar = new G4Material("ArgonGas",z=18., a= 39.95*g/mole, density=1.782*mg/cm3);

   Silicon = new G4Material("Silicon", z=14., a= 28.09*g/mole, density= 2.33*g/cm3);

   Scinti = new G4Material("Scintillator", density= 1.032*g/cm3, nel=2);
   Scinti->AddElement(C, 9);
   Scinti->AddElement(H, 10);

   a = 183.84*g/mole;
   elW = new G4Element("elW","W", z=74., a);
   density = 19.3*g/cm3;
   Tungsten = new G4Material("Tungsten",density,1);
   Tungsten->AddElement(elW,1);


   ////////////////////
   //---- Helium ----//
   ////////////////////

   //--For the target (5 atm)
   G4int ncomponents;
   G4double fractionmass;
   G4double temperature, pressure;

   a = 4.003*g/mole;
   G4double a_noUnit = 4.003;
   elHe  = new G4Element("Helium","He" , 2., a);
   pressure    = 3.0*atmosphere;
   G4double pre_noUnit = 3.0;
   temperature = 298.*kelvin;
   G4double tempe_noUnit = 298.;
   density     = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; //0.164*kg/m3 at 1 atm;
   He_Target = new G4Material("He", density, ncomponents=1,
         kStateGas,temperature,pressure);
   He_Target->AddElement(elHe, fractionmass=1.);
   //--

   //--For the clear space around the target (1 atm)
   pressure    = 1.0*atmosphere;
   pre_noUnit = 1.0;
   temperature = 298.*kelvin;
   tempe_noUnit = 298.;
   G4double d_He = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; //0.164*kg/m3 at 1 atm;
   He_ClearS = new G4Material("He_ClearS", d_He, ncomponents=1,
         kStateGas,temperature,pressure);
   He_ClearS->AddElement(elHe, fractionmass=1.);
   //--
   G4cout << "He density is " << d_He/(kg/m3) << "." << G4endl;

   ///////////////////
   //---- Xenon ----//
   ///////////////////
   a = 131.293*g/mole;
   a_noUnit = 131.293;
   elXe  = new G4Element("Xenon","Xe" , 54., a);
   pressure    = 0.3*atmosphere;
   pre_noUnit = 0.3;
   temperature = 298.*kelvin;
   tempe_noUnit = 298.;
   density     = (a_noUnit*pre_noUnit)/(0.0821*tempe_noUnit)*kg/m3; //5.37*kg/m3 at 1 atm;
   Xe_varPT = new G4Material("Xe", density, ncomponents=1,
         kStateGas,temperature,pressure);
   Xe_varPT->AddElement(elXe, fractionmass=1.);
   G4cout << "Xe density is " << density/(kg/m3) << "." << G4endl;

   ////////////////////
   //-- iso-Butane --//
   ////////////////////
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
   G4cout << "Isobutane density is " << density/(kg/m3) << "." << G4endl;


   ///////////////////
   //-- Deuterium --//
   ///////////////////
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

   /////////////////
   //---- Air ----//
   /////////////////

   density = 1.290*mg/cm3;
   Air = new G4Material("Air  ",density,2);
   Air->AddElement(N, 70*perCent);
   Air->AddElement(O, 30*perCent);


   ////////////////
   //-- Vacuum --//
   ////////////////
   //   Vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole,universe_mean_density,
   //                           kStateGas, 2.73*kelvin, 3.e-18*pascal);



   ////////////////
   //-- Kapton --//
   ////////////////
   density = 1.42*g/cm3;
   Kapton = new G4Material("Kapton",density, nel=4);
   Kapton->AddElement(H, fractionmass = 0.0273);
   Kapton->AddElement(C, fractionmass = 0.7213);
   Kapton->AddElement(N, fractionmass = 0.0765);
   Kapton->AddElement(O, fractionmass = 0.1749);



   ////////////////////////
   //-- Carbon dioxyde --//
   ////////////////////////


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

   G4cout << "C02 density is " << d_CO2/(kg/m3) << "." << G4endl;


   //////////////////////
   //-- Gas mixtures --//
   //////////////////////

   // 90% He - 10% CO2
   density = (d_He*90./100. + d_CO2*10./100.);
   He10CO2 = new G4Material("He10CO2"  , density, ncomponents=2);
   He10CO2->AddElement(elHe, fractionmass = 0.9);
   He10CO2->AddMaterial(CO2,    fractionmass = 0.1);

   G4cout << "He10C02 density is " << density/(kg/m3) << "." << G4endl;

   // 90% He - 10% iC4H10
   density = (d_He*90./100. + d_iC4H10*10./100.);
   HeiC4H10 = new G4Material("HeC4H10"  , density, ncomponents=2);
   HeiC4H10->AddElement(elHe,   fractionmass = 0.9);
   HeiC4H10->AddMaterial(isobutane, fractionmass = 0.1);

   G4cout << "HeiC4H10 density is " << density/(kg/m3) << "." << G4endl;
}
//______________________________________________________________________________

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
   //std::cout << "============================================================" << std::endl;
   //std::cout << "B1DetectorConstruction::Construct()" << std::endl;

   DefineMaterials();

   // Get nist material manager

   G4NistManager* nist = G4NistManager::Instance();

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
   if(!world_log)   world_log = new G4LogicalVolume( world_solid, world_mat, "world_log" );
   if(!world_phys)  world_phys  = new G4PVPlacement( 0, G4ThreeVector(), world_log, "world_phys", 0, false, 0, checkOverlaps );

   G4Colour            world_color {red, green, blue, alpha }; 
   G4VisAttributes   * world_vis   = new G4VisAttributes(world_color);
   //(*world_vis) = G4VisAttributes::GetInvisible();
   world_vis->SetForceWireframe(true);
   world_log->SetVisAttributes(world_vis);

   // ------------------------------------------------------------------------
   // Drift chamber
   // ------------------------------------------------------------------------
   //double d_layer_sl1 = 0.3861*cm;
   //double h_hex       = d_layer_sl1 * TMath::Sqrt(3.0);
   //double l_hex_side  = d_layer_sl1 * 2.0;

   //std::vector<G4ExtrudedSolid::ZSection> zsections;
   //zsections.push_back({ -1.0*cm, G4TwoVector{0.0,0.0},1.0 });
   //zsections.push_back({  1.0*cm, G4TwoVector{0.0,0.0},1.0 });

   //std::vector<G4TwoVector> polygon;
   //polygon.push_back({ -d_layer_sl1,  h_hex} );
   //polygon.push_back({  d_layer_sl1,  h_hex} );
   //polygon.push_back({        h_hex,    0.0} );
   //polygon.push_back({  d_layer_sl1, -h_hex} );
   //polygon.push_back({ -d_layer_sl1, -h_hex} );
   //polygon.push_back({       -h_hex,    0.0} );
   //G4Box * subtraction_box = new G4Box("a_hex_solid",0.01*cm,0.01*cm,0.01*cm);

   fDriftChamber = SimulationManager::GetInstance()->GetDriftDetectorGeometry();
   fDriftChamber->BuildLogicalVolumes();

   // Sectors
   for(int i = 1; i<=6; i++ ) {

      // Region I
      fDriftChamber->PlacePhysicalVolume( world_log, i, 1);

      // Region II
      fDriftChamber->PlacePhysicalVolume( world_log, i, 2);

      // Region III
      fDriftChamber->PlacePhysicalVolume( world_log, i, 3);
   }

   // ------------------------------------------------------------------------
   // Recoil Chamber
   // ------------------------------------------------------------------------

   fRecoilChamber = SimulationManager::GetInstance()->GetRecoilDetectorGeometry();
   fRecoilChamber->He10CO2   = He10CO2;
   fRecoilChamber->HeiC4H10  = HeiC4H10;
   fRecoilChamber->Tungsten  = Tungsten; 
   fRecoilChamber->Mylar     = Mylar;
   fRecoilChamber->PlacePhysicalVolume( world_log);

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

   G4NistManager* man = nist;//G4NistManager::Instance();

   // Way to access a material
   Mylar = man->FindOrBuildMaterial("G4_MYLAR");


   ///////////////////////////////////////////////////////////////////////////
   //------------------------- Magnetic field ------------------------------//
   ///////////////////////////////////////////////////////////////////////////

   // create a field
   //G4UniformMagField* magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue)); // create a field
   C12MagneticField * magField = new C12MagneticField(true,true);

   // set it as the default field
   G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager(); // set it as the default field
   fieldMgr->SetDetectorField(magField);
   fieldMgr->CreateChordFinder(magField); // create the objects which calculate the trajectory
   // change the accuracy of volume intersection
   //fieldMgr->GetChordFinder()->SetDeltaChord(acc_inter);
   double pos[4] = {0.0,0.0,0.0,0.0};
   G4ThreeVector myfield = magField->GetFieldValue(pos);
   G4cout << "magnetic field at origin : ( " << myfield[0] << " , " << myfield[1] << " , " << G4BestUnit(myfield[2],"Magnetic flux density") << ")" << G4endl;;

   // Relative accuracy values:
   G4double minEps= 1.0e-6;  //   Minimum & value for smallest steps
   G4double maxEps= 1.0e-5;  //   Maximum & value for largest steps

   fieldMgr->SetMinimumEpsilonStep( minEps );
   fieldMgr->SetMaximumEpsilonStep( maxEps );
   fieldMgr->SetDeltaOneStep( 0.1e-4 * mm );  // 0.5 micrometer

   G4cout << "EpsilonStep: set min= " << minEps << " max= " << maxEps << G4endl;




   ///////////////////////////////////////////////////////////////////////////
   //------------------------------ Target ---------------------------------//
   ///////////////////////////////////////////////////////////////////////////

   //--Creating the geometry
   G4ThreeVector target_pos {target_posX, target_posY, target_posZ};

   G4Tubs* target = new G4Tubs("target_solid", innerRadiusOfTheTarget,outerRadiusOfTheTarget, hightOfTheTarget, startAngleOfTheTarget,spanningAngleOfTheTarget);

   G4Material * target_mat = Deuterium;
   G4LogicalVolume* logicTarget = new G4LogicalVolume(target, target_mat,"target_log");

   new G4PVPlacement(0, target_pos,   logicTarget,  "Target_phys", world_log,    
         false,           //no boolean operation
         0,               //copy number
         checkOverlaps);  //overlaps checking

   // Definition of visualisation attributes
   // Instantiation of a set of visualization attributes with cyan colour
   G4VisAttributes * TargetVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.));
   // Set the forced wireframe style 
   TargetVisAtt->SetForceWireframe(true);
   // Assignment of the visualization attributes to the logical volume
   logicTarget->SetVisAttributes(TargetVisAtt);
   //--





   ///////////////////////////////////////////////////////////////////////////
   //---------------------------- Kapton layer -----------------------------//
   ///////////////////////////////////////////////////////////////////////////

   //--Creating the geometry
   G4ThreeVector kapton_pos = G4ThreeVector(target_posX, target_posY, target_posZ);

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
   //--





   ///////////////////////////////////////////////////////////////////////////
   //---------------------------- Clear Space ------------------------------//
   ///////////////////////////////////////////////////////////////////////////

   //--Creating the geometry
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
   //--



   ///////////////////////////////////////////////////////////////////////////
   //--------------------- Mylar layer around clear space-------------------//
   ///////////////////////////////////////////////////////////////////////////

   //--Creating the geometry
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
   //--




   ///////////////////////////////////////////////////////////////////////////
   //------------------------ Outisde Mylar layer --------------------------//
   ///////////////////////////////////////////////////////////////////////////

   //--Creating the geometry
   G4ThreeVector okapton_pos = G4ThreeVector(target_posX, target_posY, target_posZ);

   G4Tubs* okapton_cyl = new G4Tubs("OKaptonCylinder", innerRadiusOfTheOKapton, outerRadiusOfTheOKapton, hightOfTheOKapton, startAngleOfTheOKapton, spanningAngleOfTheOKapton);

   G4LogicalVolume* logicOKapton =                         
      new G4LogicalVolume(okapton_cyl,          //its solid
            Mylar,                //its material
            "OKaptonCylinder");   //its name

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
   //--



   ///////////////////////////////////////////////////////////////////////////
   //--------------------------- Sci Detector ------------------------------//
   ///////////////////////////////////////////////////////////////////////////

   //--Creating the geometry
   G4ThreeVector SiDetector_pos = G4ThreeVector(SiDetector_posX, SiDetector_posY, SiDetector_posZ);

   G4Tubs* SiDetector = new G4Tubs("SiDetector", innerRadiusOfTheSiDetector,outerRadiusOfTheSiDetector, hightOfTheSiDetector, startAngleOfTheSiDetector, spanningAngleOfTheSiDetector);

   G4LogicalVolume* logicSiDetector =                         
      new G4LogicalVolume(SiDetector,          //its solid
            Scinti,              //its material
            "SiDetector");       //its name

   new G4PVPlacement(0,           //no rotation
         SiDetector_pos,          //at position
         logicSiDetector,         //its logical volume
         "SiDetector",            //its name
         world_log,               //its mother  volume
         false,                   //no boolean operation
         0,                       //copy number
         checkOverlaps);          //overlaps checking

   // Definition of visualisation attributes
   // Instantiation of a set of visualization attributes with cyan colour
   G4VisAttributes * SiDetectorVisAtt = new G4VisAttributes(G4Colour(0.5,0.1,0.2));
   // Set the forced wireframe style 
   SiDetectorVisAtt->SetForceWireframe(true);
   // Assignment of the visualization attributes to the logical volume
   logicSiDetector->SetVisAttributes(SiDetectorVisAtt);
   //--

   //--Making it a sensitive detector
   //G4String SiDetectorSDname = "/mydet/SiDetector";
   //SiDetectorSD * siDetectorSD = new SiDetectorSD(SiDetectorSDname);
   //SDman->AddNewDetector(siDetectorSD);
   //logicSiDetector->SetSensitiveDetector(siDetectorSD);
   //--



   //fScoringVolume = scoring2_log;
   fHasBeenBuilt = true;

   return world_phys;
}
//___________________________________________________________________

void B1DetectorConstruction::SetCollimatorMaterial(G4String materialName)
{
   fCollimatorMatName = materialName;
   if(fHasBeenBuilt) Rebuild();
}
//______________________________________________________________________________

void     B1DetectorConstruction::SetRadiatorCollimatorGap(G4double l)
{   
   radiator_collimator_gap = l; 
   if(fHasBeenBuilt) Rebuild();
}
//______________________________________________________________________________

void     B1DetectorConstruction::SetInnerCollimatorOD(G4double l)
{   
   collimator_OD       = l;
   outer_collimator_ID = l;
   if(fHasBeenBuilt) Rebuild();
}
//______________________________________________________________________________

void     B1DetectorConstruction::SetCollimatorLength(G4double l)
{   
   collimator_length = l;
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

