#include "ScintTileDetectorGeometry.hh"
#include <cmath>
#include <algorithm>

#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4GenericTrap.hh"
#include "G4Box.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4UserLimits.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTrap.hh"
#include "G4Torus.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"

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
#include "G4ExtrudedSolid.hh"

#include "BeamTestSD.hh"
#include "SiPMSD.hh"


ScintTileDetectorGeometry::ScintTileDetectorGeometry()
{
   using namespace CLHEP;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();

   double hc_over_eV  = 1.239842e-6*m*eV;

   fScintPhotonEnergies = {
      2.48989*eV, 2.53049*eV, 2.58855*eV, 2.63802*eV, 2.70399*eV, 
      2.73661*eV, 2.7762*eV, 2.82651*eV, 2.86553*eV, 2.90224*eV, 
      2.93292*eV, 2.95713*eV, 2.97453*eV, 2.98844*eV, 2.99188*eV, 
      3.00986*eV, 3.02077*eV, 3.0318*eV, 3.04672*eV, 3.06565*eV, 
      3.06947*eV, 3.07719*eV, 3.08495*eV, 3.10056*eV, 3.12048*eV, 
      3.12048*eV, 3.12867*eV, 3.15315*eV, 3.17815*eV, 3.20767*eV, 
      3.25013*eV, 3.26726*eV, 3.35135*eV
   };

   fScintLightOutput = {
      0.0290686, 0.0480032, 0.0665378, 0.0986732, 0.163811, 0.229749, 
      0.308887, 0.407827, 0.460364, 0.512968, 0.572371, 0.65191, 0.744916, 
      0.864658, 0.911261, 0.957597, 0.984066, 0.9972, 0.990266, 0.963264, 
      0.95653, 0.929729, 0.902927, 0.85599, 0.748983, 0.748983, 0.668845, 
      0.515101, 0.328022, 0.147543, 0.0468698, 0.0466031, 0.0
   };

   for(int i = 0; i < fScintWavelengths.size(); i++ ){
      fScintWavelengths[i]      = hc_over_eV/fScintPhotonEnergies[i];
      fScintLightOutput_slow[i] = 0.0;
      fScintRefIndex[i]         = 1.58;
   }

   fScint_pos  = {0.0, 0.0, 0.0};

   //fSensitiveDetector = new FTScintSensitiveDetector("FTScint",6*6*6*112);
   //SDMan->AddNewDetector(fSensitiveDetector);

}
//______________________________________________________________________________

ScintTileDetectorGeometry::~ScintTileDetectorGeometry()
{
   // Scintillator material

}
//______________________________________________________________________________

void ScintTileDetectorGeometry::BuildMaterials()
{
   using std::vector;
   int     natoms           = 0;
   int     ncomponents      = 0;
   double  A                = 0.0;
   double  Z                = 0.0;
   double  density          = 0.0;
   double  pressure         = 0.0;
   double  temperature      = 0.0;
   double  a                = 0.0;

   // ----------------------------------------------------------
   // Relative light output for different particles
   // refscint is the reference scintillation yield for the relative measuremnts.
   // It appears to be the scintillation yield from a 160 MeV proton.

   double refscint = 10000.0/MeV;

   vector< vector< double > > relativeLightOutput = { 
      { // Particle kinetic energy
         0*MeV, 10*MeV, 20*MeV, 30*MeV, 40*MeV, 50*MeV, 50*MeV, 70*MeV, 90*MeV, 110*MeV, 130*MeV, 150*MeV, 170*MeV, 190*MeV, 190*MeV, 
         240*MeV, 290*MeV, 340*MeV, 390*MeV, 440*MeV, 490*MeV, 540*MeV, 590*MeV, 640*MeV, 690*MeV, 740*MeV, 790*MeV, 840*MeV, 
         890*MeV, 940*MeV, 990*MeV},
      { // 1 electron
         0.0010197*refscint, 0.070038*refscint, 0.139817*refscint,
         0.210325*refscint , 0.281509*refscint, 0.353295*refscint,
         0.353295*refscint , 0.498268*refscint, 0.644225*refscint,
         0.789808*refscint , 0.932235*refscint, 1.07667*refscint ,
         1.22111*refscint  , 1.36555*refscint , 1.36555*refscint ,
         1.72665*refscint  , 2.08774*refscint , 2.44884*refscint ,
         2.80994*refscint  , 3.17103*refscint , 3.53213*refscint ,
         3.89323*refscint  , 4.25432*refscint , 4.61542*refscint ,
         4.97652*refscint  , 5.33761*refscint , 5.69871*refscint ,
         6.0598*refscint   , 6.4209*refscint  , 6.782*refscint   ,
         7.14309*refscint },
      { // 2 proton
         0*refscint, 0.0458257*refscint, 0.0964081*refscint, 
         0.149766*refscint, 0.20579*refscint, 0.264325*refscint, 0.264325*refscint, 0.388084*refscint, 0.518921*refscint, 0.654003*refscint,
         0.786399*refscint, 0.919298*refscint, 1.0522*refscint, 1.1851*refscint, 1.1851*refscint, 1.51734*refscint, 1.84959*refscint, 
         2.18184*refscint, 2.51409*refscint, 2.84633*refscint, 3.17858*refscint, 3.51083*refscint, 3.84308*refscint, 4.17533*refscint, 
         4.50757*refscint, 4.83982*refscint, 5.17207*refscint, 5.50432*refscint, 5.83656*refscint, 6.16881*refscint, 
         6.50106*refscint},
      { // 3 deuteron
         0.000228012*refscint, 0.0357616*refscint, 0.0758948*refscint, 0.120514*refscint, 0.169428*refscint, 
         0.222373*refscint, 0.222373*refscint, 0.338909*refscint, 0.466477*refscint, 0.600215*refscint, 0.725919*refscint, 
         0.853917*refscint, 0.981914*refscint, 1.10991*refscint, 1.10991*refscint, 1.42991*refscint, 1.7499*refscint, 2.06989*refscint, 
         2.38989*refscint, 2.70988*refscint, 3.02988*refscint, 3.34987*refscint, 3.66987*refscint, 3.98986*refscint, 4.30986*refscint, 
         4.62985*refscint, 4.94984*refscint, 5.26984*refscint, 5.58983*refscint, 5.90983*refscint, 6.22982*refscint},
      { // 4 triton
         0.00327668*refscint, 0.0257885*refscint, 0.0518588*refscint, 0.081431*refscint, 0.114411*refscint, 0.150668*refscint, 0.150668*refscint, 
         0.232295*refscint, 0.324506*refscint, 0.424895*refscint, 0.530595*refscint, 0.636537*refscint, 0.742478*refscint, 0.84842*refscint,
         0.84842*refscint, 1.11327*refscint, 1.37813*refscint, 1.64298*refscint, 1.90784*refscint, 2.17269*refscint, 2.43754*refscint, 
         2.7024*refscint, 2.96725*refscint, 3.23211*refscint, 3.49696*refscint, 3.76181*refscint, 4.02667*refscint, 4.29152*refscint, 
         4.55638*refscint, 4.82123*refscint, 5.08609*refscint},
      { // 5 He3
         0.00217241*refscint, 0.0331461*refscint, 0.0682954*refscint, 
         0.107539*refscint, 0.150743*refscint, 0.197719*refscint, 0.197719*refscint, 0.301959*refscint, 0.417672*refscint, 
         0.541409*refscint, 0.667533*refscint, 0.792031*refscint, 0.916529*refscint, 1.04103*refscint, 1.04103*refscint, 1.35227*refscint, 
         1.66352*refscint, 1.97476*refscint, 2.28601*refscint, 2.59725*refscint, 2.9085*refscint, 3.21974*refscint, 3.53099*refscint, 
         3.84223*refscint, 4.15348*refscint, 4.46472*refscint, 4.77597*refscint, 5.08721*refscint, 5.39846*refscint, 5.7097*refscint, 
         6.02094*refscint},
      { // 6 alpha
         0.00438096*refscint, 0.0184309*refscint, 0.0354221*refscint, 0.0553226*refscint, 0.078079*refscint, 
         0.103617*refscint, 0.103617*refscint, 0.16263*refscint, 0.23134*refscint, 0.308381*refscint, 0.393657*refscint, 0.481043*refscint, 
         0.568428*refscint, 0.655813*refscint, 0.655813*refscint, 0.874276*refscint, 1.09274*refscint, 1.3112*refscint, 1.52967*refscint, 
         1.74813*refscint, 1.96659*refscint, 2.18506*refscint, 2.40352*refscint, 2.62198*refscint, 2.84045*refscint, 3.05891*refscint, 
         3.27737*refscint, 3.49584*refscint, 3.7143*refscint, 3.93276*refscint, 4.15123*refscint}
   };

   // ------------------------------------------------
   // 
   a = 1.01*g/mole;
   G4Element* elH  = new G4Element("Hydrogen", "H",  1.0, a);

   a = 12.01*g/mole;
   G4Element* elC  = new G4Element("Carbon"  , "C",  6.0, a);

   density   = 1.032*g/cm3;
   fScint_mat = new G4Material("Scintillator", density, 2);
   fScint_mat->AddElement(elC, 9);
   fScint_mat->AddElement(elH, 10);

   // ------------------------------------------------
   // 
   G4MaterialPropertiesTable* Scnt_MPT = new G4MaterialPropertiesTable();
   Scnt_MPT->AddProperty("FASTCOMPONENT", fScintPhotonEnergies.data(), fScintLightOutput.data()     , fScintPhotonEnergies.size());
   Scnt_MPT->AddProperty("SLOWCOMPONENT", fScintPhotonEnergies.data(), fScintLightOutput_slow.data(), fScintPhotonEnergies.size());
   Scnt_MPT->AddProperty("RINDEX"       , fScintPhotonEnergies.data(), fScintRefIndex.data()        , fScintPhotonEnergies.size());
   Scnt_MPT->AddConstProperty("SCINTILLATIONYIELD", refscint  );
   Scnt_MPT->AddConstProperty("RESOLUTIONSCALE"   ,   5.0     );
   Scnt_MPT->AddConstProperty("FASTTIMECONSTANT"  ,   0.6*ns  );
   Scnt_MPT->AddConstProperty("SLOWTIMECONSTANT"  , 100.*ns   );
   Scnt_MPT->AddConstProperty("YIELDRATIO"        ,   1.0    );

   Scnt_MPT->AddProperty("ELECTRONSCINTILLATIONYIELD", relativeLightOutput[0].data(), relativeLightOutput[1].data(), relativeLightOutput[0].size());
   Scnt_MPT->AddProperty("PROTONSCINTILLATIONYIELD"  , relativeLightOutput[0].data(), relativeLightOutput[2].data(), relativeLightOutput[0].size());
   Scnt_MPT->AddProperty("DEUTERONSCINTILLATIONYIELD", relativeLightOutput[0].data(), relativeLightOutput[3].data(), relativeLightOutput[0].size());
   Scnt_MPT->AddProperty("TRITONSCINTILLATIONYIELD"  , relativeLightOutput[0].data(), relativeLightOutput[4].data(), relativeLightOutput[0].size());
   Scnt_MPT->AddProperty("HE3SCINTILLATIONYIELD"     , relativeLightOutput[0].data(), relativeLightOutput[5].data(), relativeLightOutput[0].size());
   Scnt_MPT->AddProperty("ALPHASCINTILLATIONYIELD"   , relativeLightOutput[0].data(), relativeLightOutput[6].data(), relativeLightOutput[0].size());

   // The relative strength of the fast component as a fraction of total scintillation yield is given by the YIELDRATIO. 
   fScint_mat->SetMaterialPropertiesTable(Scnt_MPT);

   // -----------------------------------
   // Wavelength shifting fiber material
   density   = 1.032*g/cm3;
   fWLSFiber_mat = new G4Material("WLSFiber_mat", density, 2);
   fWLSFiber_mat->AddElement(elC, 9);
   fWLSFiber_mat->AddElement(elH, 10);

   const G4int nEntries = 96;
   std::array<double, nEntries> PhotonEnergy = {
      1.9*eV,1.92*eV,1.94*eV,1.96*eV,1.98*eV,
      2.*eV,2.02*eV,2.04*eV,2.06*eV,2.08*eV,
      2.1*eV,2.12*eV,2.14*eV,2.16*eV,2.18*eV,
      2.2*eV,2.22*eV,2.24*eV,2.26*eV,2.28*eV,
      2.3*eV,2.32*eV,2.34*eV,2.36*eV,2.38*eV,
      2.4*eV,2.42*eV,2.44*eV,2.46*eV,2.48*eV,
      2.5*eV,2.52*eV,2.54*eV,2.56*eV,2.58*eV,
      2.6*eV,2.62*eV,2.64*eV,2.66*eV,2.68*eV,
      2.7*eV,2.72*eV,2.74*eV,2.76*eV,2.78*eV,
      2.8*eV,2.82*eV,2.84*eV,2.86*eV,2.88*eV,
      2.9*eV,2.92*eV,2.94*eV,2.96*eV,2.98*eV,
      3.*eV,3.02*eV,3.04*eV,3.06*eV,3.08*eV,
      3.1*eV,3.12*eV,3.14*eV,3.16*eV,3.18*eV,
      3.2*eV,3.22*eV,3.24*eV,3.26*eV,3.28*eV,
      3.3*eV,3.32*eV,3.34*eV,3.36*eV,3.38*eV,
      3.4*eV,3.42*eV,3.44*eV,3.46*eV,3.48*eV,
      3.5*eV,3.52*eV,3.54*eV,3.56*eV,3.58*eV,
      3.6*eV,3.62*eV,3.64*eV,3.66*eV,3.68*eV,
      3.7*eV,3.72*eV,3.74*eV,3.76*eV,3.78*eV,
      3.8*eV
   };

   std::array<double, nEntries> RIndexFiber = {
      1.59,1.59,1.59,1.59,1.59, 1.59,1.59,1.59,1.59,1.59,
      1.59,1.59,1.59,1.59,1.59, 1.59,1.59,1.59,1.59,1.59,
      1.59,1.59,1.59,1.59,1.59, 1.59,1.59,1.59,1.59,1.59,
      1.59,1.59,1.59,1.59,1.59, 1.59,1.59,1.59,1.59,1.59,
      1.59,1.59,1.59,1.59,1.59, 1.59,1.59,1.59,1.59,1.59,
      1.59,1.59,1.59,1.59,1.59, 1.59,1.59,1.59,1.59,1.59,
      1.59,1.59,1.59,1.59,1.59, 1.59,1.59,1.59,1.59,1.59,
      1.59,1.59,1.59,1.59,1.59, 1.59,1.59,1.59,1.59,1.59,
      1.59,1.59,1.59,1.59,1.59, 1.59,1.59,1.59,1.59,1.59,
      1.59,1.59,1.59,1.59,1.59, 1.59
   };

   // Absorption length, AbsLen(lambda) = N*(1 - A(lambda)) 
   // where A is the absorption spectrum and N is the normalization
   // point where the absorption length was measred.
   // N = AbsLen(lambda=445nm)/(1-A(445nm))
   std::array<double, nEntries> AbsFiber = {
      15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 
      15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 
      15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 
      15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 
      15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 15.7512*m, 
      15.6116*m, 15.5711*m, 15.5557*m, 15.5093*m, 15.332*m, 14.562*m, 
      12.6202*m, 10.0419*m, 7.77181*m, 4.55662*m, 2.6489*m, 2.18268*m, 
      2.54234*m, 3.13265*m, 3.57652*m, 3.30194*m, 2.29932*m, 1.10024*m, 
      0.488114*m, 0.207704*m, 1.42615*m, 2.63382*m, 3.8237*m, 4.7017*m, 
      5.39579*m, 5.99016*m, 6.42452*m, 6.85888*m, 7.4273*m, 8.08595*m, 
      8.74459*m, 9.42673*m, 10.1186*m, 10.8106*m, 11.4386*m, 12.0395*m, 
      12.6184*m, 13.0177*m, 13.417*m, 13.8164*m, 14.2157*m, 14.4387*m, 
      14.6114*m, 14.7842*m, 14.9569*m, 15.0293*m, 15.0754*m, 15.1216*m, 
      15.1677*m, 15.2139*m, 15.26*m, 15.3061*m, 15.3478*m, 15.3849*m, 
      15.4219*m, 15.459*m, 15.4961*m, 15.5332*m, 15.5703*m, 15.6073*m, 
      15.6444*m, 15.6815*m, 15.7186*m, 15.7512*m, 15.7512*m, 15.7512*m
   };
   //std::for_each(AbsFiber.begin(), AbsFiber.end(), [](double &a){ a*=0.5;});

   std::array<double, nEntries> EmissionFiber = {
      0,0.00721253,0.0189329,0.0226868,0.0263995,
      0.02979,0.0321353,0.0344806,0.0368259,0.0392154,
      0.045909,0.0526026,0.0592962,0.0755797,0.0937883,
      0.117489,0.14813,0.178771,0.215699,0.253016,
      0.286764,0.31805,0.356039,0.443485,0.5372,
      0.630915,0.727269,0.783258,0.762788,0.733277,
      0.70717,0.728632,0.796122,0.890933,0.972089,
      1.00397,0.880085,0.686637,0.471271,0.301391,
      0.193752,0.113452,0.0696161,0.0551339,0.0433935,
      0.0316532,0.0199128,0.00938565,0.00509222,0.000798779,
      0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
      0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
      0,0,0,0,0, 0
   };

   G4MaterialPropertiesTable* MPTFiber = new G4MaterialPropertiesTable();

   MPTFiber->AddProperty("RINDEX"      , PhotonEnergy.data(), RIndexFiber.data()  , nEntries);
   MPTFiber->AddProperty("WLSABSLENGTH", PhotonEnergy.data(), AbsFiber.data()     , nEntries);
   MPTFiber->AddProperty("WLSCOMPONENT", PhotonEnergy.data(), EmissionFiber.data(), nEntries);
   MPTFiber->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);

   fWLSFiber_mat->SetMaterialPropertiesTable(MPTFiber);

}
//______________________________________________________________________________

void ScintTileDetectorGeometry::BuildLogicalVolumes()
{
   using namespace CLHEP;

   // ------------------------------------------------------------------------
   // 
   BuildMaterials();

   // ------------------------------------------------------------------------
   // 
   fScintWidth = CLHEP::twopi*(fRadius-fScintThickness)/double(fNScint) - 0.02*CLHEP::mm;

   std::cout << " fScintWidth = " << fScintWidth/CLHEP::mm << std::endl;

   // ------------------------------------------------------------------------
   // 

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
   double  a                = 0.0;

   G4NistManager* nist = G4NistManager::Instance();

   // ------------------------------------------------------------------------
   // Parameterization
   // ------------------------------------------------------------------------

   // Here we subtract the diamter so that the fiber leaves the face instead
   // of exactly the corner.
   double clearence_dist = 1.5*fWLSFiberDiameter;
   double fiber_angle    = atan((fScintLength-clearence_dist)/fScintThickness);
   double scint_diagonal = sqrt((fScintLength-clearence_dist)*(fScintLength-clearence_dist) + fScintThickness*fScintThickness);

   fWLSFiberEmbedAngle   = fiber_angle;

   fWLSFiber_rot = G4RotationMatrix::IDENTITY;
   fWLSFiber_rot.rotateX(fiber_angle);

   fCurvedFiber_rot = G4RotationMatrix::IDENTITY;
   fCurvedFiber_rot.rotateY(90.0*CLHEP::degree);

   fCurvedWLSF_rot = G4RotationMatrix::IDENTITY;
   fCurvedWLSF_rot.rotateY(90.0*CLHEP::degree);
   fCurvedWLSF_rot.rotateX(-90.0*CLHEP::degree);

   // The position of the center of the wls fiber tube  the extra 1.5*diameter is so that the corner doersn't interfear
   fWLSFiber_pos = {0.0,0.0,-fWLSFiberLength/2.0 + scint_diagonal/2.0 - 1.5*fWLSFiberDiameter };
   fWLSFiber_pos  *= fWLSFiber_rot;

   // Curved fiber (torus)
   double  fiber_curve_angle = 90.0*CLHEP::degree - fiber_angle;

   G4ThreeVector fiber_pos_2 = {0.0, 0.0, -fWLSFiberLength/2.0 };
   fiber_pos_2              *= fWLSFiber_rot;
   G4ThreeVector fiber_pos_5 = {0.0, 0.0, -fWLSFiberLength/2.0 };

   fCurvedFiber_pos = { 0*cm, 0.0, fFiberCurveRadius };
   fCurvedFiber_pos.rotateX( -fiber_curve_angle );
   fCurvedFiber_pos += (fWLSFiber_pos + fiber_pos_2);

   fCurvedWLSF_pos = { 0*cm, fFiberCurveRadius,  -fWLSFiberLength/2.0};
   //fCurvedWLSF_pos.rotateX( -fiber_curve_angle );
   //fCurvedWLSF_pos += ( fiber_pos_5);

   fSiPM_pos = { 0.0, 
      fFiberCurveRadius*(1.0-cos(fiber_curve_angle-fSiPMThickness/2.0/fFiberCurveRadius)),
      -fWLSFiberLength/2.0-fFiberCurveRadius*sin(fiber_curve_angle-fSiPMThickness/2.0/fFiberCurveRadius) };


   G4ThreeVector fiber_pos_3   = { 0*cm, 0.0, fFiberCurveRadius };
   fiber_pos_3      -= fCurvedFiber_pos;


   G4ThreeVector fiber_pos_4 = {0.0, 0.0, -fWLSFiberLength/2.0 };
   fiber_pos_4       *= fWLSFiber_rot;
   fiber_pos_4       += fiber_pos_3;

   // ------------------------------------------------------------------------
   // scintillation radiator target centered at origin
   // ------------------------------------------------------------------------
   red       = 256.0/256.0;
   green     = 1.0/256.0;
   blue      = 1.0/256.0;
   alpha     = 0.4;

   if(fScint_phys) delete fScint_phys;
   if(fScint_log)   delete fScint_log;
   if(fScint_solid) delete fScint_solid;

   G4VSolid * embeded_fiber_solid = new G4Tubs("embeded_fiber_solid", 0.0,fWLSFiberDiameter/2.0, fWLSFiberLength/2.0, 0.0, CLHEP::twopi );
   G4VSolid * scint_box_solid     = new G4Box("scint_box_solid", fScintWidth/2.0, fScintLength/2.0, fScintThickness/2.0 );
   G4VSolid * scint_with_fibers_solid = scint_box_solid;

   // --------------------------------
   // Curved fiber
   fCurvedFiber_solid = new G4Torus("fCurvedFiber_solid", 0.0, fWLSFiberDiameter/2.0,          fFiberCurveRadius ,-fiber_curve_angle , fiber_curve_angle);

   fCurvedFiber_mat  = fWLSFiber_mat;
   //fCurvedFiber_log  = new G4LogicalVolume(fCurvedFiber_solid, fCurvedFiber_mat,"fCurvedFiber_log");
   //fCurvedFiber_phys = new G4PVPlacement(G4Transform3D(fCurvedFiber_rot, fCurvedFiber_pos), 
   //      fCurvedFiber_log, "fCurvedFiber_phys", mother,false,0,checkOverlaps);                                  

   // 4 WLS fibers
   for(int j=0; j<2; j++) {

      G4ThreeVector fiber_extra_pos = {-fScintWidth/2.0 + (fScintWidth/3.0)*double(j+1),-1.5*fWLSFiberDiameter,0.0 };
      //G4ThreeVector curved_fiber_extra_pos = {-fScintWidth/2.0 + (fScintWidth/5.0)*double(1),-1.5*fWLSFiberDiameter,0.0 };

      // union with the WLS tube
      std::string  temp_name =  "scint_with_fibers_solid" + std::to_string(j);
      scint_with_fibers_solid = new G4UnionSolid(temp_name.data(), scint_with_fibers_solid, 
            embeded_fiber_solid, G4Transform3D( fWLSFiber_rot, fWLSFiber_pos + fiber_extra_pos )  );

      // Union with the torus part
      scint_with_fibers_solid = new G4UnionSolid(temp_name.data(), scint_with_fibers_solid, 
            fCurvedFiber_solid, G4Transform3D( fCurvedFiber_rot, fCurvedFiber_pos + fiber_extra_pos )  );
   }

   fScint_solid   = scint_with_fibers_solid;
   fScint_log     = new G4LogicalVolume(fScint_solid, fScint_mat,"fScint_log");
   //fScint_phys  = new G4PVPlacement(0,fScint_pos, fScint_log, "fScint_phys",world_log,false,0,checkOverlaps);                                  

   G4Colour            fScint_color {red, green, blue, alpha };   // Gray 
   G4VisAttributes   * fScint_vis   = new G4VisAttributes(fScint_color);
   fScint_log->SetVisAttributes(fScint_vis);

   G4UserLimits * fScint_limits = new G4UserLimits(0.5*mm);
   fScint_log->SetUserLimits(fScint_limits);

   // ------------------------------------------------------------------------
   // embedded WLS fiber 
   // ------------------------------------------------------------------------
   red       = 0.0/256.0;
   green     = 200.0/256.0;
   blue      = 10.0/256.0;
   alpha     = 0.4;

   if(fWLSFiber_log) delete fWLSFiber_log;
   if(fWLSFiber_solid) delete fWLSFiber_solid;

   fWLSFiber_solid    = new G4Tubs( "fWLSFiber_solid"  , 0.0, fWLSFiberDiameter/2.0 /*-0.005*mm*/, fWLSFiberLength/2.0, 0.0, CLHEP::twopi );
   fCurvedWLSF_solid  = new G4Torus("fCurvedWLSF_solid", 0.0, fWLSFiberDiameter/2.0 /*-0.005*mm*/, fFiberCurveRadius ,0 , fiber_curve_angle);

   fWLSFiber_solid = new G4UnionSolid("fWLSFiber_solid", fWLSFiber_solid, fCurvedWLSF_solid, G4Transform3D( fCurvedWLSF_rot, fCurvedWLSF_pos )  );

   fWLSFiber_log   = new G4LogicalVolume(fWLSFiber_solid, fWLSFiber_mat,"fWLSFiber_log");

   G4Colour            fWLSFiber_color {red, green, blue, alpha };   // Gray 
   G4VisAttributes   * fWLSFiber_vis   = new G4VisAttributes(fWLSFiber_color);
   fWLSFiber_log->SetVisAttributes(fWLSFiber_vis);

   G4UserLimits * fWLSFiber_limits = new G4UserLimits(0.5*mm);
   fWLSFiber_log->SetUserLimits(fWLSFiber_limits);

   // --------------------------------
   // SiPM
   red       = 0.0/256.0;
   green     = 0.0/256.0;
   blue      = 230.0/256.0;
   alpha     = 0.5;

   fSiPM_solid = new G4Tubs("fSiPM_solid", 0.0, fWLSFiberDiameter/2.0-0.01*mm, fSiPMThickness/16.0, 0.0, CLHEP::twopi );
   //fSiPM_pos   = { 0.0, 0.0, -1.0*(fWLSFiberLength/2.0 - 1.0*mm/2.0) };

   fSiPM_mat  = fWLSFiber_mat;
   fSiPM_rot.rotateX(fiber_curve_angle);
   fSiPM_log  = new G4LogicalVolume(fSiPM_solid, fSiPM_mat,"fSiPM_log");
   fSiPM_phys = new G4PVPlacement(G4Transform3D(fSiPM_rot,fSiPM_pos), fSiPM_log, "fSiPM_phys", fWLSFiber_log,false,0,checkOverlaps);                                  

   G4Colour            fSiPM_color {red, green, blue, alpha };
   G4VisAttributes   * fSiPM_vis   = new G4VisAttributes(fSiPM_color);
   fSiPM_log->SetVisAttributes(fSiPM_vis);

   //G4UserLimits * fWLSFiber_limits = new G4UserLimits(0.5*mm);
   //fWLSFiber_log->SetUserLimits(fWLSFiber_limits);

   if(!fSiPM_det) fSiPM_det = new SiPMSD("SiPM");
   fSiPM_det->SetCopyNoParent(2);
   //SetSensitiveDetector("fPhotonDet1_log",fPhotonDet1_det);
   G4SDManager* SDMan = G4SDManager::GetSDMpointer();
   SDMan->AddNewDetector(fSiPM_det);
   fSiPM_log->SetSensitiveDetector(fSiPM_det);

   //
   // 2 WLS fibers
   for(int j=0; j<2; j++){

      G4ThreeVector fiber_extra_pos = {-fScintWidth/2.0 + (fScintWidth/3.0)*double(j+1),-1.5*fWLSFiberDiameter,0.0 };
      std::string   temp_name       =  "fWLSFiber_phys" + std::to_string(j);

      fWLSFibers_phys[j]  = new G4PVPlacement( G4Transform3D( fWLSFiber_rot, fWLSFiber_pos +fiber_extra_pos ), fWLSFiber_log,temp_name.data(),fScint_log,false,j,checkOverlaps);                                  
   }

   red       = 0.0/256.0;
   green     = 0.0/256.0;
   blue      = 230.0/256.0;
   alpha     = 0.5;

   //double  fiber_curve_angle = 90.0*CLHEP::degree - fiber_angle;
   //fCurvedFiber_solid = new G4Torus("fCurvedFiber_solid", 0.0, fWLSFiberDiameter/2.0, fFiberCurveRadius ,-fiber_curve_angle , fiber_curve_angle);

   //G4ThreeVector fiber_pos_2 = {0.0, 0.0, -fWLSFiberLength/2.0 };
   ////G4ThreeVector fiber_pos_2 = {0.0,0.0,-fWLSFiberLength - scint_diagonal + 2.0*1.5*fWLSFiberDiameter };
   //fiber_pos_2       *= fiber_rot;

   //fCurvedFiber_pos   = { 0*cm, 0.0, fFiberCurveRadius };
   //fCurvedFiber_pos.rotateX(-fiber_curve_angle);

   //G4ThreeVector curved_fiber_extra_pos = {-fScintWidth/2.0 + (fScintWidth/5.0)*double(1),-1.5*fWLSFiberDiameter,0.0 };
   //fCurvedFiber_pos += (fiber_pos + fiber_pos_2 + curved_fiber_extra_pos);

   //fCurvedFiber_mat  = fWLSFiber_mat;
   //fCurvedFiber_log  = new G4LogicalVolume(fCurvedFiber_solid, fCurvedFiber_mat,"fCurvedFiber_log");

}
//______________________________________________________________________________

G4VPhysicalVolume * ScintTileDetectorGeometry::PlacePhysicalVolume(
      G4LogicalVolume * mother, 
      G4VPhysicalVolume * mother_phys,
      G4VPhysicalVolume * adjacent_phys  )
{
   bool checkOverlaps = false;

   //fScint_phys  = new G4PVPlacement(0,fScint_pos, fScint_log, "fScint_phys",mother,false,0,checkOverlaps);                                  
   // scintillator surface
   const    G4int NUM           = 2;
   G4double pp2[NUM]            = {1.6*eV, 8.0*eV};
   G4double specularlobe2[NUM]  = {0.3, 0.3};
   G4double specularspike2[NUM] = {0.2, 0.2};
   G4double backscatter2[NUM]   = {0.1, 0.1};
   G4double rindex2[NUM]        = {1.58, 1.58};
   G4double rindex3[NUM]        = {1.60, 1.60};
   G4double reflectivity2[NUM]  = {0.88, 0.88};
   G4double reflectivity3[NUM]  = {0.99, 0.99};
   G4double efficiency2[NUM]    = {0.0, 0.0};

   // Scintillator internal reflective surface 
   G4OpticalSurface * OpSurface = new G4OpticalSurface("scintillator fiber");
   OpSurface->SetType(dielectric_metal);
   OpSurface->SetFinish(ground);
   OpSurface->SetModel(glisur);

   G4MaterialPropertiesTable *OpSurfaceProperty = new G4MaterialPropertiesTable();
   OpSurfaceProperty->AddProperty("REFLECTIVITY",pp2,reflectivity2,NUM);
   OpSurfaceProperty->AddProperty("EFFICIENCY",  pp2,efficiency2,NUM);
   OpSurfaceProperty->AddProperty("RINDEX",      pp2,rindex2,NUM);
   OpSurface->SetMaterialPropertiesTable(OpSurfaceProperty);

   // WLS fiber - Scintillator surface 
   G4OpticalSurface * OpSurface2 = new G4OpticalSurface("scintillator WLS fiber");
   OpSurface2->SetType(dielectric_dielectric);
   OpSurface2->SetModel(glisur);
   //OpSurface2->SetFinish(polishedlumirrorglue);
   OpSurface2->SetFinish(polished);
   G4MaterialPropertiesTable *OpSurfaceProperty2 = new G4MaterialPropertiesTable();
   OpSurfaceProperty2->AddProperty("RINDEX",      pp2,rindex3,NUM);
   OpSurface2->SetMaterialPropertiesTable(OpSurfaceProperty2);

   // Scintillator internal reflective surface 
   G4OpticalSurface * OpSurface3 = new G4OpticalSurface("fiber transport");
   OpSurface3->SetType(dielectric_metal);
   OpSurface3->SetFinish(ground);
   OpSurface3->SetModel(glisur);

   G4MaterialPropertiesTable * OpSurfaceProperty3 = new G4MaterialPropertiesTable();
   OpSurfaceProperty3->AddProperty("REFLECTIVITY", pp2, reflectivity3, NUM);
   OpSurfaceProperty3->AddProperty("EFFICIENCY"  , pp2, efficiency2  , NUM);
   OpSurfaceProperty3->AddProperty("RINDEX"      , pp2, rindex2      , NUM);
   OpSurface3->SetMaterialPropertiesTable(OpSurfaceProperty3);

   // ---------------------------------------------------------------

   G4ThreeVector s1_pos = {-fScintThickness/2.0 + 0.01*CLHEP::mm + fRadius,0.0, fScintLength*(double(fNScintPerStrip)/2.0-0.5)};
   //G4ThreeVector step( 0.0, fScintLength + fScintGap, 0.0 );
   double delta_phi = 360.0*CLHEP::degree/double(fNScint);

   for(int i = 0; i<fNScint*fNScintPerStrip; i++) {

      fScint_positions[i] = s1_pos;

      //fScint_positions[i] = fScint_pos + double(i)*step;

      fScint_rotations[i] = G4RotationMatrix::IDENTITY;
      fScint_rotations[i].rotateX(-90.0*CLHEP::degree);
      fScint_rotations[i].rotateZ(90.0*CLHEP::degree + delta_phi*double(i));


      G4String name       = "fScint_physicals_" + std::to_string(i);
      //std::cout << name <<std::endl;
      fScint_physicals[i] = new G4PVPlacement(
            G4Transform3D(fScint_rotations[i],fScint_positions[i]),
            fScint_log, name,mother,false,i,checkOverlaps); 

      name = "fScint_borders_" + std::to_string(i);
      fScint_borders[i]  = new G4LogicalBorderSurface(name, fScint_physicals[i], mother_phys, OpSurface );

      for(int j=0; j<2; j++){

         std::string name = "WLSFiberIntoScint_border_" + std::to_string(i);
         new G4LogicalBorderSurface(name, fWLSFibers_phys[j],      fScint_physicals[i], OpSurface2 );

         name = "ScintIntoWLSFiber_border_" + std::to_string(i);
         new G4LogicalBorderSurface(name, fScint_physicals[i], fWLSFibers_phys[j],      OpSurface2 );

         name = "FiberTransport_border_" + std::to_string(i);
         new G4LogicalBorderSurface(name, fWLSFibers_phys[j], mother_phys ,    OpSurface3 );

      }

      // Rotate to next position
      s1_pos.rotateZ(delta_phi);
      if( (i+1)%fNScint == 0 ) {
         s1_pos += G4ThreeVector(0,0,-fScintLength-0.005*CLHEP::mm);
      }
   }

   return fScint_physicals[0];
}
//______________________________________________________________________________


