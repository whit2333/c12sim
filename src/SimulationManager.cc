#include "SimulationManager.h"
#include "SimulationMessenger.h"

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4ScoringManager.hh"
#include "G4VSDFilter.hh"
#include "G4PSPassageCellCurrent.hh"
#include "G4PSFlatSurfaceFlux.hh"
#include "G4PSTrackLength.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include <sys/file.h>
#include "RecoilChamberDetectorGeometry.h"
#include "DriftChamberDetectorGeometry.h"
#include "HTCCDetectorGeometry.h"


SimulationManager* SimulationManager::fgSimulationManager = 0;
//______________________________________________________________________________

SimulationManager::SimulationManager () {

   fSimulationMessenger = new SimulationMessenger ( this );
   fEvent               = new clas12::hits::CLAS12HitsEvent();
   fOutputFile          = 0;
   fOutputTree          = 0;
   fOutputDirectoryName = "data/rootfiles";
   fOutputFileName      = "clas12sim";
   fOutputTreeName      = "clasdigi_hits";
   fRecoilChamberGeo    = nullptr;
   fDriftChamberGeo     = nullptr;
   fHTCCGeo             = nullptr;
   fInputFileName       = "data/lundfiles/eg_LH2_full_4.lund";
}
//______________________________________________________________________________

SimulationManager::~SimulationManager()
{
   delete fSimulationMessenger;
}
//______________________________________________________________________________

SimulationManager* SimulationManager::GetInstance (  )
{
   if ( fgSimulationManager == 0 )
   {
      fgSimulationManager = new SimulationManager ( ) ;
   }
   return fgSimulationManager;
}
//______________________________________________________________________________

int SimulationManager::ReadRunNumber(const char * fname )
{
   
   std::string  lock_filename = Form(".lock_file_%s",fname);
   int fd = tryGetLock( lock_filename.c_str()  );
   while( fd == -1 ) {
      fd = tryGetLock(lock_filename.c_str());
      //std::cout << "file locked\n";
      std::this_thread::sleep_for (std::chrono::seconds(1));
   }
   int run_number = 0;
   //std::this_thread::sleep_for (std::chrono::seconds(20));

   std::ifstream input_file;
   input_file.open ( fname );
   input_file >> run_number ;
   input_file.close(); 
   run_number++;
   std::ofstream output_file;
   output_file.open ( fname ,std::ios::trunc); // this incremtents a number so that it acts like a normal DAQ
   output_file << run_number ;
   output_file.close();
   releaseLock(fd, lock_filename.c_str());
   SetRunNumber(run_number);
   return(run_number);
}
//______________________________________________________________________________

std::string SimulationManager::OutputFileName() const {
   std::string filename = fOutputDirectoryName;
   filename += "/";
   filename += fOutputFileName;
   filename += std::to_string(GetRunNumber());
   filename += ".root";
   return filename;
}
//______________________________________________________________________________

std::string SimulationManager::OutputTreeName() const {
   std::string filename = fOutputTreeName;
   return filename;
}
//______________________________________________________________________________

std::string SimulationManager::InputFileName() const {
   std::string filename = fInputFileName;
   return filename;
}
//______________________________________________________________________________

int SimulationManager::tryGetLock( char const *lockName )
{
   mode_t m = umask( 0 );
   int fd = open( lockName, O_RDWR|O_CREAT, 0666 );
   umask( m );
   if( fd >= 0 && flock( fd, LOCK_EX | LOCK_NB ) < 0 )
   {
      close( fd );
      fd = -1;
   }
   return fd;
}
//______________________________________________________________________________

void SimulationManager::releaseLock( int fd, char const *lockName )
{
   if( fd < 0 )
      return;
   remove( lockName );
   close( fd );
}
//______________________________________________________________________________

DriftChamberDetectorGeometry * SimulationManager::GetDriftDetectorGeometry(){
   if(!fDriftChamberGeo) {
      fDriftChamberGeo = new DriftChamberDetectorGeometry();
   }
   return fDriftChamberGeo;
}
//______________________________________________________________________________

RecoilChamberDetectorGeometry * SimulationManager::GetRecoilDetectorGeometry(){
   if(!fRecoilChamberGeo) {
      fRecoilChamberGeo = new RecoilChamberDetectorGeometry();
   }
   return fRecoilChamberGeo;
}
//______________________________________________________________________________

HTCCDetectorGeometry * SimulationManager::GetHTCCDetectorGeometry(){
   if(!fHTCCGeo) {
      fHTCCGeo = new HTCCDetectorGeometry();
   }
   return fHTCCGeo;
}
//______________________________________________________________________________

//int SimulationManager::InitializeNewRun(int run) {
//   fRunNumber = RetrieveRunNumber(run);
//   return(fRunNumber);
//}
////______________________________________________________________________________
//
//void SimulationManager::Dispose()
//{
//   if ( fgSimulationManager != 0 )
//   {
//      delete fgSimulationManager;
//      fgSimulationManager = 0;
//   }
//}
////_________________________________________________________________//
//
//void SimulationManager::showPlot ( int arg )
//{
//   /*   plotVis = arg;*/
//}
////_________________________________________________________________//
//
//void SimulationManager::write()
//{
//
//}
////_________________________________________________________________//
//
//void SimulationManager::SetDetectorVerbosity(const char * detName, int level) {
//   if( !strcmp(detName,"GasCherenkov") ) {
//      std::cout << "Setting Gas Cherenkov verbosity level to " << level << "\n";
//      fGasCherenkovVerbosity = level;
//   } else if( !strcmp(detName , "Bigcal") ) {
//      std::cout << "Setting Bigcal verbosity level to " << level << "\n";
//      fBigcalVerbosity = level;
//   } else {
//      std::cout << "No such detector, " << detName << "\n";
//   }
//}
////_________________________________________________________________//
//
//int SimulationManager::GetDetectorVerbosity(const char * detName) {
//
//   if( !strcmp(detName,"GasCherenkov") ) {
//      return(fGasCherenkovVerbosity);
//   } else if( !strcmp(detName , "Bigcal") ) {
//      return(1);
//      //    return(fBigcalVerbosity);
//   }  else {
//      std::cout << "No such detector, " << detName << "\n";
//   }
//   return(-1);
//}
////______________________________________________________________________________
//int SimulationManager::RetrieveRunNumber(int run) {
//   // If the run argument is supplied the run is simply set to that number.
//   // otherwise the run file is used (and incrmeneted).
//   if(run == 0){
//      std::ifstream input_file;
//      input_file.open ( "run.txt" );
//      input_file >> fRunNumber ;
//      input_file.close(); 
//      IncrementRunNumber();
//   } else {
//      fRunNumber = run;
//   }
//   return(fRunNumber);
//}
////_________________________________________________________________//
//
//int SimulationManager::IncrementRunNumber()   
//{
//   fRunNumber++;
//   std::ofstream output_file;
//   output_file.open ( "run.txt" ,std::ios::trunc); // this incremtents a number so that it acts like a normal DAQ
//   output_file << fRunNumber ;
//   output_file.close();
//   return fRunNumber;
//}
////_________________________________________________________________//
//
//int SimulationManager::Reset()  {
//
//   // if(fDetectorTree) delete fDetectorTree;
//   // fDetectorTree=0;
//   // if(fRootFile) delete fRootFile;
//   // fRootFile=0;
//
//   //    if(betaEvent)   betaEvent->FreeMemory();
//   //    if(mcEvent)    mcEvent->FreeMemory();
//   //   if(betaEvent) delete betaEvent;
//   //   if(hmsEvent) delete hmsEvent;
//   //   if(beamEvent) delete beamEvent;
//   //   if(mcEvent) delete mcEvent;
//   return(0);
//}
////_________________________________________________________________//
//
//int SimulationManager::InitScoring()  {
//
//   //G4SDManager * sensitiveDetManager = G4SDManager::GetSDMpointer();
//   G4String filterName, particleName;
//
//
//   // Scoring and sensitive volumes
//
//   fTrackerScoring = 
//      new G4MultiFunctionalDetector("tracker");
//   fCherenkovScoring = 
//      new G4MultiFunctionalDetector("cherenkov");
//   fBigcalScoring = 
//      new G4MultiFunctionalDetector("bigcal");
//   fHodoscopeScoring = 
//      new G4MultiFunctionalDetector("hodoscope");
//   //
//
//
//   // Sensitive volume filters
//   protonFilter = new G4SDParticleFilter(filterName="protonFilter");//,"proton");
//   //    protonFilter->add("proton");
//   electronFilter = new G4SDParticleFilter(filterName="electronFilter");
//   //    electronFilter->add("e-");
//
//   chargeFilter = new G4SDChargedFilter("chargeFilter");
//   opticalPhotonFilter = new G4SDParticleFilter(filterName="opticalPhotonFilter");
//   //    opticalPhotonFilter->add("opticalphoton");
//   bigcalEnergyFilter = new G4SDKineticEnergyFilter("bigcalEnergyThreshFilter");
//
//   // DefineScoringFilters();
//
//   // scoring primitives
//   protonSurfFlux   = new G4PSFlatSurfaceFlux("pSurfFlux",1);
//   protonSurfFlux->SetFilter(protonFilter);
//
//   electronSurfFlux   = new G4PSPassageCellCurrent("eHit");
//   electronSurfFlux->SetFilter(electronFilter);
//
//   chargeSurfFlux   = new G4PSPassageCellCurrent("chargeHit");
//   chargeSurfFlux->SetFilter(chargeFilter);
//
//   photonSurfFlux   = new G4PSPassageCellCurrent("photonSurf");
//   photonSurfFlux->SetFilter(opticalPhotonFilter);
//
//   electronTracklength   = new G4PSTrackLength("eTrackLength");
//   electronTracklength->SetFilter(electronFilter);
//
//   calEnergyDeposit = new G4PSEnergyDeposit("bcEnergyDeposit");
//   calEnergyDeposit->SetFilter(bigcalEnergyFilter);
//
//   // Register Scoring volume with primitive(s)
//   fTrackerScoring->RegisterPrimitive(chargeSurfFlux);
//   fHodoscopeScoring->RegisterPrimitive(chargeSurfFlux);
//   fCherenkovScoring->RegisterPrimitive(photonSurfFlux);
//   calEnergyDeposit->SetMultiFunctionalDetector(fBigcalScoring);
//   fBigcalScoring->RegisterPrimitive(calEnergyDeposit);
//
//   fBigcalScoring->RegisterPrimitive(protonSurfFlux);
//
//   // Below does not work!
//   //  fBigcalDetector->RegisterPrimitive((G4VPrimitiveScorer*)calEnergyDeposit);
//
//   /*
//      sensitiveDetManager->AddNewDetector(fTrackerDetector);
//      sensitiveDetManager->AddNewDetector(fCherenkovDetector);
//      sensitiveDetManager->AddNewDetector(fBigcalDetector);
//      sensitiveDetManager->AddNewDetector(fHodoscopeDetector);*/
//
//   return(0);
//}
////______________________________________________________________________________
//int SimulationManager::DefineScoringFilters() {
//   // Sensitive volume filters
//   protonFilter->add("proton");
//
//   electronFilter->add("e-");
//
//   opticalPhotonFilter->add("opticalphoton");
//
//   bigcalEnergyFilter->SetKineticEnergy(0.010*GeV,6.0*GeV);
//
//   return(0);
//}
////______________________________________________________________________________
//int SimulationManager::AddDetectors(int runNumber) {
//   // Creates the detectors
//
//   fBETADetectorPackage = new BETADetectorPackage("BETADetectorPackage",runNumber);
//
//   fCherenkovDetector = fBETADetectorPackage->fCherenkovDetector;//new GasCherenkovDetector(runNumber);
//   fCherenkovDetector->SetEventAddress(fEvents->BETA->fGasCherenkovEvent);
//
//   fBigcalDetector = fBETADetectorPackage->fBigcalDetector;//= new BigcalDetector(runNumber);
//   fBigcalDetector->SetEventAddress(fEvents->BETA->fBigcalEvent);
//
//   fHodoscopeDetector = fBETADetectorPackage->fHodoscopeDetector;//= new LuciteHodoscopeDetector(runNumber);
//   fHodoscopeDetector->SetEventAddress(fEvents->BETA->fLuciteHodoscopeEvent);
//
//   fTrackerDetector = fBETADetectorPackage->fTrackerDetector;//= new ForwardTrackerDetector(runNumber);
//   fTrackerDetector->SetEventAddress(fEvents->BETA->fForwardTrackerEvent);
//
//   return(0);
//}
////______________________________________________________________________________
//void SimulationManager::PrintSummary(){
//   std::cout << " fSimulateCherenkovOptics " << fSimulateCherenkovOptics << std::endl;
//   std::cout << " fSimulateTrackerOptics   " << fSimulateTrackerOptics   << std::endl;
//   std::cout << " fSimulateHodoscopeOptics " << fSimulateHodoscopeOptics << std::endl;
//   std::cout << " fSimulateBigcalOptics    " << fSimulateBigcalOptics    << std::endl;
//   std::cout << " fSimulateTrigger         " << fSimulateTrigger         << std::endl;
//}
////______________________________________________________________________________


