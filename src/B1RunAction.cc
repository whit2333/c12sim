#include "B1RunAction.hh"

#include "B1PrimaryGeneratorAction.hh"
#include "SimplePrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1Run.hh"
#include "B1Analysis.hh"
#include "dollar.hpp"
#include "SiPMSD.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include <string>
#include <sstream>
#include <fstream>
#include "DAQManager.h"
#include "SimulationManager.h"
#include "Crate.h"
#include "Module.h"
#include "B1EventAction.hh"
#include "TSystem.h"
#include "G4TrackingManager.hh"


using ss = std::stringstream;

//______________________________________________________________________________

B1RunAction::B1RunAction(G4int rn) : G4UserRunAction(),
   fRunNumber(rn), fRun(nullptr)
{ 
   //fOutputFile               = 0;

   SimulationManager * simManager = SimulationManager::GetInstance();
   simManager->SetRunNumber(fRunNumber);
}
//______________________________________________________________________________

B1RunAction::~B1RunAction()
{ }
//______________________________________________________________________________

G4Run* B1RunAction::GenerateRun()
{ $
   //std::cout << " start B1RunAction::GenerateRun\n";
   std::cout << " Sensitive Detectors:\n";
   G4SDManager::GetSDMpointer()->ListTree();

   G4RunManager      * runManager = G4RunManager::GetRunManager();
   SimulationManager * simManager = SimulationManager::GetInstance();
   //simManager->SetRunNumber(fRunNumber);

   std::cout << " - Creating Run Number " << fRunNumber << "   " << std::endl;
   std::string filename = simManager->OutputFileName();
   gSystem->mkdir(simManager->fOutputDirectoryName.c_str(),true);

   simManager->fOutputFile = new TFile(filename.c_str(), "RECREATE" );
   simManager->fOutputFile->cd();

   simManager->fOutputTree = new TTree(simManager->fOutputTreeName.c_str(),"Clas12 simulation output");

   std::cout << " tree created\n";

   simManager->fOutputTree->Branch(
         "HitsEvent",
         "clas12::hits::CLAS12HitsEvent",
         &(simManager->fEvent)   );
   std::cout << " b1 created\n";
   simManager->fOutputTree->Branch(
         "RHEvent",
         "clas12::hits::RecoilScintEvent",
         &(simManager->fEvent->fRHEvent)   );
   std::cout << " b2 created\n";
   /*simManager->fOutputTree->Branch(
         "TrajectoryVerticies",
         &(simManager->fTrajectoryVerticies) );*/


   // ---------------------------------------------------------
   // Create the EG tree if it is a certain Primary generator

   const G4VUserPrimaryGeneratorAction * gen_action = runManager->GetUserPrimaryGeneratorAction();
   if( const B1PrimaryGeneratorAction * prigen = dynamic_cast<const B1PrimaryGeneratorAction*>(gen_action) ) {
      simManager->fOutputTree->Branch(
            "PrimaryEvent",
            "clas12::sim::ThrownEvent",
            &(prigen->fThrownEvent)   );
   }
   else if( const SimplePrimaryGeneratorAction * prigen = dynamic_cast<const SimplePrimaryGeneratorAction*>(gen_action) ) {
      simManager->fOutputTree->Branch(
            "PrimaryEvent",
            "clas12::sim::ThrownEvent",
            &(prigen->fThrownEvent)   );
   }


   if( fRunConf ) delete fRunConf;
   fRunConf = new clas12::sim::RunConfiguration(fRunNumber);
   fRunConf->fOutputTreeName = simManager->OutputTreeName();
   fRunConf->fOutputFileName = simManager->OutputFileName();
   fRunConf->fInputFileName  = simManager->InputFileName();

   fRun = new B1Run(fRunNumber);
   //fRun = new G4Run();

   //fOutputFile = new TFile(Form("c12_sim_output_%d.root",fRunNumber), "UPDATE");
   //fOutputFile->cd();


   //B1Run * run = new B1Run(fRunNumber);

   return fRun; 
}
//______________________________________________________________________________

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
   G4cout <<"=================== RUN #" << fRunNumber << " ===================" << G4endl;

   //inform the runManager to save random number seed
   G4RunManager::GetRunManager()->SetRandomNumberStore(false);

   // Get analysis manager
   //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

   // Open an output file
   ss file_name;
   file_name << "c12sim_output_" << fRunNumber;
   //analysisManager->OpenFile(file_name.str().c_str());

   SimulationManager * simManager = SimulationManager::GetInstance();
   simManager->PrintConfig();
   fRunConf->Print();
   fRunConf->fInputFileName = simManager->OutputFileName();

}
//______________________________________________________________________________

void B1RunAction::EndOfRunAction(const G4Run* run)
{
   G4int nofEvents = run->GetNumberOfEvent();
   if (nofEvents == 0) return;

   const B1Run* b1Run = static_cast<const B1Run*>(run);

   // Run conditions
   //  note: There is no primary generator action object for "master"
   //        run manager for multi-threaded mode.
   const B1PrimaryGeneratorAction* generatorAction = static_cast<const B1PrimaryGeneratorAction*>
      (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

   //G4String runCondition;
   //if (generatorAction)
   //{
   //   const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
   //   runCondition += particleGun->GetParticleDefinition()->GetParticleName();
   //   runCondition += " of ";
   //   G4double particleEnergy = particleGun->GetParticleEnergy();
   //   runCondition += G4BestUnit(particleEnergy,"Energy");
   //}

   //if (IsMaster()) {
   //   G4cout
   //      << G4endl
   //      << "--------------------End of Global Run-----------------------";
   //   //fOutputFile->Write();
   //   //fOutputFile->Close();
   //   // Save histograms
   //   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   //   analysisManager->Write();
   //   analysisManager->CloseFile();
   //}
   //else {
   //   G4cout
   //      << G4endl
   //      << "--------------------End of Local Run------------------------";
   //}
   //}

   using namespace clas12::DAQ;

   G4cout << "\e[1m=================== \e[41mEND RUN #" << fRunNumber << "\e[49m===================\e[0m" << G4endl;

   SimulationManager * simManager = SimulationManager::GetInstance();

   fRunConf->fNSimulated = simManager->fEvent->fEventNumber + 1;
   fRunConf->Print();
   simManager->PrintConfig();

   std::ofstream log_out(std::string("data/log/run")+std::to_string(fRunNumber)+std::string(".log"));
   simManager->PrintConfig(log_out);
   fRunConf->Print(log_out);


   //DAQManager      * daq_manager   = clas12::DAQ::DAQManager::GetManager();
   //Crate           * crate         = daq_manager->GetCrate(0);
   //Module<Scaler>  * scaler_module = crate->GetCrateModule<Scaler>(2);
   //scaler_module->Print("v");

   auto SD_man = G4SDManager::GetSDMpointer();
   SiPMSD* pm     = dynamic_cast<SiPMSD*>(SD_man->FindSensitiveDetector("/PM3"));
   if(pm) {
      std::cout << "Tile Scintillator neutron dose:" << std::endl;
      std::cout << "  neutron  count: " << pm->fNeutronCount << std::endl;
      std::cout << "          energy: " << pm->fNeutronEnergy << std::endl;
      std::cout << "      energy dep: " << pm->fNeutronEnergyDep << std::endl;
   }


   simManager->fOutputTree->Write();
   simManager->fOutputFile->WriteObject(fRunConf,"Run_Config");
   simManager->fOutputFile->Write();
   simManager->fOutputFile->Close();

   // Get the event action and reset it
   // This is needed so that if another run is executed, the tree pointer  is
   // reloaded. This is achieved by setting it to zero with MEventAction::Reset()
   G4RunManager * runManager   = G4RunManager::GetRunManager();
   B1EventAction * evAction  = (B1EventAction*)(runManager->GetUserEventAction());
   if(evAction ) evAction->Reset();

   simManager->SetRunNumber(simManager->GetRunNumber()+1);

}
//______________________________________________________________________________

