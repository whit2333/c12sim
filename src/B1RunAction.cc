#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1Run.hh"
#include "B1Analysis.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include <string>
#include <sstream>
#include "DAQManager.h"
#include "SimulationManager.h"
#include "Crate.h"
#include "Module.h"
#include "B1EventAction.hh"
#include "TSystem.h"


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
{
   //std::cout << " start B1RunAction::GenerateRun\n";
   std::cout << " Sensitive Detectors:\n";
   G4SDManager::GetSDMpointer()->ListTree();

   G4RunManager      * runManager = G4RunManager::GetRunManager();
   SimulationManager * simManager = SimulationManager::GetInstance();
   //simManager->SetRunNumber(fRunNumber);

   G4cout << " - Creating Run Number " << fRunNumber << "   " << G4endl;

   std::string filename = simManager->OutputFileName();

   //G4cout << " - file " << filename << "   " << G4endl;
   //G4cout << " - Creating dir " << simManager->fOutputDirectoryName << "   " << G4endl;
   gSystem->mkdir(simManager->fOutputDirectoryName.c_str(),true);

   simManager->fOutputFile = new TFile(filename.c_str(), "RECREATE" );
   simManager->fOutputFile->cd();

   simManager->fOutputTree = new TTree(simManager->fOutputTreeName.c_str(),"Clas12 simulation output");
   simManager->fOutputTree->Branch(
         "HitsEvent",
         "clas12::hits::CLAS12HitsEvent",
         &(simManager->fEvent)   );

   fRunConf = new clas12::sim::RunConfiguration(fRunNumber);
   fRunConf->fInputTreeName = simManager->OutputTreeName();
   fRunConf->fInputFileName = simManager->OutputFileName();

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

   using namespace clas12::DAQ;

   G4cout << "\e[1m=================== \e[41mEND RUN #" << fRunNumber << "\e[49m===================\e[0m" << G4endl;

   SimulationManager * simManager = SimulationManager::GetInstance();

   fRunConf->fNSimulated = simManager->fEvent->fEventNumber;
   fRunConf->Print();

   //DAQManager      * daq_manager   = clas12::DAQ::DAQManager::GetManager();
   //Crate           * crate         = daq_manager->GetCrate(0);
   //Module<Scaler>  * scaler_module = crate->GetCrateModule<Scaler>(2);
   //scaler_module->Print("v");

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

