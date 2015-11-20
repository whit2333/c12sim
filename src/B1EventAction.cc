#include "B1EventAction.hh"
#include "B1Run.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"


B1EventAction::B1EventAction() : G4UserEventAction(), fEdep(0.)
{ 
   event_number = 0;
   fTree = 0;
   SimulationManager * simManager = SimulationManager::GetInstance();
   fCLAS12HitsEvent =  simManager->fEvent;
}
//______________________________________________________________________________

B1EventAction::~B1EventAction()
{}
//______________________________________________________________________________

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
   fCLAS12HitsEvent->Clear();
   fCLAS12HitsEvent->fEventNumber            = event_number;
   fCLAS12HitsEvent->fHTCCEvent.fEventNumber = event_number;

   if(!fTree) {
      SimulationManager * simManager = SimulationManager::GetInstance();
      fTree = simManager->fOutputTree;
   }

   if(event_number%100 == 0 )
   {
      SimulationManager * simManager = SimulationManager::GetInstance();
      std::cout  << " Begin of event " << event_number 
         << " Run Number : " <<  simManager->GetRunNumber()
         << std::endl;
      std::cout << " Random Number: " << G4UniformRand() << std::endl;
   }
}
//______________________________________________________________________________

void B1EventAction::EndOfEventAction(const G4Event*)
{   

   fTree->Fill();
   //if(evtN%10 == 0 )
   //// Increase event number. Notice: this is different than evt->GetEventID()
   //evtN++;
   //return;


   if(event_number%10 == 0 )
      std::cout <<  " End of Event " << event_number << " Routine..." << std::endl;

   // Increase event number. Notice: this is different than evt->GetEventID()
   event_number++;
  // Called after G4Run::RecordEvent

  // accumulate statistics in B1Run
  //B1Run* run = static_cast<B1Run*>( G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  //run->AddEdep(fEdep);
  //if( run->GetNumberOfEvent()%1000 == 0 ) std::cout << "EndOfEventAction\n";
}
//______________________________________________________________________________

