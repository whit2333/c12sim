#include "B1EventAction.hh"
#include "B1Run.hh"

#include "dollar.hpp"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryPoint.hh"
#include "TParticle.h"


B1EventAction::B1EventAction() : G4UserEventAction(), fEdep(0.)
{ 
   event_number = 0;
   fTree        = nullptr;
   //fEGTree      = nullptr;
   SimulationManager * simManager = SimulationManager::GetInstance();
   fCLAS12HitsEvent =  simManager->fEvent;
   fTrajectoryVerticies = simManager->fTrajectoryVerticies;
}
//______________________________________________________________________________

B1EventAction::~B1EventAction()
{}
//______________________________________________________________________________

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
   // Clear the event
   fEdep = 0.;
   fCLAS12HitsEvent->Clear();

   // Set the run and event numbers
   fCLAS12HitsEvent->fEventNumber            = event_number;
   fCLAS12HitsEvent->fHTCCEvent.fEventNumber = event_number;
   fTrajectoryVerticies->Clear();

   if(!fTree) {
      SimulationManager * simManager = SimulationManager::GetInstance();
      fTree   = simManager->fOutputTree;
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

void B1EventAction::EndOfEventAction(const G4Event * event)
{ $ 

   using namespace CLHEP;
   G4TrajectoryContainer * traj_container = event->GetTrajectoryContainer();

   if( traj_container ) { $
      int n_traj = traj_container->entries();
      //std::cout << "Event has " << n_traj << " trajectores\n";

      for(int i = 0; i<n_traj; i++ ) { $

         G4VTrajectory * traj = (*traj_container)[i];
         TParticle     * part =  new( (*fTrajectoryVerticies)[i] ) TParticle();
         //TParticle(Int_t pdg, Int_t status, Int_t mother1, Int_t mother2, Int_t daughter1, Int_t daughter2, 
         //Double_t px, Double_t py, Double_t pz, Double_t etot, Double_t vx, Double_t vy, Double_t vz, Double_t time)
         part->SetMomentum(
               traj->GetInitialMomentum().x()/GeV,
               traj->GetInitialMomentum().y()/GeV,
               traj->GetInitialMomentum().z()/GeV,
               0);
         part->SetProductionVertex(
               traj->GetPoint(0)->GetPosition().x()/cm,
               traj->GetPoint(0)->GetPosition().y()/cm,
               traj->GetPoint(0)->GetPosition().z()/cm,
               0);
         part->SetFirstMother( traj->GetTrackID()     );
         part->SetLastMother(  traj->GetParentID()    );
         part->SetPdgCode(     traj->GetPDGEncoding() );
      }
   }

   fTree->Fill();

   //if(evtN%10 == 0 )
   //// Increase event number. Notice: this is different than evt->GetEventID()
   //evtN++;
   //return;

   if(event_number%1000 == 0 ) {
      std::cout <<  " End of Event " << event_number << " " << std::endl;
   }

   // Increase event number. Notice: this is different than evt->GetEventID()
   event_number++;
  // Called after G4Run::RecordEvent

  // accumulate statistics in B1Run
  //B1Run* run = static_cast<B1Run*>( G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  //run->AddEdep(fEdep);
  //if( run->GetNumberOfEvent()%1000 == 0 ) std::cout << "EndOfEventAction\n";
}
//______________________________________________________________________________

