#include "RecoilScintSensitiveDetector.hh"
#include "RecoilScintHit.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "SimulationManager.h"
#include "G4Run.hh"
#include "RecoilScintEvent.h"
#include "RecoilScintHit.h"
#include "PhotonHit.h"
#include "ParticleHit.h"

#include "dollar.hpp"


//______________________________________________________________________________

RecoilScintSensitiveDetector::RecoilScintSensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
   G4String HCname;
   collectionName.insert( HCname = "particleHits" );
   HCID = -1;

   fOpticalPhoton = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

   SimulationManager * sim_manager  = SimulationManager::GetInstance();
   fRecoilScintEvent = &(sim_manager->fEvent->fRHEvent);
}
//______________________________________________________________________________

RecoilScintSensitiveDetector::~RecoilScintSensitiveDetector()
{}
//______________________________________________________________________________

void RecoilScintSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
   fRecoilScintEvent->Clear();

   hitsCollection = new RecoilScintHitsCollection( SensitiveDetectorName, collectionName[0] ); 
   if(HCID<0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID,hitsCollection);

}
//______________________________________________________________________________

G4bool RecoilScintSensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{ $
   using namespace CLHEP;

   G4StepPoint        * preStep   = aStep->GetPreStepPoint();
   G4TouchableHistory * touchable = (G4TouchableHistory*)(preStep->GetTouchable());
   int pdgcode = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

   // check if the particle type is Optical Photon
   bool is_OpticalPhoton = (aStep->GetTrack()->GetDefinition() == fOpticalPhoton);

   if( preStep->GetStepStatus() == fGeomBoundary ) {

         int channel      = touchable->GetReplicaNumber();
         int trk_id      = aStep->GetTrack()->GetTrackID();

         G4ThreeVector mom     = preStep->GetMomentum();
         G4ThreeVector pos_global = preStep->GetPosition();

         G4TouchableHistory * touchable  = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
         G4ThreeVector pos = touchable->GetHistory()->GetTopTransform().TransformPoint(pos_global);

         double        total_energy = aStep->GetPreStepPoint()->GetTotalEnergy()/GeV;
         double        aTime        = aStep->GetTrack()->GetGlobalTime()/CLHEP::ns;
         TLorentzVector global_4vec(pos_global.x()/cm,pos_global.y()/cm,pos_global.z()/cm,aTime);

         clas12::hits::ParticleHit * pHit = fRecoilScintEvent->AddParticleHit(channel);

         pHit->fPDGCode            = pdgcode;
         pHit->fTrackID            = trk_id;
         pHit->fPosition           = TLorentzVector(pos.x()/cm, pos.y()/cm, pos.z()/cm, aTime );
         pHit->fGlobalPosition     = global_4vec;
         pHit->fMomentum           = TLorentzVector(mom.x()/GeV, mom.y()/GeV, mom.z()/GeV, total_energy);

   }

   if(!is_OpticalPhoton  ) {

      int     channel  = touchable->GetReplicaNumber();
      double  e_dep    = aStep->GetTotalEnergyDeposit()/GeV;
      //fRecoilScintEvent->fScintChannelHits[channel].fChannel = channel; // derp
      //fRecoilScintEvent->fScintChannelHits[channel].fSteps++;
      //fRecoilScintEvent->fScintChannelHits[channel].fEDep += e_dep;
   }

   return true;
}
//______________________________________________________________________________

void RecoilScintSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
   int event_number = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent();
   int run_number   = SimulationManager::GetInstance()->GetRunNumber();

   fRecoilScintEvent->fEventNumber = event_number;
   fRecoilScintEvent->fRunNumber   = run_number;
}
//______________________________________________________________________________

void RecoilScintSensitiveDetector::clear()
{;} 
//______________________________________________________________________________

void RecoilScintSensitiveDetector::DrawAll()
{;} 
//______________________________________________________________________________

void RecoilScintSensitiveDetector::PrintAll()
{;} 
//______________________________________________________________________________

