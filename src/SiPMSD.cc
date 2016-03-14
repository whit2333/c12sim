#include "SiPMSD.hh"
#include "SiPMHit.hh"
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
#include <map>


SiPMSD::SiPMSD(G4String name) : G4VSensitiveDetector(name)
{
   G4String HCname;
   collectionName.insert( HCname = "SiPM_Hits" );
   fHCID = -1;

   fOpticalPhoton = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
}
//______________________________________________________________________________

SiPMSD::~SiPMSD()
{}
//______________________________________________________________________________

void SiPMSD::Initialize(G4HCofThisEvent* HCE)
{
   if(!fRecoilScintEvent) {
      SimulationManager * simManager = SimulationManager::GetInstance();
      fRecoilScintEvent    = new clas12::hits::RecoilScintEvent();
      simManager->fOutputTree->Branch(
            Form("scint_%s_%s",SensitiveDetectorName.data(),collectionName[0].data()),
            "clas12::hits::RecoilScintEvent",
            &fRecoilScintEvent   );
   }
   fRecoilScintEvent->Clear();

   fHitsCollection = new SiPMHitsCollection( SensitiveDetectorName, collectionName[0] ); 
   if(fHCID<0) {
      fHCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(fHCID,fHitsCollection);
}
//______________________________________________________________________________

G4bool SiPMSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

   using namespace CLHEP;
   //G4StepPoint* preStep = aStep->GetPreStepPoint();
   //G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());

   G4StepPoint        * preStep   = aStep->GetPreStepPoint();
   G4TouchableHistory * touchable = (G4TouchableHistory*)(preStep->GetTouchable());
   int pdgcode = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

   // check if the particle type is Optical Photon
   bool is_OpticalPhoton = (aStep->GetTrack()->GetDefinition() == fOpticalPhoton);

   if( (pdgcode == 0) && is_OpticalPhoton ) {
      if( preStep->GetStepStatus() == fGeomBoundary ) {

         //std::cout << " pdgcode " << pdgcode << " \n";
         int event_number = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent();
         //int run_number   = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
         int run_number   = SimulationManager::GetInstance()->GetRunNumber();
         int channel      = touchable->GetReplicaNumber(GetCopyNoParent());

         G4ThreeVector momentum     = aStep->GetPreStepPoint()->GetMomentum();
         G4ThreeVector pos          = aStep->GetPreStepPoint()->GetPosition();

         double        total_energy   = aStep->GetPreStepPoint()->GetTotalEnergy()/MeV;
         double        kinetic_energy = aStep->GetPreStepPoint()->GetKineticEnergy()/MeV;

         double        lambda       = (hbarc*twopi/(total_energy*MeV))/nm;

         double        time         = aStep->GetTrack()->GetGlobalTime()/CLHEP::ns;

         double        pz           = momentum.z()/MeV;
         double        px           = momentum.x()/MeV;
         double        py           = momentum.y()/MeV;

         double        p            = momentum.mag()/MeV;
         double        theta_p      = momentum.theta();
         double        phi_p        = momentum.phi();

         double        theta        = pos.theta();
         double        phi          = pos.phi();

         double        z           = pos.z()/cm;
         double        x           = pos.x()/cm;
         double        y           = pos.y()/cm;

         //std::cout << "channel(" << GetCopyNoParent() << ") = " << channel << std::endl ;
         //std::cout << "touchable->GetReplicaNumber(2)  = " << touchable->GetReplicaNumber(2) << std::endl ;
         //std::cout << "touchable->GetReplicaNumber(4)  = " << touchable->GetReplicaNumber(4) << std::endl ;

         // -------------------------------------------------

         SiPMHit* newHit = new SiPMHit();
         newHit->fChannel =   channel ;
         newHit->SetStripNo(  channel );
         newHit->SetPosition( pos );
         newHit->SetMomentum( momentum );
         newHit->SetEnergy(   aStep->GetPreStepPoint()->GetTotalEnergy() );
         newHit->SetParticle( aStep->GetTrack()->GetDefinition() );
         newHit->fLambda = lambda;
         newHit->fTime   = time;

         fHitsCollection->insert( newHit );
         //hitsCollection->insert( newHit );

         clas12::hits::PhotonHit * pHit = fRecoilScintEvent->AddPhotonHit(channel);
         pHit->fChannel  = channel;
         pHit->fTime     = time;
         pHit->fLambda   = lambda;
         pHit->fEnergy   = total_energy;
         pHit->fPosition = {pos.x(),pos.y(),pos.z()};
         pHit->fMomentum = {px,py,pz};

         // -------------------------------------------------

         aStep->GetTrack()->SetTrackStatus(fStopAndKill);
      }
   }
   return true;
}
//______________________________________________________________________________

void SiPMSD::EndOfEvent(G4HCofThisEvent*)
{
   int event_number = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent();
   //int run_number   = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
   int run_number   = SimulationManager::GetInstance()->GetRunNumber();

   fRecoilScintEvent->fEventNumber = event_number;
   fRecoilScintEvent->fRunNumber   = run_number;

   std::map<int,double> average_lambda ;
   std::map<int,double> average_time   ;
   std::map<int,int>    channel_count ;

   int nhits = fHitsCollection->GetSize();
   for( int i = 0; i < nhits; i++ ) {
      SiPMHit * aHit = (*fHitsCollection)[i];
      if(average_time.count(aHit->fChannel) == 0 ) {
         average_time[aHit->fChannel]   =  0.0;
         average_lambda[aHit->fChannel] =  0.0;
         channel_count[aHit->fChannel]  =  0;
      }

      average_time[aHit->fChannel]   += aHit->fTime;
      average_lambda[aHit->fChannel] += aHit->fLambda;
      channel_count[aHit->fChannel]++  ;
   }

   for(auto chan : average_time ) {
      average_time[chan.first]   *= (1.0/double(channel_count[chan.first]));
      average_lambda[chan.first] *= (1.0/double(channel_count[chan.first]));

      clas12::hits::RecoilScintHit * rHit = fRecoilScintEvent->AddHit(chan.first);
      rHit->fTime   = average_time[chan.first]  ;
      rHit->fLambda = average_lambda[chan.first];
      rHit->fChannel = chan.first;
   }

}
//______________________________________________________________________________

void SiPMSD::clear()
{;} 
//______________________________________________________________________________

void SiPMSD::DrawAll()
{;} 
//______________________________________________________________________________

void SiPMSD::PrintAll()
{;} 
//______________________________________________________________________________

