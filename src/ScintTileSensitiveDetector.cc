#include "ScintTileSensitiveDetector.hh"
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


ScintTileSensitiveDetector::ScintTileSensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
   G4String HCname;
   collectionName.insert( HCname = "particleHits" );
   //collectionName.insert( HCname = "channelHits" );
   HCID = -1;

   fOpticalPhoton = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

   //SimulationManager * sim_manager  = SimulationManager::GetInstance();
   //fRecoilScintEvent = &(sim_manager->fEvent->fRecoilScintEvent);
}
//______________________________________________________________________________

ScintTileSensitiveDetector::~ScintTileSensitiveDetector()
{}
//______________________________________________________________________________

void ScintTileSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
   if(!fRecoilScintEvent) {
      SimulationManager * simManager = SimulationManager::GetInstance();
      fRecoilScintEvent    = new clas12::hits::RecoilScintEvent();
      simManager->fOutputTree->Branch(
            Form("%s_%s",SensitiveDetectorName.data(),collectionName[0].data()),
            "clas12::hits::RecoilScintEvent",
            &fRecoilScintEvent   );
   }
   fRecoilScintEvent->Clear();

   hitsCollection = new RecoilScintHitsCollection( SensitiveDetectorName, collectionName[0] ); 
   if(HCID<0) {
      HCID = GetCollectionID(0);
   }
   HCE->AddHitsCollection(HCID,hitsCollection);

   //fHitsByChannel = new ScintHitsCollection( SensitiveDetectorName, collectionName[1] ); 
   //if(fHCIDByChannel<0) {
   //   fHCIDByChannel = GetCollectionID(1);
   //}
   //HCE->AddHitsCollection(fHCIDByChannel,fHitsByChannel);

}
//______________________________________________________________________________

G4bool ScintTileSensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
   using namespace CLHEP;
   //G4StepPoint* preStep = aStep->GetPreStepPoint();
   //G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());

   G4StepPoint        * preStep   = aStep->GetPreStepPoint();
   G4TouchableHistory * touchable = (G4TouchableHistory*)(preStep->GetTouchable());
   int pdgcode = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

   // check if the particle type is Optical Photon
   bool is_OpticalPhoton = (aStep->GetTrack()->GetDefinition() == fOpticalPhoton);

   if(!is_OpticalPhoton  ) 
   {

      int     channel  = touchable->GetReplicaNumber();
      double  e_dep    = aStep->GetTotalEnergyDeposit()/GeV;
      fRecoilScintEvent->fScintChannelHits[channel].fChannel = channel; // derp
      fRecoilScintEvent->fScintChannelHits[channel].fSteps++;
      fRecoilScintEvent->fScintChannelHits[channel].fEDep += e_dep;

      if( preStep->GetStepStatus() == fGeomBoundary ) {
         G4ThreeVector mom        = preStep->GetMomentum();
         G4ThreeVector pos_global = preStep->GetPosition();

         G4TouchableHistory * touchable  = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
         G4ThreeVector pos = touchable->GetHistory()->GetTopTransform().TransformPoint(pos_global);

         double        total_energy = aStep->GetPreStepPoint()->GetTotalEnergy()/GeV;
         double        aTime        = aStep->GetTrack()->GetGlobalTime()/CLHEP::ns;
         TLorentzVector global_4vec(pos_global.x()/cm,pos_global.y()/cm,pos_global.z()/cm,aTime);

         clas12::hits::ParticleHit * pHit = fRecoilScintEvent->AddParticleHit(channel);
         pHit->fPDGCode            = pdgcode;
         pHit->fPosition           = TLorentzVector(pos.x()/cm, pos.y()/cm, pos.z()/cm, aTime );
         pHit->fGlobalPosition     = global_4vec;
         pHit->fMomentum           = TLorentzVector(mom.x()/GeV, mom.y()/GeV, mom.z()/GeV, total_energy);
      }

   }

   //if( (pdgcode == 0) && is_OpticalPhoton ) {
   //   if( preStep->GetStepStatus() == fGeomBoundary ) {

   //      //std::cout << " pdgcode " << pdgcode << " \n";
   //      int event_number = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent();
   //      //int run_number   = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
   //      int run_number   = SimulationManager::GetInstance()->GetRunNumber();
   //      int channel      = touchable->GetReplicaNumber(1);

   //      G4ThreeVector momentum     = aStep->GetPreStepPoint()->GetMomentum();
   //      G4ThreeVector pos          = aStep->GetPreStepPoint()->GetPosition();

   //      double        total_energy   = aStep->GetPreStepPoint()->GetTotalEnergy()/MeV;
   //      double        kinetic_energy = aStep->GetPreStepPoint()->GetKineticEnergy()/MeV;

   //      double        lambda       = (hbarc*twopi/(total_energy*MeV))/nm;

   //      double        time         = aStep->GetTrack()->GetGlobalTime()/CLHEP::ns;

   //      double        pz           = momentum.z()/MeV;
   //      double        px           = momentum.x()/MeV;
   //      double        py           = momentum.y()/MeV;

   //      double        p            = momentum.mag()/MeV;
   //      double        theta_p      = momentum.theta();
   //      double        phi_p        = momentum.phi();

   //      double        theta        = pos.theta();
   //      double        phi          = pos.phi();

   //      double        z           = pos.z()/cm;
   //      double        x           = pos.x()/cm;
   //      double        y           = pos.y()/cm;
   //      // -------------------------------------------------
   //      RecoilScintHit* newHit = new RecoilScintHit();
   //      newHit->SetStripNo(  touchable->GetReplicaNumber(0) );
   //      newHit->SetPosition( pos );
   //      newHit->SetMomentum( momentum );
   //      newHit->SetEnergy(   aStep->GetPreStepPoint()->GetTotalEnergy() );
   //      newHit->SetParticle( aStep->GetTrack()->GetDefinition() );
   //      newHit->fLambda = lambda;
   //      newHit->fTime   = time;

   //      theHits[channel]->insert( newHit );
   //      //hitsCollection->insert( newHit );

   //      clas12::hits::PhotonHit * pHit = fRecoilScintEvent->AddPhotonHit(channel);
   //      pHit->fChannel  = channel;
   //      pHit->fTime     = time;
   //      pHit->fLambda   = lambda;
   //      pHit->fEnergy   = total_energy;
   //      pHit->fPosition = {pos.x(),pos.y(),pos.z()};
   //      pHit->fMomentum = {px,py,pz};
   //      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
   //   }
   //}
   return true;
}
//______________________________________________________________________________

void ScintTileSensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
   int event_number = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent();
   //int run_number   = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
   int run_number   = SimulationManager::GetInstance()->GetRunNumber();

   fRecoilScintEvent->fEventNumber = event_number;
   fRecoilScintEvent->fRunNumber   = run_number;

   //for( int i_col = 0; i_col<5 ; i_col++ ) {

   //   double average_lambda = 0.0;
   //   double average_time   = 0.0;


   //   int nhits = theHits[i_col]->GetSize();
   //   //int nhits = hitsCollection->GetSize();
   //   for( int i = 0; i < nhits; i++ ) {
   //      RecoilScintHit * aHit = (*(theHits[i_col]))[i];
   //      average_time   += aHit->fTime;
   //      average_lambda += aHit->fLambda;
   //   }

   //   if( nhits > 0 ) {
   //      average_time   *= (1.0/double(nhits));
   //      average_lambda *= (1.0/double(nhits));
   //   }

   //   //man->FillH1( fTimeAvg_hID,   average_time   );
   //   //man->FillH1( fLambdaAvg_hID, average_lambda );
   //   //man->FillH2( fTimeVsLambdaAvg_hID, average_lambda, average_time );

   //   if( nhits > 0 ) {
   //      clas12::hits::RecoilScintHit * rHit = fRecoilScintEvent->AddHit(i_col);
   //      rHit->fTime   = average_time;
   //      rHit->fLambda = average_lambda;
   //      rHit->fChannel = i_col;
   //   }
   //}
   
}
//______________________________________________________________________________

void ScintTileSensitiveDetector::clear()
{;} 
//______________________________________________________________________________

void ScintTileSensitiveDetector::DrawAll()
{;} 
//______________________________________________________________________________

void ScintTileSensitiveDetector::PrintAll()
{;} 
//______________________________________________________________________________

