#include "ECSensitiveDetector.h"
#include "G4LogicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4OpticalPhoton.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "globals.hh"
//#include "RunAction.h"
#include "G4RunManager.hh"

#include "CLAS12DAQ.h"
#include "DAQManager.h"
#include "Scaler.h"
#include "SimulationManager.h"
#include "ParticleHit.h"


//______________________________________________________________________________

ECSensitiveDetector::ECSensitiveDetector(G4String name, G4int Nchan)
   : G4VSensitiveDetector(name), fNumberOfChannels(Nchan),
   fCountAllPhotons(true), fSavePhotonPositions(false)
{
   SimulationManager * sim_manager  = SimulationManager::GetInstance();
   fECHitsEvent = &(sim_manager->fEvent->fECEvent);
}
//______________________________________________________________________________
ECSensitiveDetector::~ECSensitiveDetector()
{
}
//______________________________________________________________________________
void ECSensitiveDetector::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent)
{
   //fCrate->Clear();
   //fTDCModule->Clear();
   //fDiscModule->Clear();
   //fTrigDiscModule->Clear();
   //fTrigDiscModule->GetChannel(0).Print();
   //fTrigDiscModule->GetChannel(0).Clear();
   //fTrigDiscModule->GetChannel(0).Print();

}
//______________________________________________________________________________
G4bool ECSensitiveDetector::ProcessHits ( G4Step* aStep, G4TouchableHistory* ) {

   using namespace clas12::hits;
   using namespace CLHEP;

   bool first_step = aStep->IsFirstStepInVolume();

   if( first_step ) {

      const G4ThreeVector& pos_global = aStep->GetPreStepPoint()->GetPosition();
      double time = aStep->GetPreStepPoint()->GetGlobalTime()/ns;

      G4TouchableHistory * touchable  = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
      G4ThreeVector pos = touchable->GetHistory()->GetTopTransform().TransformPoint(pos_global);

      // grouping is the placement of the sector/region
      int grouping    = touchable->GetReplicaNumber(0);
      int sector      = grouping/3+1;
      int region      = grouping%3+1;
      int pdg         = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
      int trk_id      = aStep->GetTrack()->GetTrackID();
      TLorentzVector global_4vec(pos_global.x()/cm,pos_global.y()/cm,pos_global.z()/cm,time);

      if( first_step ) {

         G4ThreeVector  mom  = aStep->GetTrack()->GetMomentum();
         double         Etot = aStep->GetTrack()->GetTotalEnergy();
         if(Etot/MeV > 1.0) { 
            ParticleHit * part_hit = fECHitsEvent->AddRegionHit();
            part_hit->fPDGCode            = pdg;
            part_hit->fTrackID            = trk_id;
            part_hit->fPosition           = TLorentzVector(pos.x()/cm, pos.y()/cm, pos.z()/cm, time );
            part_hit->fGlobalPosition     = global_4vec;
            part_hit->fMomentum           = TLorentzVector(mom.x()/GeV, mom.y()/GeV, mom.z()/GeV, Etot/GeV);
            //part_hit->fDCWire.fSector     = sector;
            //part_hit->fDCWire.fRegion     = region;
         }
      }

   }
   return true;
}
//______________________________________________________________________________

void ECSensitiveDetector::EndOfEvent ( G4HCofThisEvent* ) {

   //clas12::hits::ADCHit * adc_hit = 0;
   //for(int i = 0; i<fNumberOfChannels; i++) {
   //   int adc_value = fADCModule->GetChannel(i).Readout();
   //   adc_hit = fHTCCHitsEvent->AddADCHit(i,adc_value);
   //   //adc_hit->Print();
   //}

   //clas12::DAQ::DAQManager& daq_manager = clas12::DAQ::DAQManager::GetManager();
   //for(int i = 0; i< daq_manager.fScalers.size(); i++ ) {
   //   daq_manager.fScalers[i].Print();
   //}

}
//______________________________________________________________________________
G4double ECSensitiveDetector::QE ( G4double photonEnergy ) const {
   return 0.20;
}
//______________________________________________________________________________

