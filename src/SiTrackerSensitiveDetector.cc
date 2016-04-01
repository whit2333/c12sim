#include "SiTrackerSensitiveDetector.h"
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
//#include "B1RunAction.h"
#include "G4RunManager.hh"

#include "CLAS12DAQ.h"
#include "DAQManager.h"
#include "Scaler.h"
#include "SimulationManager.h"
#include "SVTEvent.h"
#include "ParticleHit.h"


//______________________________________________________________________________

SiTrackerSensitiveDetector::SiTrackerSensitiveDetector(G4String name, G4int Nchan)
   : G4VSensitiveDetector(name), fNumberOfChannels(Nchan),
   fCountAllPhotons(true), fSavePhotonPositions(false)
{
   SimulationManager * sim_manager  = SimulationManager::GetInstance();
   fSVTEvent = &(sim_manager->fEvent->fSVTEvent);
}
//______________________________________________________________________________

SiTrackerSensitiveDetector::~SiTrackerSensitiveDetector()
{
}
//______________________________________________________________________________

void SiTrackerSensitiveDetector::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent)
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

G4bool SiTrackerSensitiveDetector::ProcessHits ( G4Step* aStep, G4TouchableHistory* )
{

   using namespace clas12::hits;
   using namespace CLHEP;
   //G4Track * theTrack = aStep->GetTrack();
   double step_length = aStep->GetStepLength()/cm;

   bool ion_pair   = false;//does_step_create_ion_pair( step_length );
   bool first_step = false; //aStep->IsFirstStepInVolume();
   if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary)  {
      // First step in volume
      // For some reason IsFirstStepInVolume doesn't work for parallel geometry
      first_step = true;
   }

   if( ion_pair || first_step ) {

      const G4ThreeVector& pos_global = aStep->GetPreStepPoint()->GetPosition();
      double time = aStep->GetPreStepPoint()->GetGlobalTime()/ns;

      G4TouchableHistory * touchable  = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
      G4ThreeVector pos = touchable->GetHistory()->GetTopTransform().TransformPoint(pos_global);


      // layer/superlayer/wire comes from channel
      int Ti   = touchable->GetReplicaNumber(0);
      int Zi   = touchable->GetReplicaNumber(1);
      int lay  = touchable->GetReplicaNumber(2);
      int lad  = touchable->GetReplicaNumber(3);
      int wire = 0;//touchable->GetReplicaNumber(4);
      int channel  = fSVTGeometry.GetChannelNumber(lad,Ti,Zi);

      //std::cout << "Ti  : " << Ti  << std::endl;
      ////std::cout << touchable->GetVolume(0)->GetName() << std::endl;
      //std::cout << "Zi  : " << Zi  << std::endl;
      ////std::cout << touchable->GetVolume(1)->GetName() << std::endl;
      //std::cout << "lad : " << lad << std::endl;
      ////std::cout << touchable->GetVolume(2)->GetName() << std::endl;
      //std::cout << "lay: " << lay << std::endl;
      ////std::cout << touchable->GetVolume(3)->GetName() << std::endl;
      //std::cout << "wire: " << wire << std::endl;

      // grouping is the placement of the sector/region
      int grouping    = 0;
      int sector      = 0;
      int region      = 0;
      int super_layer = 0;
      int pdg         = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
      TLorentzVector global_4vec(pos_global.x()/cm,pos_global.y()/cm,pos_global.z()/cm,time);


      //if(ion_pair) {
      //   DriftChamberIonPairHit * ahit = fSVTEvent->AddIonPairHit( pos.x()/cm, pos.y()/cm, pos.z()/cm, time );
      //   ahit->fGlobalPosition         = global_4vec;
      //   ahit->fStepLength             = step_length;
      //   ahit->fChannel                = channel;
      //   ahit->fDCWire.fSector         = sector;
      //   ahit->fDCWire.fRegion         = region;
      //   ahit->fDCWire.fSuperLayer     = super_layer;
      //   ahit->fDCWire.fLayer          = layer;
      //   ahit->fDCWire.fWire           = wire;
      //   ahit->fPDGCode                = pdg;
      //}

      if( first_step ) {

         G4ThreeVector  mom  = aStep->GetTrack()->GetMomentum();
         double         Etot = aStep->GetTrack()->GetTotalEnergy();
         if(Etot/MeV > 10.0) { 
            ParticleHit * part_hit = fSVTEvent->AddParticleHit();
            part_hit->fPDGCode            = pdg;
            part_hit->fChannel            = channel;
            part_hit->fPosition           = TLorentzVector(pos.x()/cm, pos.y()/cm, pos.z()/cm, time );
            part_hit->fGlobalPosition     = global_4vec;
            part_hit->fMomentum           = TLorentzVector(mom.x()/GeV, mom.y()/GeV, mom.z()/GeV, Etot/GeV);
         }
      }

   }
   return true;
}
//______________________________________________________________________________

void SiTrackerSensitiveDetector::EndOfEvent ( G4HCofThisEvent* ) {

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
G4double SiTrackerSensitiveDetector::QE ( G4double photonEnergy ) const {
   return 0.20;
}
//______________________________________________________________________________

