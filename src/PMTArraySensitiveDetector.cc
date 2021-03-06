#include "PMTArraySensitiveDetector.h"
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
#include "G4RunManager.hh"

#include "CLAS12DAQ.h"
#include "DAQManager.h"
#include "Scaler.h"
#include "SimulationManager.h"

//G4Allocator<MyHit> MyHitAllocator;

//______________________________________________________________________________
PMTArraySensitiveDetector::PMTArraySensitiveDetector(G4String name, G4int Nchan)
   : G4VSensitiveDetector(name), fNumberOfChannels(Nchan),
   fCountAllPhotons(true), fSavePhotonPositions(false)
{
   using namespace clas12::DAQ;
   DAQManager * daq_manager = DAQManager::GetManager();

   fCrate = BuildCrate();
   daq_manager->AddCrate(fCrate);

   fDiscModule       = fCrate->GetCrateModule<Discriminator>(0);
   if( !fDiscModule ) {
      std::cout << "Error : No DiscModule \n";
   }

   fTDCModule        = fCrate->GetCrateModule<TDC>(1);
   if( !fTDCModule ) {
      std::cout << "Error : No fTDCModule \n";
   }

   fScalerModule = fCrate->GetCrateModule<Scaler>(2);
   if( !fScalerModule ) {
      std::cout << "Error : No fTDCModule \n";
   }

   fADCModule        = fCrate->GetCrateModule<ADC>(3);
   if( !fADCModule ) {
      std::cout << "Error : No fADCModule \n";
   }

   fTrigDiscModule   = fCrate->GetCrateModule<Discriminator>(8);
   if( !fTrigDiscModule ) {
      std::cout << "Error : No fTrigDiscModule \n";
   }

   SimulationManager * sim_manager  = SimulationManager::GetInstance();
   fHTCCHitsEvent = &(sim_manager->fEvent->fHTCCEvent);

}
//______________________________________________________________________________
PMTArraySensitiveDetector::~PMTArraySensitiveDetector() {
}
//______________________________________________________________________________
void PMTArraySensitiveDetector::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent)
{
   fCrate->Clear();
   //fTDCModule->Clear();
   //fDiscModule->Clear();
   //fTrigDiscModule->Clear();
   //fTrigDiscModule->GetChannel(0).Print();
   //fTrigDiscModule->GetChannel(0).Clear();
   //fTrigDiscModule->GetChannel(0).Print();

}
//______________________________________________________________________________
G4bool PMTArraySensitiveDetector::ProcessHits ( G4Step* aStep, G4TouchableHistory* ) {

   G4Track * theTrack = aStep->GetTrack();

   if( theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() &&
        aStep->GetPreStepPoint()->GetStepStatus()== fGeomBoundary ) 
   {
      // The photon has just entered the boundary

      G4int   copyNo     = theTrack->GetVolume()->GetCopyNo();
      double  track_time = theTrack->GetGlobalTime()/CLHEP::ns;

      clas12::hits::TDCHit * tdc_hit = 0;

      if(copyNo < fNumberOfChannels)
      {
         //std::cout << track_time << "\n";

         bool disc_fired      = fDiscModule->GetChannel(copyNo).Count(track_time);
         bool trig_disc_fired = fTrigDiscModule->GetChannel(0).Count(track_time);

         bool latched  = fDiscModule->GetChannel(copyNo).fLatched;
         bool latched2 = fTrigDiscModule->GetChannel(0).fLatched;

         fADCModule->GetChannel(copyNo).Count();

         if( (trig_disc_fired && latched ) || (disc_fired && latched2 ) )
         {
            int tdc_value = fTDCModule->GetChannel(copyNo).Readout();
            tdc_hit = fHTCCHitsEvent->AddTDCHit(copyNo,tdc_value,track_time);
            //tdc_hit->Print();
         }
      }

      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
   }
   return true;
}
//______________________________________________________________________________
void PMTArraySensitiveDetector::EndOfEvent ( G4HCofThisEvent* ) {

   clas12::hits::ADCHit * adc_hit = 0;
   for(int i = 0; i<fNumberOfChannels; i++) {
      int adc_value = fADCModule->GetChannel(i).Readout();
      adc_hit = fHTCCHitsEvent->AddADCHit(i,adc_value);
      //adc_hit->Print();
   }

   //clas12::DAQ::DAQManager& daq_manager = clas12::DAQ::DAQManager::GetManager();
   //for(int i = 0; i< daq_manager.fScalers.size(); i++ ) {
   //   daq_manager.fScalers[i].Print();
   //}

}
//______________________________________________________________________________
G4double PMTArraySensitiveDetector::QE ( G4double photonEnergy ) const {
   return 0.20;
}
//______________________________________________________________________________

