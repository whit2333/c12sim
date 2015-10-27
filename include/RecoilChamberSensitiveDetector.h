#ifndef RecoilChamberSensitiveDetector_HH
#define RecoilChamberSensitiveDetector_HH 1

#include "G4VSensitiveDetector.hh"
#include "TClonesArray.h"

#include "ADC.h"
#include "TDC.h"
#include "Discriminator.h"
#include "Scaler.h"
#include "ADCHit.h"
#include "TDCHit.h"

#include "TTree.h"
#include "HTCCHitsEvent.h"
#include "RCHitsEvent.h"
#include "Crate.h"
#include "Module.h"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;


class RecoilChamberSensitiveDetector : public G4VSensitiveDetector {

   protected:

      int    fNumberOfChannels;
      bool   fCountAllPhotons;
      bool   fSavePhotonPositions;


   public:

      clas12::hits::RCHitsEvent    * fRCHitsEvent;

      RecoilChamberSensitiveDetector( G4String name, G4int Nchan = 12);
      virtual ~RecoilChamberSensitiveDetector();

      virtual void SetNumberOfChannels(int n) { fNumberOfChannels = n; };
      virtual int  GetNumberOfChannels() const { return(fNumberOfChannels); };

      void    Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);
      G4bool  ProcessHits(G4Step* aStep,G4TouchableHistory* history);
      void    EndOfEvent(G4HCofThisEvent* hitsCollectionOfThisEvent);


      G4double QE(G4double photonEnergy) const ;

};

#endif

