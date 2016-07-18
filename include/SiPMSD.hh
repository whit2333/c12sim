#ifndef SiPMSD_h
#define SiPMSD_h 1

#include "G4VSensitiveDetector.hh"
#include "SiPMHit.hh"
#include <array>
#include "RecoilScintEvent.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class G4ParticleDefinition;

class SiPMSD : public G4VSensitiveDetector
{
   private:

      G4int                fCopyNoParent     = 0;
      G4int                fNChannels        = 600;
      G4int                fHCID             = -1;
      SiPMHitsCollection * fHitsCollection   = nullptr;
      bool                 fRecordAllPhotons = true;

      std::map<int,clas12::hits::PhotonCounterHit>  * fPhotonCounterHits0 = nullptr;
      //std::map<clas12::hits::PhotonCounterHit>  * fPhotonCounterHits1 = nullptr;
      //std::map<clas12::hits::PhotonCounterHit>  * fPhotonCounterHitsTile = nullptr;

   public:
      clas12::hits::RecoilScintEvent * fRHEvent = nullptr;

      G4ParticleDefinition * fOpticalPhoton = nullptr;

   public:
      SiPMSD(G4String name);
      virtual ~SiPMSD();
      void SetGroupNumber(int i);

      void  SetCopyNoParent(int i) { fCopyNoParent = i; } 
      G4int GetCopyNoParent() const {return fCopyNoParent; }

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();


};




#endif

