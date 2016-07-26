#ifndef SiPMSD_h
#define SiPMSD_h 1

#include "G4VSensitiveDetector.hh"
#include "SiPMHit.hh"
#include "RecoilScintEvent.h"
#include <array>
#include <map>

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

      std::map<int,clas12::hits::PhotonCounterHit>  * fPhotonCounterHits = nullptr;

   public:
      clas12::hits::RecoilScintEvent * fRHEvent = nullptr;

      G4ParticleDefinition * fOpticalPhoton = nullptr;
      G4ParticleDefinition * fNeutron       = nullptr;

      G4int                fNeutronCount     = 0;
      G4double             fNeutronEnergyDep = 0.0; // Accumulated E-dep from neutrons
      G4double             fNeutronEnergy    = 0.0; // Accumulated neutron energy flux

   public:
      SiPMSD(G4String name);
      virtual ~SiPMSD();

      void SetGroup(int i);

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

