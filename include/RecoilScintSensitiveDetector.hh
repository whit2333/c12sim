#ifndef RecoilScintSensitiveDetector_h
#define RecoilScintSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "RecoilScintHit.hh"
#include <array>
#include "RecoilScintEvent.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class G4ParticleDefinition;

class RecoilScintSensitiveDetector : public G4VSensitiveDetector
{
   private:
      G4int HCID;
      RecoilScintHitsCollection *hitsCollection;

      std::array<RecoilScintHitsCollection*,5> theHits;
      std::array<int,5>                     theHCIDs;

   public:

      clas12::hits::RecoilScintEvent * fRecoilScintEvent = nullptr;


      G4ParticleDefinition * fOpticalPhoton;

   public:
      RecoilScintSensitiveDetector(G4String name);
      ~RecoilScintSensitiveDetector();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

};

#endif

