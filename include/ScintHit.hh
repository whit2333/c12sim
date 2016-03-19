#ifndef ScintHit_h
#define ScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "B1Analysis.hh"
#include "ScintChannelHit.h"

class ScintHit : public G4VHit {

   public:

      clas12::hits::ScintChannelHit   fChannelHit;

   public:

      ScintHit();
      virtual ~ScintHit();
      G4int operator==(const ScintHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();


   public:
      inline void SetChannel(G4int chan)
      { fChannelHit.fChannel=chan; }
      inline G4int GetChannel()
      { return fChannelHit.fChannel; }
      inline void AddEnergy(G4double ene)
      { fChannelHit.fEDep += ene; }
      inline G4double GetEnergyDeposited()
      { return fChannelHit.fEDep; }

};

typedef G4THitsMap<ScintHit> ScintHitsCollection;

extern G4Allocator<ScintHit> ScintHitAllocator;

inline void* ScintHit::operator new(size_t)
{
   void *aHit;
   aHit = (void *) ScintHitAllocator.MallocSingle();
   return aHit;
}

inline void ScintHit::operator delete(void *aHit)
{
   ScintHitAllocator.FreeSingle((ScintHit*) aHit);
}

#endif

