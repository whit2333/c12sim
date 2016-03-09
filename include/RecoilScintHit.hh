#ifndef RecoilScintHit_h
#define RecoilScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "B1Analysis.hh"

class RecoilScintHit : public G4VHit {

   public:

      G4double               fLambda;
      G4double               fTime;
      G4ThreeVector          fPosition;
      G4ThreeVector          fMomentum;
      G4int                  stripNo;
      G4double               energy;
      G4ParticleDefinition * particle;

   public:

      RecoilScintHit();
      virtual ~RecoilScintHit();
      G4int operator==(const RecoilScintHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();


   public:
      inline void SetStripNo(G4int strip)
      { stripNo=strip; }
      inline G4int GetStripNo()
      { return stripNo; }
      inline void SetPosition(G4ThreeVector pos)
      { fPosition=pos; }
      inline G4ThreeVector GetPosition()
      { return fPosition; }
      inline void SetMomentum(G4ThreeVector mom)
      { fMomentum = mom; }
      inline G4ThreeVector GetMomentum()
      { return fMomentum; }
      inline void SetEnergy(G4double ene)
      { energy = ene; }
      inline G4double GetEnergy()
      { return energy; }
      inline void SetParticle(G4ParticleDefinition* pdef)
      { particle = pdef; }
      inline G4ParticleDefinition* GetParticle()
      { return particle; }

};

typedef G4THitsCollection<RecoilScintHit> RecoilScintHitsCollection;

extern G4Allocator<RecoilScintHit> RecoilScintHitAllocator;

inline void* RecoilScintHit::operator new(size_t)
{
   void *aHit;
   aHit = (void *) RecoilScintHitAllocator.MallocSingle();
   return aHit;
}

inline void RecoilScintHit::operator delete(void *aHit)
{
   RecoilScintHitAllocator.FreeSingle((RecoilScintHit*) aHit);
}

#endif

