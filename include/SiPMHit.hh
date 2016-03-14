#ifndef SiPMHit_h
#define SiPMHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "B1Analysis.hh"

class SiPMHit : public G4VHit {

   public:

      G4double               fLambda;
      G4double               fTime;
      G4ThreeVector          fPosition;
      G4ThreeVector          fMomentum;
      G4int                  fChannel;
      G4int                  stripNo;
      G4double               energy;
      G4ParticleDefinition * particle;

   public:

      SiPMHit();
      virtual ~SiPMHit();
      G4int operator==(const SiPMHit &right) const;

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

typedef G4THitsCollection<SiPMHit> SiPMHitsCollection;

extern G4Allocator<SiPMHit> SiPMHitAllocator;

inline void* SiPMHit::operator new(size_t)
{
   void *aHit;
   aHit = (void *) SiPMHitAllocator.MallocSingle();
   return aHit;
}

inline void SiPMHit::operator delete(void *aHit)
{
   SiPMHitAllocator.FreeSingle((SiPMHit*) aHit);
}

#endif

