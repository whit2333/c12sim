#ifndef HTCCSensitiveDetector_HH
#define HTCCSensitiveDetector_HH 1

#include "G4VSensitiveDetector.hh"
#include "TClonesArray.h"

#include "ADC.h"
#include "TDC.h"
#include "Discriminator.h"
#include "Scaler.h"
#include "ADCHit.h"
#include "TDCHit.h"

#include "TTree.h"
//#include "G4VHit.hh"
//#include "G4THitsCollection.hh"
//#include "G4Allocator.hh"
#include "HTCCHitsEvent.h"
#include "Crate.h"
#include "Module.h"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

//class MyHit : public G4VHit {
//   private:
//      int fChannel;
//
//   public:
//      MyHit(int i) : fChannel(i) { }
//      virtual ~MyHit(){ }
//
//      virtual void Print() { }
//      virtual void Draw() { }
//
//      inline void *operator new(size_t);
//      inline void operator delete(void *aHit);
//
//};
//
//typedef G4THitsCollection<MyHit> MyHitsCollection;
//
//extern G4Allocator<MyHit> MyHitAllocator;
//
//inline void* MyHit::operator new(size_t)
//{
//  void* aHit;
//  aHit = (void*)MyHitAllocator.MallocSingle();
//  return aHit;
//}
//
//inline void MyHit::operator delete(void* aHit)
//{
//  MyHitAllocator.FreeSingle((MyHit*) aHit);
//}


/** A concrete G4VSensitiveDetector class for an array of photomultiplier tubes
 * 
 *  The number of Channels should be set only once.
 *
 * \ingroup Detectors
 */
class HTCCSensitiveDetector : public G4VSensitiveDetector {

   protected:

      int    fNumberOfChannels;
      bool   fCountAllPhotons;
      bool   fSavePhotonPositions;

      //G4LogicalVolume * fPMTLogicalVolume;

      //G4int                      fHCID;
      //BETAG4PMTHitsCollection  * fHitsCollection;

   public:

      clas12::hits::HTCCHitsEvent    * fHTCCHitsEvent;
      clas12::DAQ::Crate             * fCrate;

      clas12::DAQ::Module<clas12::DAQ::Discriminator>  *  fDiscModule;
      clas12::DAQ::Module<clas12::DAQ::Discriminator>  *  fTrigDiscModule;
      clas12::DAQ::Module<clas12::DAQ::TDC>            *  fTDCModule;
      clas12::DAQ::Module<clas12::DAQ::ADC>            *  fADCModule;
      clas12::DAQ::Module<clas12::DAQ::Scaler>         *  fScalerModule;
      //clas12::DAQ::Module<FlashADC>            fFlashADCModule;

      TTree * fTree;
      //G4int collectionID;
      //MyHitsCollection * hitsCollection;

      HTCCSensitiveDetector( G4String name, G4int Nchan = 12);
      virtual ~HTCCSensitiveDetector();

      virtual void SetNumberOfChannels(int n) { fNumberOfChannels = n; };
      virtual int  GetNumberOfChannels() const { return(fNumberOfChannels); };

      void    Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);
      G4bool  ProcessHits(G4Step* aStep,G4TouchableHistory* history);
      void    EndOfEvent(G4HCofThisEvent* hitsCollectionOfThisEvent);


      G4double QE(G4double photonEnergy) const ;

};

#endif

