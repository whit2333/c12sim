#ifndef BeamOnTargetPrimaryGeneratorAction_h
#define BeamOnTargetPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <fstream>
#include "ThrownEvent.h"

class G4ParticleGun;
class G4Event;
class G4Box;
class G4ParticleTable;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued 
/// in front of the phantom across 80% of the (X,Y) phantom size.

class BeamOnTargetPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  private:
    G4ParticleGun *  fParticleGun; // pointer a to G4 gun class
    G4Box         *  fEnvelopeBox;

    std::ifstream    fInputLundFile;
    clas12::sim::ThrownEvent fThrownEvent;
    G4ParticleTable* particleTable;
  public:
    BeamOnTargetPrimaryGeneratorAction();    
    virtual ~BeamOnTargetPrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         
  
    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
};

#endif


