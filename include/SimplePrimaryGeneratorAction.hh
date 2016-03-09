#ifndef SimplePrimaryGeneratorAction_h
#define SimplePrimaryGeneratorAction_h 1

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

class SimplePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  private:
    G4ParticleGun *  fParticleGun; // pointer a to G4 gun class
    G4Box         *  fEnvelopeBox;

    double        fP_min     = 0.001*CLHEP::MeV;
    double        fP_max     = 1.000*CLHEP::GeV;
    double        fTheta_min = 0.0;
    double        fTheta_max = 100.0*CLHEP::degree;
    double        fPhi_min   =-180.0*CLHEP::degree;
    double        fPhi_max   = 180.0*CLHEP::degree;

    double        fdP        = 1.0;
    double        fdTheta    = 1.0;
    double        fdPhi      = 1.0;

    double        fDelta_z   = 5.0*CLHEP::cm;

    clas12::sim::ThrownEvent fThrownEvent;
    G4ParticleTable* particleTable;

  public:
    SimplePrimaryGeneratorAction();    
    virtual ~SimplePrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         
  
    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
};

#endif


