#include "BeamOnTargetPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "SimulationManager.h"
#include "ThrownEvent.h"
#include "TParticle.h"


BeamOnTargetPrimaryGeneratorAction::BeamOnTargetPrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0), fEnvelopeBox(0)
{
   G4int n_particle = 1;
   fParticleGun  = new G4ParticleGun(n_particle);

   // default particle kinematic
   particleTable = G4ParticleTable::GetParticleTable();
   G4String  particleName;
   G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
   fParticleGun->SetParticleDefinition(particle);
   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
   fParticleGun->SetParticleEnergy(11.0*GeV); // kinetic energy (not total)

   //std::cout << " LUND FILE " << SimulationManager::GetInstance()->InputFileName() << "\n";
   fInputLundFile.open(SimulationManager::GetInstance()->InputFileName().c_str());
}
//______________________________________________________________________________

BeamOnTargetPrimaryGeneratorAction::~BeamOnTargetPrimaryGeneratorAction()
{
   delete fParticleGun;
}
//______________________________________________________________________________

void BeamOnTargetPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   // Beam on target  

   G4double envSizeXY = 2.0*mm;
   G4double envSizeZ  = 0.0*cm;

   //double phi    = 2.*CLHEP::pi* G4UniformRand();
   //double theta  = 1.0*deg+60.0*deg*G4UniformRand();
   ////double cosTheta = -1. + 2. * G4UniformRand();
   //double cosTheta = TMath::Cos(theta);//1. * G4UniformRand();
   //double sinTheta = TMath::Sin(theta);//sqrt(1. - cosTheta * cosTheta);
   //double ux= sinTheta * cos(phi);
   //double uy= sinTheta * sin(phi);
   //double uz =cosTheta;

   //if( G4UniformRand() < 1.0 ) {
   //   G4ParticleDefinition* particle = particleTable->FindParticle("proton");
   //   fParticleGun->SetParticleDefinition(particle);
   //} else {
   //   G4ParticleDefinition* particle = particleTable->FindParticle("e-");
   //   fParticleGun->SetParticleDefinition(particle);
   //}

   //fParticleGun->SetParticleEnergy((3.5 + 0.1*G4UniformRand())*GeV); // kinetic energy (not total)
   //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));

   G4double size = 1.0; 
   G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
   G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
   G4double z0 = -10.0*cm;

   fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
   fParticleGun->GeneratePrimaryVertex(anEvent);
}
//______________________________________________________________________________

