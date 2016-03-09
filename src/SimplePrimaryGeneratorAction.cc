#include "SimplePrimaryGeneratorAction.hh"

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


SimplePrimaryGeneratorAction::SimplePrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0), fEnvelopeBox(0)
{
   G4int n_particle = 1;

   fThrownEvent.Clear();
   fThrownEvent.AddParticle();
   TParticle * part = fThrownEvent.GetParticle(0);
   int pdgcode = 2212;//
   part->SetPdgCode(pdgcode);

   fParticleGun  = new G4ParticleGun(n_particle);
   // default particle kinematic
   particleTable = G4ParticleTable::GetParticleTable();
   G4String  particleName;

   //G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
   G4ParticleDefinition* particle = particleTable->FindParticle(pdgcode);
   fParticleGun->SetParticleDefinition(particle);
   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
   fParticleGun->SetParticleEnergy(11.0*GeV); // kinetic energy (not total)

   fdP     = fP_max-fP_min;
   fdTheta = fTheta_max-fTheta_min;
   fdPhi   = fPhi_max-fPhi_min;
}
//______________________________________________________________________________

SimplePrimaryGeneratorAction::~SimplePrimaryGeneratorAction()
{
   delete fParticleGun;
}
//______________________________________________________________________________

void SimplePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   fdP     = fP_max-fP_min;
   fdTheta = fTheta_max-fTheta_min;
   fdPhi   = fPhi_max-fPhi_min;

   int npart = fThrownEvent.GetNParticles();
   for(int i = 0; i<npart; i++) {
      TParticle * part = fThrownEvent.GetParticle(i);
      int pdgcode = part->GetPdgCode();
      G4ParticleDefinition* particle = particleTable->FindParticle(pdgcode);
      fParticleGun->SetParticleDefinition(particle);

      double p_rand   = fP_min     + fdP*G4UniformRand();
      double th_rand  = fTheta_min + fdTheta*G4UniformRand();
      double phi_rand = fPhi_min   + fdPhi*G4UniformRand();

      double mass     = part->GetMass();
      double E_tot    = TMath::Sqrt(p_rand*p_rand + mass*mass );

      TLorentzVector P4vec = {1,0,0,1}; 
      P4vec.SetRho(p_rand);
      P4vec.SetTheta(th_rand);
      P4vec.SetPhi(phi_rand);
      P4vec.SetE(E_tot);
      part->SetMomentum(P4vec);
      part->SetProductionVertex(TLorentzVector(0.0,0.0,fDelta_z*(-0.5 + G4UniformRand()), 0.0));

      double KE = part->Energy() - mass;
      double ux = part->Px();
      double uy = part->Py();
      double uz = part->Pz();
      double x0 = part->Vx()*cm;
      double y0 = part->Vy()*cm;
      double z0 = part->Vz()*cm;
      fParticleGun->SetParticleEnergy( KE*GeV );
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
      fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
      fParticleGun->GeneratePrimaryVertex(anEvent);
   } 
   //fThrownEvent.Print();

   //} else {
   //   // Beam on target  

   //   G4double envSizeXY = 2.0*mm;
   //   G4double envSizeZ  = 0.0*cm;

   //   //double phi    = 2.*CLHEP::pi* G4UniformRand();
   //   //double theta  = 1.0*deg+60.0*deg*G4UniformRand();
   //   ////double cosTheta = -1. + 2. * G4UniformRand();
   //   //double cosTheta = TMath::Cos(theta);//1. * G4UniformRand();
   //   //double sinTheta = TMath::Sin(theta);//sqrt(1. - cosTheta * cosTheta);
   //   //double ux= sinTheta * cos(phi);
   //   //double uy= sinTheta * sin(phi);
   //   //double uz =cosTheta;

   //   //if( G4UniformRand() < 1.0 ) {
   //   //   G4ParticleDefinition* particle = particleTable->FindParticle("proton");
   //   //   fParticleGun->SetParticleDefinition(particle);
   //   //} else {
   //   //   G4ParticleDefinition* particle = particleTable->FindParticle("e-");
   //   //   fParticleGun->SetParticleDefinition(particle);
   //   //}

   //   //fParticleGun->SetParticleEnergy((3.5 + 0.1*G4UniformRand())*GeV); // kinetic energy (not total)
   //   //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));

   //   G4double size = 1.0; 
   //   G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
   //   G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
   //   G4double z0 = -10.0*cm;

   //   fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
   //   fParticleGun->GeneratePrimaryVertex(anEvent);
   //}
}
//______________________________________________________________________________

