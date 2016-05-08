#include "B1PrimaryGeneratorAction.hh"

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

#include "dollar.hpp"

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction(int ev_start) : 
   G4VUserPrimaryGeneratorAction(),
   fParticleGun(0), fEnvelopeBox(0), fStartEvent(ev_start)
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
   if( !(fInputLundFile.good()) ) {
      std::cout << "Error: LUND file, " << SimulationManager::GetInstance()->InputFileName() << ", not found.\nAborting.\n";
      std::exit(EXIT_FAILURE);
   }

   // Move to starting event
   for(int i = 0; i<fStartEvent; i++){
      fThrownEvent.ReadLundEvent(fInputLundFile);
      if( fInputLundFile.eof() ) {
         std::cout << "Error: LUND file, " 
            << SimulationManager::GetInstance()->InputFileName() 
            << ", contains only " << i  
            << " events and starting event is set to " << fStartEvent << ".\nAborting.\n";
         std::exit(EXIT_FAILURE);
      }
   }
}
//______________________________________________________________________________

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
   delete fParticleGun;
}
//______________________________________________________________________________

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ $

   fThrownEvent.ReadLundEvent(fInputLundFile);
   if( fInputLundFile.eof() ) {
      std::cout << "Warning: Reached end of LUND file, " 
         << SimulationManager::GetInstance()->InputFileName() 
         << ". Rewinding to beginning.\n";
      fInputLundFile.clear();                 // clear fail and eof bits
      fInputLundFile.seekg(0, std::ios::beg); // back to the start!
      fThrownEvent.ReadLundEvent(fInputLundFile);
   }

   int npart = fThrownEvent.GetNParticles();
   //std::cout << " npart = " << npart << std::endl;
   for(int i = 0; i<npart; i++) {
      TParticle * part = fThrownEvent.GetParticle(i);
      int pdgcode = part->GetPdgCode();
      G4ParticleDefinition* particle = particleTable->FindParticle(pdgcode);
      fParticleGun->SetParticleDefinition(particle);
      double KE =  part->Energy() - part->GetMass();
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
   
}
//______________________________________________________________________________

