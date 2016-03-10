#include "SimplePrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "SimulationManager.h"
#include "ThrownEvent.h"
#include "TParticle.h"


SimplePrimaryGeneratorAction::SimplePrimaryGeneratorAction(int pdgcode) : 
   G4VUserPrimaryGeneratorAction()
{
   G4int n_particle = 1;

   fThrownEvent.Clear();
   fThrownEvent.AddParticle();
   TParticle * part = fThrownEvent.GetParticle(0);
   part->SetPdgCode(pdgcode);

   fParticleGun  = new G4GeneralParticleSource();
   // default particle kinematic
   particleTable = G4ParticleTable::GetParticleTable();
   G4String  particleName;

   //G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
   G4ParticleDefinition* particle = particleTable->FindParticle(pdgcode);
   fParticleGun->SetParticleDefinition(particle);
   //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
   //fParticleGun->SetParticleEnergy(11.0*GeV); // kinetic energy (not total)

   fdP     = fP_max-fP_min;
   fdTheta = fTheta_max-fTheta_min;
   fdPhi   = fPhi_max-fPhi_min;

   std::cout << " SIMPLE PRIMARY EVENT GENERATOR ACTION!!!!!!!!!! \n";

   G4UImanager* UImanager = G4UImanager::GetUIpointer();
   UImanager->ApplyCommand("/gps/pos/type Volume");
   UImanager->ApplyCommand("/gps/pos/shape Cylinder");
   UImanager->ApplyCommand("/gps/pos/centre 0 0 0 cm");
   UImanager->ApplyCommand("/gps/pos/halfz 2.5 cm ");
   UImanager->ApplyCommand("/gps/pos/radius 0.1 mm ");

   UImanager->ApplyCommand("/gps/ang/type iso");
   UImanager->ApplyCommand("/gps/ang/mintheta 0   deg");
   UImanager->ApplyCommand("/gps/ang/maxtheta 100 deg");
   UImanager->ApplyCommand("/gps/ang/minphi  -180 deg");
   UImanager->ApplyCommand("/gps/ang/maxphi   180 deg");

   //UImanager->ApplyCommand("/gps/ene/type   Exp");
   //UImanager->ApplyCommand("/gps/ene/ezero 100 GeV");
   UImanager->ApplyCommand("/gps/ene/type   Lin");
   UImanager->ApplyCommand("/gps/ene/intercept 1.0");
   UImanager->ApplyCommand("/gps/ene/gradient  0.0");

   UImanager->ApplyCommand("/gps/ene/min  1 MeV");
   UImanager->ApplyCommand("/gps/ene/max  1 GeV");
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
      //int pdgcode = part->GetPdgCode();
      //G4ParticleDefinition* particle = particleTable->FindParticle(pdgcode);
      //fParticleGun->SetParticleDefinition(particle);

      fParticleGun->GeneratePrimaryVertex(anEvent);
      int pdgcode = fParticleGun->GetParticleDefinition()->GetPDGEncoding();

      double mass          = part->GetMass();
      double E_kin         = fParticleGun->GetParticleEnergy()/GeV;
      double E_tot         = (E_kin+mass);
      double p_rand        = TMath::Sqrt( E_tot*E_tot - mass*mass );
      G4ThreeVector pvec   = fParticleGun->GetParticleMomentumDirection();
      pvec.setMag(p_rand); 
      G4ThreeVector pos   = fParticleGun->GetParticlePosition();

      //double p_rand   = fP_min     + fdP*G4UniformRand();
      //double th_rand  = fTheta_min + fdTheta*G4UniformRand();
      //double phi_rand = fPhi_min   + fdPhi*G4UniformRand();

      TLorentzVector P4vec = {pvec.x(),pvec.y(),pvec.z(),E_tot}; 
      part->SetPdgCode(pdgcode);
      part->SetMomentum(P4vec);
      part->SetProductionVertex(TLorentzVector(pos.x()/cm,pos.y()/cm,pos.z()/cm, 0.0));

      //std::cout << pdgcode << std::endl;
      //std::cout << " E_kin = " << E_kin/MeV << std::endl;
      //std::cout << " theta = " << pvec.getTheta()/CLHEP::degree << std::endl;
      //std::cout << " phi   = " << pvec.getPhi()/CLHEP::degree << std::endl;
      //std::cout << " mass  = " << mass << std::endl;

      //double KE = part->Energy() - mass;
      //double ux = part->Px();
      //double uy = part->Py();
      //double uz = part->Pz();
      //double x0 = part->Vx()*cm;
      //double y0 = part->Vy()*cm;
      //double z0 = part->Vz()*cm;
      //fParticleGun->SetParticleEnergy( KE*GeV );
      //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
      //fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
   } 
   //fThrownEvent.Print();

}
//______________________________________________________________________________

