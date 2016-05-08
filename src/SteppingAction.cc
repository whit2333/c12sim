#include "SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"
#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction() : G4UserSteppingAction()
{ }
//______________________________________________________________________________

SteppingAction::~SteppingAction()
{ }
//______________________________________________________________________________

void SteppingAction::UserSteppingAction(const G4Step* step)
{
   using namespace CLHEP;
   G4Track * track = step->GetTrack();
   double track_length = track->GetTrackLength();
   //std::cout << " Track Length = " << track_length/cm << std::endl;
   if( track->GetVolume()->GetLogicalVolume() ==  fRecoilHodoscope_log ) {
      if( track_length > twopi*90.0*cm ){
         //std::cout << " killing track\n";
         track->SetTrackStatus(fStopAndKill);
      }
   }
   if( track_length > 9.0*m ){
      //std::cout << " killing track\n";
      track->SetTrackStatus(fStopAndKill);
   }
}
//______________________________________________________________________________

