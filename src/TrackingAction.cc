#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack) 
{
  fpTrackingManager->SetStoreTrajectory(true);
  fpTrackingManager->SetTrajectory(new G4Trajectory(aTrack));
}
//______________________________________________________________________________

void TrackingAction::PostUserTrackingAction(const G4Track* /*aTrack*/)
{;}
//______________________________________________________________________________


