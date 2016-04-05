#include "B1ActionInitialization.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "SimplePrimaryGeneratorAction.hh"
#include "B1RunAction.hh"
#include "B1EventAction.hh"
#include "B1SteppingAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"


B1ActionInitialization::B1ActionInitialization(G4int rn) : 
   G4VUserActionInitialization(), fRunNumber(rn)
{ }
//______________________________________________________________________________

B1ActionInitialization::~B1ActionInitialization()
{ }
//______________________________________________________________________________

void B1ActionInitialization::BuildForMaster() const
{
   
   SetUserAction(new B1RunAction(fRunNumber));
}
//______________________________________________________________________________

void B1ActionInitialization::Build() const
{

   SetUserAction(new TrackingAction);
   //SetUserAction(new B1PrimaryGeneratorAction);
   //SetUserAction(new SimplePrimaryGeneratorAction);
   //
   SetUserAction(new SteppingAction);

   SetUserAction(new B1RunAction(fRunNumber));

   B1EventAction* eventAction = new B1EventAction;
   SetUserAction(eventAction);

   //SetUserAction(new B1SteppingAction(eventAction));
}  
//______________________________________________________________________________

