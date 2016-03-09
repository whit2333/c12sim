#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
//#include "B1PrimaryGeneratorAction.h"

#include "TTree.h"
#include "CLAS12HitsEvent.h"
#include "SimulationManager.h"


/// Event action class
///

class B1EventAction : public G4UserEventAction
{
   public:
      B1EventAction();
      virtual ~B1EventAction();

      //TTree                         * fEGTree;
      TTree                         * fTree;
      clas12::hits::CLAS12HitsEvent * fCLAS12HitsEvent;
      TClonesArray                  * fTrajectoryVerticies;

      virtual void BeginOfEventAction(const G4Event* event);
      virtual void EndOfEventAction(const G4Event* event);

      void AddEdep(G4double edep) { fEdep += edep; }

      void Reset(){
         fTree = 0;
      }

   private:
      G4double  fEdep;
      int event_number;
      
};


#endif

    
