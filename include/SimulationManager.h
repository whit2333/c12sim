#ifndef SimulationManager_h
#define SimulationManager_h 1

#include <fstream>
#include <string>

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "CLAS12HitsEvent.h"

#include "TFile.h"
#include "TTree.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>

class SimulationMessenger;

/** The simulation manager singleton class.
 *
 */
class SimulationManager {

   public:
      SimulationMessenger           * fSimulationMessenger;
      clas12::hits::CLAS12HitsEvent * fEvent;
      std::string                     fOutputDirectoryName;
      TFile                         * fOutputFile;
      std::string                     fOutputFileName;
      TTree                         * fOutputTree;
      std::string                     fOutputTreeName;

   private: 
      SimulationManager( );
      static SimulationManager* fgSimulationManager;

   public:
      ~SimulationManager();

      static SimulationManager* GetInstance();

      static void Dispose();

      void SetRunNumber(int i) { fEvent->fRunNumber = i; }
      int  GetRunNumber() const { return fEvent->fRunNumber; }

      std::string OutputFileName() const ;
      std::string OutputTreeName() const ;

      //G4double GetBeamEnergy() const {return(fBeamEnergy);}
      //void     SetBeamEnergy(G4double en){fBeamEnergy=en;}
      //G4double GetTargetAngle() const {   return(fTargetAngle);}
      //void     SetTargetAngle(G4double a){fTargetAngle=a;}

      int GetEventNumber() const { return(fEvent->fEventNumber); };


      int ReadRunNumber(const char * fname = "run_number.txt");

      ///**
      // * Increments the run number in memory and in file/database
      // * returns the run number
      // * \todo Make the source of run number a database, not a text file. 
      // */
      //int IncrementRunNumber();

      ///**
      // * Gets the run number from a file/database
      // * \todo Make the source of run number a database, not a text file. 
      // */
      //int RetrieveRunNumber(int run = 0);

      ///**
      // * Allocate Event and Hit memory
      // */
      //int AllocateTreeMemory();

      ///**
      // * Free Event and Hit memory
      // */
      //int Reset();

      ///** Creates the InSANEDetector Base classes and sets event addresses
      // *  used in pedestal and noise simulation
      // */
      //int AddDetectors(int runNumber = 0) ;

      ///**
      // * Sets up detector's scoring
      // * This defines what should be detected but not which volume it is associated with
      // * this should be done within the det const class
      // */
      //int InitScoring();

      ///**
      // *  Defines the scoring filters
      // *  This must come during detector construction...??
      // */ 
      //int DefineScoringFilters();

      ///**
      // *  Only increments some counter in text file.
      // */
      //int InitializeNewRun(int run = 0);

   private:

      /*! Try to get lock. Return its file descriptor or -1 if failed.
       *
       *  @param lockName Name of file used as lock (i.e. '/var/lock/myLock').
       *  @return File descriptor of lock file, or -1 if failed.
       */
      int tryGetLock( char const *lockName )
      {
         mode_t m = umask( 0 );
         int fd = open( lockName, O_RDWR|O_CREAT, 0666 );
         umask( m );
         if( fd >= 0 && flock( fd, LOCK_EX | LOCK_NB ) < 0 )
         {
            close( fd );
            fd = -1;
         }
         return fd;
      }

      /*! Release the lock obtained with tryGetLock( lockName ).
       *
       *  @param fd File descriptor of lock returned by tryGetLock( lockName ).
       *  @param lockName Name of file used as lock (i.e. '/var/lock/myLock').
       */
      void releaseLock( int fd, char const *lockName )
      {
         if( fd < 0 )
            return;
         remove( lockName );
         close( fd );
      }



};

#endif

