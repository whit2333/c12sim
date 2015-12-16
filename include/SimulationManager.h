#ifndef SimulationManager_h
#define SimulationManager_h 1

#include <fstream>
#include <string>

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "CLAS12HitsEvent.h"
#include "TriggerEvent.h"

#include "TFile.h"
#include "TTree.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>

class RecoilChamberDetectorGeometry;
class DriftChamberDetectorGeometry;
class HTCCDetectorGeometry;
class SolenoidDetectorGeometry;
class TorusDetectorGeometry;
class BeamlineDetectorGeometry;

class SimulationMessenger;

/** The simulation manager singleton class.
 */
class SimulationManager {

   private:
      RecoilChamberDetectorGeometry * fRecoilChamberGeo = nullptr;
      DriftChamberDetectorGeometry  * fDriftChamberGeo  = nullptr;
      HTCCDetectorGeometry          * fHTCCGeo          = nullptr;
      SolenoidDetectorGeometry      * fSolenoidGeo      = nullptr;
      TorusDetectorGeometry         * fTorusGeo         = nullptr;
      BeamlineDetectorGeometry      * fBeamlineGeo      = nullptr;

   public:
      DriftChamberDetectorGeometry  * GetDriftDetectorGeometry();
      RecoilChamberDetectorGeometry * GetRecoilDetectorGeometry();
      HTCCDetectorGeometry          * GetHTCCDetectorGeometry();
      SolenoidDetectorGeometry      * GetSolenoidDetectorGeometry();
      TorusDetectorGeometry         * GetTorusDetectorGeometry();
      BeamlineDetectorGeometry      * GetBeamlineDetectorGeometry();

   public:
      SimulationMessenger           * fSimulationMessenger;
      TFile                         * fOutputFile;
      TTree                         * fOutputTree;
      std::string                     fOutputDirectoryName;
      std::string                     fOutputFileName;
      std::string                     fOutputTreeName;
      std::string                     fInputFileName;

      clas12::hits::CLAS12HitsEvent * fEvent;
      clas12::hits::TriggerEvent    * fTriggerEvent;

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
      std::string InputFileName() const ;

      void SetOutputDirectoryName( std::string  s ) { fOutputDirectoryName = s;}
      void SetOutputFileName     ( std::string  s ) { fOutputFileName      = s;}
      void SetOutputTreeName     ( std::string  s ) { fOutputTreeName      = s;}
      void SetInputFileName      ( std::string  s ) { fInputFileName       = s;}

      int GetEventNumber() const { return(fEvent->fEventNumber); };
      int ReadRunNumber(const char * fname = "run_number.txt");


   private:

      /*! Try to get lock. Return its file descriptor or -1 if failed.
       *
       *  @param lockName Name of file used as lock (i.e. '/var/lock/myLock').
       *  @return File descriptor of lock file, or -1 if failed.
       */
      int tryGetLock( char const *lockName );
      

      /*! Release the lock obtained with tryGetLock( lockName ).
       *
       *  @param fd File descriptor of lock returned by tryGetLock( lockName ).
       *  @param lockName Name of file used as lock (i.e. '/var/lock/myLock').
       */
      void releaseLock( int fd, char const *lockName );

};

#endif

