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
#include "TClonesArray.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <iostream>

class RecoilChamberDetectorGeometry;
class RecoilHodoDetectorGeometry;
class RecoilHodoDetectorGeometry2;
class RecoilHodoDetectorGeometry3;
class DriftChamberDetectorGeometry;
class HTCCDetectorGeometry;
class SolenoidDetectorGeometry;
class ECDetectorGeometry;
class TorusDetectorGeometry;
class BeamlineDetectorGeometry;
class SimulationMessenger;
class MicromegasVertexTrackerDetectorGeometry;
class SiliconVertexTrackerDetectorGeometry;
class SiVertexTrackerDetectorGeometry;

using MVTDetectorGeometry = MicromegasVertexTrackerDetectorGeometry;
using SVTDetectorGeometry = SiliconVertexTrackerDetectorGeometry;


/** The simulation manager singleton class.
 */
class SimulationManager {

   private:
      RecoilChamberDetectorGeometry * fRecoilChamberGeo = nullptr;
      RecoilHodoDetectorGeometry    * fRecoilHodoGeo    = nullptr;
      RecoilHodoDetectorGeometry2   * fRecoilHodoGeo2   = nullptr;
      RecoilHodoDetectorGeometry3   * fRecoilHodoGeo3   = nullptr;
      DriftChamberDetectorGeometry  * fDriftChamberGeo  = nullptr;
      HTCCDetectorGeometry          * fHTCCGeo          = nullptr;
      SolenoidDetectorGeometry      * fSolenoidGeo      = nullptr;
      TorusDetectorGeometry         * fTorusGeo         = nullptr;
      BeamlineDetectorGeometry      * fBeamlineGeo      = nullptr;
      MVTDetectorGeometry           * fMVTGeo           = nullptr;
      SVTDetectorGeometry           * fSVTGeo           = nullptr;
      ECDetectorGeometry            * fECGeo            = nullptr;
      SiVertexTrackerDetectorGeometry * fSiVertexTrackerGeo           = nullptr;


   public:
      DriftChamberDetectorGeometry  * GetDriftDetectorGeometry();
      RecoilHodoDetectorGeometry    * GetRecoilHodoDetectorGeometry();
      RecoilHodoDetectorGeometry2   * GetRecoilHodoDetectorGeometry2();
      RecoilHodoDetectorGeometry3   * GetRecoilHodoDetectorGeometry3();
      RecoilChamberDetectorGeometry * GetRecoilDetectorGeometry();
      HTCCDetectorGeometry          * GetHTCCDetectorGeometry();
      SolenoidDetectorGeometry      * GetSolenoidDetectorGeometry();
      TorusDetectorGeometry         * GetTorusDetectorGeometry();
      BeamlineDetectorGeometry      * GetBeamlineDetectorGeometry();
      MVTDetectorGeometry           * GetMVTDetectorGeometry();
      SVTDetectorGeometry           * GetSVTDetectorGeometry();
      ECDetectorGeometry            * GetECDetectorGeometry();
      SiVertexTrackerDetectorGeometry           * GetSiVertexTrackerDetectorGeometry();

   public:
      SimulationMessenger           * fSimulationMessenger;
      TFile                         * fOutputFile;
      TTree                         * fOutputTree;
      TTree                         * fEGTree;
      std::string                     fOutputDirectoryName;
      std::string                     fOutputFileName;
      std::string                     fOutputTreeName;
      std::string                     fInputFileName;

      clas12::hits::CLAS12HitsEvent * fEvent;
      clas12::hits::TriggerEvent    * fTriggerEvent;
      TClonesArray                  * fTrajectoryVerticies;

      double                          fToroidFieldScale = 0.0;
      double                          fSolenoidFieldScale = 0.0;

   private: 
      SimulationManager( );
      static SimulationManager* fgSimulationManager;

   public:
      ~SimulationManager();

      static SimulationManager* GetInstance();

      //static void Dispose();

      void SetRunNumber(int i);// { fEvent->fRunNumber = i; }
      int  GetRunNumber() const { return fEvent->fRunNumber; }

      std::string OutputFileName() const ;
      std::string OutputTreeName() const ;
      std::string InputFileName() const ;

      void SetOutputDirectoryName( std::string  s ) { fOutputDirectoryName = s;}
      void SetOutputFileName     ( std::string  s ) { fOutputFileName      = s;}
      void SetOutputTreeName     ( std::string  s ) { fOutputTreeName      = s;}
      void SetInputFileName      ( std::string  s ) { fInputFileName       = s;}

      int  GetEventNumber() const { return(fEvent->fEventNumber); };
      int  ReadRunNumber(const char * fname = "run_number.txt");

      void PrintConfig(std::ostream& s = std::cout) const;


      double    GetToroidFieldScale()   const { return fToroidFieldScale ;}
      double    GetSolenoidFieldScale() const { return fSolenoidFieldScale ;}
      void      SetToroidFieldScale(  double v) { fToroidFieldScale   = v;}
      void      SetSolenoidFieldScale(double v) { fSolenoidFieldScale = v;}

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

