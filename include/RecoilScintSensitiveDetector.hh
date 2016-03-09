#ifndef RecoilScintSensitiveDetector_h
#define RecoilScintSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "RecoilScintHit.hh"
#include <array>
#include "RecoilScintEvent.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class G4ParticleDefinition;

class RecoilScintSensitiveDetector : public G4VSensitiveDetector
{
   public:
      G4AnalysisManager * fAnalysisManager;

      clas12::hits::RecoilScintEvent * fRecoilScintEvent = nullptr;

     int    fNtuple_ID;
     int    fEvNum_ntID;
     int    fRunNum_ntID;
     int    fChannel_ntID;
     int    fLambda_ntID;
     int    fTime_ntID;
     int    fTotalEnergy_ntID;
     int    fKineticEnergy_ntID;
     int    fMomentum_ntID;
     int    fTheta_ntID;
     int    fPhi_ntID;
     int    fThetaP_ntID;
     int    fPhiP_ntID;
     int    fPID_ntID;
     int    fPx_ntID;
     int    fPy_ntID;
     int    fPz_ntID;
     int    fx_ntID;
     int    fy_ntID;
     int    fz_ntID;

     int    fTimeAvg_hID;
     int    fLambdaAvg_hID;
     int    fTimeVsLambdaAvg_hID;

     //int    fNtuple2_ID;
     //int    fEvNum2_ntID;
     //int    fRunNum2_ntID;
     //int    fChannel2_ntID;
     //int    fLambda2_ntID;
     //int    fTime2_ntID;

     G4ParticleDefinition * fOpticalPhoton;

  public:
      RecoilScintSensitiveDetector(G4String name);
      ~RecoilScintSensitiveDetector();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      G4int HCID;
      RecoilScintHitsCollection *hitsCollection;

      std::array<RecoilScintHitsCollection*,5> theHits;
      std::array<int,5>                     theHCIDs;

  public:
};




#endif

