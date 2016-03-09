#ifndef ECDetectorGeometry_HH
#define ECDetectorGeometry_HH

#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Region.hh"
#include "DriftChamberSensitiveDetector.h"
#include "ECSensitiveDetector.h"

#include "DCWire.h"
#include "G4TwoVector.hh"

class ECDetectorGeometry {

   private:

      // Following CLAS12 Note 2014-008

      double fL1 = 721.723*CLHEP::cm; // nominal distance from target to P (point on face where line is normal)
      G4ThreeVector fL1_pos = { 0.0, fL1*TMath::Sin(25.0*CLHEP::degree), fL1*TMath::Cos(25.0*CLHEP::degree) };

      double fLPO = 95.088*CLHEP::cm; // distance from P to 
      G4ThreeVector fS_pos = { 0.0, -fLPO*TMath::Cos(25.0*CLHEP::degree), fLPO*TMath::Sin(25.0*CLHEP::degree) };

      int    fNStrips = 36;
      int    fNViewLayers = 13;

      double fScintThickness     = 1.0*CLHEP::cm;
      double fPbThickness        = 0.2381*CLHEP::cm;
      double fLayerHalfThickness = 6.19*CLHEP::mm;
      double fTheta_0            = 62.88*CLHEP::degree;
      double fTanTheta_0         = TMath::Tan(fTheta_0);

      std::array<double,10>  fa_params = {{
         0.0,  // dummy to keep the indicies the same (starts at 1) 
         0.0856*CLHEP::mm,
         1864.6*CLHEP::mm,
         4.45635*CLHEP::mm,
         4.3708*CLHEP::mm,
         103.66*CLHEP::mm,
         0.2476*CLHEP::mm,
         94.701*CLHEP::mm,
         0.2256*CLHEP::mm,
         94.926*CLHEP::mm
      }};

      double fd2_param      = 36.4*CLHEP::mm;
      double fd2prime_param = 25.4*CLHEP::mm;

      double y_center(int L) const { return( fa_params.at(1)*double(L-1) ); } 
      double dy_center(int L) const { return( fa_params.at(2) + fa_params.at(3)*double(L-1) ); } 

      double y_pos(int L, double x) const { return( y_center(L) - dy_center(L) + fTanTheta_0*x ); }
      double y_neg(int L, double x) const { return( y_center(L) - dy_center(L) - fTanTheta_0*x ); }
      double y_outer(int L) const { return( y_center(L) + dy_center(L) ); }
      //______________________________________________________________

      double Y_U(int L,int U,double x) const {
         return( -fa_params.at(2) - fa_params.at(4)*double(L-1)  + w_U(L)*double(U-1));
      } 
      double Y_V(int L,int V,double x) const {
         return( y_center(L) - dy_center(L) + w_V(L)*(36+1-V)*TMath::Sqrt( 1.0+fTanTheta_0*fTanTheta_0) + fTanTheta_0*x );
      } 
      double Y_W(int L,int W,double x) const {
         return( y_center(L) - dy_center(L) + w_W(L)*(36+1-W)*TMath::Sqrt( 1.0+fTanTheta_0*fTanTheta_0) - fTanTheta_0*x );
      } 

      double w_U(int L) const { return( fa_params.at(5) + fa_params.at(6)*(L-1) ); } 
      double w_V(int L) const { return( fa_params.at(7) + fa_params.at(8)*(L-2) ); } 
      double w_W(int L) const { return( fa_params.at(9) + fa_params.at(8)*(L-3) ); } 
      //______________________________________________________________

      double d_1(int L) const {
         double d2 = fd2_param ;
         if( L > 15 ) d2 = fd2prime_param;
         return( d2 + (w_V(L)/2.0)*((fTanTheta_0*fTanTheta_0 - 1.0)/(fTanTheta_0)) );
      }

      //______________________________________________________________
      double x_corner_U1(int L, int U) const {
         return( (y_center(L) - dy_center(L) + fa_params.at(2) + fa_params.at(4)*double(L-1) - double(U-1)*w_U(L))/fTanTheta_0 );
      }
      double x_corner_U2(int L, int U) const {
         return( -1.0*x_corner_U1(L,U) );
      }
      double y_corner_U1(int L, int U) const {
         return( y_neg(L, x_corner_U1(L,U) ) );
      }
      double y_corner_U2(int L, int U) const {
         return( y_pos(L, x_corner_U2(L,U) ) );
      }
      //______________________________________________________________

      double x_corner_V1(int L, int V) const {
         return( (-1.0*y_outer(L) + y_center(L) - dy_center(L) + 
                double(fNStrips+1-V)*w_V(L)*TMath::Sqrt(1.0+fTanTheta_0*fTheta_0))/fTanTheta_0 );
      }
      double x_corner_V2(int L, int V) const {
         return( (double(fNStrips+1-V)*w_V(L)*TMath::Sqrt(1.0+fTanTheta_0*fTheta_0))/(2.0*fTanTheta_0) );
      }
      double y_corner_V1(int L, int V) const {
         return( y_outer(L) );
      }
      double y_corner_V2(int L, int V) const {
         return( y_pos(L,x_corner_V2(L,V)) );
      }
      //______________________________________________________________

      double x_corner_W1(int L, int W) const {
         return( (y_outer(L) - 1.0*y_center(L) + dy_center(L) + 
                double(fNStrips+1-W)*w_W(L)*TMath::Sqrt(1.0+fTanTheta_0*fTheta_0))/fTanTheta_0 );
      }
      double x_corner_W2(int L, int W) const {
         return( -1.0*(double(fNStrips+1-W)*w_W(L)*TMath::Sqrt(1.0+fTanTheta_0*fTheta_0))/(2.0*fTanTheta_0) );
      }
      double y_corner_W1(int L, int W) const {
         return( y_outer(L) );
      }
      double y_corner_W2(int L, int W) const {
         return( y_neg(L,x_corner_W2(L,W)) );
      }
      //______________________________________________________________

      //  4 _________ 3
      //    \ front /
      //     \     /
      //      \   /
      //       \_/
      //      1   2
      std::vector<G4TwoVector> fContainer_points;

      G4VSolid*        fContainer_solid;
      G4LogicalVolume* fContainer_log;

      //std::array<G4LogicalVolume*,3> fEmptyRegions_log;
      //std::array<G4VSolid*,3>        fClippingRegions_solid;

      std::array<G4VSolid*,3>        fEndPlates_solid;
      std::array<G4LogicalVolume*,3> fEndPlates_log;
      std::array<G4ThreeVector,3>    fLeftEndPlates_pos;
      std::array<G4ThreeVector,3>    fRightEndPlates_pos;

      std::array< G4RotationMatrix,6>   fLeftEndPlates_rot;
      std::array< G4RotationMatrix,6>   fRightEndPlates_rot;

      G4Region * fRegion_g4Region = nullptr;

   public:

      DriftChamberSensitiveDetector * fSensitiveDetector;
      ECSensitiveDetector       * fSensitiveDetector2;

      //G4VSolid * fRegion1_solid;
      //G4VSolid * fRegion2_solid;
      //G4VSolid * fRegion3_solid;

      //G4VSolid * fClippingRegion1_solid;
      //G4VSolid * fClippingRegion2_solid;
      //G4VSolid * fClippingRegion3_solid;

      //G4LogicalVolume * fRegion1_log;
      //G4LogicalVolume * fRegion2_log;
      //G4LogicalVolume * fRegion3_log;

      //G4LogicalVolume * fEmptyRegion1_log;
      //G4LogicalVolume * fEmptyRegion2_log;
      //G4LogicalVolume * fEmptyRegion3_log;


      G4Material * fGasMaterial;

   public:
      ECDetectorGeometry();
      ~ECDetectorGeometry();

      void BuildLogicalVolumes();
      
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother, int sec, int region);
      G4VPhysicalVolume * PlaceParallelPhysicalVolume(G4LogicalVolume * mother, int sec, int region);

};

#endif

