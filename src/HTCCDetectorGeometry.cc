#include "HTCCDetectorGeometry.h"

#include "G4UserLimits.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Paraboloid.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4IntersectionSolid.hh"
#include "G4GenericTrap.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4IntersectionSolid.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4UserLimits.hh"
#include "HTCCSensitiveDetector.h"

HTCCDetectorGeometry::HTCCDetectorGeometry()
{
   using namespace clas12::geo;
   using namespace CLHEP;

   G4SDManager* SDMan = G4SDManager::GetSDMpointer();
   fSensitiveDetector = new HTCCSensitiveDetector("HTCC",4*2*6);
   SDMan->AddNewDetector(fSensitiveDetector);

   G4NistManager* nist = G4NistManager::Instance();
   G4Material * default_mat   = nist->FindOrBuildMaterial("G4_AIR");

   //htccBigGasVolume  | root  |volume containing cherenkov gas  | 0 0 0  | 0 0 0  |ee99ff5   | Polycone  |0*deg 360*deg 4*counts 0*mm 15.8*mm 91.5*mm 150*mm 1742*mm 2300*mm 2300*mm 1589*mm -275*mm 181*mm 1046*mm 1740*mm  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  | no  | no 
   int zplanes = 4; // number of planes in z directions
   double rInner[] = { 0*mm, 15.8*mm, 91.5*mm, 150*mm };
   double rOuter[] = { 1742*mm, 2300*mm, 2300*mm, 1589*mm };
   double zPlane[] = { -275*mm, 181*mm, 1046*mm, 1740*mm };

   htccBigGasVolume_solid = new G4Polycone("htccBigGasVolume_solid",            ///< name
         0,  ///< Initial Phi starting angle
         360*deg,  ///< Total Phi angle
         zplanes,        ///< Number of z planes
         zPlane,         ///< z coordinate of corners
         rInner,         ///< Tangent distance to inner surface
         rOuter);        ///< Tangent distance to outer surface

   //______________________________________________________________________________
   //   htccEntryDishVolume  |                root  |        HTCC entry dish volume  |                                             0 0 0  |                                   0 0 0  | ee99ff   |            Polycone  |0.0*deg 360.0*deg 17*counts 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 1416.05*mm 1410.081*mm 1404.493*mm 1398.905*mm 1393.317*mm 1387.729*mm 1387.729*mm 1363.4466*mm 1339.1388*mm 1305.8648*mm 1272.5908*mm 1239.3168*mm 1206.0174*mm 1158.24*mm 1105.408*mm 996.0356*mm 945.896*mm -276.4122*mm -178.6222*mm -87.1822*mm 4.2578*mm 95.6978*mm 187.1378*mm 278.5778*mm 370.0178*mm 461.4578*mm 552.8978*mm 644.3378*mm 735.7778*mm 827.2178*mm 918.6578*mm 991.378*mm 1080.278*mm 1107.71*mm  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 
   zplanes = 17; // number of planes in z directions
   double rInner_2[] = { 0*mm, 0*mm, 0*mm, 0*mm, 0*mm, 0*mm, 0*mm, 0*mm, 0*mm, 0*mm,
      0*mm, 0*mm, 0*mm, 0*mm, 0*mm, 0*mm, 0*mm};
   double rOuter_2[] = {1416.05*mm, 1410.081*mm, 1404.493*mm, 1398.905*mm, 1393.317*mm,
      1387.729*mm, 1387.729*mm, 1363.4466*mm, 1339.1388*mm, 1305.8648*mm,
      1272.5908*mm, 1239.3168*mm, 1206.0174*mm, 1158.24*mm, 1105.408*mm,
      996.0356*mm, 945.896*mm};
   double zPlane_2[] = {-276.4122*mm, -178.6222*mm, -87.1822*mm, 4.2578*mm, 95.6978*mm,
      187.1378*mm, 278.5778*mm, 370.0178*mm, 461.4578*mm, 552.8978*mm,
      644.3378*mm, 735.7778*mm, 827.2178*mm, 918.6578*mm, 991.378*mm,
      1080.278*mm, 1107.71*mm };
   htccEntryDishVolume_solid = new G4Polycone("htccEntryDishVolume_solid",
         0,         ///< Initial Phi starting angle
         360*deg,   ///< Total Phi angle
         zplanes,        ///< Number of z planes
         zPlane_2,         ///< z coordinate of corners
         rInner_2,         ///< Tangent distance to inner surface
         rOuter_2);        ///< Tangent distance to outer surface

   //______________________________________________________________________________
   //   htccEntryConeVolume  |                root  |        HTCC entry cone volume  |                                             0 0 0  |                                   0 0 0  | ee99ff   |            Polycone  |0.0*deg 360.0*deg 9*counts 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 257.505*mm 323.952*mm 390.373*mm 456.819*mm 525.831*mm 599.872*mm 673.913*mm 747.979*mm 827.151*mm 380*mm 470.17*mm 561.61*mm 653.05*mm 744.49*mm 835.93*mm 927.37*mm 1018.81*mm 1116.6*mm  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 
   zplanes = 9; // number of planes in z directions
   double rInner_3[] = {0*mm, 0*mm, 0*mm,
      0*mm, 0*mm, 0*mm,
      0*mm, 0*mm, 0*mm};
   double rOuter_3[] = { 257.505*mm, 323.952*mm, 390.373*mm,
      456.819*mm, 525.831*mm, 599.872*mm,
      673.913*mm, 747.979*mm, 827.151*mm};

   double zPlane_3[] = {
      380*mm, 470.17*mm, 561.61*mm,
      653.05*mm, 744.49*mm, 835.93*mm,
      927.37*mm, 1018.81*mm, 1116.6*mm};
   htccEntryConeVolume_solid = new G4Polycone("htccEntryDishVolume_solid",
         0,         ///< Initial Phi starting angle
         360*deg,   ///< Total Phi angle
         zplanes,        ///< Number of z planes
         zPlane_3,         ///< z coordinate of corners
         rInner_3,         ///< Tangent distance to inner surface
         rOuter_3);        ///< Tangent distance to outer surface

   //______________________________________________________________________________
   //   htccEntryDishCone  |                root  | subtraction entry dish - cone  |                                             0 0 0  |                                   0 0 0  | ee99ff   |Operation:@ htccEntryDishVolume - htccEntryConeVolume  |                                                           0  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 

   htccEntryDishCone_solid = new G4SubtractionSolid("htccEntryDishCone_solid",htccEntryDishVolume_solid , htccEntryConeVolume_solid);

   //   htcc  |                root  |           gas volume for htcc  |                                             0 0 0  |                                   0 0 0  |0000ff3   |Operation:@ htccBigGasVolume - htccEntryDishCone  |                                                           0  |             htccGas  |                  no  |     1   |     1   |     1   |   1   |   1   |                  no  |                  no  |                                      no 

   htcc_solid = new G4SubtractionSolid("htccEntryDishCone_solid",htccBigGasVolume_solid , htccEntryDishCone_solid);

   // wedge for sector placement
   G4VSolid * temp = new G4Tubs("temp_wedge",0,3.0*m, 4.0*m, 60.0*deg, 60.0*deg);
   sector_wedge_solid = new G4IntersectionSolid("sector_wedge_solid",htcc_solid , temp);

}
//______________________________________________________________________________

HTCCDetectorGeometry::~HTCCDetectorGeometry()
{
}
//______________________________________________________________________________
void HTCCDetectorGeometry::BuildMirrors()
{
   using namespace CLHEP;

   // ------------------------------------------------------------------------
   //HTCC_OuterCutCylinder  |                htcc  |            Outer cut cylinder  |                                 0*mm 0*mm 1500*mm  |                                   0 0 0  | 0000ff   |                Tube  |                    |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 

   HTCC_OuterCutCylinder_solid = new G4Tubs("HTCC_OuterCutCylinder_solid", 1165.352*mm, 1500.0*mm, 500*mm, 0*deg, 360*deg);
   G4ThreeVector HTCC_OuterCutCylinder_trans(0*mm, 0*mm, 1500*mm);
   G4RotationMatrix HTCC_OuterCutCylinder_rot;

   // ------------------------------------------------------------------------
   //HTCC_InnerCutCone  |                htcc  |                Inner cut cone  |                                 0*mm 0*mm 1000*mm  |                                   0 0 0  | 0000ff   |                Cons  |   0.0*mm 0.0*mm 0.0*mm 174.97733*mm 1000.0*mm 0*deg 360*deg  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 
   HTCC_InnerCutCone_solid = new G4Cons("HTCC_InnerCutCone_solid", 0.0*mm, 0.0*mm, 0.0*mm, 174.97733*mm, 1000.0*mm, 0*deg, 360*deg );
   G4ThreeVector HTCC_InnerCutCone_trans(0*mm, 0*mm, 1000*mm);
   G4RotationMatrix HTCC_InnerCutCone_rot;

   // ------------------------------------------------------------------------
   //Barrel_sect0half0  |                htcc  |Barrel defining mirror back surface  |                                             0 0 0  |ordered: yzx 0*rad   0.2617993878*rad    1.392074611*rad  | ee99ff   |            Polycone  |0.0*deg 360.0*deg 21*counts 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 1435.3392*mm 1461.162*mm 1483.4626*mm 1502.4172*mm 1518.1665*mm 1530.8215*mm 1540.4676*mm 1547.168*mm 1550.9654*mm 1551.8836*mm 1549.9283*mm 1545.0874*mm 1537.3303*mm 1526.6071*mm 1512.8469*mm 1495.9556*mm 1475.8119*mm 1452.2631*mm 1425.1181*mm 1394.1386*mm 1359.0258*mm 360.01139*mm 423.09863*mm 486.18587*mm 549.27312*mm 612.36036*mm 675.4476*mm 738.53484*mm 801.62208*mm 864.70932*mm 927.79656*mm 990.8838*mm 1053.971*mm 1117.0583*mm 1180.1455*mm 1243.2328*mm 1306.32*mm 1369.4073*mm 1432.4945*mm 1495.5817*mm 1558.669*mm 1621.7562*mm  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 
   int zplanes = 21; // number of planes in z directions
   double rInner[] = {0*mm, 0*mm, 0*mm, 0*mm, 0*mm,
      0*mm, 0*mm, 0*mm, 0*mm, 0*mm,
      0*mm, 0*mm, 0*mm, 0*mm, 0*mm,
      0*mm, 0*mm, 0*mm, 0*mm, 0*mm,
      0*mm};
   double rOuter[] = { 1435.3392*mm, 1461.162*mm, 1483.4626*mm, 1502.4172*mm, 1518.1665*mm,
      1530.8215*mm, 1540.4676*mm, 1547.168*mm, 1550.9654*mm, 1551.8836*mm,
      1549.9283*mm, 1545.0874*mm, 1537.3303*mm, 1526.6071*mm, 1512.8469*mm,
      1495.9556*mm, 1475.8119*mm, 1452.2631*mm, 1425.1181*mm, 1394.1386*mm,
      1359.0258*mm};

   double zPlane[] = { 360.01139*mm, 423.09863*mm, 486.18587*mm, 549.27312*mm, 612.36036*mm,
      675.4476*mm, 738.53484*mm, 801.62208*mm, 864.70932*mm, 927.79656*mm,
      990.8838*mm, 1053.971*mm, 1117.0583*mm, 1180.1455*mm, 1243.2328*mm,
      1306.32*mm, 1369.4073*mm, 1432.4945*mm, 1495.5817*mm, 1558.669*mm,
      1621.7562*mm };
   Barrel_sect0half0_solid = new G4Polycone("Barrel_sect0half0_solid",
         0,         ///< Initial Phi starting angle
         360*deg,   ///< Total Phi angle
         zplanes,        ///< Number of z planes
         zPlane,         ///< z coordinate of corners
         rInner,         ///< Tangent distance to inner surface
         rOuter);        ///< Tangent distance to outer surface
   G4ThreeVector Barrel_sect0half0_trans(0,0,0);
   G4RotationMatrix Barrel_sect0half0_rot;
   Barrel_sect0half0_rot.rotateY(   0.0*rad);
   Barrel_sect0half0_rot.rotateZ(   0.2617993878*rad );
   Barrel_sect0half0_rot.rotateX(   1.392074611*rad);


   // ------------------------------------------------------------------------
   //phicut_sect0half0 |   htcc |   Half-sector phi cut |      0*mm 0*mm 1500*mm |                  0 0 0 | 443388  |        Tube |0.0*mm 1500.0*mm 500.0*mm  1.047197551*rad  0.5225987756*rad |      Component |         no |   1  |   1  |   1  |  1  |  0  |         no |         no |                   no 
   phicut_sect0half0_solid = new G4Tubs("phicut_sect0half0_solid", 0.0*mm,1500.0*mm, 500.0*mm,   1.047197551*rad,  0.5225987756*rad );
   G4ThreeVector    phicut_sect0half0_trans(0,0,1500*mm);
   G4RotationMatrix phicut_sect0half0_rot;

   // ------------------------------------------------------------------------
   // Mirror 4 
   //Mirror_sect0mirr0half0  |                htcc  |Ellipsoid defining mirror surface  |        208.7414292*mm 779.0336195*mm -31.02335*mm  |ordered: yzx 0*rad   0.2617993878*rad    1.609243305*rad  | ff8080   |           Ellipsoid  |               1728.673*mm 1728.673*mm 1907.810*mm 0*mm 0*mm  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 
   Mirror_sect0mirr0half0_solid = new G4Ellipsoid("Mirror_sect0mirr0half0_solid",1728.673*mm, 1728.673*mm, 1907.810*mm, 0*mm, 0*mm); 
   //Mirror_sect0mirr0half0_solid = new G4Sphere("Mirror_sect0mirr0half0_solid",0, 1907.810*mm, 0*deg, 360*deg, 0.0*deg, 180.0*deg); 
   G4ThreeVector    Mirror_sect0mirr0half0_trans(208.7414292*mm, 779.0336195*mm, -31.02335*mm);
   G4RotationMatrix Mirror_sect0mirr0half0_rot;
   Mirror_sect0mirr0half0_rot.rotateY( 0*rad  );
   Mirror_sect0mirr0half0_rot.rotateZ( 0.2617993878*rad );
   Mirror_sect0mirr0half0_rot.rotateX( 1.609243305*rad);

   // ------------------------------------------------------------------------
   //BarrelEllipseCut_sect0mirr0half0  |                htcc  |subtraction of barrel and ellipse  |                                    0*mm 0*mm 0*mm  |ordered: yzx 0*rad   0.2617993878*rad    1.392074611*rad  | 00ff00   |Operation:@ Barrel_sect0half0 - Mirror_sect0mirr0half0  |                                                           0  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 
   G4RotationMatrix &r1   = Barrel_sect0half0_rot;
   G4RotationMatrix &r2   = Mirror_sect0mirr0half0_rot;
   G4RotationMatrix r_tot = r2*(r1.inverse());
   G4ThreeVector    &t1   = Barrel_sect0half0_trans;
   G4ThreeVector    &t2   = Mirror_sect0mirr0half0_trans;
   G4ThreeVector    t_tot = (t2-t1);
   t_tot *= r1;
   BarrelEllipseCut_sect0mirr0half0_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr0half0_solid",
         Barrel_sect0half0_solid,
         Mirror_sect0mirr0half0_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector    BarrelEllipseCut_sect0mirr0half0_trans(0,0,0);
   G4RotationMatrix BarrelEllipseCut_sect0mirr0half0_rot;
   BarrelEllipseCut_sect0mirr0half0_rot.rotateY( 0*rad  );
   BarrelEllipseCut_sect0mirr0half0_rot.rotateZ( 0.2617993878*rad );
   BarrelEllipseCut_sect0mirr0half0_rot.rotateX( 1.392074611*rad);


   // ------------------------------------------------------------------------
   //Boxcut12_sect0mirr0half0  |                htcc  |Box defining dividing plane 1/2  |     -184.7943978*mm 242.3741479*mm 978.1545571*mm  |ordered: yzx 0*rad   0.2617980211*rad   0.7951867457*rad  | ff0000   |                 Box  |                                     1000*mm 1000*mm 1000*mm  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 
   Boxcut12_sect0mirr0half0_solid = new G4Box("Boxcut12_sect0mirr0half0_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut12_sect0mirr0half0_trans(-184.7943978*mm, 242.3741479*mm, 978.1545571*mm);
   G4RotationMatrix Boxcut12_sect0mirr0half0_rot;
   Boxcut12_sect0mirr0half0_rot.rotateY( 0*rad  );
   Boxcut12_sect0mirr0half0_rot.rotateZ( 0.2617980211*rad    );
   Boxcut12_sect0mirr0half0_rot.rotateX( 0.7951867457*rad);

   // ------------------------------------------------------------------------
   //MirrorBoxCut_sect0mirr0half0  |                htcc  |subtraction of box and (barrel - ellipse)  |                                    0*mm 0*mm 0*mm  |ordered: yzx 0*rad   0.2617993878*rad    1.392074611*rad  | 00ff00   |Operation:@ BarrelEllipseCut_sect0mirr0half0 - Boxcut12_sect0mirr0half0  |                                                           0  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 
   r1   = BarrelEllipseCut_sect0mirr0half0_rot;
   r2   = Boxcut12_sect0mirr0half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = BarrelEllipseCut_sect0mirr0half0_trans;
   t2   = Boxcut12_sect0mirr0half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_sect0mirr0half0_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr0half0_solid",
         BarrelEllipseCut_sect0mirr0half0_solid, 
         Boxcut12_sect0mirr0half0_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_sect0mirr0half0_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_sect0mirr0half0_rot;
   MirrorBoxCut_sect0mirr0half0_rot.rotateY( 0*rad  );
   MirrorBoxCut_sect0mirr0half0_rot.rotateZ( 0.2617993878*rad    );
   MirrorBoxCut_sect0mirr0half0_rot.rotateX( 1.392074611*rad);

   // -----------------------------------------------------------------------
   //MirrorCylinderCut_sect0mirr0half0  |                htcc  |subtraction of cylinder and (box - (barrel - ellipse))  | 0*mm 0*mm 0*mm  |ordered: yzx 0*rad   0.2617993878*rad    1.392074611*rad  | 00ff00   |Operation:@ MirrorBoxCut_sect0mirr0half0 - HTCC_OuterCutCylinder  |                                                           0  |           Component  |                  no  |     1   |     1   |     1   |   1   |   0   |                  no  |                  no  |                                      no 
   r1   = MirrorBoxCut_sect0mirr0half0_rot;
   r2   = HTCC_OuterCutCylinder_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorBoxCut_sect0mirr0half0_trans;
   t2   = HTCC_OuterCutCylinder_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorCylinderCut_sect0mirr0half0_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr0half0_solid",
         MirrorBoxCut_sect0mirr0half0_solid,
         HTCC_OuterCutCylinder_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorCylinderCut_sect0mirr0half0_trans(0,0,0);
   G4RotationMatrix MirrorCylinderCut_sect0mirr0half0_rot;
   MirrorCylinderCut_sect0mirr0half0_rot.rotateY( 0*rad    );
   MirrorCylinderCut_sect0mirr0half0_rot.rotateZ( 0.2617993878*rad );
   MirrorCylinderCut_sect0mirr0half0_rot.rotateX( 1.392074611*rad);

   // ------------------------------------------------------------------------
   //mirror_4_sector2_2  | htcc  |  htcc mirror 4, sector2 right  | 0*mm 0*mm 0*mm  |ordered: yzx 0*rad   0.2617993878*rad    1.392074611*rad  | ff8080   |Operation:@ MirrorCylinderCut_sect0mirr0half0 * phicut_sect0half0  |                                                           0  |          rohacell31  |                  no  |     1   |     1   |     1   |   1   |   1   | mirror: htcc_AlMgF2  |              mirror  |                             id manual 0 
   r1   = MirrorCylinderCut_sect0mirr0half0_rot;
   r2   = phicut_sect0half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorCylinderCut_sect0mirr0half0_trans;
   t2   = phicut_sect0half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   mirror_4_sector2_2_solid = new G4IntersectionSolid(
         "mirror_4_sector2_2_solid",
         MirrorCylinderCut_sect0mirr0half0_solid , 
         phicut_sect0half0_solid,
         &r_tot,
         t_tot
         );
   mirror_4_sector2_2_rot =  G4RotationMatrix::IDENTITY ;
   mirror_4_sector2_2_rot.rotateY( 0*rad    );
   mirror_4_sector2_2_rot.rotateZ( 0.2617993878*rad );
   mirror_4_sector2_2_rot.rotateX( 1.392074611*rad);

   // ------------------------------------------------------------------------
   // Mirror 3 
   //Mirror_sect0mirr1half0 | htcc |Ellipsoid defining mirror surface | 231.4737648*mm 863.8718507*mm 81.66635*mm |ordered: yzx 0*rad 0.2617993878*rad 1.479734815*rad | 8080ff | Ellipsoid | 1612.998*mm 1612.998*mm 1846.155*mm 0*mm 0*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Mirror_sect0mirr1half0_solid = new G4Ellipsoid("Mirror_sect0mirr1half0_solid",1612.998*mm,1612.998*mm,1846.155*mm,0*mm,0*mm); 
   G4ThreeVector    Mirror_sect0mirr1half0_trans(231.4737648*mm,863.8718507*mm,81.66635*mm);
   G4RotationMatrix Mirror_sect0mirr1half0_rot;
   Mirror_sect0mirr1half0_rot.rotateY( 0*rad  );
   Mirror_sect0mirr1half0_rot.rotateZ( 0.2617993878*rad );
   Mirror_sect0mirr1half0_rot.rotateX( 1.479734815*rad);
   // ------------------------------------------------------------------------
   //BarrelEllipseCut_sect0mirr1half0 | htcc |subtraction of barrel and ellipse |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | 00ff00 |Operation:@ Barrel_sect0half0 - Mirror_sect0mirr1half0 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = Barrel_sect0half0_rot ;
   r2   = Mirror_sect0mirr1half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = Barrel_sect0half0_trans;
   t2   = Mirror_sect0mirr1half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   BarrelEllipseCut_sect0mirr1half0_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr1half0_solid",
         Barrel_sect0half0_solid,
         Mirror_sect0mirr1half0_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector    BarrelEllipseCut_sect0mirr1half0_trans(0,0,0);
   G4RotationMatrix BarrelEllipseCut_sect0mirr1half0_rot;
   BarrelEllipseCut_sect0mirr1half0_rot.rotateY( 0*rad  );
   BarrelEllipseCut_sect0mirr1half0_rot.rotateZ( 0.2617993878*rad );
   BarrelEllipseCut_sect0mirr1half0_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //Boxcut_up_sect0mirr1half0 | htcc | upper plane cut | 184.7943978*mm 1621.705852*mm 2378.357443*mm |ordered: yzx 0*rad 0.2617980211*rad 0.7951867457*rad | 00ff00 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut_up_sect0mirr1half0_solid = new G4Box("Boxcut_up_sect0mirr1half0_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut_up_sect0mirr1half0_trans(184.7943978*mm,1621.705852*mm,2378.357443*mm);
   G4RotationMatrix Boxcut_up_sect0mirr1half0_rot;
   Boxcut_up_sect0mirr1half0_rot.rotateY( 0*rad  );
   Boxcut_up_sect0mirr1half0_rot.rotateZ( 0.2617980211*rad    );
   Boxcut_up_sect0mirr1half0_rot.rotateX( 0.7951867457*rad);
   // ------------------------------------------------------------------------
   //Boxcut_down_sect0mirr1half0 | htcc | lower plane cut | -139.5367254*mm 119.9598692*mm 821.4432351*mm |ordered: yzx 0*rad 0.2618003565*rad 0.5693996995*rad | 00ff00 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut_down_sect0mirr1half0_solid = new G4Box("Boxcut_down_sect0mirr1half0_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut_down_sect0mirr1half0_trans(-139.5367254*mm,119.9598692*mm,821.4432351*mm);
   G4RotationMatrix Boxcut_down_sect0mirr1half0_rot;
   Boxcut_down_sect0mirr1half0_rot.rotateY( 0*rad  );
   Boxcut_down_sect0mirr1half0_rot.rotateZ( 0.2617980211*rad    );
   Boxcut_down_sect0mirr1half0_rot.rotateX( 0.5693996995*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_up_sect0mirr1half0 | htcc |subtraction of upper box from (barrel - ellipse) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | 0000ff |Operation:@ BarrelEllipseCut_sect0mirr1half0 - Boxcut_up_sect0mirr1half0 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = BarrelEllipseCut_sect0mirr1half0_rot;
   r2   = Boxcut_up_sect0mirr1half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = BarrelEllipseCut_sect0mirr1half0_trans;
   t2   = Boxcut_up_sect0mirr1half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_up_sect0mirr1half0_solid = new G4SubtractionSolid(
         "MirrorBoxCut_up_sect0mirr1half0_solid",
         BarrelEllipseCut_sect0mirr1half0_solid, 
         Boxcut_up_sect0mirr1half0_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_up_sect0mirr1half0_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_up_sect0mirr1half0_rot;
   MirrorBoxCut_up_sect0mirr1half0_rot.rotateY( 0*rad  );
   MirrorBoxCut_up_sect0mirr1half0_rot.rotateZ( 0.2617993878*rad    );
   MirrorBoxCut_up_sect0mirr1half0_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_down_sect0mirr1half0 | htcc |subtraction of lower box from (barrel - ellipse - upper box) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | 0000ff |Operation:@ MirrorBoxCut_up_sect0mirr1half0 - Boxcut_down_sect0mirr1half0 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = MirrorBoxCut_up_sect0mirr1half0_rot;
   r2   = Boxcut_down_sect0mirr1half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorBoxCut_up_sect0mirr1half0_trans;
   t2   = Boxcut_down_sect0mirr1half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_down_sect0mirr1half0_solid = new G4SubtractionSolid(
         "MirrorBoxCut_down_sect0mirr1half0_solid",
         MirrorBoxCut_up_sect0mirr1half0_solid, 
         Boxcut_down_sect0mirr1half0_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_down_sect0mirr1half0_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_down_sect0mirr1half0_rot;
   MirrorBoxCut_down_sect0mirr1half0_rot.rotateY( 0*rad  );
   MirrorBoxCut_down_sect0mirr1half0_rot.rotateZ( 0.2617993878*rad    );
   MirrorBoxCut_down_sect0mirr1half0_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //mirror_3_sector2_2 | htcc | htcc mirror 3, sector2 right |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | 8080ff |Operation:@ MirrorBoxCut_down_sect0mirr1half0 * phicut_sect0half0 |  0 | rohacell31 | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 1 
   r1   = MirrorBoxCut_down_sect0mirr1half0_rot;
   r2   = phicut_sect0half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorBoxCut_down_sect0mirr1half0_trans;
   t2   = phicut_sect0half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   mirror_3_sector2_2_solid = new G4IntersectionSolid(
         "mirror_3_sector2_2_solid",
         MirrorBoxCut_down_sect0mirr1half0_solid , 
         phicut_sect0half0_solid,
         &r_tot,
         t_tot
         );
   mirror_3_sector2_2_trans =  {0,0,0};
   mirror_3_sector2_2_rot =  G4RotationMatrix::IDENTITY ;
   mirror_3_sector2_2_rot.rotateY( 0*rad    );
   mirror_3_sector2_2_rot.rotateZ( 0.2617993878*rad );
   mirror_3_sector2_2_rot.rotateX( 1.392074611*rad);

   // ------------------------------------------------------------------------
   // Mirror 2
   //Mirror_sect0mirr2half0 | htcc |Ellipsoid defining mirror surface | 246.1895557*mm 918.7919301*mm 212.74075*mm |ordered: yzx 0*rad 0.2617993878*rad 1.35076309*rad | ff8080 | Ellipsoid | 1497.604*mm 1497.604*mm 1786.859*mm 0*mm 0*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Mirror_sect0mirr2half0_solid = new G4Ellipsoid("Mirror_sect0mirr2half0_solid",1497.604*mm,1497.604*mm,1786.859*mm,0*mm,0*mm); 
   G4ThreeVector    Mirror_sect0mirr2half0_trans(246.1895557*mm, 918.7919301*mm, 212.74075*mm);
   G4RotationMatrix Mirror_sect0mirr2half0_rot;
   Mirror_sect0mirr2half0_rot.rotateY( 0*rad  );
   Mirror_sect0mirr2half0_rot.rotateZ( 0.2617993878*rad );
   Mirror_sect0mirr2half0_rot.rotateX(  1.35076309*rad);
   // ------------------------------------------------------------------------
   //BarrelEllipseCut_sect0mirr2half0 | htcc |subtraction of barrel and ellipse |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | 00ff00 |Operation:@ Barrel_sect0half0 - Mirror_sect0mirr2half0 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = Barrel_sect0half0_rot ;
   r2   = Mirror_sect0mirr2half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = Barrel_sect0half0_trans;
   t2   = Mirror_sect0mirr2half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   BarrelEllipseCut_sect0mirr2half0_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr2half0_solid",
         Barrel_sect0half0_solid,
         Mirror_sect0mirr2half0_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector    BarrelEllipseCut_sect0mirr2half0_trans(0,0,0);
   G4RotationMatrix BarrelEllipseCut_sect0mirr2half0_rot;
   BarrelEllipseCut_sect0mirr2half0_rot.rotateY( 0*rad  );
   BarrelEllipseCut_sect0mirr2half0_rot.rotateZ( 0.2617993878*rad );
   BarrelEllipseCut_sect0mirr2half0_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //Boxcut_up_sect0mirr2half0 | htcc | upper plane cut | 139.5367254*mm 1161.472131*mm 2505.792765*mm |ordered: yzx 0*rad 0.2618003565*rad 0.5693996995*rad | 00ff00 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut_up_sect0mirr2half0_solid = new G4Box("Boxcut_up_sect0mirr2half0_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut_up_sect0mirr2half0_trans(139.5367254*mm,1161.472131*mm,2505.792765*mm);
   G4RotationMatrix Boxcut_up_sect0mirr2half0_rot;
   Boxcut_up_sect0mirr2half0_rot.rotateY( 0*rad  );
   Boxcut_up_sect0mirr2half0_rot.rotateZ( 0.2618003565*rad );
   Boxcut_up_sect0mirr2half0_rot.rotateX( 0.5693996995*rad);
   // ------------------------------------------------------------------------
   //Boxcut_down_sect0mirr2half0 | htcc | lower plane cut | -87.7977486*mm 42.96431111*mm 652.2335457*mm |ordered: yzx 0*rad 0.2617993647*rad 0.3460923731*rad | 00ff00 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut_down_sect0mirr2half0_solid = new G4Box("Boxcut_down_sect0mirr2half0_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut_down_sect0mirr2half0_trans(-87.7977486*mm,42.96431111*mm,652.2335457*mm);
   G4RotationMatrix Boxcut_down_sect0mirr2half0_rot;
   Boxcut_down_sect0mirr2half0_rot.rotateY( 0*rad  );
   Boxcut_down_sect0mirr2half0_rot.rotateZ( 0.2617993647*rad );
   Boxcut_down_sect0mirr2half0_rot.rotateX( 0.3460923731*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_up_sect0mirr2half0 | htcc |subtraction of upper box from (barrel - ellipse) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | 0000ff |Operation:@ BarrelEllipseCut_sect0mirr2half0 - Boxcut_up_sect0mirr2half0 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = BarrelEllipseCut_sect0mirr2half0_rot;
   r2   = Boxcut_up_sect0mirr2half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = BarrelEllipseCut_sect0mirr2half0_trans;
   t2   = Boxcut_up_sect0mirr2half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_up_sect0mirr2half0_solid = new G4SubtractionSolid(
         "MirrorBoxCut_up_sect0mirr2half0_solid",
         BarrelEllipseCut_sect0mirr2half0_solid, 
         Boxcut_up_sect0mirr2half0_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_up_sect0mirr2half0_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_up_sect0mirr2half0_rot;
   MirrorBoxCut_up_sect0mirr2half0_rot.rotateY( 0*rad  );
   MirrorBoxCut_up_sect0mirr2half0_rot.rotateZ( 0.2617993878*rad    );
   MirrorBoxCut_up_sect0mirr2half0_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_down_sect0mirr2half0 | htcc |subtraction of lower box from (barrel - ellipse - upper box) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | 0000ff |Operation:@ MirrorBoxCut_up_sect0mirr2half0 - Boxcut_down_sect0mirr2half0 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1    = MirrorBoxCut_up_sect0mirr2half0_rot;
   r2    = Boxcut_down_sect0mirr2half0_rot;
   r_tot = r2*(r1.inverse());
   t1    = MirrorBoxCut_up_sect0mirr2half0_trans;
   t2    = Boxcut_down_sect0mirr2half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_down_sect0mirr2half0_solid = new G4SubtractionSolid(
         "MirrorBoxCut_down_sect0mirr2half0_solid",
         MirrorBoxCut_up_sect0mirr2half0_solid, 
         Boxcut_down_sect0mirr2half0_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_down_sect0mirr2half0_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_down_sect0mirr2half0_rot;
   MirrorBoxCut_down_sect0mirr2half0_rot.rotateY( 0*rad  );
   MirrorBoxCut_down_sect0mirr2half0_rot.rotateZ( 0.2617993878*rad    );
   MirrorBoxCut_down_sect0mirr2half0_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //mirror_2_sector2_2 | htcc | htcc mirror 2, sector2 right |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | ff8080 |Operation:@ MirrorBoxCut_down_sect0mirr2half0 * phicut_sect0half0 |  0 | rohacell31 | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 2 
   r1   = MirrorBoxCut_down_sect0mirr2half0_rot;
   r2   = phicut_sect0half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorBoxCut_down_sect0mirr2half0_trans;
   t2   = phicut_sect0half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   mirror_2_sector2_2_solid = new G4IntersectionSolid(
         "mirror_2_sector2_2_solid",
         MirrorBoxCut_down_sect0mirr2half0_solid , 
         phicut_sect0half0_solid,
         &r_tot,
         t_tot
         );
   mirror_2_sector2_2_trans =  {0,0,0};
   mirror_2_sector2_2_rot =  G4RotationMatrix::IDENTITY ;
   mirror_2_sector2_2_rot.rotateY( 0*rad    );
   mirror_2_sector2_2_rot.rotateZ( 0.2617993878*rad );
   mirror_2_sector2_2_rot.rotateX( 1.392074611*rad);

   // ------------------------------------------------------------------------
   // Mirror 1
   //Mirror_sect0mirr3half0 | htcc |Ellipsoid defining mirror surface | 251.9535205*mm 940.3033397*mm 353.86945*mm |ordered: yzx 0*rad 0.2617993878*rad 1.222135079*rad | 8080ff | Ellipsoid | 1383.621*mm 1383.621*mm 1728.375*mm 0*mm 0*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Mirror_sect0mirr3half0_solid = new G4Ellipsoid("Mirror_sect0mirr3half0_solid",1497.604*mm,1497.604*mm,1786.859*mm,0*mm,0*mm); 
   G4ThreeVector    Mirror_sect0mirr3half0_trans(246.1895557*mm, 918.7919301*mm, 212.74075*mm);
   G4RotationMatrix Mirror_sect0mirr3half0_rot;
   Mirror_sect0mirr3half0_rot.rotateY( 0*rad  );
   Mirror_sect0mirr3half0_rot.rotateZ( 0.2617993878*rad );
   Mirror_sect0mirr3half0_rot.rotateX( 1.222135079*rad);
   // ------------------------------------------------------------------------
   //BarrelEllipseCut_sect0mirr3half0 | htcc |subtraction of barrel and ellipse |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | 00ff00 |Operation:@ Barrel_sect0half0 - Mirror_sect0mirr3half0 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = Barrel_sect0half0_rot ;
   r2   = Mirror_sect0mirr3half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = Barrel_sect0half0_trans;
   t2   = Mirror_sect0mirr3half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   BarrelEllipseCut_sect0mirr3half0_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr3half0_solid",
         Barrel_sect0half0_solid,
         Mirror_sect0mirr3half0_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector    BarrelEllipseCut_sect0mirr3half0_trans(0,0,0);
   G4RotationMatrix BarrelEllipseCut_sect0mirr3half0_rot;
   BarrelEllipseCut_sect0mirr3half0_rot.rotateY( 0*rad  );
   BarrelEllipseCut_sect0mirr3half0_rot.rotateZ( 0.2617993878*rad );
   BarrelEllipseCut_sect0mirr3half0_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //Boxcut_up_sect0mirr3half0 | htcc | upper plane cut | 87.7977486*mm 698.2956889*mm 2533.544454*mm |ordered: yzx 0*rad 0.2617993647*rad 0.3460923731*rad | ff0000 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut_up_sect0mirr3half0_solid = new G4Box("Boxcut_up_sect0mirr3half0_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut_up_sect0mirr3half0_trans(87.7977486*mm,698.2956889*mm,2533.544454*mm);
   G4RotationMatrix Boxcut_up_sect0mirr3half0_rot;
   Boxcut_up_sect0mirr3half0_rot.rotateY( 0*rad  );
   Boxcut_up_sect0mirr3half0_rot.rotateZ( 0.2618003565*rad );
   Boxcut_up_sect0mirr3half0_rot.rotateX( 0.3460923731*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_up_sect0mirr3half0 | htcc |subtraction of upper box from (barrel - ellipse) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | ffffff |Operation:@ BarrelEllipseCut_sect0mirr3half0 - Boxcut_up_sect0mirr3half0 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = BarrelEllipseCut_sect0mirr3half0_rot;
   r2   = Boxcut_up_sect0mirr3half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = BarrelEllipseCut_sect0mirr3half0_trans;
   t2   = Boxcut_up_sect0mirr3half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_up_sect0mirr3half0_solid = new G4SubtractionSolid(
         "MirrorBoxCut_up_sect0mirr3half0_solid",
         BarrelEllipseCut_sect0mirr3half0_solid, 
         Boxcut_up_sect0mirr3half0_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_up_sect0mirr3half0_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_up_sect0mirr3half0_rot;
   MirrorBoxCut_up_sect0mirr3half0_rot.rotateY( 0*rad  );
   MirrorBoxCut_up_sect0mirr3half0_rot.rotateZ( 0.2617993878*rad    );
   MirrorBoxCut_up_sect0mirr3half0_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //MirrorConeCut_sect0mirr3half0 | htcc |subtraction of lower cone from (barrel - ellipse - box) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | ffffff |Operation:@ MirrorBoxCut_up_sect0mirr3half0 - HTCC_InnerCutCone |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = MirrorBoxCut_up_sect0mirr3half0_rot;
   r2   = HTCC_InnerCutCone_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorBoxCut_up_sect0mirr3half0_trans;
   t2   = HTCC_InnerCutCone_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorConeCut_sect0mirr3half0_solid = new G4SubtractionSolid(
         "MirrorConeCut_sect0mirr3half0_solid",
         MirrorBoxCut_up_sect0mirr3half0_solid,
         HTCC_InnerCutCone_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorConeCut_sect0mirr3half0_trans(0,0,0);
   G4RotationMatrix MirrorConeCut_sect0mirr3half0_rot;
   MirrorConeCut_sect0mirr3half0_rot.rotateY( 0*rad  );
   MirrorConeCut_sect0mirr3half0_rot.rotateZ( 0.2617993878*rad    );
   MirrorConeCut_sect0mirr3half0_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //mirror_1_sector2_2 | htcc | htcc mirror 1, sector2 right |  0*mm 0*mm 0*mm |ordered: yzx 0*rad 0.2617993878*rad 1.392074611*rad | 8080ff |Operation:@ MirrorConeCut_sect0mirr3half0 * phicut_sect0half0 |  0 | rohacell31 | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 3 
   r1   = MirrorConeCut_sect0mirr3half0_rot;
   r2   = phicut_sect0half0_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorConeCut_sect0mirr3half0_trans;
   t2   = phicut_sect0half0_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   mirror_1_sector2_2_solid = new G4IntersectionSolid(
         "mirror_1_sector2_2_solid",
         MirrorConeCut_sect0mirr3half0_solid,
         phicut_sect0half0_solid,
         &r_tot,
         t_tot
         );
   mirror_1_sector2_2_trans =  {0,0,0};
   mirror_1_sector2_2_rot =  G4RotationMatrix::IDENTITY ;
   mirror_1_sector2_2_rot.rotateY( 0*rad    );
   mirror_1_sector2_2_rot.rotateZ( 0.2617993878*rad );
   mirror_1_sector2_2_rot.rotateX( 1.392074611*rad);

   // ------------------------------------------------------------------------
   // Other side
   //Barrel_sect0half1 | htcc |Barrel defining mirror back surface |  0 0 0 |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | ee99ff | Polycone |0.0*deg 360.0*deg 21*counts 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 0*mm 1435.3392*mm 1461.162*mm 1483.4626*mm 1502.4172*mm 1518.1665*mm 1530.8215*mm 1540.4676*mm 1547.168*mm 1550.9654*mm 1551.8836*mm 1549.9283*mm 1545.0874*mm 1537.3303*mm 1526.6071*mm 1512.8469*mm 1495.9556*mm 1475.8119*mm 1452.2631*mm 1425.1181*mm 1394.1386*mm 1359.0258*mm 360.01139*mm 423.09863*mm 486.18587*mm 549.27312*mm 612.36036*mm 675.4476*mm 738.53484*mm 801.62208*mm 864.70932*mm 927.79656*mm 990.8838*mm 1053.971*mm 1117.0583*mm 1180.1455*mm 1243.2328*mm 1306.32*mm 1369.4073*mm 1432.4945*mm 1495.5817*mm 1558.669*mm 1621.7562*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   // These solids are identical
   Barrel_sect0half1_solid = Barrel_sect0half0_solid;
   G4ThreeVector Barrel_sect0half1_trans(0,0,0);
   G4RotationMatrix Barrel_sect0half1_rot;
   Barrel_sect0half1_rot.rotateY(   0.0*rad);
   Barrel_sect0half1_rot.rotateZ(  -0.2617993878*rad );
   Barrel_sect0half1_rot.rotateX(   1.392074611*rad);
   // ------------------------------------------------------------------------
   //phicut_sect0half1 | htcc | Half-sector phi cut |  0*mm 0*mm 1500*mm |  0 0 0 | 443388 | Tube |0.0*mm 1500.0*mm 500.0*mm 1.570796327*rad 0.5225987756*rad | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   phicut_sect0half1_solid = new G4Tubs("phicut_sect0half1_solid", 0.0*mm, 1500.0*mm, 500.0*mm, 1.570796327*rad, 0.5225987756*rad );
   G4ThreeVector    phicut_sect0half1_trans(0,0,1500*mm);
   G4RotationMatrix phicut_sect0half1_rot;
   // ------------------------------------------------------------------------
   // Mirror 5
   //Mirror_sect0mirr0half1 | htcc |Ellipsoid defining mirror surface | -208.7414292*mm 779.0336195*mm -31.02335*mm |ordered: yzx 0*rad -0.2617993878*rad 1.609243305*rad | 80ff80 | Ellipsoid | 1728.673*mm 1728.673*mm 1907.810*mm 0*mm 0*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Mirror_sect0mirr0half1_solid = new G4Ellipsoid("Mirror_sect0mirr0half1_solid",1728.673*mm,1728.673*mm,1907.810*mm,0*mm,0*mm); 
   G4ThreeVector    Mirror_sect0mirr0half1_trans(-208.7414292*mm,779.0336195*mm,-31.02335*mm);
   G4RotationMatrix Mirror_sect0mirr0half1_rot;
   Mirror_sect0mirr0half1_rot.rotateY( 0*rad  );
   Mirror_sect0mirr0half1_rot.rotateZ( -0.2617993878*rad );
   Mirror_sect0mirr0half1_rot.rotateX(  1.609243305*rad );
   // ------------------------------------------------------------------------
   //BarrelEllipseCut_sect0mirr0half1 | htcc |subtraction of barrel and ellipse |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 00ff00 |Operation:@ Barrel_sect0half1 - Mirror_sect0mirr0half1 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = Barrel_sect0half1_rot ;
   r2   = Mirror_sect0mirr0half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = Barrel_sect0half1_trans ;
   t2   = Mirror_sect0mirr0half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   BarrelEllipseCut_sect0mirr0half1_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr0half1_solid",
         Barrel_sect0half1_solid,
         Mirror_sect0mirr0half1_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector    BarrelEllipseCut_sect0mirr0half1_trans(0,0,0);
   G4RotationMatrix BarrelEllipseCut_sect0mirr0half1_rot;
   BarrelEllipseCut_sect0mirr0half1_rot.rotateY( 0*rad  );
   BarrelEllipseCut_sect0mirr0half1_rot.rotateZ(-0.2617993878*rad );
   BarrelEllipseCut_sect0mirr0half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //Boxcut12_sect0mirr0half1 | htcc |Box defining dividing plane 1/2 | 184.7943978*mm 242.3741479*mm 978.1545571*mm |ordered: yzx 0*rad -0.2617980211*rad 0.7951867457*rad | ff0000 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut12_sect0mirr0half1_solid = new G4Box("Boxcut12_sect0mirr0half1_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut12_sect0mirr0half1_trans( 184.7943978*mm, 242.3741479*mm, 978.1545571*mm);
   G4RotationMatrix Boxcut12_sect0mirr0half1_rot;
   Boxcut12_sect0mirr0half1_rot.rotateY( 0*rad  );
   Boxcut12_sect0mirr0half1_rot.rotateZ( -0.2617980211*rad    );
   Boxcut12_sect0mirr0half1_rot.rotateX( 0.7951867457*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_sect0mirr0half1 | htcc |subtraction of box and (barrel - ellipse) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 00ff00 |Operation:@ BarrelEllipseCut_sect0mirr0half1 - Boxcut12_sect0mirr0half1 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = BarrelEllipseCut_sect0mirr0half1_rot;
   r2   = Boxcut12_sect0mirr0half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = BarrelEllipseCut_sect0mirr0half1_trans;
   t2   = Boxcut12_sect0mirr0half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_sect0mirr0half1_solid = new G4SubtractionSolid(
         "MirrorBoxCut_sect0mirr0half1_solid",
         BarrelEllipseCut_sect0mirr0half1_solid, 
         Boxcut12_sect0mirr0half1_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_sect0mirr0half1_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_sect0mirr0half1_rot;
   MirrorBoxCut_sect0mirr0half1_rot.rotateY( 0*rad  );
   MirrorBoxCut_sect0mirr0half1_rot.rotateZ(-0.2617993878*rad    );
   MirrorBoxCut_sect0mirr0half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //MirrorCylinderCut_sect0mirr0half1 | htcc |subtraction of cylinder and (box - (barrel - ellipse)) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 00ff00 |Operation:@ MirrorBoxCut_sect0mirr0half1 - HTCC_OuterCutCylinder |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = MirrorBoxCut_sect0mirr0half1_rot;
   r2   = HTCC_OuterCutCylinder_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorBoxCut_sect0mirr0half1_trans;
   t2   = HTCC_OuterCutCylinder_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorCylinderCut_sect0mirr0half1_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr0half1_solid",
         MirrorBoxCut_sect0mirr0half1_solid,
         HTCC_OuterCutCylinder_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorCylinderCut_sect0mirr0half1_trans(0,0,0);
   G4RotationMatrix MirrorCylinderCut_sect0mirr0half1_rot;
   MirrorCylinderCut_sect0mirr0half1_rot.rotateY( 0*rad    );
   MirrorCylinderCut_sect0mirr0half1_rot.rotateZ(-0.2617993878*rad );
   MirrorCylinderCut_sect0mirr0half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //mirror_4_sector3_1 | htcc | htcc mirror 4, sector3 left |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 80ff80 |Operation:@ MirrorCylinderCut_sect0mirr0half1 * phicut_sect0half1 |  0 | rohacell31 | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 10 
   r1   = MirrorCylinderCut_sect0mirr0half1_rot;
   r2   = phicut_sect0half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorCylinderCut_sect0mirr0half1_trans;
   t2   = phicut_sect0half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   mirror_4_sector3_1_solid = new G4IntersectionSolid(
         "mirror_4_sector3_1_solid",
         MirrorCylinderCut_sect0mirr0half1_solid , 
         phicut_sect0half1_solid,
         &r_tot,
         t_tot
         );
   mirror_4_sector3_1_rot =  G4RotationMatrix::IDENTITY ;
   mirror_4_sector3_1_rot.rotateY( 0*rad    );
   mirror_4_sector3_1_rot.rotateZ( -0.2617993878*rad );
   mirror_4_sector3_1_rot.rotateX( 1.392074611*rad);
   // Mirror 6
   // ------------------------------------------------------------------------
   //Mirror_sect0mirr1half1 | htcc |Ellipsoid defining mirror surface | -231.4737648*mm 863.8718507*mm 81.66635*mm |ordered: yzx 0*rad -0.2617993878*rad 1.479734815*rad | f0f0f0 | Ellipsoid | 1612.998*mm 1612.998*mm 1846.155*mm 0*mm 0*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Mirror_sect0mirr1half1_solid = new G4Ellipsoid("Mirror_sect0mirr1half0_solid",1612.998*mm,1612.998*mm,1846.155*mm,0*mm,0*mm); 
   G4ThreeVector    Mirror_sect0mirr1half1_trans(-231.4737648*mm,863.8718507*mm,81.66635*mm);
   G4RotationMatrix Mirror_sect0mirr1half1_rot;
   Mirror_sect0mirr1half1_rot.rotateY( 0*rad  );
   Mirror_sect0mirr1half1_rot.rotateZ(-0.2617993878*rad );
   Mirror_sect0mirr1half1_rot.rotateX( 1.479734815*rad);
   // ------------------------------------------------------------------------
   //BarrelEllipseCut_sect0mirr1half1 | htcc |subtraction of barrel and ellipse |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 00ff00 |Operation:@ Barrel_sect0half1 - Mirror_sect0mirr1half1 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = Barrel_sect0half1_rot ;
   r2   = Mirror_sect0mirr1half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = Barrel_sect0half1_trans;
   t2   = Mirror_sect0mirr1half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   BarrelEllipseCut_sect0mirr1half1_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr1half1_solid",
         Barrel_sect0half1_solid,
         Mirror_sect0mirr1half1_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector    BarrelEllipseCut_sect0mirr1half1_trans(0,0,0);
   G4RotationMatrix BarrelEllipseCut_sect0mirr1half1_rot;
   BarrelEllipseCut_sect0mirr1half1_rot.rotateY( 0*rad  );
   BarrelEllipseCut_sect0mirr1half1_rot.rotateZ(-0.2617993878*rad );
   BarrelEllipseCut_sect0mirr1half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //Boxcut_up_sect0mirr1half1 | htcc | upper plane cut | -184.7943978*mm 1621.705852*mm 2378.357443*mm |ordered: yzx 0*rad -0.2617980211*rad 0.7951867457*rad | 00ff00 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut_up_sect0mirr1half1_solid = new G4Box("Boxcut_up_sect0mirr1half0_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut_up_sect0mirr1half1_trans(-184.7943978*mm,1621.705852*mm,2378.357443*mm);
   G4RotationMatrix Boxcut_up_sect0mirr1half1_rot;
   Boxcut_up_sect0mirr1half1_rot.rotateY( 0*rad  );
   Boxcut_up_sect0mirr1half1_rot.rotateZ(-0.2617980211*rad    );
   Boxcut_up_sect0mirr1half1_rot.rotateX( 0.7951867457*rad);
   // ------------------------------------------------------------------------
   //Boxcut_down_sect0mirr1half1 | htcc | lower plane cut | 139.5367254*mm 119.9598692*mm 821.4432351*mm |ordered: yzx 0*rad -0.2618003565*rad 0.5693996995*rad | 00ff00 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut_down_sect0mirr1half1_solid = new G4Box("Boxcut_down_sect0mirr1half0_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut_down_sect0mirr1half1_trans( 139.5367254*mm,119.9598692*mm,821.4432351*mm);
   G4RotationMatrix Boxcut_down_sect0mirr1half1_rot;
   Boxcut_down_sect0mirr1half1_rot.rotateY( 0*rad  );
   Boxcut_down_sect0mirr1half1_rot.rotateZ(-0.2617980211*rad    );
   Boxcut_down_sect0mirr1half1_rot.rotateX( 0.5693996995*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_up_sect0mirr1half1 | htcc |subtraction of upper box from (barrel - ellipse) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 0000ff |Operation:@ BarrelEllipseCut_sect0mirr1half1 - Boxcut_up_sect0mirr1half1 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = BarrelEllipseCut_sect0mirr1half1_rot;
   r2   = Boxcut_up_sect0mirr1half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = BarrelEllipseCut_sect0mirr1half1_trans;
   t2   = Boxcut_up_sect0mirr1half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_up_sect0mirr1half1_solid = new G4SubtractionSolid(
         "MirrorBoxCut_up_sect0mirr1half1_solid",
         BarrelEllipseCut_sect0mirr1half1_solid, 
         Boxcut_up_sect0mirr1half1_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_up_sect0mirr1half1_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_up_sect0mirr1half1_rot;
   MirrorBoxCut_up_sect0mirr1half1_rot.rotateY( 0*rad  );
   MirrorBoxCut_up_sect0mirr1half1_rot.rotateZ(-0.2617993878*rad    );
   MirrorBoxCut_up_sect0mirr1half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_down_sect0mirr1half1 | htcc |subtraction of lower box from (barrel - ellipse - upper box) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 0000ff |Operation:@ MirrorBoxCut_up_sect0mirr1half1 - Boxcut_down_sect0mirr1half1 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = MirrorBoxCut_up_sect0mirr1half1_rot;
   r2   = Boxcut_down_sect0mirr1half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorBoxCut_up_sect0mirr1half1_trans;
   t2   = Boxcut_down_sect0mirr1half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_down_sect0mirr1half1_solid = new G4SubtractionSolid(
         "MirrorBoxCut_down_sect0mirr1half1_solid",
         MirrorBoxCut_up_sect0mirr1half1_solid, 
         Boxcut_down_sect0mirr1half1_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_down_sect0mirr1half1_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_down_sect0mirr1half1_rot;
   MirrorBoxCut_down_sect0mirr1half1_rot.rotateY( 0*rad  );
   MirrorBoxCut_down_sect0mirr1half1_rot.rotateZ(-0.2617993878*rad    );
   MirrorBoxCut_down_sect0mirr1half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //mirror_3_sector3_1 | htcc | htcc mirror 3, sector3 left |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | f0f0f0 |Operation:@ MirrorBoxCut_down_sect0mirr1half1 * phicut_sect0half1 |  0 | rohacell31 | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 11 
   r1   = MirrorBoxCut_down_sect0mirr1half1_rot;
   r2   = phicut_sect0half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorBoxCut_down_sect0mirr1half1_trans;
   t2   = phicut_sect0half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   mirror_3_sector3_1_solid = new G4IntersectionSolid(
         "mirror_3_sector3_1_solid",
         MirrorBoxCut_down_sect0mirr1half1_solid , 
         phicut_sect0half1_solid,
         &r_tot,
         t_tot
         );
   mirror_3_sector3_1_trans =  {0,0,0};
   mirror_3_sector3_1_rot =  G4RotationMatrix::IDENTITY ;
   mirror_3_sector3_1_rot.rotateY( 0*rad    );
   mirror_3_sector3_1_rot.rotateZ(-0.2617993878*rad );
   mirror_3_sector3_1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //Mirror_sect0mirr2half1 | htcc |Ellipsoid defining mirror surface | -246.1895557*mm 918.7919301*mm 212.74075*mm |ordered: yzx 0*rad -0.2617993878*rad 1.35076309*rad | 80ff80 | Ellipsoid | 1497.604*mm 1497.604*mm 1786.859*mm 0*mm 0*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Mirror_sect0mirr2half1_solid = new G4Ellipsoid("Mirror_sect0mirr2half0_solid",1497.604*mm,1497.604*mm,1786.859*mm,0*mm,0*mm); 
   G4ThreeVector    Mirror_sect0mirr2half1_trans(-246.1895557*mm, 918.7919301*mm, 212.74075*mm);
   G4RotationMatrix Mirror_sect0mirr2half1_rot;
   Mirror_sect0mirr2half1_rot.rotateY( 0*rad  );
   Mirror_sect0mirr2half1_rot.rotateZ(-0.2617993878*rad );
   Mirror_sect0mirr2half1_rot.rotateX(  1.35076309*rad);
   // ------------------------------------------------------------------------
   //BarrelEllipseCut_sect0mirr2half1 | htcc |subtraction of barrel and ellipse |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 00ff00 |Operation:@ Barrel_sect0half1 - Mirror_sect0mirr2half1 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = Barrel_sect0half1_rot ;
   r2   = Mirror_sect0mirr2half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = Barrel_sect0half1_trans;
   t2   = Mirror_sect0mirr2half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   BarrelEllipseCut_sect0mirr2half1_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr2half1_solid",
         Barrel_sect0half1_solid,
         Mirror_sect0mirr2half1_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector    BarrelEllipseCut_sect0mirr2half1_trans(0,0,0);
   G4RotationMatrix BarrelEllipseCut_sect0mirr2half1_rot;
   BarrelEllipseCut_sect0mirr2half1_rot.rotateY( 0*rad  );
   BarrelEllipseCut_sect0mirr2half1_rot.rotateZ(-0.2617993878*rad );
   BarrelEllipseCut_sect0mirr2half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //Boxcut_up_sect0mirr2half1 | htcc | upper plane cut | -139.5367254*mm 1161.472131*mm 2505.792765*mm |ordered: yzx 0*rad -0.2618003565*rad 0.5693996995*rad | 00ff00 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut_up_sect0mirr2half1_solid = new G4Box("Boxcut_up_sect0mirr2half0_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut_up_sect0mirr2half1_trans(-139.5367254*mm,1161.472131*mm,2505.792765*mm);
   G4RotationMatrix Boxcut_up_sect0mirr2half1_rot;
   Boxcut_up_sect0mirr2half1_rot.rotateY( 0*rad  );
   Boxcut_up_sect0mirr2half1_rot.rotateZ(-0.2618003565*rad );
   Boxcut_up_sect0mirr2half1_rot.rotateX( 0.5693996995*rad);
   // ------------------------------------------------------------------------
   //Boxcut_down_sect0mirr2half1 | htcc | lower plane cut | 87.7977486*mm 42.96431111*mm 652.2335457*mm |ordered: yzx 0*rad -0.2617993647*rad 0.3460923731*rad | 00ff00 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut_down_sect0mirr2half1_solid = new G4Box("Boxcut_down_sect0mirr2half0_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut_down_sect0mirr2half1_trans(87.7977486*mm,42.96431111*mm,652.2335457*mm);
   G4RotationMatrix Boxcut_down_sect0mirr2half1_rot;
   Boxcut_down_sect0mirr2half1_rot.rotateY( 0*rad  );
   Boxcut_down_sect0mirr2half1_rot.rotateZ(-0.2617993647*rad );
   Boxcut_down_sect0mirr2half1_rot.rotateX( 0.3460923731*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_up_sect0mirr2half1 | htcc |subtraction of upper box from (barrel - ellipse) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 0000ff |Operation:@ BarrelEllipseCut_sect0mirr2half1 - Boxcut_up_sect0mirr2half1 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = BarrelEllipseCut_sect0mirr2half1_rot;
   r2   = Boxcut_up_sect0mirr2half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = BarrelEllipseCut_sect0mirr2half1_trans;
   t2   = Boxcut_up_sect0mirr2half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_up_sect0mirr2half1_solid = new G4SubtractionSolid(
         "MirrorBoxCut_up_sect0mirr2half1_solid",
         BarrelEllipseCut_sect0mirr2half1_solid, 
         Boxcut_up_sect0mirr2half1_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_up_sect0mirr2half1_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_up_sect0mirr2half1_rot;
   MirrorBoxCut_up_sect0mirr2half1_rot.rotateY( 0*rad  );
   MirrorBoxCut_up_sect0mirr2half1_rot.rotateZ(-0.2617993878*rad    );
   MirrorBoxCut_up_sect0mirr2half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_down_sect0mirr2half1 | htcc |subtraction of lower box from (barrel - ellipse - upper box) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 0000ff |Operation:@ MirrorBoxCut_up_sect0mirr2half1 - Boxcut_down_sect0mirr2half1 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1    = MirrorBoxCut_up_sect0mirr2half1_rot;
   r2    = Boxcut_down_sect0mirr2half1_rot;
   r_tot = r2*(r1.inverse());
   t1    = MirrorBoxCut_up_sect0mirr2half1_trans;
   t2    = Boxcut_down_sect0mirr2half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_down_sect0mirr2half1_solid = new G4SubtractionSolid(
         "MirrorBoxCut_down_sect0mirr2half1_solid",
         MirrorBoxCut_up_sect0mirr2half1_solid, 
         Boxcut_down_sect0mirr2half1_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_down_sect0mirr2half1_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_down_sect0mirr2half1_rot;
   MirrorBoxCut_down_sect0mirr2half1_rot.rotateY( 0*rad  );
   MirrorBoxCut_down_sect0mirr2half1_rot.rotateZ(-0.2617993878*rad    );
   MirrorBoxCut_down_sect0mirr2half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //mirror_2_sector3_1 | htcc | htcc mirror 2, sector3 left |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 80ff80 |Operation:@ MirrorBoxCut_down_sect0mirr2half1 * phicut_sect0half1 |  0 | rohacell31 | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 12 
   r1   = MirrorBoxCut_down_sect0mirr2half1_rot;
   r2   = phicut_sect0half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorBoxCut_down_sect0mirr2half1_trans;
   t2   = phicut_sect0half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   mirror_2_sector3_1_solid = new G4IntersectionSolid(
         "mirror_2_sector3_1_solid",
         MirrorBoxCut_down_sect0mirr2half1_solid , 
         phicut_sect0half1_solid,
         &r_tot,
         t_tot
         );
   mirror_2_sector3_1_trans =  {0,0,0};
   mirror_2_sector3_1_rot =  G4RotationMatrix::IDENTITY ;
   mirror_2_sector3_1_rot.rotateY( 0*rad    );
   mirror_2_sector3_1_rot.rotateZ(-0.2617993878*rad );
   mirror_2_sector3_1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //Mirror_sect0mirr3half1 | htcc |Ellipsoid defining mirror surface | -251.9535205*mm 940.3033397*mm 353.86945*mm |ordered: yzx 0*rad -0.2617993878*rad 1.222135079*rad | f0f0f0 | Ellipsoid | 1383.621*mm 1383.621*mm 1728.375*mm 0*mm 0*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Mirror_sect0mirr3half1_solid = new G4Ellipsoid("Mirror_sect0mirr3half1_solid",1497.604*mm,1497.604*mm,1786.859*mm,0*mm,0*mm); 
   G4ThreeVector    Mirror_sect0mirr3half1_trans(-246.1895557*mm, 918.7919301*mm, 212.74075*mm);
   G4RotationMatrix Mirror_sect0mirr3half1_rot;
   Mirror_sect0mirr3half1_rot.rotateY( 0*rad  );
   Mirror_sect0mirr3half1_rot.rotateZ(-0.2617993878*rad );
   Mirror_sect0mirr3half1_rot.rotateX( 1.222135079*rad);
   // ------------------------------------------------------------------------
   //BarrelEllipseCut_sect0mirr3half1 | htcc |subtraction of barrel and ellipse |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | 00ff00 |Operation:@ Barrel_sect0half1 - Mirror_sect0mirr3half1 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = Barrel_sect0half1_rot ;
   r2   = Mirror_sect0mirr3half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = Barrel_sect0half1_trans;
   t2   = Mirror_sect0mirr3half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   BarrelEllipseCut_sect0mirr3half1_solid = new G4SubtractionSolid(
         "BarrelEllipseCut_sect0mirr3half1_solid",
         Barrel_sect0half1_solid,
         Mirror_sect0mirr3half1_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector    BarrelEllipseCut_sect0mirr3half1_trans(0,0,0);
   G4RotationMatrix BarrelEllipseCut_sect0mirr3half1_rot;
   BarrelEllipseCut_sect0mirr3half1_rot.rotateY( 0*rad  );
   BarrelEllipseCut_sect0mirr3half1_rot.rotateZ(-0.2617993878*rad );
   BarrelEllipseCut_sect0mirr3half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //Boxcut_up_sect0mirr3half1 | htcc | upper plane cut | -87.7977486*mm 698.2956889*mm 2533.544454*mm |ordered: yzx 0*rad -0.2617993647*rad 0.3460923731*rad | ff0000 | Box |  1000*mm 1000*mm 1000*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   Boxcut_up_sect0mirr3half1_solid = new G4Box("Boxcut_up_sect0mirr3half1_solid",1000*mm, 1000*mm, 1000*mm);
   G4ThreeVector    Boxcut_up_sect0mirr3half1_trans(-87.7977486*mm,698.2956889*mm,2533.544454*mm);
   G4RotationMatrix Boxcut_up_sect0mirr3half1_rot;
   Boxcut_up_sect0mirr3half1_rot.rotateY( 0*rad  );
   Boxcut_up_sect0mirr3half1_rot.rotateZ(-0.2618003565*rad );
   Boxcut_up_sect0mirr3half1_rot.rotateX( 0.3460923731*rad);
   // ------------------------------------------------------------------------
   //MirrorBoxCut_up_sect0mirr3half1 | htcc |subtraction of upper box from (barrel - ellipse) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | ffffff |Operation:@ BarrelEllipseCut_sect0mirr3half1 - Boxcut_up_sect0mirr3half1 |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = BarrelEllipseCut_sect0mirr3half1_rot;
   r2   = Boxcut_up_sect0mirr3half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = BarrelEllipseCut_sect0mirr3half1_trans;
   t2   = Boxcut_up_sect0mirr3half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorBoxCut_up_sect0mirr3half1_solid = new G4SubtractionSolid(
         "MirrorBoxCut_up_sect0mirr3half1_solid",
         BarrelEllipseCut_sect0mirr3half1_solid, 
         Boxcut_up_sect0mirr3half1_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorBoxCut_up_sect0mirr3half1_trans(0,0,0);
   G4RotationMatrix MirrorBoxCut_up_sect0mirr3half1_rot;
   MirrorBoxCut_up_sect0mirr3half1_rot.rotateY( 0*rad  );
   MirrorBoxCut_up_sect0mirr3half1_rot.rotateZ(-0.2617993878*rad    );
   MirrorBoxCut_up_sect0mirr3half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //MirrorConeCut_sect0mirr3half1 | htcc |subtraction of lower cone from (barrel - ellipse - box) |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | ffffff |Operation:@ MirrorBoxCut_up_sect0mirr3half1 - HTCC_InnerCutCone |  0 | Component | no | 1 | 1 | 1 | 1 | 0 | no | no |  no 
   r1   = MirrorBoxCut_up_sect0mirr3half1_rot;
   r2   = HTCC_InnerCutCone_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorBoxCut_up_sect0mirr3half1_trans;
   t2   = HTCC_InnerCutCone_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   MirrorConeCut_sect0mirr3half1_solid = new G4SubtractionSolid(
         "MirrorConeCut_sect0mirr3half1_solid",
         MirrorBoxCut_up_sect0mirr3half1_solid,
         HTCC_InnerCutCone_solid,
         &r_tot,
         t_tot
         );
   G4ThreeVector MirrorConeCut_sect0mirr3half1_trans(0,0,0);
   G4RotationMatrix MirrorConeCut_sect0mirr3half1_rot;
   MirrorConeCut_sect0mirr3half1_rot.rotateY( 0*rad  );
   MirrorConeCut_sect0mirr3half1_rot.rotateZ(-0.2617993878*rad    );
   MirrorConeCut_sect0mirr3half1_rot.rotateX( 1.392074611*rad);
   // ------------------------------------------------------------------------
   //mirror_1_sector3_1 | htcc | htcc mirror 1, sector3 left |  0*mm 0*mm 0*mm |ordered: yzx 0*rad -0.2617993878*rad 1.392074611*rad | f0f0f0 |Operation:@ MirrorConeCut_sect0mirr3half1 * phicut_sect0half1 |  0 | rohacell31 | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 13 
   r1   = MirrorConeCut_sect0mirr3half1_rot;
   r2   = phicut_sect0half1_rot;
   r_tot = r2*(r1.inverse());
   t1   = MirrorConeCut_sect0mirr3half1_trans;
   t2   = phicut_sect0half1_trans;
   t_tot = (t2-t1);
   t_tot *= r1;
   mirror_1_sector3_1_solid = new G4IntersectionSolid(
         "mirror_1_sector3_1_solid",
         MirrorConeCut_sect0mirr3half1_solid,
         phicut_sect0half1_solid,
         &r_tot,
         t_tot
         );
   mirror_1_sector3_1_trans =  {0,0,0};
   mirror_1_sector3_1_rot =  G4RotationMatrix::IDENTITY ;
   mirror_1_sector3_1_rot.rotateY( 0*rad    );
   mirror_1_sector3_1_rot.rotateZ(-0.2617993878*rad );
   mirror_1_sector3_1_rot.rotateX( 1.392074611*rad);

   // ------------------------------------------------------------------------

//pmt_4_sector3_1 | htcc | htcc pmt 4, sector3 left | -417.60646*mm 1558.5278*mm -63.469128*mm | 162.04762*deg 4.7212992*deg 0*deg | 80ff80 | Tube | 0*mm 55*mm 1.5*mm 0*deg 360*deg | HTCCPMTquartz | no | 1 | 1 | 1 | 1 | 1 | htcc | htcc | sector manual 3 ring manual 4 half manual 1
   pmt_4_sector3_1_solid = new G4Tubs("pmt_4_sector3_1_solid", 0*mm,55*mm,1.5*mm,0*deg,360*deg );
   pmt_4_sector3_1_trans = {-417.60646*mm,1558.5278*mm,-63.469128*mm};
   pmt_4_sector3_1_rot = G4RotationMatrix::IDENTITY;
   pmt_4_sector3_1_rot.rotateX( 162.04762*deg    );
   pmt_4_sector3_1_rot.rotateY( 4.7212992*deg );
   pmt_4_sector3_1_rot.rotateZ( 0*deg);
   //pmt_4_sector2_2 | htcc | htcc pmt 4, sector2 right | 417.60646*mm 1558.5278*mm -63.469128*mm | 162.04762*deg -4.7212992*deg 0*deg | ff8080 | Tube | 0*mm 55*mm 1.5*mm 0*deg 360*deg | HTCCPMTquartz | no | 1 | 1 | 1 | 1 | 1 | htcc | htcc | sector manual 2 ring manual 4 half manual 2
   pmt_4_sector2_2_solid = new G4Tubs("pmt_4_sector2_2_solid", 0*mm,55*mm,1.5*mm,0*deg,360*deg );
   pmt_4_sector2_2_trans = {417.60646*mm,1558.5278*mm,-63.469128*mm};
   pmt_4_sector2_2_rot = G4RotationMatrix::IDENTITY;
   pmt_4_sector2_2_rot.rotateX( 162.04762*deg    );
   pmt_4_sector2_2_rot.rotateY( -4.7212992*deg );
   pmt_4_sector2_2_rot.rotateZ( 0*deg);

   // ----------------------------------------------------------------
   // Winston cone
   //wc_4_sector3_1inner | htcc | htcc wc 4, sector3 leftinner | 0 0 0 | 0 0 0 | 999999 | Paraboloid | 125.25*mm 51.338961*mm 76.723114*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_4_sector3_1inner_solid = new G4Paraboloid("wc_4_sector3_1inner_solid",125.25*mm,51.338961*mm,76.723114*mm);
   G4ThreeVector    wc_4_sector3_1inner_trans = {0,0,0};
   G4RotationMatrix wc_4_sector3_1inner_rot   = G4RotationMatrix::IDENTITY;

   //wc_4_sector3_1outer | htcc | htcc wc 4, sector3 leftouter | -409.64307*mm 1528.8081*mm 28.25812*mm | -17.952376*deg -4.7212992*deg 0*deg | 999999 | Paraboloid | 95.25*mm 56*mm 75.1426*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_4_sector3_1outer_solid = new G4Paraboloid("wc_4_sector3_1outer_solid",95.25*mm,56*mm,75.1426*mm);
   G4ThreeVector    wc_4_sector3_1outer_trans = {-409.64307*mm,1528.8081*mm,28.25812*mm};
   G4RotationMatrix wc_4_sector3_1outer_rot   = G4RotationMatrix::IDENTITY;
   wc_4_sector3_1outer_rot.rotateX( -17.952376*deg    );
   wc_4_sector3_1outer_rot.rotateY(  -4.7212992*deg );
   wc_4_sector3_1outer_rot.rotateZ( 0*deg);

   //wc_4_sector3_1 | htcc | htcc wc 4, sector3 left | -409.64307*mm 1528.8081*mm 28.25812*mm | -17.952376*deg -4.7212992*deg 0*deg | 80ff802 | Operation: wc_4_sector3_1outer - wc_4_sector3_1inner | 0 | G4_Al | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 1
   wc_4_sector3_1_solid = new G4SubtractionSolid(
         "wc_4_sector3_1_solid",
         wc_4_sector3_1outer_solid,
         wc_4_sector3_1inner_solid,
         0,
         wc_4_sector3_1inner_trans
         );
   wc_4_sector3_1_trans = {-409.64307*mm,1528.8081*mm,28.25812*mm};
   wc_4_sector3_1_rot = G4RotationMatrix::IDENTITY;
   wc_4_sector3_1_rot.rotateX( -17.952376*deg    );
   wc_4_sector3_1_rot.rotateY(  -4.7212992*deg );
   wc_4_sector3_1_rot.rotateZ( 0*deg);

   // ----------------------------------------------------------------
   //wc_4_sector2_2inner | htcc | htcc wc 4, sector2 rightinner | 0 0 0 | 0 0 0 | 999999 | Paraboloid | 125.25*mm 51.338961*mm 76.723114*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_4_sector2_2inner_solid = new G4Paraboloid("wc_4_sector2_2inner_solid",125.25*mm,51.338961*mm,76.723114*mm);
   G4ThreeVector    wc_4_sector2_2inner_trans = {0,0,0};
   G4RotationMatrix wc_4_sector2_2inner_rot   = G4RotationMatrix::IDENTITY;

   //wc_4_sector2_2outer | htcc | htcc wc 4, sector2 rightouter | 409.64307*mm 1528.8081*mm 28.25812*mm | -17.952376*deg 4.7212992*deg 0*deg | 999999 | Paraboloid | 95.25*mm 56*mm 75.1426*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_4_sector2_2outer_solid = new G4Paraboloid("wc_4_sector2_2outer_solid",95.25*mm,56*mm,75.1426*mm);
   G4ThreeVector    wc_4_sector2_2outer_trans = {-409.64307*mm,1528.8081*mm,28.25812*mm};
   G4RotationMatrix wc_4_sector2_2outer_rot   = G4RotationMatrix::IDENTITY;
   wc_4_sector2_2outer_rot.rotateX( -17.952376*deg    );
   wc_4_sector2_2outer_rot.rotateY(  4.7212992*deg );
   wc_4_sector2_2outer_rot.rotateZ( 0*deg);

   //wc_4_sector2_2 | htcc | htcc wc 4, sector2 right | 409.64307*mm 1528.8081*mm 28.25812*mm | -17.952376*deg 4.7212992*deg 0*deg | ff80802 | Operation: wc_4_sector2_2outer - wc_4_sector2_2inner | 0 | G4_Al | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 1
   wc_4_sector2_2_solid = new G4SubtractionSolid(
         "wc_4_sector2_2_solid",
         wc_4_sector2_2outer_solid,
         wc_4_sector2_2inner_solid,
         0,
         wc_4_sector2_2inner_trans
         );
   wc_4_sector2_2_trans = {409.64307*mm,1528.8081*mm,28.25812*mm};
   wc_4_sector2_2_rot = G4RotationMatrix::IDENTITY;
   wc_4_sector2_2_rot.rotateX( -17.952376*deg    );
   wc_4_sector2_2_rot.rotateY(  4.7212992*deg );
   wc_4_sector2_2_rot.rotateZ( 0*deg);

   //---------------------------------------------------------
   // row 3
   //pmt_3_sector3_1 | htcc | htcc pmt 3, sector3 left | -463.16694*mm 1728.5611*mm 162.09428*mm | 146.59016*deg 8.392807*deg 0*deg | f0f0f0 | Tube | 0*mm 55*mm 1.5*mm 0*deg 360*deg | HTCCPMTquartz | no | 1 | 1 | 1 | 1 | 1 | htcc | htcc | sector manual 3 ring manual 3 half manual 1
   pmt_3_sector3_1_solid = pmt_4_sector3_1_solid;
   pmt_3_sector3_1_trans = {-463.16694*mm, 1728.5611*mm, 162.09428*mm};
   pmt_3_sector3_1_rot = G4RotationMatrix::IDENTITY;
   pmt_3_sector3_1_rot.rotateX( 146.59016*deg    );
   pmt_3_sector3_1_rot.rotateY( 8.392807*deg );
   pmt_3_sector3_1_rot.rotateZ( 0*deg);
   //pmt_3_sector2_2 | htcc | htcc pmt 3, sector2 right | 463.16694*mm 1728.5611*mm 162.09428*mm | 146.59016*deg -8.392807*deg 0*deg | 8080ff | Tube | 0*mm 55*mm 1.5*mm 0*deg 360*deg | HTCCPMTquartz | no | 1 | 1 | 1 | 1 | 1 | htcc | htcc | sector manual 2 ring manual 3 half manual 2
   pmt_3_sector2_2_solid = pmt_4_sector3_1_solid;
   pmt_3_sector2_2_trans = {463.16694*mm, 1728.5611*mm, 162.09428*mm};
   pmt_3_sector2_2_rot = G4RotationMatrix::IDENTITY;
   pmt_3_sector2_2_rot.rotateX( 146.59016*deg    );
   pmt_3_sector2_2_rot.rotateY(-8.392807*deg );
   pmt_3_sector2_2_rot.rotateZ( 0*deg);

   // -------------------------------------
   //wc_3_sector3_1inner | htcc | htcc wc 3, sector3 leftinner | 0 0 0 | 0 0 0 | 999999 | Paraboloid | 125.25*mm 51.338961*mm 76.723114*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_3_sector3_1inner_solid = wc_4_sector3_1inner_solid ;
   G4ThreeVector    wc_3_sector3_1inner_trans = {0,0,0};
   G4RotationMatrix wc_3_sector3_1inner_rot   = G4RotationMatrix::IDENTITY;

   //wc_3_sector3_1outer | htcc | htcc wc 3, sector3 leftouter | -449.04543*mm 1675.8588*mm 241.99168*mm | -33.409841*deg -8.392807*deg 0*deg | 999999 | Paraboloid | 95.25*mm 56*mm 75.1426*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_3_sector3_1outer_solid = wc_4_sector3_1outer_solid;
   G4ThreeVector    wc_3_sector3_1outer_trans = {-449.04543*mm,1675.8588*mm,241.99168*mm};
   G4RotationMatrix wc_3_sector3_1outer_rot   = G4RotationMatrix::IDENTITY;
   wc_3_sector3_1outer_rot.rotateX( -33.409841*deg    );
   wc_3_sector3_1outer_rot.rotateY(   -8.392807*deg );
   wc_3_sector3_1outer_rot.rotateZ( 0*deg);

   //wc_3_sector3_1 | htcc | htcc wc 3, sector3 left | -449.04543*mm 1675.8588*mm 241.99168*mm | -33.409841*deg -8.392807*deg 0*deg | f0f0f02 | Operation: wc_3_sector3_1outer - wc_3_sector3_1inner | 0 | G4_Al | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 1
   wc_3_sector3_1_solid = new G4SubtractionSolid(
         "wc_3_sector3_1_solid",
         wc_3_sector3_1outer_solid,
         wc_3_sector3_1inner_solid,
         0,
         wc_3_sector3_1inner_trans
         );
   wc_3_sector3_1_trans = {-449.04543*mm,1675.8588*mm,241.99168*mm};
   wc_3_sector3_1_rot = G4RotationMatrix::IDENTITY;
   wc_3_sector3_1_rot.rotateX( -33.409841*deg    );
   wc_3_sector3_1_rot.rotateY(   -8.392807*deg );
   wc_3_sector3_1_rot.rotateZ( 0*deg);

   // ----------------------------------------------------------------
   //wc_3_sector2_2inner | htcc | htcc wc 3, sector2 rightinner | 0 0 0 | 0 0 0 | 999999 | Paraboloid | 125.25*mm 51.338961*mm 76.723114*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_3_sector2_2inner_solid = wc_4_sector2_2inner_solid;
   G4ThreeVector    wc_3_sector2_2inner_trans = {0,0,0};
   G4RotationMatrix wc_3_sector2_2inner_rot   = G4RotationMatrix::IDENTITY;

   //wc_3_sector2_2outer | htcc | htcc wc 3, sector2 rightouter | 449.04543*mm 1675.8588*mm 241.99168*mm | -33.409841*deg 8.392807*deg 0*deg | 999999 | Paraboloid | 95.25*mm 56*mm 75.1426*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_3_sector2_2outer_solid = wc_4_sector2_2outer_solid;
   G4ThreeVector    wc_3_sector2_2outer_trans = {449.04543*mm,1675.8588*mm,241.99168*mm};
   G4RotationMatrix wc_3_sector2_2outer_rot   = G4RotationMatrix::IDENTITY;
   wc_3_sector2_2outer_rot.rotateX( -33.409841*deg    );
   wc_3_sector2_2outer_rot.rotateY(   8.392807*deg );
   wc_3_sector2_2outer_rot.rotateZ( 0*deg);

   //wc_3_sector2_2 | htcc | htcc wc 3, sector2 right | 449.04543*mm 1675.8588*mm 241.99168*mm | -33.409841*deg 8.392807*deg 0*deg | 8080ff2 | Operation: wc_3_sector2_2outer - wc_3_sector2_2inner | 0 | G4_Al | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 1
   wc_3_sector2_2_solid = new G4SubtractionSolid(
         "wc_3_sector2_2_solid",
         wc_3_sector2_2outer_solid,
         wc_3_sector2_2inner_solid,
         0,
         wc_3_sector2_2inner_trans
         );
   wc_3_sector2_2_trans = {449.04543*mm,1675.8588*mm,241.99168*mm};
   wc_3_sector2_2_rot = G4RotationMatrix::IDENTITY;
   wc_3_sector2_2_rot.rotateX( -33.409841*deg    );
   wc_3_sector2_2_rot.rotateY(   8.392807*deg );
   wc_3_sector2_2_rot.rotateZ( 0*deg);

   //pmt_2_sector3_1 | htcc | htcc pmt 2, sector3 left | -492.67437*mm 1838.6863*mm 424.50853*mm | 131.44772*deg 11.356395*deg 0*deg | 80ff80 | Tube | 0*mm 55*mm 1.5*mm 0*deg 360*deg | HTCCPMTquartz | no | 1 | 1 | 1 | 1 | 1 | htcc | htcc | sector manual 3 ring manual 2 half manual 1
   pmt_2_sector3_1_solid = pmt_4_sector3_1_solid;
   pmt_2_sector3_1_trans = {-492.67437*mm, 1838.6863*mm, 424.50853*mm };
   pmt_2_sector3_1_rot = G4RotationMatrix::IDENTITY;
   pmt_2_sector3_1_rot.rotateX( 131.44772*deg   );
   pmt_2_sector3_1_rot.rotateY( 11.356395*deg );
   pmt_2_sector3_1_rot.rotateZ( 0*deg);
   //pmt_2_sector2_2 | htcc | htcc pmt 2, sector2 right | 492.67437*mm 1838.6863*mm 424.50853*mm | 131.44772*deg -11.356395*deg 0*deg | ff8080 | Tube | 0*mm 55*mm 1.5*mm 0*deg 360*deg | HTCCPMTquartz | no | 1 | 1 | 1 | 1 | 1 | htcc | htcc | sector manual 2 ring manual 2 half manual 2
   pmt_2_sector2_2_solid = pmt_4_sector3_1_solid;
   pmt_2_sector2_2_trans = {492.67437*mm, 1838.6863*mm, 424.50853*mm};
   pmt_2_sector2_2_rot = G4RotationMatrix::IDENTITY;
   pmt_2_sector2_2_rot.rotateX( 131.44772*deg   );
   pmt_2_sector2_2_rot.rotateY(-11.356395*deg  );
   pmt_2_sector2_2_rot.rotateZ( 0*deg);

   // ------------------------------
   //wc_2_sector3_1inner | htcc | htcc wc 2, sector3 leftinner | 0 0 0 | 0 0 0 | 999999 | Paraboloid | 125.25*mm 51.338961*mm 76.723114*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_2_sector3_1inner_solid = wc_4_sector3_1inner_solid;
   G4ThreeVector    wc_2_sector3_1inner_trans = {0,0,0};
   G4RotationMatrix wc_2_sector3_1inner_rot   = G4RotationMatrix::IDENTITY;

   //wc_2_sector3_1outer | htcc | htcc wc 2, sector3 leftouter | -473.6232*mm 1767.5862*mm 487.29704*mm | -48.552284*deg -11.356395*deg 0*deg | 999999 | Paraboloid | 95.25*mm 56*mm 75.1426*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_2_sector3_1outer_solid = wc_4_sector3_1outer_solid;
   G4ThreeVector    wc_2_sector3_1outer_trans = {-473.6232*mm,1767.5862*mm,487.29704*mm };
   G4RotationMatrix wc_2_sector3_1outer_rot   = G4RotationMatrix::IDENTITY;
   wc_2_sector3_1outer_rot.rotateX( -48.552284*deg    );
   wc_2_sector3_1outer_rot.rotateY(  -11.356395*deg  );
   wc_2_sector3_1outer_rot.rotateZ( 0*deg);

   //wc_2_sector3_1 | htcc | htcc wc 2, sector3 left | -473.6232*mm 1767.5862*mm 487.29704*mm | -48.552284*deg -11.356395*deg 0*deg | 80ff802 | Operation: wc_2_sector3_1outer - wc_2_sector3_1inner | 0 | G4_Al | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 1
   wc_2_sector3_1_solid = new G4SubtractionSolid(
         "wc_2_sector3_1_solid",
         wc_2_sector3_1outer_solid,
         wc_2_sector3_1inner_solid,
         0,
         wc_2_sector3_1inner_trans
         );
   wc_2_sector3_1_trans = {-473.6232*mm,1767.5862*mm,487.29704*mm};
   wc_2_sector3_1_rot = G4RotationMatrix::IDENTITY;
   wc_2_sector3_1_rot.rotateX( -48.552284*deg    );
   wc_2_sector3_1_rot.rotateY( -11.356395*deg  );
   wc_2_sector3_1_rot.rotateZ( 0*deg);

   //wc_2_sector2_2inner | htcc | htcc wc 2, sector2 rightinner | 0 0 0 | 0 0 0 | 999999 | Paraboloid | 125.25*mm 51.338961*mm 76.723114*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_2_sector2_2inner_solid = wc_4_sector2_2inner_solid;
   G4ThreeVector    wc_2_sector2_2inner_trans = {0,0,0};
   G4RotationMatrix wc_2_sector2_2inner_rot   = G4RotationMatrix::IDENTITY;

   //wc_2_sector2_2outer | htcc | htcc wc 2, sector2 rightouter | 473.6232*mm 1767.5862*mm 487.29704*mm | -48.552284*deg 11.356395*deg 0*deg | 999999 | Paraboloid | 95.25*mm 56*mm 75.1426*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_2_sector2_2outer_solid = wc_4_sector2_2outer_solid;
   G4ThreeVector    wc_2_sector2_2outer_trans = {473.6232*mm,1767.5862*mm,487.29704*mm };
   G4RotationMatrix wc_2_sector2_2outer_rot   = G4RotationMatrix::IDENTITY;
   wc_2_sector2_2outer_rot.rotateX(   -48.552284*deg);
   wc_2_sector2_2outer_rot.rotateY(    11.356395*deg );
   wc_2_sector2_2outer_rot.rotateZ( 0*deg);

   //wc_2_sector2_2 | htcc | htcc wc 2, sector2 right | 473.6232*mm 1767.5862*mm 487.29704*mm | -48.552284*deg 11.356395*deg 0*deg | ff80802 | Operation: wc_2_sector2_2outer - wc_2_sector2_2inner | 0 | G4_Al | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 1
   wc_2_sector2_2_solid = new G4SubtractionSolid(
         "wc_2_sector2_2_solid",
         wc_2_sector2_2outer_solid,
         wc_2_sector2_2inner_solid,
         0,
         wc_2_sector2_2inner_trans
         );
   wc_2_sector2_2_trans = {473.6232*mm,1767.5862*mm,487.29704*mm};
   wc_2_sector2_2_rot = G4RotationMatrix::IDENTITY;
   wc_2_sector2_2_rot.rotateX(   -48.552284*deg);
   wc_2_sector2_2_rot.rotateY(    11.356395*deg );
   wc_2_sector2_2_rot.rotateZ( 0*deg);

   //pmt_1_sector3_1 | htcc | htcc pmt 1, sector3 left | -504.25588*mm 1881.909*mm 707.08098*mm | 116.81116*deg 13.449397*deg 0*deg | f0f0f0 | Tube | 0*mm 55*mm 1.5*mm 0*deg 360*deg | HTCCPMTquartz | no | 1 | 1 | 1 | 1 | 1 | htcc | htcc | sector manual 3 ring manual 1 half manual 1
   pmt_1_sector3_1_solid = pmt_4_sector3_1_solid;
   pmt_1_sector3_1_trans = {-504.25588*mm,1881.909*mm,707.08098*mm };
   pmt_1_sector3_1_rot = G4RotationMatrix::IDENTITY;
   pmt_1_sector3_1_rot.rotateX( 116.81116*deg  );
   pmt_1_sector3_1_rot.rotateY( 13.449397*deg   );
   pmt_1_sector3_1_rot.rotateZ( 0*deg);
   //pmt_1_sector2_2 | htcc | htcc pmt 1, sector2 right | 504.25588*mm 1881.909*mm 707.08098*mm | 116.81116*deg -13.449397*deg 0*deg | 8080ff | Tube | 0*mm 55*mm 1.5*mm 0*deg 360*deg | HTCCPMTquartz | no | 1 | 1 | 1 | 1 | 1 | htcc | htcc | sector manual 2 ring manual 1 half manual 2
   pmt_1_sector2_2_solid = pmt_4_sector3_1_solid;
   pmt_1_sector2_2_trans = {504.25588*mm,1881.909*mm,707.08098*mm};
   pmt_1_sector2_2_rot = G4RotationMatrix::IDENTITY;
   pmt_1_sector2_2_rot.rotateX( 116.81116*deg  );
   pmt_1_sector2_2_rot.rotateY( -13.449397*deg   );
   pmt_1_sector2_2_rot.rotateZ( 0*deg);

   // ------------------------------
   //wc_1_sector3_1inner | htcc | htcc wc 1, sector3 leftinner | 0 0 0 | 0 0 0 | 999999 | Paraboloid | 125.25*mm 51.338961*mm 76.723114*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_1_sector3_1inner_solid = wc_4_sector3_1inner_solid;
   G4ThreeVector    wc_1_sector3_1inner_trans = {0,0,0};
   G4RotationMatrix wc_1_sector3_1inner_rot   = G4RotationMatrix::IDENTITY;

   //wc_1_sector3_1outer | htcc | htcc wc 1, sector3 leftouter | -481.75313*mm 1797.9279*mm 749.52343*mm | -63.188841*deg -13.449397*deg 0*deg | 999999 | Paraboloid | 95.25*mm 56*mm 75.1426*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_1_sector3_1outer_solid = wc_4_sector3_1outer_solid;
   G4ThreeVector    wc_1_sector3_1outer_trans = {-481.75313*mm,1797.9279*mm,749.52343*mm };
   G4RotationMatrix wc_1_sector3_1outer_rot   = G4RotationMatrix::IDENTITY;
   wc_1_sector3_1outer_rot.rotateX( -63.188841*deg     );
   wc_1_sector3_1outer_rot.rotateY( -13.449397*deg );
   wc_1_sector3_1outer_rot.rotateZ(  0*deg);

   //wc_1_sector3_1 | htcc | htcc wc 1, sector3 left | -481.75313*mm 1797.9279*mm 749.52343*mm | -63.188841*deg -13.449397*deg 0*deg | f0f0f02 | Operation: wc_1_sector3_1outer - wc_1_sector3_1inner | 0 | G4_Al | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 1
   wc_1_sector3_1_solid = new G4SubtractionSolid(
         "wc_1_sector3_1_solid",
         wc_1_sector3_1outer_solid,
         wc_1_sector3_1inner_solid,
         0,
         wc_1_sector3_1inner_trans
         );
   wc_1_sector3_1_trans = {-481.75313*mm,1797.9279*mm,749.52343*mm};
   wc_1_sector3_1_rot = G4RotationMatrix::IDENTITY;
   wc_1_sector3_1_rot.rotateX( -63.188841*deg     );
   wc_1_sector3_1_rot.rotateY( -13.449397*deg );
   wc_1_sector3_1_rot.rotateZ(  0*deg);

   //wc_1_sector2_2inner | htcc | htcc wc 1, sector2 rightinner | 0 0 0 | 0 0 0 | 999999 | Paraboloid | 125.25*mm 51.338961*mm 76.723114*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_1_sector2_2inner_solid = wc_4_sector2_2inner_solid;
   G4ThreeVector    wc_1_sector2_2inner_trans = {0,0,0};
   G4RotationMatrix wc_1_sector2_2inner_rot   = G4RotationMatrix::IDENTITY;

   //wc_1_sector2_2outer | htcc | htcc wc 1, sector2 rightouter | 481.75313*mm 1797.9279*mm 749.52343*mm | -63.188841*deg 13.449397*deg 0*deg | 999999 | Paraboloid | 95.25*mm 56*mm 75.1426*mm | Component | no | 1 | 1 | 1 | 1 | 0 | no | no | no
   G4VSolid *       wc_1_sector2_2outer_solid = wc_4_sector2_2outer_solid;
   G4ThreeVector    wc_1_sector2_2outer_trans = {481.75313*mm,1797.9279*mm,749.52343*mm};
   G4RotationMatrix wc_1_sector2_2outer_rot   = G4RotationMatrix::IDENTITY;
   wc_1_sector2_2outer_rot.rotateX(  -63.188841*deg    );
   wc_1_sector2_2outer_rot.rotateY(   8.392807*deg );
   wc_1_sector2_2outer_rot.rotateZ( 0*deg);

   //wc_1_sector2_2 | htcc | htcc wc 1, sector2 right | 481.75313*mm 1797.9279*mm 749.52343*mm | -63.188841*deg 13.449397*deg 0*deg | 8080ff2 | Operation: wc_1_sector2_2outer - wc_1_sector2_2inner | 0 | G4_Al | no | 1 | 1 | 1 | 1 | 1 | mirror: htcc_AlMgF2 | mirror | id manual 1
   wc_1_sector2_2_solid = new G4SubtractionSolid(
         "wc_1_sector2_2_solid",
         wc_1_sector2_2outer_solid,
         wc_1_sector2_2inner_solid,
         0,
         wc_1_sector2_2inner_trans
         );
   wc_1_sector2_2_trans = {481.75313*mm,1797.9279*mm,749.52343*mm};
   wc_1_sector2_2_rot = G4RotationMatrix::IDENTITY;
   wc_1_sector2_2_rot.rotateX(  -63.188841*deg    );
   wc_1_sector2_2_rot.rotateY(    13.449397*deg );
   wc_1_sector2_2_rot.rotateZ( 0*deg);

}
//______________________________________________________________________________

void HTCCDetectorGeometry::BuildLogicalVolumes()
{

   using namespace CLHEP;
   bool check_overlaps = false;

   BuildMirrors();

   G4NistManager* nist = G4NistManager::Instance();
   G4Material * default_mat   = nist->FindOrBuildMaterial("G4_AIR");

   fGasMaterial = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");

   const G4int NUMENTRIES = 7; //32;

   G4double ppckov[NUMENTRIES] = {
      1.9074494*eV, 1.9372533*eV, 1.9680033*eV, 1.9997453*eV, 2.0325280*eV,
      6.1992105*eV, 6.5254848*eV
   };
   //   2.0664035*eV, 2.1014273*eV, 2.1376588*eV, 2.1751616*eV, 2.2140038*eV,
   //   2.2542584*eV, 2.2960039*eV, 2.3393247*eV, 2.3843117*eV, 2.4310630*eV,
   //   2.4796842*eV, 2.5302900*eV, 2.5830044*eV, 2.6379619*eV, 2.6953089*eV,
   //   2.7552047*eV, 2.8178230*eV, 2.8833537*eV, 2.9520050*eV, 3.0240051*eV,
   //   3.0996053*eV, 3.1790823*eV, 3.2627424*eV, 3.3509246*eV, 3.4440059*eV,
   //   3.5424060*eV, 3.6465944*eV, 3.7570973*eV, 3.8745066*eV, 3.9994907*eV,
   //   4.1328070*eV, 4.2753176*eV, 4.4280075*eV, 4.5920078*eV, 4.7686235*eV,
   //   4.9593684*eV, 5.1660088*eV, 5.3906179*eV, 5.6356459*eV, 5.9040100*eV,
   //   6.1992105*eV, 6.5254848*eV
   //};
   G4double rindex[NUMENTRIES] = {
      1.0004473, 1.0004475, 1.0004477, 1.000448, 1.0004483,
      1.0005262, 1.0005378
   };
   //   1.0004486, 1.0004489, 1.0004492, 1.0004495, 1.0004498,
   //   1.0004502, 1.0004506, 1.000451, 1.0004514, 1.0004518,
   //   1.0004523, 1.0004528, 1.0004534, 1.0004539, 1.0004545,
   //   1.0004552, 1.0004559, 1.0004566, 1.0004574, 1.0004583,
   //   1.0004592, 1.0004602, 1.0004613, 1.0004625, 1.0004638,
   //   1.0004652, 1.0004668, 1.0004685, 1.0004704, 1.0004724,
   //   1.0004748, 1.0004773, 1.0004803, 1.0004835, 1.0004873,
   //   1.0004915, 1.0004964, 1.0005021, 1.0005088, 1.0005167,
   //   1.0005262, 1.0005378
   //};
   G4double absorption[NUMENTRIES] = {
      1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m,
      802.8323273*m, 0.7465970*m
   };
   //   1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m,
   //   1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m,
   //   1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m,
   //   1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m,
   //   1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m,
   //   1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m, 1000.0000000*m,
   //   82.8323273*m, 0.7465970*m
   //};

   // ---------------------------------------------

   G4MaterialPropertiesTable *MPT = new G4MaterialPropertiesTable();
   if(fUseIndexOfRefraction){
      //MPT->AddConstProperty("SCINTILLATIONYIELD",100./MeV);
      MPT->AddProperty("RINDEX",ppckov,rindex,NUMENTRIES);
      MPT->AddProperty("ABSLENGTH",ppckov,absorption,NUMENTRIES);
   }
   fGasMaterial->SetMaterialPropertiesTable(MPT);

   // ---------------------------------------------

   OpSurface = new G4OpticalSurface("mirror");
   G4double sigma_alpha = 0.1;
   OpSurface->SetType(dielectric_metal);
   OpSurface->SetFinish(polished);
   OpSurface->SetModel(unified);
   //OpSurface -> SetSigmaAlpha(sigma_alpha);

   // ---------------------------------------------

   const G4int NUM = 2;
   G4double pp[NUM] = {1.9074494*eV,6.5254848*eV};
   G4double specularlobe[NUM]  = {0.3, 0.3};
   G4double specularspike[NUM] = {0.2, 0.2};
   G4double backscatter[NUM]   = {0.0, 0.0};
   G4double rindex2[NUM]       = {1.05, 1.40};
   G4double reflectivity[NUM]  = {0.9, 0.9}; // made up
   G4double efficiency[NUM]    = {0.0, 0.0};

   G4MaterialPropertiesTable* SMPT = new G4MaterialPropertiesTable();
   //SMPT -> AddProperty("RINDEX",pp,rindex2,NUM);
   //SMPT -> AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM);
   //SMPT -> AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM);
   //SMPT -> AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM);
   SMPT -> AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
   //SMPT -> AddProperty("EFFICIENCY",pp,efficiency,NUM);
   OpSurface->SetMaterialPropertiesTable(SMPT);

   // ---------------------------------------------

   //htccBigGasVolume_log    = new G4LogicalVolume(htccBigGasVolume_solid, default_mat, "htccBigGasVolume_log");
   //htccEntryDishVolume_log = new G4LogicalVolume(htccEntryDishVolume_solid, default_mat, "htccEntryDishVolume_log");
   //htccEntryConeVolume_log = new G4LogicalVolume(htccEntryConeVolume_solid, default_mat, "htccEntryConeVolume_log");
   //htccEntryDishCone_log = new G4LogicalVolume(htccEntryDishCone_solid, default_mat, "htccEntryDishCone_log");

   G4VisAttributes * htcc_vis = new G4VisAttributes(G4VisAttributes::GetInvisible());//;
   htcc_vis->SetDaughtersInvisible(false);
   htcc_vis->SetForceSolid(false);
   htcc_vis->SetForceWireframe(true);

   htcc_log               = new G4LogicalVolume(htcc_solid, fGasMaterial, "htcc_log");
   htcc_log->SetVisAttributes(htcc_vis);//G4VisAttributes::GetInvisible());

   G4VisAttributes * vs3 = new G4VisAttributes(G4Colour(0.8,0.1,0.1));//G4VisAttributes::GetInvisible());//;
   vs3->SetDaughtersInvisible(false);
   vs3->SetForceSolid(false);
   vs3->SetForceWireframe(true);

   for(int i = 0; i<8; i++) {
      fMirror_vis[i] = new G4VisAttributes(G4Colour(0.0,0.1+float(i)/10.0,0.8-float(i)/10.0, 0.5 ));//G4VisAttributes::GetInvisible());//;
      fMirror_vis[i]->SetForceSolid(true);
      //vs3->SetForceWireframe(true);
   }


   sector_wedge_log = new G4LogicalVolume(sector_wedge_solid, fGasMaterial, "sector_wedge_log");
   sector_wedge_log->SetVisAttributes(vs3);

   //-----------------------------------------------------------------------------

   G4Material * mirror_mat   = nist->FindOrBuildMaterial ( "G4_Pyrex_Glass" );
   G4double rindex3[NUM]     = {1.0, 1.0};
   G4double absorption3[NUM] = {1000.0*m, 1000.0*m};
   G4MaterialPropertiesTable * mirror_MPT = new G4MaterialPropertiesTable();
   mirror_MPT->AddProperty("RINDEX",   pp,rindex3,NUM);
   //mirror_MPT->AddProperty("ABSLENGTH",pp,absorption3,NUM);
   //mirror_mat->SetMaterialPropertiesTable(mirror_MPT);

   G4Material * pmt_mat   = nist->FindOrBuildMaterial ( "G4_SILICON_DIOXIDE" );
   G4double rindex4[NUM]     = {1.5, 1.5};
   G4double absorption4[NUM] = {1000.0*m, 1000.0*m};
   G4MaterialPropertiesTable * pmt_MPT = new G4MaterialPropertiesTable();
   pmt_MPT->AddProperty("RINDEX",   pp,rindex4,NUM);
   //pmt_MPT->AddProperty("ABSLENGTH",pp,absorption4,NUM);
   pmt_mat->SetMaterialPropertiesTable(pmt_MPT);

   // ---------------------------------------------

   mirror_4_sector2_2_log = new G4LogicalVolume(mirror_4_sector2_2_solid, mirror_mat, "mirror_4_sector2_2_log");
   mirror_3_sector2_2_log = new G4LogicalVolume(mirror_3_sector2_2_solid, mirror_mat, "mirror_3_sector2_2_log");
   mirror_2_sector2_2_log = new G4LogicalVolume(mirror_2_sector2_2_solid, mirror_mat, "mirror_2_sector2_2_log");
   mirror_1_sector2_2_log = new G4LogicalVolume(mirror_1_sector2_2_solid, mirror_mat, "mirror_1_sector2_2_log");
   mirror_4_sector3_1_log = new G4LogicalVolume(mirror_4_sector3_1_solid, mirror_mat, "mirror_4_sector3_1_log");
   mirror_3_sector3_1_log = new G4LogicalVolume(mirror_3_sector3_1_solid, mirror_mat, "mirror_3_sector3_1_log");
   mirror_2_sector3_1_log = new G4LogicalVolume(mirror_2_sector3_1_solid, mirror_mat, "mirror_2_sector3_1_log");
   mirror_1_sector3_1_log = new G4LogicalVolume(mirror_1_sector3_1_solid, mirror_mat, "mirror_1_sector3_1_log");

   fMirror_log[0] = mirror_4_sector2_2_log ;
   fMirror_log[1] = mirror_3_sector2_2_log ;
   fMirror_log[2] = mirror_2_sector2_2_log ;
   fMirror_log[3] = mirror_1_sector2_2_log ;
   fMirror_log[4] = mirror_4_sector3_1_log ;
   fMirror_log[5] = mirror_3_sector3_1_log ;
   fMirror_log[6] = mirror_2_sector3_1_log ;
   fMirror_log[7] = mirror_1_sector3_1_log ;

   //-----------------------------------------------------------------------------

   pmt_4_sector3_1_log = new G4LogicalVolume(pmt_4_sector3_1_solid , pmt_mat , "pmt_4_sector3_1_log");
   pmt_4_sector2_2_log = new G4LogicalVolume(pmt_4_sector2_2_solid , pmt_mat , "pmt_4_sector2_2_log");
   wc_4_sector3_1_log  = new G4LogicalVolume(wc_4_sector3_1_solid  , mirror_mat  , "wc_4_sector3_1_log");
   wc_4_sector2_2_log  = new G4LogicalVolume(wc_4_sector2_2_solid  , mirror_mat  , "wc_4_sector2_2_log");

   pmt_3_sector3_1_log = new G4LogicalVolume(pmt_3_sector3_1_solid , pmt_mat , "pmt_3_sector3_1_log");
   pmt_3_sector2_2_log = new G4LogicalVolume(pmt_3_sector2_2_solid , pmt_mat , "pmt_3_sector2_2_log");
   wc_3_sector3_1_log  = new G4LogicalVolume(wc_3_sector3_1_solid  , mirror_mat  , "wc_3_sector3_1_log");
   wc_3_sector2_2_log  = new G4LogicalVolume(wc_3_sector2_2_solid  , mirror_mat  , "wc_3_sector2_2_log");

   pmt_2_sector3_1_log = new G4LogicalVolume(pmt_2_sector3_1_solid , pmt_mat , "pmt_2_sector3_1_log");
   pmt_2_sector2_2_log = new G4LogicalVolume(pmt_2_sector2_2_solid , pmt_mat , "pmt_2_sector2_2_log");
   wc_2_sector3_1_log  = new G4LogicalVolume(wc_2_sector3_1_solid  , mirror_mat  , "wc_2_sector3_1_log");
   wc_2_sector2_2_log  = new G4LogicalVolume(wc_2_sector2_2_solid  , mirror_mat  , "wc_2_sector2_2_log");

   pmt_1_sector3_1_log = new G4LogicalVolume(pmt_1_sector3_1_solid , pmt_mat , "pmt_1_sector3_1_log");
   pmt_1_sector2_2_log = new G4LogicalVolume(pmt_1_sector2_2_solid , pmt_mat , "pmt_1_sector2_2_log");
   wc_1_sector3_1_log  = new G4LogicalVolume(wc_1_sector3_1_solid  , mirror_mat  , "wc_1_sector3_1_log");
   wc_1_sector2_2_log  = new G4LogicalVolume(wc_1_sector2_2_solid  , mirror_mat  , "wc_1_sector2_2_log");

   //BarrelEllipseCut_sect0mirr0half0_log = new G4LogicalVolume(BarrelEllipseCut_sect0mirr0half0_solid, default_mat, "BarrelEllipseCut_sect0mirr0half0_log");
   pmt_4_sector3_1_log->SetSensitiveDetector(fSensitiveDetector);
   pmt_4_sector2_2_log->SetSensitiveDetector(fSensitiveDetector);
   pmt_3_sector3_1_log->SetSensitiveDetector(fSensitiveDetector);
   pmt_3_sector2_2_log->SetSensitiveDetector(fSensitiveDetector);
   pmt_2_sector3_1_log->SetSensitiveDetector(fSensitiveDetector);
   pmt_2_sector2_2_log->SetSensitiveDetector(fSensitiveDetector);
   pmt_1_sector3_1_log->SetSensitiveDetector(fSensitiveDetector);
   pmt_1_sector2_2_log->SetSensitiveDetector(fSensitiveDetector);

   fPMT_log[0] = pmt_4_sector2_2_log ;
   fPMT_log[1] = pmt_3_sector2_2_log ;
   fPMT_log[2] = pmt_2_sector2_2_log ;
   fPMT_log[3] = pmt_1_sector2_2_log ;
   fPMT_log[4] = pmt_4_sector3_1_log ;
   fPMT_log[5] = pmt_3_sector3_1_log ;
   fPMT_log[6] = pmt_2_sector3_1_log ;
   fPMT_log[7] = pmt_1_sector3_1_log ;

   fWC_log[0] = wc_4_sector2_2_log ;
   fWC_log[1] = wc_3_sector2_2_log ;
   fWC_log[2] = wc_2_sector2_2_log ;
   fWC_log[3] = wc_1_sector2_2_log ;
   fWC_log[4] = wc_4_sector3_1_log ;
   fWC_log[5] = wc_3_sector3_1_log ;
   fWC_log[6] = wc_2_sector3_1_log ;
   fWC_log[7] = wc_1_sector3_1_log ;

   // Visual attributes
   for(int i = 0; i< 8; i++) {
    fMirror_log[i]->SetVisAttributes( fMirror_vis[i] );
    fPMT_log[i]->SetVisAttributes( fMirror_vis[i] );
    fWC_log[i]->SetVisAttributes( fMirror_vis[i] );
   }


}
//______________________________________________________________________________

G4VPhysicalVolume * HTCCDetectorGeometry::PlacePhysicalVolume(G4LogicalVolume * mother, int sec, int region )
{
   std::cout << " HTCCDetectorGeometry::PlacePhysicalVolume \n";
   using namespace CLHEP;
   using namespace clas12::geo;
   int index    = region-1;
   int grouping = (sec-1)*3 + (region-1);

   if(!htcc_phys) {
      htcc_phys = new G4PVPlacement(
            0,
            G4ThreeVector(0,0,0),
            htcc_log,//htcc_log,          // its logical volume
            "htcc_phys", // its name
            mother,                       // its mother (logical) volume
            false,                        // no boolean operations
            0,                     // its copy number
            true);                        // check for overlaps
   }

   // --------------------------------------------
   G4RotationMatrix * sector_rot = new G4RotationMatrix();
   sector_rot->rotateZ( 60.0*CLHEP::deg*(sec-1) );
   G4VPhysicalVolume * sector_phys = new G4PVPlacement(
         sector_rot, 
         G4ThreeVector(0,0,0),
         sector_wedge_log,          // its logical volume
         Form("sector_wedge_phys%d",sec), // its name
         htcc_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         sec,                     // its copy number
         true);                        // check for overlaps
   //G4UserLimits * scoring_limits = new G4UserLimits(1.0*mm);
   //sector_wedge_log->SetUserLimits(scoring_limits);

   PlaceMirrors(sector_phys);

   return htcc_phys;
}
//______________________________________________________________________________

G4VPhysicalVolume * HTCCDetectorGeometry::PlaceParallelPhysicalVolume(G4LogicalVolume * mother, int sec, int region )
{
   using namespace clas12::geo;
   int index    = region-1;
   int grouping = (sec-1)*3 + (region-1);

   //G4VPhysicalVolume * phys = new G4PVPlacement(
   //      G4Transform3D(
   //         RegionRotation(sec,region),
   //         RegionTranslation(sec, region)
   //      ),
   //      fRegions_log[index],          // its logical volume
   //      Form("region%d_phys",region), // its name
   //      mother,                       // its mother (logical) volume
   //      false,                        // no boolean operations
   //      grouping,                     // its copy number
   //      false);                        // check for overlaps

   return nullptr;
}
//______________________________________________________________________________
void HTCCDetectorGeometry::PlaceMirrors(G4VPhysicalVolume * mother_phys)
{
   using namespace CLHEP;
   using namespace clas12::geo;

   // ----------------------------------------------
   // Mirrors
   bool check_mirror_overlaps = false;
   G4VPhysicalVolume * m1_phys = new G4PVPlacement(
         &mirror_4_sector2_2_rot, 
         mirror_4_sector2_2_trans,
         mirror_4_sector2_2_log,          // its logical volume
         "mirror_4_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps
   G4VPhysicalVolume * m2_phys = new G4PVPlacement(
         &mirror_3_sector2_2_rot,
         mirror_3_sector2_2_trans,
         mirror_3_sector2_2_log,          // its logical volume
         "mirror_3_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps
   G4VPhysicalVolume * m3_phys = new G4PVPlacement(
         &mirror_2_sector2_2_rot,
         mirror_2_sector2_2_trans,
         mirror_2_sector2_2_log,          // its logical volume
         "mirror_2_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps
   G4VPhysicalVolume * m4_phys = new G4PVPlacement(
         &mirror_1_sector2_2_rot, 
         mirror_1_sector2_2_trans,
         mirror_1_sector2_2_log,          // its logical volume
         "mirror_1_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps
   // other side
   G4VPhysicalVolume * m5_phys = new G4PVPlacement(
        &mirror_4_sector3_1_rot, 
         mirror_4_sector3_1_trans,
         mirror_4_sector3_1_log,          // its logical volume
        "mirror_4_sector3_1_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps
   G4VPhysicalVolume * m6_phys = new G4PVPlacement(
        &mirror_3_sector3_1_rot, 
         mirror_3_sector3_1_trans,
         mirror_3_sector3_1_log,          // its logical volume
        "mirror_3_sector3_1_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps
   G4VPhysicalVolume * m7_phys = new G4PVPlacement(
        &mirror_2_sector3_1_rot, 
         mirror_2_sector3_1_trans,
         mirror_2_sector3_1_log,          // its logical volume
        "mirror_2_sector3_1_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps
   G4VPhysicalVolume * m8_phys = new G4PVPlacement(
        &mirror_1_sector3_1_rot, 
         mirror_1_sector3_1_trans,
         mirror_1_sector3_1_log,          // its logical volume
        "mirror_1_sector3_1_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps

   // ---------------------------------
   //

   //const G4int num = 2;
   //G4double Ephoton[num]          = {1.0*eV, 10.0*eV};
   //G4double MSpecularLobe[num]    = {0.0, 0.0};
   //G4double MSpecularSpike[num]   = {0.0, 0.0};
   //G4double MBackscatter[num]     = {0.0, 0.0};
   //G4double MEfficiency[num]      = {0.0, 0.0};
   //G4double coatingReflectivity[num] = {0.95, 0.95};

   //G4OpticalSurface* OpMirrorSurface = new G4OpticalSurface ( "SilverSurface" );
   //OpMirrorSurface->SetType(   dielectric_metal );
   //OpMirrorSurface->SetFinish( polishedbackpainted );
   //OpMirrorSurface->SetModel(  unified );

   //G4MaterialPropertiesTable* reflectiveCoatingSurfacePTable = new G4MaterialPropertiesTable();

   //reflectiveCoatingSurfacePTable->AddProperty ( "REFLECTIVITY",          Ephoton, coatingReflectivity, num );
   //reflectiveCoatingSurfacePTable->AddProperty ( "SPECULARLOBECONSTANT",  Ephoton, MSpecularLobe,    num );
   //reflectiveCoatingSurfacePTable->AddProperty ( "SPECULARSPIKECONSTANT", Ephoton, MSpecularSpike,   num );
   //reflectiveCoatingSurfacePTable->AddProperty ( "BACKSCATTERCONSTANT",   Ephoton, MBackscatter,     num );
   //OpMirrorSurface->SetMaterialPropertiesTable ( reflectiveCoatingSurfacePTable );

   //G4LogicalBorderSurface* Surface1 = new G4LogicalBorderSurface("mirror_1", mother_phys, m1_phys, OpSurface);
   //G4LogicalBorderSurface* Surface2 = new G4LogicalBorderSurface("mirror_2", mother_phys, m2_phys, OpSurface);
   //G4LogicalBorderSurface* Surface3 = new G4LogicalBorderSurface("mirror_3", mother_phys, m3_phys, OpSurface);
   //G4LogicalBorderSurface* Surface4 = new G4LogicalBorderSurface("mirror_4", mother_phys, m4_phys, OpSurface);

   //G4LogicalBorderSurface* Surface5 = new G4LogicalBorderSurface("mirror_5", mother_phys, m5_phys, OpSurface);
   //G4LogicalBorderSurface* Surface6 = new G4LogicalBorderSurface("mirror_6", mother_phys, m6_phys, OpSurface);
   //G4LogicalBorderSurface* Surface7 = new G4LogicalBorderSurface("mirror_7", mother_phys, m7_phys, OpSurface);
   //G4LogicalBorderSurface* Surface8 = new G4LogicalBorderSurface("mirror_8", mother_phys, m8_phys, OpSurface);

   G4LogicalSkinSurface * surf1 = new G4LogicalSkinSurface("mirror_1", mirror_4_sector2_2_log , OpSurface);
   G4LogicalSkinSurface * surf2 = new G4LogicalSkinSurface("mirror_2", mirror_2_sector2_2_log , OpSurface);
   G4LogicalSkinSurface * surf3 = new G4LogicalSkinSurface("mirror_3", mirror_3_sector2_2_log , OpSurface);
   G4LogicalSkinSurface * surf4 = new G4LogicalSkinSurface("mirror_4", mirror_1_sector2_2_log , OpSurface);
   G4LogicalSkinSurface * surf5 = new G4LogicalSkinSurface("mirror_5", mirror_4_sector3_1_log , OpSurface);
   G4LogicalSkinSurface * surf6 = new G4LogicalSkinSurface("mirror_6", mirror_3_sector3_1_log , OpSurface);
   G4LogicalSkinSurface * surf7 = new G4LogicalSkinSurface("mirror_7", mirror_2_sector3_1_log , OpSurface);
   G4LogicalSkinSurface * surf8 = new G4LogicalSkinSurface("mirror_8", mirror_1_sector3_1_log , OpSurface);

   // ---------------------------------
   // pmt
   new G4PVPlacement(
        &pmt_4_sector3_1_rot, 
         pmt_4_sector3_1_trans,
         pmt_4_sector3_1_log,   
        "pmt_4_sector3_1_phys", 
         sector_wedge_log,      
         false,                 
         0,                     
         check_mirror_overlaps);
   // winston cone
   new G4PVPlacement(
        &wc_4_sector3_1_rot, 
         wc_4_sector3_1_trans,
         wc_4_sector3_1_log,   
        "wc_4_sector3_1_phys", 
         sector_wedge_log,     
         false,                
         0,                    
         check_mirror_overlaps);
   // pmt
   new G4PVPlacement(
        &pmt_4_sector2_2_rot, 
         pmt_4_sector2_2_trans,
         pmt_4_sector2_2_log,    
        "pmt_4_sector2_2_phys",  
         sector_wedge_log,       
         false,                  
         0+4,                      
         check_mirror_overlaps); 
   // winston cone
   new G4PVPlacement(
        &wc_4_sector2_2_rot, 
         wc_4_sector2_2_trans,
         wc_4_sector2_2_log,          // its logical volume
        "wc_4_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         0+4,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps

   // --------------------------------------
   // pmt
   new G4PVPlacement(
        &pmt_3_sector3_1_rot, 
         pmt_3_sector3_1_trans,
         pmt_3_sector3_1_log,          // its logical volume
        "pmt_3_sector3_1_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         1,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps
   // winston cone
   new G4PVPlacement(
        &wc_3_sector3_1_rot, 
         wc_3_sector3_1_trans,
         wc_3_sector3_1_log,          // its logical volume
        "wc_3_sector3_1_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         1,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps
   // pmt
   new G4PVPlacement(
        &pmt_3_sector2_2_rot, 
         pmt_3_sector2_2_trans,
         pmt_3_sector2_2_log,          // its logical volume
        "pmt_3_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         1+4,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps
   // winston cone
   new G4PVPlacement(
        &wc_3_sector2_2_rot, 
         wc_3_sector2_2_trans,
         wc_3_sector2_2_log,          // its logical volume
        "wc_3_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         1+4,                     // its copy number
         check_mirror_overlaps);                        // check for overlaps

   // --------------------------------------
   // pmt
   new G4PVPlacement(
        &pmt_2_sector3_1_rot, 
         pmt_2_sector3_1_trans,
         pmt_2_sector3_1_log,          // its logical volume
        "pmt_2_sector3_1_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         2,                     // its copy number
         false);                        // check for overlaps
   // winston cone
   new G4PVPlacement(
        &wc_2_sector3_1_rot, 
         wc_2_sector3_1_trans,
         wc_2_sector3_1_log,          // its logical volume
        "wc_2_sector3_1_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         2,                     // its copy number
         false);                        // check for overlaps
   // pmt
   new G4PVPlacement(
        &pmt_2_sector2_2_rot, 
         pmt_2_sector2_2_trans,
         pmt_2_sector2_2_log,          // its logical volume
        "pmt_2_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         2+4,                     // its copy number
         false);                        // check for overlaps
   // winston cone
   new G4PVPlacement(
        &wc_2_sector2_2_rot, 
         wc_2_sector2_2_trans,
         wc_2_sector2_2_log,          // its logical volume
        "wc_2_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         2+4,                     // its copy number
         false);                        // check for overlaps

   // --------------------------------------
   // pmt
   new G4PVPlacement(
        &pmt_1_sector3_1_rot, 
         pmt_1_sector3_1_trans,
         pmt_1_sector3_1_log,          // its logical volume
        "pmt_1_sector3_1_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         3,                     // its copy number
         false);                        // check for overlaps
   // winston cone
   new G4PVPlacement(
        &wc_1_sector3_1_rot, 
         wc_1_sector3_1_trans,
         wc_1_sector3_1_log,          // its logical volume
        "wc_1_sector3_1_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         3,                     // its copy number
         false);                        // check for overlaps
   // pmt
   new G4PVPlacement(
        &pmt_1_sector2_2_rot, 
         pmt_1_sector2_2_trans,
         pmt_1_sector2_2_log,          // its logical volume
        "pmt_1_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         3+4,                     // its copy number
         false);                        // check for overlaps
   // winston cone
   new G4PVPlacement(
        &wc_1_sector2_2_rot, 
         wc_1_sector2_2_trans,
         wc_1_sector2_2_log,          // its logical volume
        "wc_1_sector2_2_phys", // its name
         sector_wedge_log,                       // its mother (logical) volume
         false,                        // no boolean operations
         3+4,                     // its copy number
         false);                        // check for overlaps

   G4LogicalSkinSurface * wc1 = new G4LogicalSkinSurface("wc_1", wc_4_sector2_2_log , OpSurface);
   G4LogicalSkinSurface * wc2 = new G4LogicalSkinSurface("wc_2", wc_2_sector2_2_log , OpSurface);
   G4LogicalSkinSurface * wc3 = new G4LogicalSkinSurface("wc_3", wc_3_sector2_2_log , OpSurface);
   G4LogicalSkinSurface * wc4 = new G4LogicalSkinSurface("wc_4", wc_1_sector2_2_log , OpSurface);
   G4LogicalSkinSurface * wc5 = new G4LogicalSkinSurface("wc_5", wc_4_sector3_1_log , OpSurface);
   G4LogicalSkinSurface * wc6 = new G4LogicalSkinSurface("wc_6", wc_3_sector3_1_log , OpSurface);
   G4LogicalSkinSurface * wc7 = new G4LogicalSkinSurface("wc_7", wc_2_sector3_1_log , OpSurface);
   G4LogicalSkinSurface * wc8 = new G4LogicalSkinSurface("wc_8", wc_1_sector3_1_log , OpSurface);
         
   const G4int num = 2;
   G4double Ephoton[num]          = {1.0*eV, 10.0*eV};
   G4double MSpecularLobe[num]    = {0.0, 0.0};
   G4double MSpecularSpike[num]   = {0.0, 0.0};
   G4double MBackscatter[num]     = {0.0, 0.0};
   G4double MEfficiency[num]      = {1.0, 1.0};
   G4double pmtReflectivity[num] = {0.0, 0.0};

   G4OpticalSurface* OpPMTSurface = new G4OpticalSurface ( "pmt_surface" );
   OpPMTSurface->SetType(   dielectric_metal );
   OpPMTSurface->SetFinish( polished );
   OpPMTSurface->SetModel(  unified );

   G4MaterialPropertiesTable* pmtSurfacePTable = new G4MaterialPropertiesTable();

   pmtSurfacePTable->AddProperty( "REFLECTIVITY",          Ephoton, pmtReflectivity, num );
   pmtSurfacePTable->AddProperty( "SPECULARLOBECONSTANT",  Ephoton, MSpecularLobe,    num );
   pmtSurfacePTable->AddProperty( "SPECULARSPIKECONSTANT", Ephoton, MSpecularSpike,   num );
   pmtSurfacePTable->AddProperty( "BACKSCATTERCONSTANT",   Ephoton, MBackscatter,     num );
   pmtSurfacePTable->AddProperty( "EFFICIENCY",            Ephoton, MEfficiency,num);
   OpPMTSurface->SetMaterialPropertiesTable ( pmtSurfacePTable );

   G4LogicalSkinSurface * pmt1 = new G4LogicalSkinSurface("pmt_1", pmt_4_sector2_2_log , OpPMTSurface);
   G4LogicalSkinSurface * pmt2 = new G4LogicalSkinSurface("pmt_2", pmt_2_sector2_2_log , OpPMTSurface);
   G4LogicalSkinSurface * pmt3 = new G4LogicalSkinSurface("pmt_3", pmt_3_sector2_2_log , OpPMTSurface);
   G4LogicalSkinSurface * pmt4 = new G4LogicalSkinSurface("pmt_4", pmt_1_sector2_2_log , OpPMTSurface);
   G4LogicalSkinSurface * pmt5 = new G4LogicalSkinSurface("pmt_5", pmt_4_sector3_1_log , OpPMTSurface);
   G4LogicalSkinSurface * pmt6 = new G4LogicalSkinSurface("pmt_6", pmt_3_sector3_1_log , OpPMTSurface);
   G4LogicalSkinSurface * pmt7 = new G4LogicalSkinSurface("pmt_7", pmt_2_sector3_1_log , OpPMTSurface);
   G4LogicalSkinSurface * pmt8 = new G4LogicalSkinSurface("pmt_8", pmt_1_sector3_1_log , OpPMTSurface);
}
//______________________________________________________________________________


