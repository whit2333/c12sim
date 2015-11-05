#ifndef HTCCDetectorGeometry_HH
#define HTCCDetectorGeometry_HH

#include <array>
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
//#include "HTCCSensitiveDetector.h"
#include "SensitiveRegionDetector.h"

class HTCCDetectorGeometry {

   private:

      G4VSolid * fHexWireVolume_solid;

      std::array<G4LogicalVolume*,3> fRegions_log;
      std::array<G4LogicalVolume*,3> fEmptyRegions_log;
      std::array<G4VSolid*,3>        fRegions_solid;

   protected:
      G4VSolid        * HTCC_OuterCutCylinder_solid        ;
      G4VSolid        * HTCC_InnerCutCone_solid            ;
      G4VSolid        * Barrel_sect0half0_solid           ;
      G4VSolid        * phicut_sect0half0_solid           ;
      G4VSolid        * Mirror_sect0mirr0half0_solid      ;
      G4VSolid        * BarrelEllipseCut_sect0mirr0half0_solid ;
      G4VSolid        * Boxcut12_sect0mirr0half0_solid     ;
      G4VSolid        * MirrorBoxCut_sect0mirr0half0_solid ;
      G4VSolid        * MirrorCylinderCut_sect0mirr0half0_solid ;
      G4VSolid        * mirror_4_sector2_2_solid           ;
      G4LogicalVolume * mirror_4_sector2_2_log          ;
      G4LogicalVolume * BarrelEllipseCut_sect0mirr0half0_log;

      G4VSolid        * Mirror_sect0mirr1half0_solid            ;
      G4VSolid        * BarrelEllipseCut_sect0mirr1half0_solid  ;
      G4VSolid        * Boxcut_up_sect0mirr1half0_solid         ;
      G4VSolid        * Boxcut_down_sect0mirr1half0_solid       ;
      G4VSolid        * MirrorBoxCut_up_sect0mirr1half0_solid   ;
      G4VSolid        * MirrorBoxCut_down_sect0mirr1half0_solid ;
      G4VSolid        * mirror_3_sector2_2_solid                ;
      G4LogicalVolume * mirror_3_sector2_2_log          ;
      G4VSolid        * Mirror_sect0mirr2half0_solid            ;
      G4VSolid        * BarrelEllipseCut_sect0mirr2half0_solid  ;
      G4VSolid        * Boxcut_up_sect0mirr2half0_solid         ;
      G4VSolid        * Boxcut_down_sect0mirr2half0_solid       ;
      G4VSolid        * MirrorBoxCut_up_sect0mirr2half0_solid   ;
      G4VSolid        * MirrorBoxCut_down_sect0mirr2half0_solid ;
      G4VSolid        * mirror_2_sector2_2_solid                ;
      G4LogicalVolume * mirror_2_sector2_2_log          ;
      G4VSolid        * Mirror_sect0mirr3half0_solid            ;
      G4VSolid        * BarrelEllipseCut_sect0mirr3half0_solid  ;
      G4VSolid        * Boxcut_up_sect0mirr3half0_solid         ;
      G4VSolid        * MirrorBoxCut_up_sect0mirr3half0_solid   ;
      G4VSolid        * MirrorConeCut_sect0mirr3half0_solid     ;
      G4VSolid        * mirror_1_sector2_2_solid                ;
      G4LogicalVolume * mirror_1_sector2_2_log          ;
      G4VSolid        * Barrel_sect0half1_solid                 ;
      G4VSolid        * phicut_sect0half1_solid                 ;
      G4VSolid        * Mirror_sect0mirr0half1_solid            ;
      G4VSolid        * BarrelEllipseCut_sect0mirr0half1_solid  ;
      G4VSolid        * Boxcut12_sect0mirr0half1_solid          ;
      G4VSolid        * MirrorBoxCut_sect0mirr0half1_solid      ;
      G4VSolid        * MirrorCylinderCut_sect0mirr0half1_solid ;
      G4VSolid        * mirror_4_sector3_1_solid                ;
      G4LogicalVolume * mirror_4_sector3_1_log          ;
      G4VSolid        * Mirror_sect0mirr1half1_solid            ;
      G4VSolid        * BarrelEllipseCut_sect0mirr1half1_solid  ;
      G4VSolid        * Boxcut_up_sect0mirr1half1_solid         ;
      G4VSolid        * Boxcut_down_sect0mirr1half1_solid       ;
      G4VSolid        * MirrorBoxCut_up_sect0mirr1half1_solid   ;
      G4VSolid        * MirrorBoxCut_down_sect0mirr1half1_solid ;
      G4VSolid        * mirror_3_sector3_1_solid                ;
      G4LogicalVolume * mirror_3_sector3_1_log          ;
      G4VSolid        * Mirror_sect0mirr2half1_solid            ;
      G4VSolid        * BarrelEllipseCut_sect0mirr2half1_solid  ;
      G4VSolid        * Boxcut_up_sect0mirr2half1_solid         ;
      G4VSolid        * Boxcut_down_sect0mirr2half1_solid       ;
      G4VSolid        * MirrorBoxCut_up_sect0mirr2half1_solid   ;
      G4VSolid        * MirrorBoxCut_down_sect0mirr2half1_solid ;
      G4VSolid        * mirror_2_sector3_1_solid                ;
      G4LogicalVolume * mirror_2_sector3_1_log          ;
      G4VSolid        * Mirror_sect0mirr3half1_solid            ;
      G4VSolid        * BarrelEllipseCut_sect0mirr3half1_solid  ;
      G4VSolid        * Boxcut_up_sect0mirr3half1_solid         ;
      G4VSolid        * MirrorBoxCut_up_sect0mirr3half1_solid   ;
      G4VSolid        * MirrorConeCut_sect0mirr3half1_solid     ;
      G4VSolid        * mirror_1_sector3_1_solid                ;
      G4LogicalVolume * mirror_1_sector3_1_log          ;

      G4RotationMatrix mirror_4_sector2_2_rot;
      G4RotationMatrix mirror_3_sector2_2_rot;
      G4RotationMatrix mirror_2_sector2_2_rot;
      G4RotationMatrix mirror_1_sector2_2_rot;

      G4ThreeVector mirror_4_sector2_2_trans;
      G4ThreeVector mirror_3_sector2_2_trans;
      G4ThreeVector mirror_2_sector2_2_trans;
      G4ThreeVector mirror_1_sector2_2_trans;

      // -----------------
      //G4VSolid        * Barrel_sect0half1_solid           ;
      //G4VSolid        * phicut_sect0half1_solid           ;
      //G4VSolid        * Mirror_sect0mirr0half1_solid      ;
      //G4VSolid        * BarrelEllipseCut_sect0mirr0half1_solid ;
      //G4VSolid        * Boxcut12_sect0mirr0half1_solid     ;
      //G4VSolid        * MirrorBoxCut_sect0mirr0half1_solid ;
      //G4VSolid        * MirrorCylinderCut_sect0mirr0half1_solid ;
      //G4VSolid        * mirror_4_sector2_2_solid           ;
      //G4LogicalVolume * mirror_4_sector2_2_log          ;
      //G4LogicalVolume * BarrelEllipseCut_sect0mirr0half1_log;

      //G4VSolid        * Mirror_sect0mirr1half0_solid            ;
      //G4VSolid        * BarrelEllipseCut_sect0mirr1half0_solid  ;
      //G4VSolid        * Boxcut_up_sect0mirr1half0_solid         ;
      //G4VSolid        * Boxcut_down_sect0mirr1half0_solid       ;
      //G4VSolid        * MirrorBoxCut_up_sect0mirr1half0_solid   ;
      //G4VSolid        * MirrorBoxCut_down_sect0mirr1half0_solid ;
      //G4VSolid        * mirror_3_sector2_2_solid                ;
      //G4LogicalVolume * mirror_3_sector2_2_log          ;
      //G4VSolid        * Mirror_sect0mirr2half0_solid            ;
      //G4VSolid        * BarrelEllipseCut_sect0mirr2half0_solid  ;
      //G4VSolid        * Boxcut_up_sect0mirr2half0_solid         ;
      //G4VSolid        * Boxcut_down_sect0mirr2half0_solid       ;
      //G4VSolid        * MirrorBoxCut_up_sect0mirr2half0_solid   ;
      //G4VSolid        * MirrorBoxCut_down_sect0mirr2half0_solid ;
      //G4VSolid        * mirror_2_sector2_2_solid                ;
      //G4LogicalVolume * mirror_2_sector2_2_log          ;
      //G4VSolid        * Mirror_sect0mirr3half0_solid            ;
      //G4VSolid        * BarrelEllipseCut_sect0mirr3half0_solid  ;
      //G4VSolid        * Boxcut_up_sect0mirr3half0_solid         ;
      //G4VSolid        * MirrorBoxCut_up_sect0mirr3half0_solid   ;
      //G4VSolid        * MirrorConeCut_sect0mirr3half0_solid     ;
      //G4VSolid        * mirror_1_sector2_2_solid                ;
      //G4LogicalVolume * mirror_1_sector2_2_log          ;
      //G4VSolid        * Barrel_sect0half1_solid                 ;
      //G4VSolid        * phicut_sect0half1_solid                 ;
      //G4VSolid        * Mirror_sect0mirr0half1_solid            ;
      //G4VSolid        * BarrelEllipseCut_sect0mirr0half1_solid  ;
      //G4VSolid        * Boxcut12_sect0mirr0half1_solid          ;
      //G4VSolid        * MirrorBoxCut_sect0mirr0half1_solid      ;
      //G4VSolid        * MirrorCylinderCut_sect0mirr0half1_solid ;
      //G4VSolid        * mirror_4_sector3_1_solid                ;
      //G4LogicalVolume * mirror_4_sector3_1_log          ;
      //G4VSolid        * Mirror_sect0mirr1half1_solid            ;
      //G4VSolid        * BarrelEllipseCut_sect0mirr1half1_solid  ;
      //G4VSolid        * Boxcut_up_sect0mirr1half1_solid         ;
      //G4VSolid        * Boxcut_down_sect0mirr1half1_solid       ;
      //G4VSolid        * MirrorBoxCut_up_sect0mirr1half1_solid   ;
      //G4VSolid        * MirrorBoxCut_down_sect0mirr1half1_solid ;
      //G4VSolid        * mirror_3_sector3_1_solid                ;
      //G4LogicalVolume * mirror_3_sector3_1_log          ;
      //G4VSolid        * Mirror_sect0mirr2half1_solid            ;
      //G4VSolid        * BarrelEllipseCut_sect0mirr2half1_solid  ;
      //G4VSolid        * Boxcut_up_sect0mirr2half1_solid         ;
      //G4VSolid        * Boxcut_down_sect0mirr2half1_solid       ;
      //G4VSolid        * MirrorBoxCut_up_sect0mirr2half1_solid   ;
      //G4VSolid        * MirrorBoxCut_down_sect0mirr2half1_solid ;
      //G4VSolid        * mirror_2_sector3_1_solid                ;
      //G4LogicalVolume * mirror_2_sector3_1_log          ;
      //G4VSolid        * Mirror_sect0mirr3half1_solid            ;
      //G4VSolid        * BarrelEllipseCut_sect0mirr3half1_solid  ;
      //G4VSolid        * Boxcut_up_sect0mirr3half1_solid         ;
      //G4VSolid        * MirrorBoxCut_up_sect0mirr3half1_solid   ;
      //G4VSolid        * MirrorConeCut_sect0mirr3half1_solid     ;
      //G4VSolid        * mirror_1_sector3_1_solid                ;
      //G4LogicalVolume * mirror_1_sector3_1_log          ;

      G4RotationMatrix mirror_4_sector3_1_rot;
      G4RotationMatrix mirror_3_sector3_1_rot;
      G4RotationMatrix mirror_2_sector3_1_rot;
      G4RotationMatrix mirror_1_sector3_1_rot;
      G4ThreeVector mirror_4_sector3_1_trans;
      G4ThreeVector mirror_3_sector3_1_trans;
      G4ThreeVector mirror_2_sector3_1_trans;
      G4ThreeVector mirror_1_sector3_1_trans;

      G4VSolid        * pmt_4_sector3_1_solid;
      G4LogicalVolume * pmt_4_sector3_1_log;
      G4ThreeVector     pmt_4_sector3_1_trans;
      G4RotationMatrix  pmt_4_sector3_1_rot;

   public:
      void BuildMirrors();

      //HTCCSensitiveDetector * fSensitiveDetector;
      SensitiveRegionDetector       * fSensitiveDetector2;

      G4VSolid        * htccBigGasVolume_solid;
      G4LogicalVolume * htccBigGasVolume_log;
      G4VSolid        * htccEntryDishVolume_solid;
      G4LogicalVolume * htccEntryDishVolume_log;
      G4VSolid        * htccEntryConeVolume_solid;
      G4LogicalVolume * htccEntryConeVolume_log;

      G4VSolid        * htccEntryDishCone_solid;
      G4LogicalVolume * htccEntryDishCone_log;
      G4VSolid        * htcc_solid;
      G4LogicalVolume * htcc_log;

      G4VSolid * fRegion1_solid;
      G4VSolid * fRegion2_solid;
      G4VSolid * fRegion3_solid;

      G4LogicalVolume * fRegion1_log;
      G4LogicalVolume * fRegion2_log;
      G4LogicalVolume * fRegion3_log;

      G4LogicalVolume * fEmptyRegion1_log;
      G4LogicalVolume * fEmptyRegion2_log;
      G4LogicalVolume * fEmptyRegion3_log;

      G4Material * fGasMaterial;

   public:
      HTCCDetectorGeometry();
      ~HTCCDetectorGeometry();

      void BuildLogicalVolumes();
      
      G4VPhysicalVolume * PlacePhysicalVolume(G4LogicalVolume * mother, int sec, int region);
      G4VPhysicalVolume * PlaceParallelPhysicalVolume(G4LogicalVolume * mother, int sec, int region);

};

#endif

