void eve_geometry()
{
   TEveManager::Create();

   TGeoManager::Import("beamline.gdml");
   gGeoManager->DefaultColors();

   TGeoNode* node1 = gGeoManager->GetTopNode();
   TEveGeoTopNode* its = new TEveGeoTopNode(gGeoManager, node1);
   gEve->AddGlobalElement(its);

   gGeoManager = gEve->GetGeometry("beamline.gdml");
   gGeoManager->DefaultColors();

   gEve->FullRedraw3D(kTRUE);

   TGLViewer *v = gEve->GetDefaultGLViewer();

   //v->GetClipSet()->SetClipType(TGLClip::EType::kClipPlane);
   v->RefreshPadEditor(v);

   v->CurrentCamera().RotateRad(-.7, 0.5);
   v->DoDraw();

   return 0;
}
