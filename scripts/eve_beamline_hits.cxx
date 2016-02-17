void eve_beamline_hits(int nhits = 100,
      int run_number = 0
      )
{

   std::vector<int> runs = { 307, 309 };//, 303, 304, 305, 308, 306 };

   TChain * t =  new TChain("clasdigi_hits");
   for(auto i : runs) {
      t->Add(Form("data/rootfiles/clas12sim%d.root",i));
   } 

   int aNumber = run_number;//+10000;
   gStyle->SetCanvasPreferGL(true);

   if( !t ) {
      std::cout << " TREE NOT FOUND" << std::endl;
      return -111;
   }

   TEveStraightLineSet* ls = new TEveStraightLineSet();
   TEveStraightLineSet* ls2 = new TEveStraightLineSet();

   // ------------------------------------------------------

   std::vector<int> colors = { 1,2,4,6,8,9,40,41,42,43};

   TH1F * fhZ              = new TH1F("fhZ", "Z", 100,-10,600);
   TH2F * fhRvsZ           = new TH2F("fhRvsZ", "r vs Z", 100,-10,600,100,0,400);
   TH2F * fhWireNumbervsZ  = new TH2F("fhWireNumbervsZ", "wire vs Z", 100,-10,600,112,1,113);
   //fAvgNhitsVsWireAll      = new TH2F("fAvgNhitsVsWireAll", "Avg Nhits Vs Wire all",
   //      112,1,113,100,0,0.01);


   clas12::hits::CLAS12HitsEvent * event      = 0;
   TClonesArray                  * traj_verts = 0;
   t->SetBranchAddress("HitsEvent",           &event);
   t->SetBranchAddress("TrajectoryVerticies", &traj_verts);

   // --------------------------------------------------------------------------

   int ihit = 0;
   int nEvents = t->GetEntries();
   std::cout << " nEvents = " << nEvents << std::endl;

   for(int iEvent =  0 ; iEvent < nEvents ; iEvent++) {

      t->GetEntry(iEvent);

      if( iEvent%1000 == 0 ) std::cout << "event : " << iEvent << "\n";

      for(int i = 0; i < event->fDCEvent.fNParticleHits; i++) {

         clas12::hits::DriftChamberParticleHit * ahit = event->fDCEvent.GetParticleHit(i);

         if( ( ahit->fGlobalPosition.T() < 500.0     ) && // 500 ns
             ( ahit->fMomentum.E()       >  500.0e-9) ) {  // 500 eV threshold

            int trk_id  = ahit->fTrackID;
            int sec     = ahit->fDCWire.fSector;
            int region  = ahit->fDCWire.fRegion;
            int sl      = ahit->fDCWire.fSuperLayer;
            int lay     = ahit->fDCWire.fLayer;
            int wire    = ahit->fDCWire.fWire;
            int bin     = wire + (ahit->fDCWire.fLayer-1)*112;

            int trk_mother_id  = -1;

            TParticle * part = 0;
            TParticle * orig_part = 0;
            TParticle * track_part = 0;
            int ntraj = traj_verts->GetEntries();
            for(int ipart = 0; ipart<ntraj; ipart++){
               TParticle * avert = (TParticle*)((*traj_verts)[ipart]);
               if(avert->GetFirstMother() == trk_id) {
                  part = avert;
               }
               if(avert->GetSecondMother() == trk_id) {
                  orig_part = avert;
               }
            }
            //for(int ipart = 0; ipart<ntraj; ipart++){
            //   TParticle * avert = (TParticle*)((*traj_verts)[ipart]);
            //   if(avert->GetFirstMother() == trk_mother_id) {
            //      track_part = avert;
            //      break;
            //   }
            //}


               fhRvsZ->Fill(part->Vz(),part->R());
               fhZ->Fill(part->Vz());
               fhWireNumbervsZ->Fill(part->Vz(),wire);

            if(part && orig_part) {
               //if(part->Vz() > 10.0){
                  ls->AddLine( part->Vx(), part->Vy(), part->Vz(),
                        ahit->fPosition.X(), 
                        ahit->fPosition.Y(), 
                        ahit->fPosition.Z());
               //}

               //if(orig_part->Vz() > 10.0){
                  ls2->AddLine( orig_part->Vx(), orig_part->Vy(), orig_part->Vz(),
                        part->Vx(), part->Vy(), part->Vz());
                  ls2->AddMarker( part->Vx(), part->Vy(), part->Vz());
                  //}
                  ihit++;
            }
         }
         if(ihit>nhits) break;
      }
      // End of hit loop
      // ---------------------------
   }

   TEveManager::Create();

   TGeoManager::Import("beamline.gdml");
   gGeoManager->DefaultColors();

   TGeoNode* node1 = gGeoManager->GetTopNode();
   TEveGeoTopNode* its = new TEveGeoTopNode(gGeoManager, node1);
   gEve->AddGlobalElement(its);

   gGeoManager = gEve->GetGeometry("beamline.gdml");
   gGeoManager->DefaultColors();

   gEve->FullRedraw3D(kTRUE);


   //for(Int_t i = 0; i<nlines; i++)
   //{
   //   ls->AddLine( r.Uniform(-s,s), r.Uniform(-s,s), r.Uniform(-s,s),
   //         r.Uniform(-s,s), r.Uniform(-s,s), r.Uniform(-s,s));
   //   // add random number of markers
   //   Int_t nm = Int_t(nmarkers* r.Rndm());
   //   for(Int_t m = 0; m < nm; m++) {
   //      ls->AddMarker(i, r.Rndm());
   //   }
   //}

   ls->SetMarkerSize(1.5);
   ls->SetMarkerStyle(4);
   ls2->SetMarkerColor(5);
   ls->SetLineColor(2);
   ls2->SetLineColor(4);

   gEve->AddElement(ls2);
   gEve->AddElement(ls);
   gEve->Redraw3D();

   TGLViewer *v = gEve->GetDefaultGLViewer();

   //v->GetClipSet()->SetClipType(TGLClip::EType::kClipPlane);
   v->RefreshPadEditor(v);

   v->CurrentCamera().RotateRad(-.7, 0.5);
   v->DoDraw();

   return 0;
}
