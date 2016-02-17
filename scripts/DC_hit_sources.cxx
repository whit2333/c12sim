#include <map>
#include "DCWire.h"
#include "CLAS12HitsEvent.h"
#include "TH2Poly.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"


TH1F * get_hist_from_map(std::map<int,TH1F*>& hist_map, int pid){
   if(hist_map.count(pid) == 0) {
   std::cout << " pid" << pid << std::endl;
      hist_map[pid] = new TH1F(Form("fhZ_%d",pid), Form("Z - %d",pid), 300,-10,600);
   }
   return hist_map[pid];
}


int DC_hit_sources(
      int run_number = 0
      )
{

   std::vector<int> runs;// = { 300 , 307, 309 };//, 303, 304, 305, 308, 306 };
   for(int i=300; i<320;i++) {
      runs.push_back(i);
   }
   //std::vector<int> runs;// = { 300 , 307, 309 };//, 303, 304, 305, 308, 306 };
   for(int i=400; i<420;i++) {
      runs.push_back(i);
   }

   std::map<int,TH1F*> pid_hists;

   TChain * t =  new TChain("clasdigi_hits");
   for(auto i : runs) {
      t->Add(Form("data/rootfiles/clas12sim%d.root",i));
   } 

   int aNumber = run_number;//+10000;
   gStyle->SetCanvasPreferGL(true);


   //TFile * f = new TFile(Form("data/rootfiles/clas12sim%d.root",run_number),"UPDATE");
   //TTree * t = (TTree*)gROOT->FindObject("clasdigi_hits");
   if( !t ) {
      std::cout << " TREE NOT FOUND" << std::endl;
      return -111;
   }

   // ------------------------------------------------------

   std::vector<int> colors = { 1,2,4,6,8,9,40,41,42,43};

   //init_dc_hists();
   //THStack * hs0 = new THStack("hsOccupancies","Occupancies");
   //THStack * hs  = 0;
   //std::vector< std::vector<TH1F*> > fOccupancies;
   //std::vector< std::vector< std::vector<TH1F*> > > fLayerOccupancies;
   //std::vector< THStack* >           fStacks;
   //std::vector< TH1F* >              fSLAveraged;

   //std::array<TH2F*,6> fAvgNhitsVsWire;
   //std::array<TH2F*,6> fStdDevNhitsVsWire;
   //std::array<TH2F*,6> fNhitsVsWire;
   //std::array<std::array<std::array<std::array<int,112>,6>,6>,6> nWireHits;
   //std::array<std::array<std::array<std::array<int,112>,6>,6>,6> nWireHitsTotal;
   //std::array<std::array<std::array<std::array<int,112>,6>,6>,6> nWireHitsSquaredTotal;

   TH1F * fhZ              = new TH1F("fhZ", "Z", 300,-10,600);
   TH2F * fhRvsZ           = new TH2F("fhRvsZ", "r vs Z", 300,-10,600,200,0,400);
   TH2F * fhWireNumbervsZ  = new TH2F("fhWireNumbervsZ", "wire vs Z", 300,-10,600,112,1,113);
   TH1F * fhZ2              = new TH1F("fhZ2", "mother Z", 300,-10,600);
   TH2F * fhRvsZ2           = new TH2F("fhRvsZ2", "mother r vs Z", 300,-10,600,200,0,400);
   TH2F * fhWireNumbervsZ2  = new TH2F("fhWireNumbervsZ2", "mother wire vs Z", 300,-10,600,112,1,113);
   //fAvgNhitsVsWireAll      = new TH2F("fAvgNhitsVsWireAll", "Avg Nhits Vs Wire all",
   //      112,1,113,100,0,0.01);


   clas12::hits::CLAS12HitsEvent * event      = 0;
   TClonesArray                  * traj_verts = 0;
   t->SetBranchAddress("HitsEvent",           &event);
   t->SetBranchAddress("TrajectoryVerticies", &traj_verts);

   // --------------------------------------------------------------------------

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
            int pid     = ahit->fPDGCode;


            TParticle * part = 0;
            TParticle * orig_part = 0;
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
            if(part) {
               TH1F * pid_hist = get_hist_from_map(pid_hists,pid);
               pid_hist->Fill(part->Vz());
               fhRvsZ->Fill(part->Vz(),part->R());
               fhZ->Fill(part->Vz());
               fhWireNumbervsZ->Fill(part->Vz(),wire);
            }
            if(orig_part) {

               fhRvsZ2->Fill(orig_part->Vz(),orig_part->R());
               fhZ2->Fill(orig_part->Vz());
               fhWireNumbervsZ2->Fill(orig_part->Vz(),wire);
            }

            //nWireHits[sec-1][sl-1][lay-1][wire-1]++;
            //nWireHitsTotal[sec-1][sl-1][lay-1][wire-1]++;

            //double time_norm = 500.0/time_window_by_sl[sl-1];

            //DCHist * h  =  fgDCHists[sl_id];
            //double val  =  h->GetBinContent(bin) + norm_sim/(time_norm*double(norm));
            //h->SetBinContent(bin,val);


         }
      }
      // End of hit loop
      // ---------------------------

      //for(int sec = 0; sec<6; sec++) {
      //   for(int superl = 0; superl<6; superl++){
      //      for(int layer = 0; layer<6; layer++){
      //         for(int wire = 0; wire<112; wire++) {

      //            fNhitsVsWire[superl]->Fill(   wire+1, nWireHits[sec][superl][layer][wire] );

      //            nWireHitsSquaredTotal[sec][superl][layer][wire] += (nWireHits[sec][superl][layer][wire]*nWireHits[sec][superl][layer][wire]);

      //            //val = [sec][superl][layer][wire]->GetBinContent(bin) + norm_sim/(time_norm*double(norm*6));
      //            //fOccupancies[sec-1][sl-1]->SetBinContent(bin, val);

      //            // all done let's zero this wire
      //            nWireHits[sec][superl][layer][wire] = 0;
      //         }
      //      }
      //   }
      //}

   }
   // End of event loop
   // -----------------------------------------------------

   std::cout << " N events : " << nEvents << std::endl;
   gSystem->mkdir("data/results/DC_hit_sources");

   // -----------------------------------------------------

   TCanvas * c0 = 0;
   TMathText mathtex;
   mathtex.SetTextFont(43);
   mathtex.SetTextSize(20);
   mathtex.SetNDC(true);

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   fhZ->Draw("colz");
   //mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_Z_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_Z_%d.pdf",aNumber));

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   fhRvsZ->Draw("colz");
   //mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_RvsZ_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_RvsZ_%d.pdf",aNumber));

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   fhWireNumbervsZ->Draw("colz");
   //mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_wirevsZ_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_wirevsZ_%d.pdf",aNumber));

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   fhZ2->Draw("colz");
   //mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_Z2_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_Z2_%d.pdf",aNumber));

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   fhRvsZ2->Draw("colz");
   //mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_RvsZ2_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_RvsZ2_%d.pdf",aNumber));

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   fhWireNumbervsZ2->Draw("colz");
   //mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_wirevsZ2_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_wirevsZ2_%d.pdf",aNumber));

   std::cout << " N events : " << nEvents << std::endl;

   c0 = new TCanvas();
   c0->cd();
   int icol = 0;
   THStack * hs = new THStack("z for pids","");
   TLegend * leg = new TLegend(0.7,0.7,0.9,0.9);
   for(auto hp : pid_hists ) {
      hp.second->SetLineColor(colors[icol]);
      hp.second->SetLineWidth(2);
      leg->AddEntry(hp.second, Form("%d",hp.first),"l");
      hs->Add(hp.second);
      icol++;
   }
   hs->Draw("nostack");
   leg->Draw();
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_pid_Z_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_pid_Z_%d.pdf",aNumber));

   gPad->SetLogy(true);
   c0->Update();
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_pid_Z_log_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_hit_sources/DC_hit_sources_pid_Z_log_%d.pdf",aNumber));

   return 0;
}

