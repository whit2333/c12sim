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

int DC_hit_sources(
      int run_number = 0
      )
{

   int aNumber = run_number;//+10000;
   gStyle->SetCanvasPreferGL(true);


   TFile * f = new TFile(Form("data/rootfiles/clas12sim%d.root",run_number),"UPDATE");
   TTree * t = (TTree*)gROOT->FindObject("clasdigi_hits");
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

   TH2F * fhRvsZ           = new TH2F("fhRvsZ", "r vs Z", 100,0,400,100,0,300);
   TH2F * fhWireNumbervsZ  = new TH2F("fhWireNumbervsZ", "r vs Z", 100,0,400,112,1,113);
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

            TParticle * part = 0;
            int ntraj = traj_verts->GetEntries();
            for(int ipart = 0; ipart<ntraj; ipart++){
               TParticle * avert = (TParticle*)((*traj_verts)[ipart]);
               if(avert->GetSecondMother() == trk_id) {
                  part = avert;
                  break;
               }
            }
            if(part) {

               fhRvsZ->Fill(part->Vz(),part->R());

               fhWireNumbervsZ->Fill(part->Vz(),wire);
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

   std::cout << " N events : " << nEvents << std::endl;

   return 0;
}

