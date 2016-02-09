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

int DC_occupancy_curlers(
      int run_number = 30,    int norm  = 3000,
      int n_gen      = 10e6,  int n_sim = 500000
      ) {

   int aNumber = run_number;//+10000;
   gStyle->SetCanvasPreferGL(true);

   std::vector<double> time_window_by_sl = {250.0,250.0, 400.0, 400.0, 500.0, 500.0};

   double norm_sim = n_gen/n_sim;
   std::cout << " norm_sim = " << norm_sim << std::endl;

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

   std::array<TH2F*,6> fAvgNhitsVsWire;
   std::array<TH2F*,6> fStdDevNhitsVsWire;
   std::array<TH2F*,6> fNhitsVsWire;
   std::array<std::array<std::array<std::array<int,112>,6>,6>,6> nWireHits;
   std::array<std::array<std::array<std::array<int,112>,6>,6>,6> nWireHitsTotal;
   std::array<std::array<std::array<std::array<int,112>,6>,6>,6> nWireHitsSquaredTotal;

   fStdDevNhitsVsWireAll   = new TH2F("fStdDevNhitsVsWireAll", "Nhits std deviation Vs Wire all",
         112,1,113,100,0,0.5);
   fAvgNhitsVsWireAll      = new TH2F("fAvgNhitsVsWireAll", "Avg Nhits Vs Wire all",
         112,1,113,100,0,0.01);

   for(int sec = 0; sec<6; sec++) {
      fNhitsVsWire[sec]    = new TH2F(Form("fNhitsVsWire%d",sec+1), Form("Nhits Vs Wire %d",sec+1),
            112,1,113,100,1,101);
      fAvgNhitsVsWire[sec] = new TH2F(Form("fAvgNhitsVsWire%d",sec+1), Form("Avg Nhits Vs Wire %d",sec+1),
            112,1,113,100,0,0.01);
      fStdDevNhitsVsWire[sec] = new TH2F(Form("fStdDevNhitsVsWire%d",sec+1), Form("std dev of Nhits Vs Wire %d",sec+1),
         112,1,113,100,0,0.5);
      for(int superl = 0; superl<6; superl++){
         for(int layer = 0; layer<6; layer++){
            for(int wire = 0; wire<112; wire++) {
               nWireHits[sec][superl][layer][wire] = 0;
               nWireHitsTotal[sec][superl][layer][wire] = 0;
               nWireHitsSquaredTotal[sec][superl][layer][wire] = 0;
            }
         }
      }
   }

   clas12::hits::CLAS12HitsEvent * event = 0;
   t->SetBranchAddress("HitsEvent", &event);

   // --------------------------------------------------------------------------

   int nEvents = t->GetEntries();
   std::cout << " nEvents = " << nEvents << std::endl;

   for(int iEvent =  0 ; iEvent < nEvents && iEvent < n_sim; iEvent++) {

      t->GetEntry(iEvent);

      if( iEvent%1000 == 0 ) std::cout << "event : " << iEvent << "\n";

      for(int i = 0; i < event->fDCEvent.fNParticleHits; i++) {

         clas12::hits::DriftChamberParticleHit * ahit = event->fDCEvent.GetParticleHit(i);

         if( ( ahit->fGlobalPosition.T() < 500.0     ) && // 500 ns
             ( ahit->fMomentum.E()       >  500.0e-9) ) {  // 500 eV threshold

            int sec     = ahit->fDCWire.fSector;
            int region  = ahit->fDCWire.fRegion;
            int sl      = ahit->fDCWire.fSuperLayer;
            int lay     = ahit->fDCWire.fLayer;
            int wire    = ahit->fDCWire.fWire;
            int bin     = wire + (ahit->fDCWire.fLayer-1)*112;

            nWireHits[sec-1][sl-1][lay-1][wire-1]++;
            nWireHitsTotal[sec-1][sl-1][lay-1][wire-1]++;

            //double time_norm = 500.0/time_window_by_sl[sl-1];

            //DCHist * h  =  fgDCHists[sl_id];
            //double val  =  h->GetBinContent(bin) + norm_sim/(time_norm*double(norm));
            //h->SetBinContent(bin,val);


         }
      }

      for(int sec = 0; sec<6; sec++) {
         for(int superl = 0; superl<6; superl++){
            for(int layer = 0; layer<6; layer++){
               for(int wire = 0; wire<112; wire++) {

                  fNhitsVsWire[superl]->Fill(   wire+1, nWireHits[sec][superl][layer][wire] );

                  nWireHitsSquaredTotal[sec][superl][layer][wire] += (nWireHits[sec][superl][layer][wire]*nWireHits[sec][superl][layer][wire]);

                  //val = [sec][superl][layer][wire]->GetBinContent(bin) + norm_sim/(time_norm*double(norm*6));
                  //fOccupancies[sec-1][sl-1]->SetBinContent(bin, val);

                  // all done let's zero this wire
                  nWireHits[sec][superl][layer][wire] = 0;
               }
            }
         }
      }
      // End of hit loop
      // ---------------------------

   }
   // End of event loop
   // -----------------------------------------------------

   std::cout << " N events : " << nEvents << std::endl;
   gSystem->mkdir("data/results/DC_occupancy_curlers");

   // -----------------------------------------------------

   TCanvas * c0 = 0;
   TMathText mathtex;
   mathtex.SetTextFont(43);
   mathtex.SetTextSize(20);
   mathtex.SetNDC(true);

   // -----------------------------------------------------
   // Calculate the average number of hits
   TGraph * gr_mean     =  new TGraph(6*6*6*112);
   TGraph * gr_variance =  new TGraph(6*6*6*112);
   TGraph * gr_stddev   =  new TGraph(6*6*6*112);
   std::array<std::array<std::array<std::array<double,112>,6>,6>,6> nWireHitsAverage;
   std::array<std::array<std::array<std::array<double,112>,6>,6>,6> nWireHitsSquaredAverage;
   std::array<std::array<std::array<std::array<double,112>,6>,6>,6> nWireHitsVariance;
   int i_point = 0;
   for(int sec = 0; sec<6; sec++) {
      for(int superl = 0; superl<6; superl++){
         for(int layer = 0; layer<6; layer++){
            for(int wire = 0; wire<112; wire++) {

               double mean = double(nWireHitsTotal[sec][superl][layer][wire])/double(n_sim);
               double x2   = double(nWireHitsSquaredTotal[sec][superl][layer][wire])/double(n_sim);
               double squared_diff = x2 - mean*mean;
               double var = (double(n_sim)/double(n_sim-1))*squared_diff;

               nWireHitsAverage[sec][superl][layer][wire]        = mean;
               nWireHitsSquaredAverage[sec][superl][layer][wire] = x2;
               nWireHitsVariance[sec][superl][layer][wire] = var;

               gr_mean->SetPoint(    i_point,wire+1, mean);
               gr_variance->SetPoint(i_point,wire+1, var);
               gr_stddev->SetPoint(  i_point,wire+1, TMath::Sqrt(var) );

               //double val = fAvgNhitsVsWire[superl]->GetBinContent(wire+1) + nWireHitsAverage[sec][superl][layer][wire];
               fAvgNhitsVsWire[superl]->Fill( wire+1, mean );
               fStdDevNhitsVsWire[superl]->Fill( wire+1, TMath::Sqrt(var) );
               fAvgNhitsVsWireAll->Fill(      wire+1, mean );
               fStdDevNhitsVsWireAll->Fill(   wire+1, TMath::Sqrt(var) );


               i_point++;
            }
         }
      }
   }

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   gr_mean->SetMarkerStyle(20);
   gr_mean->SetMarkerSize(0.5);
   gr_mean->Draw("ap");
   gr_mean->GetXaxis()->SetTitle("wire number");
   gr_mean->GetYaxis()->SetTitle("#LT N_{hit} #GT");
   mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_avg_nhits_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_avg_nhits_%d.pdf",aNumber));

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   gr_variance->SetMarkerStyle(20);
   gr_variance->SetMarkerSize(0.5);
   gr_variance->Draw("ap");
   gr_variance->GetXaxis()->SetTitle("wire number");
   gr_variance->GetYaxis()->SetTitle("Var(N_{hit})");
   mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_var_nhits_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_var_nhits_%d.pdf",aNumber));

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   gr_stddev->SetMarkerStyle(20);
   gr_stddev->SetMarkerSize(0.5);
   gr_stddev->Draw("ap");
   gr_stddev->GetXaxis()->SetTitle("wire number");
   gr_stddev->GetYaxis()->SetTitle("N_{hit} standard deviation");
   mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_stddev_nhits_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_stddev_nhits_%d.pdf",aNumber));

   //---------------------------------------------------------
   //
   for(int sl = 1; sl<=6; sl++) {
      c0 = new TCanvas();
      c0->cd();
      fNhitsVsWire[sl-1]->Draw("colz");
      fNhitsVsWire[sl-1]->GetXaxis()->SetTitle("wire number");
      mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
      c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_%d_%d.png",sl,aNumber));
      c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_%d_%d.pdf",sl,aNumber));

      c0 = new TCanvas();
      c0->cd();
      fAvgNhitsVsWire[sl-1]->Draw("colz");
      fAvgNhitsVsWire[sl-1]->GetXaxis()->SetTitle("wire number");
      mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
      c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_avg_%d_%d.png",sl,aNumber));
      c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_avg_%d_%d.pdf",sl,aNumber));

      c0 = new TCanvas();
      c0->cd();
      fStdDevNhitsVsWire[sl-1]->Draw("colz");
      fStdDevNhitsVsWire[sl-1]->GetXaxis()->SetTitle("wire number");
      mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
      c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_StdDev_%d_%d.png",sl,aNumber));
      c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_StdDev_%d_%d.pdf",sl,aNumber));
   }

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   fAvgNhitsVsWireAll->Draw("colz");
   fAvgNhitsVsWireAll->GetXaxis()->SetTitle("wire number");
   mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_avg_all_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_avg_all_%d.pdf",aNumber));

   //---------------------------------------------------------
   //
   c0 = new TCanvas();
   c0->cd();
   fStdDevNhitsVsWireAll->Draw("colz");
   fStdDevNhitsVsWireAll->GetXaxis()->SetTitle("wire number");
   mathtex.DrawMathText(0.4,0.83,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
   c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_StdDev_all_%d.png",aNumber));
   c0->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_curlers_StdDev_all_%d.pdf",aNumber));


   std::cout << " N events : " << nEvents << std::endl;

   return 0;
}

