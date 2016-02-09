//______________________________________________________________________________
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

class DCHist : public TH2Poly {
   public: 

      void BetterHoneycomb(Double_t xstart, Double_t ystart, Double_t a, Int_t k, Int_t s)
      {
         // Bins the histogram using a honeycomb structure
         // Add the bins
         Double_t numberOfHexagonsInTheRow;
         Double_t x[6], y[6];
         Double_t xtemp, ytemp;
         double yoffset = 0;

         numberOfHexagonsInTheRow = k;
         xtemp = xstart; // Resets the temp variable
         ytemp = ystart;

         for (int kCounter = 0; kCounter < 2*numberOfHexagonsInTheRow; kCounter++) {

            ytemp = ystart;
            for (int sCounter = 0; sCounter < s; sCounter++) {

               // Go around the hexagon
               x[0] = xtemp + a/2.0;
               y[0] = ytemp + yoffset ;
               x[1] = x[0] + a;
               y[1] = y[0];

               x[2] = x[1] + a/2.0;
               y[2] = y[1] + a*TMath::Sqrt(3)/2.0;

               x[3] = x[1];
               y[3] = y[2] + a*TMath::Sqrt(3)/2.0;
               x[4] = x[0];
               y[4] = y[3];

               x[5] = x[0] - a/2.0;
               y[5] = y[2] ;

               this->AddBin(6, x, y);

               //std::cout <<  sCounter << "  " << kCounter << std::endl;

               // Go right
               ytemp += a*TMath::Sqrt(3);
            }

            xtemp += 1.5*a;

            if (kCounter%2 == 0) {yoffset = a*TMath::Sqrt(3)/2.0;    }//;
            else                 {yoffset = 0;        }//TMath::Sqrt(3)/2.0;
            // Increment the starting position
            //yloop += a*TMath::Sqrt(3);// 1.5*a;
         }
      }
  
};

std::map<clas12::geo::DCSuperLayer,DCHist*> fgDCHists;
//______________________________________________________________________________

DCHist * dc_sl_hist(const clas12::geo::DCSuperLayer& sl)
{
   if( fgDCHists.find(sl) == fgDCHists.end() ) {

      fgDCHists[sl]  = new DCHist();
      double x_start = 0.0;
      double y_start = 0.0;
      double a_side  = 0.5;
      if( (sl.fSuperLayer-1)%2 == 1 ) {
         y_start += 0.0;
         x_start += 5.0;
      }
      fgDCHists[sl]->BetterHoneycomb(x_start, y_start, a_side, 3,112);
      //gStyle->SetPalette(1);
   }
   return fgDCHists[sl];
}
//______________________________________________________________________________

void init_dc_hists(){
   for(int sector = 1; sector<=6; sector++) {
      for(int superlayer = 1; superlayer<=6; superlayer++){
      //for(int region = 1; region<=3; region++) {
         //for(int superlayer = (region-1)*2+1; superlayer<=(region-1)*2+2; superlayer++) {
            clas12::geo::DCSuperLayer sl(sector, (superlayer-1)/2 + 1, superlayer);
            dc_sl_hist(sl);
         //}
      }
   }
}
//______________________________________________________________________________

int DC_occupancy_sectors(
      int run_number = 30, int norm = 3000, 
      int n_gen = 1e6, int n_sim=500000
      ) {

   gStyle->SetCanvasPreferGL(true);

   std::vector<double> time_window_by_sl = {250.0,250.0, 400.0, 400.0, 500.0, 500.0};
   double norm_sim = n_gen/n_sim;
   std::cout << " norm_sim = " << norm_sim << std::endl;

   TFile * f = new TFile(Form("data/rootfiles/clas12sim%d.root",run_number),"UPDATE");
   TTree * t = (TTree*)gROOT->FindObject("clasdigi_hits");
   if( !t ) return -111;

   //---------------------------------------
   std::vector<int> colors = { 1,2,4,6,8,9,40,41,42,43};
   init_dc_hists();

   std::array<std::array<std::array<std::array<int,112>,6>,6>,6> nWireHits;
   for(int sec = 0; sec<6; sec++) {
      for(int superl = 0; superl<6; superl++){
         for(int layer = 0; layer<6; layer++){
            for(int wire = 0; wire<112; wire++) {
               nWireHits[sec][superl][layer][wire] = 0;
            }
         }
      }
   }

   THStack * hs0 = new THStack("hsOccupancies","Occupancies");
   THStack * hs  = 0;
   std::vector< std::vector<TH1F*> > fOccupancies;
   std::vector< std::vector< std::vector<TH1F*> > > fLayerOccupancies;
   std::vector< THStack* >           fStacks;
   std::vector< TH1F* >              fSLAveraged;
   for(int sec = 1; sec<=6; sec++) {

      hs = new THStack(Form("hsOcc_%d",sec),Form("Sector %d Occupancies",sec));
      std::vector<TH1F*> sec_occupancies;
      std::vector< std::vector<TH1F*> > sec_layer_occupancies;

      for(int superl = 1; superl<=6; superl++){

         TH1F * h1 = new TH1F(Form("Sector %d, SL %d",sec,superl),Form("Sector %d, SL %d",sec,superl),113,0,113);
         sec_occupancies.push_back(h1);
         h1->SetLineColor(colors[superl-1]);
         h1->SetLineWidth(2);
         hs0->Add(h1);
         hs->Add(h1);

         std::vector<TH1F*> sl_occupancies;

         for(int layer = 1; layer<=6; layer++){

            int sec_chan = (superl-1)*6  + layer-1;

            TH1F * h2 = new TH1F(Form("Sector_%d_SL_%d_Layer_%d",sec,superl,layer),Form("Sector %d, SL %d, Layer %d",sec,superl,layer),113,0,113);
            sl_occupancies.push_back(h2);
            h2->SetLineColor(colors[superl-1]);
            h2->SetLineWidth(2);

            if(sec ==1)  {
               TH1F * hist2 = new TH1F(
                     Form("SL_%d_Layer_%d",superl,layer),
                     Form("SL %d, Layer %d",superl,layer),
                     113,0,113);
               fSLAveraged.push_back(hist2);
            }
         }
         sec_layer_occupancies.push_back(sl_occupancies);

      }
      fLayerOccupancies.push_back(sec_layer_occupancies);
      fOccupancies.push_back(sec_occupancies);
      fStacks.push_back(hs);
   }

   clas12::hits::CLAS12HitsEvent * event = 0;
   t->SetBranchAddress("HitsEvent", &event);

   // --------------------------------------------------------------------------

   int nEvents = t->GetEntries();
   std::cout << " nEvents = " << nEvents << std::endl;

   for(int iEvent =  0 ; iEvent < nEvents && iEvent < n_sim ; iEvent++) {

      t->GetEntry(iEvent);

      //std::cout << "event : " << iEvent << "\n";

      for(int i = 0; i < event->fDCEvent.fNParticleHits; i++) {

         //std::cout << "hit : " << i << "\n";
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

            if( nWireHits[sec-1][sl-1][lay-1][wire-1] == 1 ){
               double time_norm = 500.0/time_window_by_sl[sl-1];

               clas12::geo::DCSuperLayer sl_id(sec, region, sl);

               DCHist * h  =  fgDCHists[sl_id];
               double val  =  h->GetBinContent(bin) + norm_sim/(time_norm*double(norm));
               h->SetBinContent(bin,val);

               // ----------------------
               // 1D hists
               bin = ahit->fDCWire.fWire;
               val = fOccupancies[sec-1][sl-1]->GetBinContent(bin) + norm_sim/(time_norm*double(6*norm));
               fOccupancies[sec-1][sl-1]->SetBinContent(bin, val);

               bin = ahit->fDCWire.fWire;
               val = fLayerOccupancies[sec-1][sl-1][lay-1]->GetBinContent(bin) + norm_sim/(time_norm*double(norm));
               fLayerOccupancies[sec-1][sl-1][lay-1]->SetBinContent(bin, val);

               int sec_chan = 6*(sl-1) + (lay-1);
               val = fSLAveraged[sec_chan]->GetBinContent(bin) + norm_sim/(time_norm*double(6*norm));
               fSLAveraged[sec_chan]->SetBinContent(bin, val);

            }

         }
      }

      // reset the wire counts
      for(int sec_ = 0; sec_<6; sec_++) {
         for(int superl_ = 0; superl_<6; superl_++){
            for(int layer_ = 0; layer_<6; layer_++){
               for(int wire_ = 0; wire_<112; wire_++) {
                  nWireHits[sec_][superl_][layer_][wire_] = 0;
               }
            }
         }
      }
   } 


   std::cout << " N events : " << nEvents << std::endl;
   gSystem->mkdir("data/results/DC_occupancy");

   //---------------------------------------------------------

   TMathText mathtex; 
   mathtex.SetTextFont(43);
   mathtex.SetTextSize(20);
   mathtex.SetNDC(true);

   TCanvas * c = 0;

   TLatex Tl; 
   Tl.SetTextFont(43);
   Tl.SetTextSize(20);
   Tl.SetNDC(true);

   gStyle->SetPadTopMargin    ( 0.05 );
   gStyle->SetPadBottomMargin ( 0.01 );
   gStyle->SetPadLeftMargin   ( 0.01 );
   gStyle->SetPadRightMargin  ( 0.25 );

   for(int region = 1; region<=3; region++) {
      c = new TCanvas(Form("region%d",region),Form("region%d",region),1200,900);
      c->SetBottomMargin(0.0001);
      c->SetTopMargin(0.0001);
      c->SetLeftMargin(0.0001);
      c->SetRightMargin(0.01);
      c->Divide(6,1,0.0001,0.0001);
      double min = 0;
      double max = 2.0;
      bool first = true;
      for(int sector = 1; sector<=6; sector++) {

         c->cd( sector );
         //TH2F * temp_hist = new TH2F("temphist","temp",100,-1,12,100,-1,110);
         //temp_hist->Draw("");
         THStack * stack = new THStack(
               Form("Sector%d_region%d ",sector,region),
               Form("Sector %d, Region %d",sector,region) );

         for(int superlayer = (region-1)*2+1; superlayer<=(region-1)*2+2; superlayer++) {
            clas12::geo::DCSuperLayer sl_id(sector, region, superlayer);
            if( fgDCHists.find(sl_id) != fgDCHists.end() ){
               //fgDCHists[sl_id]->Draw("colz,same");
               stack->Add(fgDCHists[sl_id]);
               //double temp =  fgDCHists[sl_id]->GetMaximum();
               //if(temp > max) max = temp;
               if(first){
                  min =  0.0;//fgDCHists[sl_id]->GetMinimum();
                  fgDCHists[sl_id]->SetMinimum(min);
                  //max =  fgDCHists[sl_id]->GetMaximum();
                  fgDCHists[sl_id]->SetMaximum(max);
                  first = false;
               } else {
                  fgDCHists[sl_id]->SetMinimum(min);
                  fgDCHists[sl_id]->SetMaximum(max);
               }

            }
         }
         stack->SetMaximum(max);
         stack->Draw("nostack,colz");
      }
      c->cd(0);
      mathtex.DrawMathText(0.4,0.98,"\\mathscr{L} = 1.3 \\times 10^{35} [cm^{-1}s^{-1}]");
      c->Update();
      c->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_sectors_region%d_%d.png",region,run_number));
      c->SaveAs(Form("data/results/DC_occupancy/DC_occupancy_sectors_region%d_%d.pdf",region,run_number));
   }


   std::cout << " N events : " << nEvents << std::endl;

   //Tl.DrawLatex(.01, .5, Form("#splitline{Sector %d}{Region %d}",i+1,region));

   return 0;
}

