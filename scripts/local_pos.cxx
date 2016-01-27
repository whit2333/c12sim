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

int local_pos(int run_number = 30, int norm = 3000)
{

   gStyle->SetCanvasPreferGL(true);

   TFile * f = new TFile(Form("data/rootfiles/clas12sim%d.root",run_number),"UPDATE");
   TTree * t = (TTree*)gROOT->FindObject("clasdigi_hits");
   if( !t ) return -111;

   //---------------------------------------
   std::vector<int> colors = { 1,2,4,6,8,9,40,41,42,43};
   init_dc_hists();

   THStack * hs0 = new THStack("hsOccupancies","Occupancies");
   THStack * hs  = 0;
   std::vector< std::vector<TH1F*> >    fOccupancies;
   std::vector< std::vector< std::vector<TH1F*> > > fLayerOccupancies;
   std::vector< THStack* >              fStacks;
   std::vector< std::vector< std::vector<TH2F*> > > fLocalPositions;

   for(int sec = 1; sec<=6; sec++) {

      hs = new THStack(Form("hsOcc_%d",sec),Form("Sector %d Occupancies",sec));
      std::vector<TH1F*> sec_occupancies;

      std::vector< std::vector<TH1F*> > sec_layer_occupancies;
      std::vector< std::vector<TH2F*> > sec_layer_localpos;

      for(int superl = 1; superl<=6; superl++){

         TH1F * h1 = new TH1F(Form("Sector %d, SL %d",sec,superl),Form("Sector %d, SL %d",sec,superl),113,0,113);
         sec_occupancies.push_back(h1);
         h1->SetLineColor(colors[superl-1]);
         h1->SetLineWidth(2);
         hs0->Add(h1);
         hs->Add(h1);

         std::vector<TH1F*> sl_occupancies;
         std::vector<TH2F*> sl_localpos;

         for(int layer = 1; layer<=6; layer++){

            TH1F * h2 = new TH1F(Form("Sector_%d_SL_%d_Layer_%d",sec,superl,layer),Form("Sector %d, SL %d, Layer %d",sec,superl,layer),113,0,113);
            sl_occupancies.push_back(h2);
            h2->SetLineColor(colors[superl-1]);
            h2->SetLineWidth(2);

            TH2F * hist= new TH2F(
                  Form("LocalPosSector_%d_SL_%d_Layer_%d",sec,superl,layer),
                  Form("Local pos Sector %d, SL %d, Layer %d",sec,superl,layer),
                  60,0,120,
                  50,-100,100
                  );
            sl_localpos.push_back(hist);

         }
         sec_layer_occupancies.push_back(sl_occupancies);
         sec_layer_localpos.push_back(sl_localpos);

      }
      fLayerOccupancies.push_back(sec_layer_occupancies);
      fLocalPositions.push_back(sec_layer_localpos);
      fOccupancies.push_back(sec_occupancies);
      fStacks.push_back(hs);
   }

   //---------------------------------------------

   clas12::hits::CLAS12HitsEvent * event = 0;
   t->SetBranchAddress("HitsEvent", &event);

   int nEvents = t->GetEntries();

   for(int iEvent =  0 ; iEvent < nEvents; iEvent++) {

      t->GetEntry(iEvent);

      //std::cout << "event : " << iEvent << "\n";

      for(int i = 0; i < event->fDCEvent.fNParticleHits; i++) {

         //std::cout << "hit : " << i << "\n";
         clas12::hits::DriftChamberParticleHit * ahit = event->fDCEvent.GetParticleHit(i);

         if( ahit->fMomentum.E() > 0.01 ) {

            int sec     = ahit->fDCWire.fSector;
            int region  = ahit->fDCWire.fRegion;
            int sl      = ahit->fDCWire.fSuperLayer;
            int lay     = ahit->fDCWire.fLayer;
            int bin     = ahit->fDCWire.fWire + (ahit->fDCWire.fLayer-1)*112;

            clas12::geo::DCSuperLayer sl_id(sec, region, sl);

            DCHist * h  =  fgDCHists[sl_id];
            double val  =  h->GetBinContent(bin) + 1.0/double(norm);
            h->SetBinContent(bin,val);

            // ----------------------
            // 1D hists
            bin = ahit->fDCWire.fWire;
            val = fOccupancies[sec-1][sl-1]->GetBinContent(bin) + 1.0/double(norm);
            fOccupancies[sec-1][sl-1]->SetBinContent(bin, val);

            bin = ahit->fDCWire.fWire;
            val = fLayerOccupancies[sec-1][sl-1][lay-1]->GetBinContent(bin) + 1.0/double(norm);
            fLayerOccupancies[sec-1][sl-1][lay-1]->SetBinContent(bin, val);

            fLocalPositions[sec-1][sl-1][lay-1]->Fill(ahit->fDCWire.fWire, ahit->fGlobalPosition.Y());
            //fLocalPositions[sec-1][sl-1][lay-1]->Fill(ahit->fMomentum.Mag()*1000,ahit->fGlobalPosition.Y() );

         }
      }

   } 


   std::cout << " N events : " << nEvents << std::endl;
   gSystem->mkdir("data/results/local_pos");

   //------------------------------------------------------------

   TCanvas * c0 = 0;
   c0 = new TCanvas();
   c0->Divide(3,2);
   TCanvas * c1 = new TCanvas();
   c1->Divide(3,2);

   for(int i = 0; i<6; i++) {
      c0->cd(i+1);
      fLocalPositions[0][4][i]->Draw("colz");

      c1->cd(i+1);
      fLocalPositions[0][5][i]->Draw("colz");
   }

   return 0;

}

