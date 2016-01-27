#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/IFunction.h"
#include "TVector3.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
//#include "Minuit2/MinimumPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"

class ParametricHelix {
   public:
      TVector3 fV = {0.0,0.0,0.0}; // Vertex position
      double fC   = 1.0;
      double fR   = 1.0;

      void SetHelixParameters(double r, double c) { fR = r; fC = c; }
      void SetVertex(const TVector3 v)  { fV = v; }

      TVector3 GetVertex()    const { return fV; }
      double   GetR()         const { return fR; }
      double   GetC()         const { return fC; }

      TVector3 operator()(double t) const {
         using namespace TMath;
         return TVector3(
               fV.x() + fR*Cos(t),
               fV.y() + fR*Sin(t),
               fV.z() + fC*t 
               );
      }

};

class HelixClosestPointFinder: public ROOT::Math::IGradientFunctionOneDim
{
   private:
      //TVector3 fV = {0.0,0.0,0.0}; // Vertex position
      mutable ParametricHelix fHelix;
      mutable TVector3 fX = {0.0,1.1,TMath::Pi()/2.0}; // Measured Position for wich we are trying to find the point on the helix

   public:
      void SetHelixParameters(double r, double c) const { fHelix.SetHelixParameters(r,c); }
      void SetVertex(const TVector3 v) const { fHelix.SetVertex(v); }

      TVector3 GetDataPoint() const { return fX; }
      TVector3 GetVertex()    const { return fHelix.GetVertex(); }
      double   GetR()         const { return fHelix.GetR(); }
      double   GetC()         const { return fHelix.GetC(); }

      double DoEval(double ti) const
      {
         using namespace TMath;
         //2 (c (c ti + z0 - z1) + r (y0 - y1) Cos[ti] - r (x0 - x1) Sin[ti])
         double r = GetR();
         double c = GetC();
         TVector3 V = GetVertex();
         TVector3 X = GetDataPoint();
         double res = 2.0*(c*(c*ti + V.z() - X.z()) 
               + r*(V.y() - X.y())*Cos(ti) - r*(V.x() - X.x())*Sin(ti));
         return res;
      }

      double DoDerivative(double ti) const
      {
         using namespace TMath;
         //2 (c^2 - r (x0 - x1) Cos[ti] - r (y0 - y1) Sin[ti])
         double r = GetR();
         double c = GetC();
         TVector3 V = GetVertex();
         TVector3 X = GetDataPoint();
         double res = 2.0*(c*c - r*(V.x() - X.x())*Cos(ti) - r*(V.y() - X.y())*Sin(ti));
         return res;
      }

      ROOT::Math::IBaseFunctionOneDim* Clone() const
      {
         return new HelixClosestPointFinder();
      }

      TVector3 GetPoint(const TVector3& p) const {
         return fHelix( FindPoint(fX) );
      }

      double FindPoint(const TVector3& p) const {
         fX = p;
         return FindPoint();
      }

      double FindPoint() const {
         // Create the Inktegrator
         ROOT::Math::Roots::Bisection brf;

         // Set parameters of the method
         brf.SetFunction( *this, 0, TMath::TwoPi() );
         brf.Solve();

         cout << brf.Root()/CLHEP::degree << endl;
         return brf.Root();
      }

}; 

class HelixFcn : public ROOT::Minuit2::FCNBase {
   protected:
      mutable HelixClosestPointFinder fHelixPointFinder;
      //mutable ParametricHelix         fHelix;
      mutable std::vector<TVector3> theMeasurements;
      mutable std::vector<TVector3> thePositions;
      mutable std::vector<double>   theMVariances;
      double  theErrorDef;

   public:
      HelixFcn(): theErrorDef(1.) {
         theMeasurements.clear();
         thePositions.clear();
         theMVariances.clear();
      }
      ~HelixFcn() {}

      void Clear() {
         theMeasurements.clear();
         thePositions.clear();
         theMVariances.clear();
      }

      virtual double Up() const {return theErrorDef;}

      virtual double operator()(const std::vector<double>& par) const {
         assert(par.size() == 3);

         fHelixPointFinder.SetHelixParameters(par[0], par[1]);
         fHelixPointFinder.SetVertex({0.0,0.0,par[2]});

         thePositions.clear();
         theMVariances.clear();
         for( auto ameas : theMeasurements) {
            thePositions.push_back( fHelixPointFinder.GetPoint(ameas) );
            theMVariances.push_back(0.001);
         }
         double chi2 = 0.;
         for(unsigned int n = 0; n < theMeasurements.size(); n++) {
            TVector3 p0 = thePositions[n];
            TVector3 p1 = theMeasurements[n];
            //p0.Print();
            //p1.Print();
            double   sig = theMVariances[n];
            double  numerator = (p0-p1).Mag2();
            chi2 += numerator/sig;
         }
            std::cout << chi2 << std::endl;
         return chi2;
      }

      void SetMeasurements(const std::vector<TVector3>& m) const {return theMeasurements = m;}
      void AddMeasurement( const TVector3& m) const {return theMeasurements.push_back(m);}

      std::vector<TVector3> measurements() const {return theMeasurements;}
      std::vector<TVector3> positions() const {return thePositions;}
      std::vector<double>   variances() const {return theMVariances;}
      void setErrorDef(double def) {theErrorDef = def;}
};

int HelixFit(int run_number = 1001)
{

   //TFile * f = new TFile(Form("data/rootfiles/clas12sim%d.root",run_number),"UPDATE");
   //TTree * t = (TTree*)gROOT->FindObject("clasdigi_hits");
   //if( !t ) return -111;

   //clas12::hits::CLAS12HitsEvent * event = 0;
   //t->SetBranchAddress("HitsEvent", &event);

   TRandom3 rand(0);
   ParametricHelix helix;

   std::vector<TVector3> test_points;
   std::vector<double>   t_points;

   // --------------------------------------------------------------------------

   //int nEvents = t->GetEntries();

   //for(int iEvent =  0 ; iEvent < nEvents; iEvent++) {

   //   t->GetEntry(iEvent);
   //   //std::cout << "event : " << iEvent << "\n";

   //   for(int i = 0; i < event->fDCEvent.fNParticleHits; i++) {

   //      //std::cout << "hit : " << i << "\n";
   //      clas12::hits::DriftChamberParticleHit * ahit = event->fDCEvent.GetParticleHit(i);

   //      if( ahit->fMomentum.E() > 0.01 ) {

   double r  = 1.0;
   double dr = 0.001;
   double c  = 2.0;
   double dc = 0.001;
   double ti = 0.0;

   double dxyz = 0.1;
   double x0 = 0.0;
   double y0 = 0.0;
   double z0 = 0.0;
   TVector3 vertex = {x0,y0,z0};

   double x_samp = x0;//rand.Uniform( x0*(1.0-dxyz), x0*(1.0+dxyz) );
   double y_samp = y0;//rand.Uniform( y0*(1.0-dxyz), y0*(1.0+dxyz) );
   double z_samp = rand.Uniform( z0-dxyz, z0+dxyz);

   helix.SetHelixParameters(r,c);
   helix.SetVertex( TVector3(x_samp,y_samp,z_samp) );


   for(int i = 0; i< 10; i++) {

      double r_samp = rand.Uniform( r*(1.0-dr), r*(1.0+dr) );
      double c_samp = rand.Uniform( c*(1.0-dc), c*(1.0+dc) );
      double t_samp = rand.Uniform( 0, TMath::Pi() );


      helix.SetHelixParameters( r_samp, c_samp );

      test_points.push_back( helix(t_samp) );
      t_points.push_back( ti );

      std::cout << "r = "<< r_samp << std::endl;
      std::cout << "c = "<< c_samp << std::endl;
      std::cout << "t = "<< t_samp << std::endl;
   }

   HelixFcn theFCN;
   theFCN.SetMeasurements(test_points);

   //{
   //   // demonstrate minimal required interface for minimization
   //   // create Minuit parameters without names
   //   // starting values for parameters
   //   std::vector<double> init_par;
   //   init_par.push_back(1.0);
   //   init_par.push_back(1.55);
   //   // starting values for initial uncertainties
   //   std::vector<double> init_err;
   //   init_err.push_back(0.1);
   //   init_err.push_back(0.1);
   //   // create minimizer (default constructor)
   //   ROOT::Minuit2::VariableMetricMinimizer theMinimizer;
   //   // minimize
   //   ROOT::Minuit2::FunctionMinimum min =
   //      theMinimizer.Minimize(theFCN, init_par, init_err);
   //   // output
   //   std::cout<<"minimum: "<<min<<std::endl;
   //}

   {
      // demonstrate standard minimization using MIGRAD
      // create Minuit parameters with names
      ROOT::Minuit2::MnUserParameters upar;
      upar.Add("r", 1.0, 0.01);
      upar.Add("c", 1.0, 0.01);
      upar.Add("z0", 0.0, 0.1);
      //upar.Add("x0", 0.0, 0.01);
      //upar.Add("y0", 0.0, 0.01);
      upar.SetLimits(0,-2.0,2.0);
      upar.SetLimits(1,0.0,5.0);
      upar.SetLimits(2,-0.2,0.2);
      upar.Fix(2);
      //upar.SetLimits(2,-1.0,1.0);
      //upar.SetLimits(3,-1.0,1.0);
      // create MIGRAD minimizer
      ROOT::Minuit2::MnMigrad migrad(theFCN, upar);
      // minimize
      ROOT::Minuit2::FunctionMinimum min = migrad();
      // output
      std::cout<<"minimum: "<<min<<std::endl;
   }


   std::cout << "r = "<< r << std::endl;
   std::cout << "c = "<< c << std::endl;
   vertex.Print();

   return 0;
}
//______________________________________________________________________________

