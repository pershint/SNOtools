// C++ Headers
#include <iostream>
#include <vector>
#include <map>
// ROOT Headers
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TFile.h>
#include <TString.h>
#include <Rtypes.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TMath.h>
// RAT Headers
#include <RAT/SunUtil.hh>
#include <RAT/DU/ReactorNuOsc.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/ReconCorrector.hh> 
// JSON Reader
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
using namespace std;

double shitZshift = -108;

// Order of operations:
// 1. Prune {fitValid, posr<6m}
// 2. Add to output tree: SunR & uDotR -- sun,dxdydz
// 3. Add header tree with: Accepted, rejected
class ReactorBro : public RAT::DU::ReactorNuOsc {
  public:
  ReactorBro()
  {
    RAT::DU::Utility::Get()->LoadDBAndBeginRun();
    this->BeginOfRun();
  }
  double OscWeight( const double nuKE, std::string coreData )
  {
    this->ReactorOsc(nuKE, coreData);
    return 1 - this->fOscProb;
  }
};

class Parser
{
  public:
    string oname = "default.root";
    vector<string> flist;
    bool isData = false;

    Parser(int argc, char** argv)
    {
      vector<string> args(argv, argv+argc);
      for(int n=1; n<argc; n++)
      {
        if( args[n] == "--out" )
        {
          this->oname = args[n+1];
          n++;
        }
        else if( args[n] == "--data" )
          this->isData = true;
        else
          flist.push_back(args[n]);
      }
      if (flist.size() == 0)
      {
        cout << "Need at least one file" << endl;
        exit(-1);
      }
    }
};

class BadChannels
{
  public:
    string iname = "badchannels.json";
    std::map<int, int> NBCMap;

    BadChannels()
    {
      // Read the input json file
      namespace pt = boost::property_tree;
      pt::ptree iroot;
      // Load the json file
      pt::read_json(iname, iroot);
      for( const auto& parent : iroot )
      {
        NBCMap[ stoi(parent.first) ] = stoi(parent.second.data());
      }
    }
};

class SmearEnergy
{
  private:
    double minE;
    double maxE;
    double smearDown;
    double smearUp;
    double bias;
    int bins;
    map<int, TF1*> smearFuncs;
    double delE;
  public:
    SmearEnergy(double minE, double maxE, int bins, 
        double smearDown, double smearUp, double bias)
    {
      this->minE = minE;
      this->maxE = maxE;
      this->bins = bins;
      this->smearDown = smearDown;
      this->smearUp = smearUp;
      this->bias = bias;
      this->delE = (maxE-minE)/bins;
      for(int i=0; i<bins; i++)
      {
        TF1* sFunc = new TF1("fitfun", fitfun, minE, maxE, 4);
        double mu = minE + i*delE + bias;
        double siglow = smearDown*mu;
        double sighi = smearUp*mu;
        sFunc->SetParameters(mu, siglow, sighi, bias);
        smearFuncs[i] = sFunc;
      }
    }
    double random(double energy)
    {
      if( energy > maxE || energy < minE )
        return -9999;
      int tbin = int((energy-minE)/delE);
      return (smearFuncs[tbin])->GetRandom();
    }
    static double fitfun(double *x, double *par)
    {
      double mu = par[0];
      double s1 = par[1];
      double s2 = par[2];
      double b = par[3];
      double val = x[0];
      double sigma = 0;
      if( val<mu )
        sigma = s1;
      else
        sigma = s2;
      return sqrt(2/TMath::Pi())/(s1+s2)*TMath::Exp(-pow((val-mu)/sigma,2)/2);
    }
};

class SmearVertex
{
  // Smear X, Y, and Z by the same? Can change later
  private:
    TF1* sFunc;
  public:
    SmearVertex(double smearSigma)
    {
      sFunc = new TF1("vFunc", "gaus(0)", 0, 10000);
      double mu=0;
      double sigma = smearSigma;
      sFunc->SetParameters(1, mu, sigma);
    }
    ~SmearVertex(){delete sFunc;}
    double random(double pos)
    {
      return sFunc->GetRandom() + pos;
    }
    double randomR(double posx, double posy, double posz)
    {
      double nx = this->random(posx);
      double ny = this->random(posy);
      double nz = this->random(posz);
      return sqrt(nx*nx + ny*ny + nz*nz);
    }
};

class SmearSun
{
  // Similar to the method used in the LETA unidoc
  private:
    double smearterm;
    TRandom3* rndm;
  public:
    SmearSun(double smearterm)
    {
      // smear term should be -1(->)1, where 1=perfect
      // and negative is for negative error bar
      // should be 1 - b1/b2 where b12 are the exponential fits
      this->smearterm = smearterm;
      rndm = new TRandom3();
    }
    double NewSunct(double sunct)
    {
      double newct = -1 + (sunct+1)*(1+smearterm);
      if(newct<-1 || newct>1)
        newct = 2*rndm->Rndm()-1;
      return newct;
    }
};

int main(int argc, char** argv)
{
  Parser args(argc, argv); 
  BadChannels* bc = new BadChannels();

  ReactorBro* rb = new ReactorBro(); 
  // Check the files in argv
  TChain* chain = new TChain("output");
  cerr << "Saving to file " << args.oname << endl;
  cerr<<"loading tchains "<< args.flist.size() - 1 <<endl;
  for(auto& it : args.flist)
    chain->Add(it.c_str());
  //for(int i=1; i<argc; i++)
  //{
  //  chain->Add(argv[i]);
  //}
  // Correction to energy
  const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();

  TFile* outfile = new TFile(args.oname.c_str(), "recreate");
  // Header tree
  TTree* header = new TTree("header", "MC Header Information");
  int accepted=0;
  double efficiency;
  int total;
  header->Branch("efficiency", &efficiency);
  header->Branch("total", &total);
  TH1F* energy_hist = new TH1F("energy_hist", "energy_hist", 50, 0, 10);
  TH1F* radius_hist = new TH1F("radius_hist", "radius_hist", 600, 0, 6000);
  TH1F* sunct_hist = new TH1F("sunct_hist", "sunct_hist", 100, -1, 1);
  TH1F* udotr_hist = new TH1F("udotr_hist", "udtor_hist", 100, -1, 1);
  TH1F* itr_hist = new TH1F("itr_hist", "itr_hist", 100, 0, 1);
  TH1F* beta14_hist = new TH1F("beta14_hist", "beta14_hist", 100, -5, 5);

  TTree* outtree = chain->CloneTree(0);
  double rmax = 6000;
  double nbcmax = 1100;

  int evIndex;
  double energy, posr;
  double posx, posy, posz;
  double dirx, diry, dirz;
  int uTDays, uTSecs, uTNSecs;
  double itr, beta14;
  ULong64_t dcFlagged;
  ULong64_t dcClean = 0x7ffe;
  TString* parentMeta = new TString(); // Do i need parentMeta1 and 2?
  double parentKE;
  int runID;

  int parentpdg1;
  bool fitValid;
  //TString* parentMeta = new TString(256); // Do i need parentMeta1 and 2?
  chain->SetBranchAddress("parentMeta1", &parentMeta);
  chain->SetBranchAddress("parentKE1", &parentKE);
  chain->SetBranchAddress("evIndex", &evIndex);
  chain->SetBranchAddress("runID", &runID);
  chain->SetBranchAddress("parentpdg1", &parentpdg1);
  chain->SetBranchAddress("energy", &energy);
  chain->SetBranchAddress("posr", &posr);
  chain->SetBranchAddress("fitValid", &fitValid);
  chain->SetBranchAddress("posx", &posx);
  chain->SetBranchAddress("posy", &posy);
  chain->SetBranchAddress("posz", &posz);
  chain->SetBranchAddress("dirx", &dirx);
  chain->SetBranchAddress("diry", &diry);
  chain->SetBranchAddress("dirz", &dirz);
  chain->SetBranchAddress("uTDays", &uTDays);
  chain->SetBranchAddress("uTSecs", &uTSecs);
  chain->SetBranchAddress("uTNSecs", &uTNSecs);
  chain->SetBranchAddress("itr", &itr);
  chain->SetBranchAddress("beta14", &beta14);
  chain->SetBranchAddress("dcFlagged", &dcFlagged);
  // New Branches
  double udotr, sunct, sunx, suny, sunz, posr3;
  int nbc;
  double oscWeight;
  outtree->Branch("posr3", &posr3);
  outtree->Branch("udotr", &udotr);
  outtree->Branch("sunct", &sunct);
  outtree->Branch("sunx", &sunx);
  outtree->Branch("suny", &suny);
  outtree->Branch("sunz", &sunz);
  outtree->Branch("oscillate", &oscWeight);
  outtree->Branch("nbc", &nbc);
  unsigned int entries = chain->GetEntries();
  total = 0;
  // Load in solar spectrum
  TFile* solarFile = new TFile("Solar/out_pee_sun.root");
  TGraph* nueSurvival = (TGraph*)solarFile->Get("pee_bs05op_b8");

  cerr<<"Begin processing: "<<entries<<" events"<<endl;
  for(int i=0; i<entries; i++)
  {
    if(i%10000==0)
      cerr<<"Proc: "<<i/double(entries)*100<<"%\r"<<flush;
    chain->GetEvent(i);
    // Oscillate reactors
    if( *parentMeta == "" )
    {
      // Check if solar neutrino, assume all nu_e are solar btw
      if( parentpdg1 == 12 ) // nu_e
        oscWeight = nueSurvival->Eval( parentKE );
      else if( parentpdg1 == 14 ) // nu_mu
        oscWeight = 1 - nueSurvival->Eval( parentKE );
      else
        oscWeight = 1;
    }
    else
    {
      oscWeight = rb->OscWeight(parentKE, static_cast<std::string>(*parentMeta));
    }
    // Trigger cut
    if( evIndex>0 )
      continue;
    total++;
    if( !fitValid )
      continue;
    if( args.isData )
      if( (dcFlagged & dcClean) != dcClean )
        continue;
    if( posr > rmax )
      continue;
    nbc = bc->NBCMap[runID];
    if( nbc >= nbcmax )
      continue;
    energy = eCorr.CorrectEnergyRSP(energy);
    accepted++;
    posz = posz + shitZshift;
    posr = sqrt( posx*posx + posy*posy + posz*posz );
    udotr = (posx*dirx+posy*diry+posz*dirz)/sqrt(posx*posx+posy*posy+posz*posz)/sqrt(dirx*dirx+diry*diry+dirz*dirz);
    TVector3 Solar = RAT::SunDirection(uTDays, uTSecs, uTNSecs);
    sunct = (Solar.X()*dirx + Solar.Y()*diry + Solar.Z()*dirz);
    sunx = Solar.X();
    suny = Solar.Y();
    sunz = Solar.Z();
    posr3 = pow((posr/6000),3);

    outtree->Fill();

    energy_hist->Fill(energy);
    radius_hist->Fill(posr);
    sunct_hist->Fill(sunct);
    udotr_hist->Fill(udotr);
    itr_hist->Fill(itr);
    beta14_hist->Fill(beta14);
  }
  efficiency = double(accepted)/total;
  header->Fill();
  cerr<<endl<<"fin"<<endl;
  outfile->Write();
}
