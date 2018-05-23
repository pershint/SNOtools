////////////////////////////////////////////////////////////////////
/// \file Sacrifice.cc
///
/// \This file takes in an ntuple and calculates the number of events
/// \Flagged as bad by either data cleaning or Beta14/ITR.  Before
/// \Checking DC/B14ITR, the event will be skipped if considered a
/// \Pathological event (i.e. one not included in the Bifurcation
/// \Analysis study.
//\ Output is a set of histograms showing how many events were flagged
/// \as a function of energy.
////////////////////////////////////////////////////////////////////
#include <TH1D.h>
#include <TCut.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TApplication.h>

#include <iostream>
#include <typeinfo>
#include <string>
#include <map>
#include <iterator>

#include "ConfigParser.hh"

using namespace std;
using namespace configuration;

int main(int argc, char** argv)
{
  //Get filenames from command line
  std::string outname = "none.txt";
  std::string infile = "none.txt";
  bool haveinput = false;
  bool haveoutput = false;
  for(int i=1; i<argc; i++)
  {
    if( argv[i] == string("-i") )
    {
      haveinput = true;
      cerr << "found it" << endl;
      continue;
    }
    if( haveinput == true )
    {
      infile = argv[i];
      haveinput = false;
      continue;
    }
    if( argv[i] == string("-o") )
    {
      haveoutput = true;
      continue;
    }
    if( haveoutput == true )
    {
      outname = argv[i];
      haveoutput = false;
      continue;
    }
  }

  int path_DC_DCmask; 
  int path_DC_trigmask;
  int cut_DCmask;
  int cut_trigmask;
  double r_cut;
  double b14_low;
  double b14_high;
  double itr_low;
  //some cut selections you could apply
  //To build this mask, see snopl.us/docs/rat/user_manual/html/node226.html
  //int prescaleonly
  try{
    configuration::CoParser configparse("../config/cuts_default.ini");
    path_DC_DCmask = configparse.getValueOfKey<int>("path_DC_DCmask");
    path_DC_trigmask = configparse.getValueOfKey<int>("path_DC_trigmask");
    cut_DCmask = configparse.getValueOfKey<int>("cut1_DCmask");
    cut_trigmask = configparse.getValueOfKey<int>("cut1_trigmask");
    r_cut = configparse.getValueOfKey<double>("r_cut");
    b14_low = configparse.getValueOfKey<double>("cut2_b14_low");
    b14_high = configparse.getValueOfKey<double>("cut2_b14_high");
    itr_low = configparse.getValueOfKey<double>("cut2_itr_low");

  }  catch (std::exception& e) {
    std::cout << "ERROR READING FROM CONFIG FILE." << std::endl;
    return EXIT_FAILURE;
  };
  //I don't even use these; it crashed on my cluster without them
  // Define our histograms
  //Define histograms to fill in
  TApplication* myapp = new TApplication("myapp",0,0);

  TCanvas* c1 = new TCanvas("c1","c1",800,1200);
  TH1D* h_DC_FlaggedEvents = new TH1D("h_DC_FlaggedEvents", "h_DC_FlaggedEvents", 58,0.0,11.6);
  TH1D* h_AllEvents = new TH1D("h_AllEvents", "h_AllEvents", 58,0.0,11.6);
  TH1D* h_DC_FracFlagged = new TH1D("h_DC_FracFlagged", "h_DC_FracFlagged", 58,0.0,11.6);
  TH1D* h_BI_FlaggedEvents = new TH1D("h_BI_FlaggedEvents", "h_BI_FlaggedEvents", 58,0.0,11.6);
  TH1D* h_BI_FracFlagged = new TH1D("h_BI_FracFlagged", "h_BI_FracFlagged", 58,0.0,11.6);
  h_AllEvents->Sumw2();
  h_DC_FracFlagged->Sumw2();
  h_DC_FlaggedEvents->Sumw2();
  h_BI_FracFlagged->Sumw2();
  h_BI_FlaggedEvents->Sumw2();
  //const char* outfile = outname;


  //open the file now
  //const string& filename = infile;
  TFile* mafile = TFile::Open(infile.c_str(),"READ");
  //Get the tree that has the entries
  TTree* T = (TTree*) mafile->Get("output");
  ULong64_t dcFlagged;   //will be filled per event; tells you what was flagged in that event
  ULong64_t dcApplied;   //should be the same for all events; tells you what was
  Bool_t fitValid;     //Uncomment if you want only the valid fits
  Int_t isN16;
  Int_t triggerWord;
  Double_t posr;
  Double_t energy;
  Double_t beta14;
  Double_t ITR;

  T->SetBranchAddress("beta14",&beta14);
  T->SetBranchAddress("itr",&ITR);
  T->SetBranchAddress("dcApplied",&dcApplied);
  T->SetBranchAddress("dcFlagged",&dcFlagged);
  T->SetBranchAddress("triggerWord",&triggerWord);
  T->SetBranchAddress("posr",&posr);
  T->SetBranchAddress("energy",&energy);
  T->SetBranchAddress("fitValid",&fitValid);
  T->SetBranchAddress("isN16",&isN16);

  //loop through entries to find clean/dirty events
  for (int entry=0; entry < T->GetEntries(); entry++)
  {
    T->GetEntry(entry);
    //First, we skip events defined as pathological
    if(posr > r_cut)
      continue;
    if(!fitValid)
      continue;
    if(!isN16)
      continue;
    //if (!(dcApplied))   //Check you actually have data cleaning run on file
    //  continue;
    if (~(dcFlagged) & path_DC_DCmask)   //true if a pathological bit is in dcFlagged
      continue;
    if (triggerWord & path_DC_trigmask)  //true if a pathological bit is in triggerWord
      continue;

    //Now, fill our histogram with all events passing pathological cuts
    h_AllEvents->Fill(energy);
    if ((~(dcFlagged) & cut_DCmask) || (triggerWord & cut_trigmask)){
      //true if an event is dirty according to DC branch or has OwlEHi
      h_DC_FlaggedEvents->Fill(energy);
    }
    //See if event is flagged by B14ITR parameters
    if ((beta14 > b14_high) | (beta14<b14_low) | (ITR < itr_low)){
        h_BI_FlaggedEvents->Fill(energy);
    }
  } //loop entries end


 
  delete T;
  mafile->Close();
  delete mafile;
  

  h_DC_FracFlagged->Divide(h_DC_FlaggedEvents,h_AllEvents,1.,1.,"b");
  h_BI_FracFlagged->Divide(h_BI_FlaggedEvents,h_AllEvents,1.,1.,"b");

  h_DC_FlaggedEvents->GetXaxis()->SetTitle("Energy(Mev)");
  h_DC_FlaggedEvents->GetYaxis()->SetTitle("Events");

  h_DC_FracFlagged->GetXaxis()->SetTitle("Energy(MeV)");
  h_DC_FracFlagged->GetYaxis()->SetTitle("Fractional Sacrifice");

  h_BI_FlaggedEvents->GetXaxis()->SetTitle("Energy(Mev)");
  h_BI_FlaggedEvents->GetYaxis()->SetTitle("Events");

  h_BI_FracFlagged->GetXaxis()->SetTitle("Energy(MeV)");
  h_BI_FracFlagged->GetYaxis()->SetTitle("Fractional Sacrifice");

  TFile* thehists = new TFile(outname.c_str(), "CREATE");
  thehists->Add(h_BI_FlaggedEvents);
  thehists->Add(h_BI_FracFlagged);
  thehists->Add(h_DC_FlaggedEvents);
  thehists->Add(h_AllEvents);
  thehists->Add(h_DC_FracFlagged);
  thehists->Write();
  thehists->Close();
} //End main
