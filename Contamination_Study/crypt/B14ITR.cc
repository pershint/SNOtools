////////////////////////////////////////////////////////////////////
/// \file ReadDCFlags_NoPlots.cc
///
/// \Program to read out data quality flags from a processing ROOT file
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

  // Define our histograms
  TApplication* myapp = new TApplication("myapp",0,0);

  TCanvas* c1 = new TCanvas("c1","c1",800,1200);
  c1->Divide(2,1);
  //Define histograms to fill in
  TH2F* h_B14ITR_fail = new TH2F("h_B14ITR_fail", "h_B14ITR_fail", 50,0.,1.,50,-0.5,2.);
  TH2F* h_B14ITR_pass = new TH2F("h_B14ITR_pass", "h_B14ITR_pass", 50,0.,1.,50,-0.5,2.);
  TH1D* h_fit_fail_E = new TH1D("h_fit_fail_E", "h_fit_fail_E",  28,5.5,11.5);
  TH1D* h_fit_E = new TH1D("h_fit_E", "h_fit_E",  28,5.5,11.5);
  TH1D* h_FracFlagged = new TH1D("h_FracFlagged", "h_FracFlagged",  28,5.5,11.5);
  TH2F* h_B14ITR = new TH2F("h_B14ITR", "h_B14ITR", 50,0.,1.,50,-0.5,2.);

  h_fit_E->Sumw2();
  h_FracFlagged->Sumw2();
  h_fit_fail_E->Sumw2();
  //some cut selections you could apply
  //To build this mask, see snopl.us/docs/rat/user_manual/html/node226.html
  //int prescaleonly
  //FIXME: don't want the location hard-coded here
  configuration::CoParser configparse("../config/cuts_default.ini")
  try{
    bool plotClean = configparse.getValueOfKey<bool>("plotClean");
    bool plotDirty = configparse.getValueOfKey<bool>("plotDirty");

    //Define cuts here
    int path_DCmask = configparse.getValueOfKey<int>("path_DCmask");
    int path_trigmask = configparse.getValueOfKey<int>("path_trigmask");
    double E_low = configparse.getValueOfKey<double>("E_low");
    double E_high = configparse.getValueOfKey<double>("E_high");
    double r_cut = configparse.getValueOfKey<double>("r_cut");
    double b14_low = configparse.getValueOfKey<double>("b14_low");
    double b14_high = configparse.getValueOfKey<double>("b14_high");
    double itr_low = configparse.getValueOfKey<double>("itr_low");
  } catch {
    std::cout << "ERROR READING FROM CONFIG FILE." << std::endl;
    return 1;
  }

  TFile* mafile = TFile::Open(infile.c_str(),"READ");
  //Get the tree that has the entries
  TTree* T = (TTree*) mafile->Get("output");
  Double_t beta14;
  Double_t ITR;
  Double_t energy;
  Double_t posr;
  Bool_t fitValid;
  Int_t isN16;
  ULong64_t dcFlagged;   //will be filled per event; tells you what was flagged in that event
  ULong64_t dcApplied;   //should be the same for all events; tells you what was
                         //looked for when Data Cleaning was run
  Int_t triggerWord;

  T->SetBranchAddress("beta14",&beta14);
  T->SetBranchAddress("itr",&ITR);
  T->SetBranchAddress("dcApplied",&dcApplied);
  T->SetBranchAddress("energy",&energy);
  T->SetBranchAddress("posr",&posr);
  T->SetBranchAddress("isN16",&isN16);
  T->SetBranchAddress("fitValid",&fitValid);
  T->SetBranchAddress("dcFlagged",&dcFlagged);
  T->SetBranchAddress("triggerWord",&triggerWord);

  for (int entry=0; entry < T->GetEntries(); entry++)
  {
    T->GetEntry(entry);
    //First, we skip events defined as pathological
    if((energy > E_high) || (energy < E_low) || (posr > r_cut))
      continue;
    if(!fitValid)
      continue;
    if(!isN16)
      continue;
    //if (!(dcApplied))   //Uncommment to Check you data cleaned the file
    //  continue;

    if (~(dcFlagged) & path_DCmask)   //true if a pathological bit is in dcFlagged
      continue;
    if (triggerWord & path_trigmask)  //true if a pathological bit is in triggerWord
      continue;

    //passed pathological cuts. Check classifiers now.
    h_B14ITR->Fill(ITR,beta14);
    h_fit_E->Fill(energy);
    if ((beta14 > b14_high) | (beta14<b14_low) | (ITR < itr_low)){
        h_B14ITR_fail->Fill(ITR,beta14);
        h_fit_fail_E->Fill(energy);
    } else{
        h_B14ITR_pass->Fill(ITR,beta14);
    }
  } //loop entries end

  delete T;
  mafile->Close();
  delete mafile;
  //Now just a bunch of axis definitions; this could be organized better...
  h_FracFlagged->Divide(h_fit_fail_E,h_fit_E,1.,1.,"b");

  h_fit_fail_E->GetXaxis()->SetTitle("Energy(Mev)");
  h_fit_fail_E->GetYaxis()->SetTitle("Events");

  h_fit_E->GetXaxis()->SetTitle("Energy(Mev)");
  h_fit_E->GetYaxis()->SetTitle("Events");

  h_fit_E->GetXaxis()->SetTitle("Energy(Mev)");
  h_fit_E->GetYaxis()->SetTitle("Events");

  h_B14ITR->SetMarkerStyle(20);
  h_B14ITR->SetMarkerColor(7);
  h_B14ITR->SetMarkerSize(0.7);
  h_B14ITR->GetXaxis()->SetTitle("ITR");
  h_B14ITR->GetYaxis()->SetTitle("Beta14");

  h_B14ITR_fail->SetMarkerStyle(20);
  h_B14ITR_fail->SetMarkerColor(2);
  h_B14ITR_fail->SetMarkerSize(0.7);
  h_B14ITR_fail->GetXaxis()->SetTitle("ITR");
  h_B14ITR_fail->GetYaxis()->SetTitle("Beta14");

  h_B14ITR_pass->SetMarkerStyle(20);
  h_B14ITR_pass->SetMarkerColor(4);
  h_B14ITR_pass->SetMarkerSize(0.7);
  h_B14ITR_pass->GetXaxis()->SetTitle("ITR");
  h_B14ITR_pass->GetYaxis()->SetTitle("Beta14");

  //Put the histograms into a final output file
  TFile* thehists = new TFile(outname.c_str(), "CREATE");
  thehists->Add(h_B14ITR);
  thehists->Add(h_B14ITR_pass);
  thehists->Add(h_B14ITR_fail);
  thehists->Add(h_fit_E);
  thehists->Add(h_fit_fail_E);
  thehists->Add(h_FracFlagged);
  thehists->Write();
  thehists->Close();
} //End main
