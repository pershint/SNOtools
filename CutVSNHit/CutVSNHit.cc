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
using namespace std;

int main(int argc, char** argv)
{
  //Define histograms to fill in
  TH1D* h_FlaggedEvents = new TH1D("h_FlaggedEvents", "h_FlaggedEvents", 200,0.,200.);
  TH1D* h_AllEvents = new TH1D("h_AllEvents", "h_AllEvents", 200,0.,200.);
  TH1D* h_FracFlagged = new TH1D("h_FracFlagged", "h_FracFlagged", 200,0.,200.);
  h_AllEvents->Sumw2();
  h_FracFlagged->Sumw2();
  h_FlaggedEvents->Sumw2();
  const char* outfile = argv[1];
  //To build this mask, see snopl.us/docs/rat/user_manual/html/node226.html
  //Or use /rat-tools/DataCleaningTools/dcflags.py
  int cut_mask = 0b1000000000000; //Flag any events with the 1 bit as "dirty"

  for (int f=2; f<argc; f++)
  {
    const string& filename = string(argv[f]);
    TFile* mafile = TFile::Open(filename.c_str(),"READ");
    //Get the tree that has the entries
    TTree* T = (TTree*) mafile->Get("output");
    ULong64_t dcFlagged;   //will be filled per event; tells you what was flagged in that event
    ULong64_t dcApplied;   //should be the same for all events; tells you what was
    Bool_t fitValid;     //Uncomment if you want only the valid fits
    Bool_t isCal;
    Int_t nhits;
    Int_t GTID;

    T->SetBranchAddress("dcApplied",&dcApplied);
    T->SetBranchAddress("dcFlagged",&dcFlagged);
    T->SetBranchAddress("nhits",&nhits);
    T->SetBranchAddress("fitValid",&fitValid);
    T->SetBranchAddress("isCal",&isCal);
    T->SetBranchAddress("eventID",&GTID);
    //T->SetBranchAddress("fitValid",&fitValid);

    for (int entry=0; entry < T->GetEntries(); entry++)
    {
      T->GetEntry(entry);
      if(nhits > 200)
        continue;
      if(!(fitValid && isCal))
        continue;
      h_AllEvents->Fill(nhits);
      if (!(dcApplied))  // Checks there was DC applied to this event
        continue;
      if (dcApplied & dcFlagged != dcFlagged){
        //Event is flagged as "clean" even though it's not in the applied mask
        //Check for bugs in DC
        continue;
      }
      if (~(dcFlagged) & cut_mask){   //true if an event was flagged cut_mask bit is in dcFlagged
        h_FlaggedEvents->Fill(nhits);
      }
    }
  delete T;
  mafile->Close();
  delete mafile;
  }

  h_FracFlagged->Divide(h_FlaggedEvents,h_AllEvents,1.,1.,"b");

  h_FlaggedEvents->GetXaxis()->SetTitle("nhit");
  h_FlaggedEvents->GetYaxis()->SetTitle("Events");

  h_FracFlagged->GetXaxis()->SetTitle("nhit");
  h_FracFlagged->GetYaxis()->SetTitle("Events");

  TFile* thehists = new TFile(outfile, "CREATE");
  thehists->Add(h_FlaggedEvents);
  thehists->Add(h_AllEvents);
  thehists->Add(h_FracFlagged);
  thehists->Write();
  thehists->Close();
} //End main
