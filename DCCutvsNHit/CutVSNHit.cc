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
  // Define our histograms

  TApplication* myapp = new TApplication("myapp",0,0);

  TCanvas* c1 = new TCanvas("c1","c1",800,1200);
  c1->Divide(2,1);
  //Define histograms to fill in
  TH1D* h_FlaggedEvents = new TH1D("h_FlaggedEvents", "h_FlaggedEvents", 200,0.,200.);
  TH1D* h_AllEvents = new TH1D("h_AllEvents", "h_AllEvents", 200,0.,200.);
  TH1D* h_FracFlagged = new TH1D("h_FracFlagged", "h_FracFlagged", 200,0.,200.);
 //some cut selections you could apply
  //To build this mask, see snopl.us/docs/rat/user_manual/html/node226.html
  //int prescaleonly
  int cut_mask = 0b100000;

  for (int f=1; f<argc; f++)
  {
    const string& filename = string(argv[f]);
    TFile* mafile = TFile::Open(filename.c_str(),"UPDATE");
    //Get the tree that has the entries
    TTree* T = (TTree*) mafile->Get("output");
    ULong64_t dcFlagged;   //will be filled per event; tells you what was flagged in that event
    ULong64_t dcApplied;   //should be the same for all events; tells you what was
    //Bool_t fitValid;     //Uncomment if you want only the valid fits
    Int_t nhits;

    T->SetBranchAddress("dcApplied",&dcApplied);
    T->SetBranchAddress("dcFlagged",&dcFlagged);
    T->SetBranchAddress("nhits",&nhits);
    //T->SetBranchAddress("fitValid",&fitValid);

    for (int entry=0; entry < T->GetEntries(); entry++)
    {
      T->GetEntry(entry);
      if(nhits > 200)
        continue;
      h_AllEvents->Fill(nhits,1.);
      if (~(dcFlagged) & cut_mask)   //true if a cut_mask bit is in dcFlagged
        h_FlaggedEvents->Fill(nhits,1.);
    }
  }

  h_FracFlagged->Divide(h_FlaggedEvents,h_AllEvents);

  h_FlaggedEvents->GetXaxis()->SetTitle("nhit");
  h_FlaggedEvents->GetYaxis()->SetTitle("Events");

  h_FracFlagged->GetXaxis()->SetTitle("nhit");
  h_FracFlagged->GetYaxis()->SetTitle("Events");

  c1->cd(1);
    h_FlaggedEvents->Draw();
  c1->cd(2);
    h_FracFlagged->Draw();
  myapp->Run();

  delete h_FlaggedEvents;
  delete h_AllEvents;
  delete h_FracFlagged;
} //End main
