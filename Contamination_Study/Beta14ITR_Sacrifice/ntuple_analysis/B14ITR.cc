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
  TH2F* h_B14ITR_dirty = new TH2F("h_B14ITR_dirty", "h_B14ITR_dirty", 50,0.,1.,50,-0.5,2.);
  TH2F* h_B14ITR_clean = new TH2F("h_B14ITR_clean", "h_B14ITR_clean", 50,0.,1.,50,-0.5,2.);
  TH1D* h_B14 = new TH1D("h_B14", "h_B14", 30, -0.5, 2.);

  //some cut selections you could apply
  //To build this mask, see snopl.us/docs/rat/user_manual/html/node226.html
  //int prescaleonly
  bool plotClean = true;
  bool plotDirty = true;
  int muonsonly = 0b10000000;
  int analysis_mask = 0b1111111111110;
  cout << analysis_mask << endl;
  //Way you could define cuts to put right in your Draw Statement
  TCut beta14_c1 = "beta14>0";

  for (int f=1; f<argc; f++)
  {
    const string& filename = string(argv[f]);
    TFile* mafile = TFile::Open(filename.c_str(),"READ");
    //Get the tree that has the entries
    TTree* T = (TTree*) mafile->Get("output");
    Double_t beta14;
    Double_t ITR;
    Double_t energy;
    Double_t radius;
    ULong64_t dcFlagged;   //will be filled per event; tells you what was flagged in that event
    ULong64_t dcApplied;   //should be the same for all events; tells you what was
                           //looked for when Data Cleaning was run
    Int_t nhits;

    T->SetBranchAddress("beta14",&beta14);
    T->SetBranchAddress("itr",&ITR);
    T->SetBranchAddress("dcApplied",&dcApplied);
    T->SetBranchAddress("energy",&energy);
    T->SetBranchAddress("posr",&radius);
    T->SetBranchAddress("dcFlagged",&dcFlagged);
    T->SetBranchAddress("nhits",&nhits);

    for (int entry=0; entry < T->GetEntries(); entry++)
    {
      T->GetEntry(entry);
      if(nhits < 40)
        continue;
      if(0.0 < beta14)   //Don't plot negative betas
        h_B14->Fill(beta14);
      if((-0.5 < beta14) & (0.0 < ITR))   //Manually applying cuts
      {
        if(((analysis_mask & dcFlagged) == analysis_mask))
          h_B14ITR_clean->Fill(ITR,beta14);
        else
          h_B14ITR_dirty->Fill(ITR,beta14);
      }
        //code to fill in beta14 histogram
      // end beta14 histogram
      

    }
  }

//  h_B14->Draw();
  h_B14ITR_dirty->SetMarkerStyle(20);
  h_B14ITR_dirty->SetMarkerColor(2);
  h_B14ITR_dirty->SetMarkerSize(0.7);
  h_B14ITR_dirty->GetXaxis()->SetTitle("ITR");
  h_B14ITR_dirty->GetYaxis()->SetTitle("Beta14");

  h_B14ITR_clean->SetMarkerStyle(20);
  h_B14ITR_clean->SetMarkerColor(4);
  h_B14ITR_clean->SetMarkerSize(0.7);
  h_B14ITR_clean->GetXaxis()->SetTitle("ITR");
  h_B14ITR_clean->GetYaxis()->SetTitle("Beta14");

  c1->cd(1);
  if(plotDirty)
    h_B14ITR_dirty->Draw();
  c1->cd(2);
  if(plotClean)
    h_B14ITR_clean->Draw("same");
  myapp->Run();

  delete h_B14;
  delete h_B14ITR_clean;
  delete h_B14ITR_dirty;
} //End main
