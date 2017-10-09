////////////////////////////////////////////////////////////////////
/// \file sacrifice.cc
///
/// \brief Thesis plots of data cleaning sacrifice.
///
////////////////////////////////////////////////////////////////////

#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>

#include <string>
#include <map>
#include <iterator>
using namespace std;

void CSHist()
{
  const string filename = "output.root";
  TFile* myfile = TFile::Open(filename.c_str(),"READ");
  int nbins = 200;    //number of bins in histograms to be loaded
  TH1D* h_CS = new TH1D("h_CS","h_CS", nbins, 0., 200.);
  TH1D* h_AllEvents = myfile->Get("h_AllEvents");
  TH1D* h_FracFlagged = myfile->Get("h_FracFlagged");
//  int ndatpoints = h_AllEvents->GetEntries();
  for (int i = 10; i<nbins; i++)
  {
    double weightedfrac = 0.0;
    int ndatpoints = 0;
    //Sum up all of the bin entries following this one. Add to h_CS
    for (int j = i; j <nbins; j++)
    {
        ndatpoints = ndatpoints + h_AllEvents->GetBinContent(j);
        weightedfrac= weightedfrac + h_FracFlagged->GetBinContent(j) * h_AllEvents->GetBinContent(j);
        if(i==10){
          cout << "FRACFLAGGED: " << h_FracFlagged->GetBinContent(j) << endl;;
          cout << "ALL EVENTS IN BIN " << j << " : " << h_AllEvents->GetBinContent(j) << endl;;
          cout << "WEIGHTEDFRAC AT STEP 1: " << weightedfrac << endl;
        }
    }
    //Divide by the total number of events
    cout << "WEIGHT FRAC FROM BIN " << i << " : " << weightedfrac;
    double binval = weightedfrac / ndatpoints;
    cout << "WEIGHT FRAC DIVIDED BY NUM DATA POINTS: " << binval << endl;
    h_CS->Fill(i, binval);
  }
  TFile* thenewhist = new TFile("CSout.root", "CREATE");
  thenewhist->Add(h_CS);
  thenewhist->Write();
  thenewhist->Close();
}
