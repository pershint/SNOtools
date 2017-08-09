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

  //Constants that we need for our analysis
  double ECut = 5.5;    //in MeV
  double RCut = 5500;   //in mm
  double b14_low = -0.12;
  double b14_high = 0.95;
  double itr_low = 0.55;

  //Masks for preliminary cut and Data Cleaning Cuts in analysis
  int pathologicalMask = 0b100000011100000;
  int dcCutMask = 0b1111100011100;
  //FIXME: Still need to add in OwlEHi cut on the dcCuts segment
  //
  //Integers used in bifurcated analysis equation
  int a = 0;
  int b = 0;
  int c = 0;
  int d = 0;

  for (int f=1; f<argc; f++)
  {
    const string& filename = string(argv[f]);
    TFile* mafile = TFile::Open(filename.c_str(),"UPDATE");
    //Get the tree that has the entries
    TTree* T = (TTree*) mafile->Get("output");
    Double_t beta14;
    Double_t ITR;
    Double_t energy;
    Double_t radius;
    ULong64_t dcFlagged;   //will be filled per event; tells you what was flagged in that event
    ULong64_t dcApplied;   //should be the same for all events; tells you what was
    cout << "SETTIG BRANCH NAMES" << endl;
    //Set the branches that will fill the doubles we define per entry
    T->SetBranchAddress("beta14",&beta14);
    T->SetBranchAddress("itr",&ITR);
    T->SetBranchAddress("dcApplied",&dcApplied);
    T->SetBranchAddress("dcFlagged",&dcFlagged);
    T->SetBranchAddress("energy",&energy);
    T->SetBranchAddress("posr",&radius);

    for (int entry=0; entry < T->GetEntries(); entry++)
    {
      T->GetEntry(entry);
      bool DC_clean = 1;
      bool Class_clean = 1;
      cout << "RADIUS:" << radius << endl; 
      //First, ignore events that are pathological
      if((energy < ECut) | (radius > RCut))
        continue;
      cout << "PASSED ENERGY AND RADIUS CUT" << endl;
      if(dcFlagged & pathologicalMask) //skip entry if dcFlagged has any PMask
        continue;

      //Now, first see if the event passes the classifiers
      if((beta14 < b14_low) | (b14_high < beta14))
        Class_clean = 0;
      if(ITR < itr_low)
        Class_clean = 0;

      //Next, see if the event is clean accoding to the defined DC branch 
      if(dcFlagged & dcCutMask) //is dirty if dcFlagged has any dcCutMask bits
        DC_clean = 0;

      //Finally, increment the a,b,c,d values according to classification
      if (DC_clean && Class_clean)
        a++;
      else if (DC_clean && !(Class_clean))
        b++;
      else if (!(DC_clean) && Class_clean)
        c++;
      else if (!(DC_clean) && !(Class_clean))
        d++;
      else
        cout << "HOW ARE WE HERE?" << endl;
    } //End entry
  } //End ntuple file

  cout << "Printing variables to see if stuff works..." << endl;
  cout << "a: " << a << endl;
  cout << "b: " << b << endl;
  cout << "c: " << c << endl;
  cout << "d: " << d << endl;
} //End main
