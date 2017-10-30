#include <TH1D.h>
#include <TCut.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TApplication.h>

#include <stdlib.h>
#include <iostream>
#include <typeinfo>
#include <string>
#include <map>
#include <iterator>
using namespace std;

int main(int argc, char** argv)
{
  //Get filenames from command line
  int cut_DCmask = 0b1111100011100;  //DC branch of contamination study
  int DC_trigmask = 0b1000000000; //OwlEHi trigger bit
  bool B14only = false;
  bool ITRonly = false;
  bool commandline_DC = false;
  bool commandline_trig = false;
  int numentries = 1; //Have the filename itself
  for(int i=1; i<argc; i++)
  {
    if( argv[i] == string("-d") )
    {
      numentries++;
      commandline_DC = true;
      cerr << "found it" << endl;
      continue;
    }
    if( commandline_DC == true )
    {
      numentries++;
      cut_DCmask = atoi(argv[i]);
      commandline_DC = false;
      continue;
    }
    if( argv[i] == string("-b") )
    {
      numentries++;
      B14only = true;
      continue;
    }
    if( argv[i] == string("-i") )
    {
      numentries++;
      ITRonly = true;
      continue;
    }
  }

 // Define our histograms

  TApplication* myapp = new TApplication("myapp",0,0);

  TCanvas* c1 = new TCanvas("c1","c1",800,1200);
  c1->Divide(2,1);
  //Define histograms to fill in
  TH2F* h_B14ITR_dirty = new TH2F("h_B14ITR_dirty", "h_B14ITR_dirty", 50,0.,1.,50,-0.5,2.);


  //Define pathological cuts here
  int path_DCmask = 0b1110000011100010;  //Pathological cuts for contamination study
  int path_trigmask = 0b1010001100000;  //ESum, PGD, and PED triggers
  double E_low = 5.5;   //MeV
  double E_high = 9.0;  //MeV
  double r_cut = 5500;  //mm
  double b14_low = -0.12;
  double b14_high = 0.95;
  double itr_low = 0.55;


  //Integers used in bifurcated analysis equation
  int a = 0;
  int b = 0;
  int c = 0;
  int d = 0;

  for (int f=numentries; f<argc; f++)
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
    Int_t trigWord;
    Bool_t fitValid;
    //Set the branches that will fill the doubles we define per entry
    T->SetBranchAddress("beta14",&beta14);
    T->SetBranchAddress("itr",&ITR);
    T->SetBranchAddress("dcApplied",&dcApplied);
    T->SetBranchAddress("dcFlagged",&dcFlagged);
    T->SetBranchAddress("triggerWord",&trigWord);
    T->SetBranchAddress("fitValid",&fitValid);
    T->SetBranchAddress("energy",&energy);
    T->SetBranchAddress("posr",&radius);

    for (int entry=0; entry < T->GetEntries(); entry++)
    {
      T->GetEntry(entry);
      bool DC_clean = 1;
      bool Class_clean = 1;
      //First, we look in our background rich regions 
      if(!fitValid)
        continue;
      if((energy < E_low) | (energy > E_high) | (radius > r_cut))
        continue;
      if(~(dcFlagged) & path_DCmask) //skip entry if dcFlagged has any PMask
        continue;
      if(trigWord & path_trigmask) //skip entry if trigger has an Esum trigger
        continue;
      //Now, first see if the event passes the classifiers
      if(B14only){
        if((beta14 < b14_low) | (b14_high < beta14))
          Class_clean = 0;
      }
      else if(ITRonly){
        if(ITR < itr_low)
          Class_clean = 0;
      }
      else{
        if((beta14 < b14_low) | (b14_high < beta14))
          Class_clean = 0;
        if(ITR < itr_low)
          Class_clean = 0;
      }

      //Next, see if the event is clean accoding to the defined DC branch 
      if(~(dcFlagged) & cut_DCmask) //is dirty if dcFlagged has any cut_DCmask bits
        DC_clean = 0;
      if((trigWord) & DC_trigmask) //is dirty if it has the OwlEHi bit
        DC_clean = 0;

      //Finally, increment the a,b,c,d values according to classification
      if (DC_clean & Class_clean)
        a++;
      else if (DC_clean & !(Class_clean))
        b++;
      else if (!(DC_clean) & Class_clean)
        c++;
      else if (!(DC_clean) & !(Class_clean))
        d++;
      else
        cout << "HOW ARE WE HERE?" << endl;
    } //End entry
  delete T;
  mafile->Close();
  delete mafile;
  //Get the tree that has the entries
  } //End ntuple file

  cout << "Data cleaning mask used in DC branch: " << cut_DCmask << endl;
  cout << "Used ITR only?" << ITRonly << endl;
  cout << "Used B14 only?" << B14only << endl;
  cout << "a: " << a << endl;
  cout << "b: " << b << endl;
  cout << "c: " << c << endl;
  cout << "d: " << d << endl;
} //End main
