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
  //These are the cut masks used if you are comparing two cuts
  int cut_DC1 = 0;  //DC branch of contamination study
  int cut_DC2 = 0;  //DC branch of contamination study
  //Get filenames from command line
  int cut_DCmask = 0b1111100011100;  //default DC branch of contamination study
  int DC_trigmask = 0b1000000000; //OwlEHi trigger bit

  //Define pathological cuts here
  int path_DCmask = 0b1110000011100010;  //default Pathological cuts study
  int path_trigmask = 0b1010001100000;  //ESum, PGD, and PED triggers
  double E_low = 5.5;                   //MeV (default 5.5)
  double E_high = 9.0;                  //MeV (default 9.0)
  double r_cut = 5500;                  //mm  (default 5500)
  double b14_low = -0.12;               //default -0.12
  double b14_high = 0.95;               //default 0.95
  double itr_low = 0.55;                //default 0.55

  //Values turned to true depending on if the commandline is used to
  //Choose what cuts are in each branch
  bool DCsonly = false;
  bool B14only = false;
  bool ITRonly = false;
  bool commandline_DC = false;
  bool commandline_twoDCs = false;
  bool commandline_trig = false;
  bool Fitsonly = false;
  int numentries = 1; //Have the filename itself

  //Command line parsing
  for(int i=1; i<argc; i++)
  {
    if( argv[i] == string("-c") )
    {
      numentries++;
      Fitsonly = true;
      continue;
    }
   if( argv[i] == string("-d") )
    {
      numentries++;
      commandline_DC = true;
      continue;
    }
    if( commandline_DC == true )
    {
      numentries++;
      cut_DCmask = atoi(argv[i]);
      commandline_DC = false;
      continue;
    }
    //use this flag to feed in two data cuts for bifur analysis
    if( argv[i] == string("-dd") )
    {
      numentries++;
      commandline_twoDCs = true;
      continue;
    }
    if( commandline_twoDCs == true )
    {
      numentries=numentries + 2;
      cut_DC1 = atoi(argv[i]);
      cut_DC2 = atoi(argv[i+1]);
      DCsonly = true;
      commandline_twoDCs = false;
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

 // Define our histograms; code crashes without this.  Thanks ROOT

  TApplication* myapp = new TApplication("myapp",0,0);
  TCanvas* c1 = new TCanvas("c1","c1",800,1200);
  c1->Divide(2,1);
  //Define histograms to fill in
  TH2F* h_B14ITR_dirty = new TH2F("h_B14ITR_dirty", "h_B14ITR_dirty", 50,0.,1.,50,-0.5,2.);


  //Integers used in bifurcated analysis equation
  int a = 0;
  int b = 0;
  int c = 0;
  int d = 0;

  //For each ntuple, we'll check if events pass or fail each bifurcation branch cut
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
      bool cut1_clean = 1;  //originally DC_clean
      bool cut2_clean = 1;  //originally Class_clean
      //First, we look in our background rich regions 
      if(!fitValid)
        continue;
      if((energy < E_low) | (energy > E_high) | (radius > r_cut))
        continue;
      if(~(dcFlagged) & path_DCmask) //skip entry if dcFlagged has any PMask
        continue;
      if(trigWord & path_trigmask) //skip entry if trigger has an Esum trigger
        continue;

      if(!DCsonly && !Fitsonly){ //We're using a DC branch and a fit branch
        //Label the event dirty if it fails
        if(B14only){
          if((beta14 < b14_low) | (b14_high < beta14))
            cut2_clean = 0;
        }
        else if(ITRonly){
          if(ITR < itr_low)
            cut2_clean = 0;
        }
        else {
          if((beta14 < b14_low) | (b14_high < beta14))
            cut2_clean = 0;
          if(ITR < itr_low)
            cut2_clean = 0;
        }

        //Next, see if the event is clean accoding to the defined DC branch 
        if(~(dcFlagged) & cut_DCmask) //is dirty if dcFlagged has any cut_DCmask bits
          cut1_clean = 0;
        //FIXME: Need a clean way to choose whether or not to use OwlEHi also
//        if((trigWord) & DC_trigmask) //is dirty if it has the OwlEHi bit
//          cut1_clean = 0;
      }

      else if(DCsonly){   //We're only using the DC cuts
        //Next, see if the event is clean accoding to the defined DC branch 
        if(~(dcFlagged) & cut_DC1) 
          cut1_clean = 0;
        if(~(dcFlagged) & cut_DC2) 
          cut2_clean = 0;
      }

     else if(Fitsonly)
     {
       if((beta14 < b14_low) | (b14_high < beta14))
         cut1_clean = 0;
       if(ITR < itr_low)
         cut2_clean = 0;
     }
    //Finally, increment the a,b,c,d values according to classification
      if (cut1_clean & cut2_clean)
        a++;
      else if (cut1_clean & !(cut2_clean))
        b++;
      else if (!(cut1_clean) & cut2_clean)
        c++;
      else if (!(cut1_clean) & !(cut2_clean))
        d++;
      else
        cout << "HOW ARE WE HERE?" << endl;
    } //End entry loop
    delete T;
    mafile->Close();
    delete mafile;
  } //End ntuple file loop

  cout << "Used DC branch and Fit branch?" << (!DCsonly & !Fitsonly) << endl;
  if(!DCsonly)
    cout << "Data cleaning mask used in DC branch: " << cut_DCmask << endl;
  else
    cout << "Data cleaning masks used: " << cut_DC1 << ":" << cut_DC2 << endl;
  cout << "Used ITR only?" << ITRonly << endl;
  cout << "Used B14 only?" << B14only << endl;
  cout << "a: " << a << endl;
  cout << "b: " << b << endl;
  cout << "c: " << c << endl;
  cout << "d: " << d << endl;
} //End main
