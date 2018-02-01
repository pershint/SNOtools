#include <TFile.h>
#include <TTree.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <typeinfo>
#include <string>
#include <map>
#include <iterator>

#include "ConfigParser.hh"

using namespace std;
using namespace configuration;


void write_to_out_file( const std::string &text )
{
    std::ofstream out_file(
        "BifurResults.out", std::ios_base::out | std::ios_base::app );
    out_file << text << std::endl;
}

int main(int argc, char** argv)
{
  //Define pathological cuts here
  int path_DCmask; 
  int path_trigmask;
  int cut1_DCmask;
  int cut1_trigmask;
  int cut2_DCmask;
  int cut2_trigmask;
  double E_low; 
  double E_high;
  double r_cut;
  double cut1_b14_low;
  double cut1_b14_high;
  double cut1_itr_low;
  double cut2_b14_low;
  double cut2_b14_high;
  double cut2_itr_low;
  bool dcs_only;
  bool fits_only;
  //some cut selections you could apply
  //To build this mask, see snopl.us/docs/rat/user_manual/html/node226.html
  //int prescaleonly
  configuration::CoParser configparse("../config/cuts_default.ini");
  try{
    dcs_only = configparse.getValueOfKey<bool>("dcs_only");
    fits_only = configparse.getValueOfKey<bool>("fits_only");
    if (dcs_only && fits_only){
      cout << "You can't use DCs only AND fits only. Choose both false or one true"
          << endl;
      exit(EXIT_FAILURE);
    }
    
    path_DCmask = configparse.getValueOfKey<int>("path_DCmask");
    path_trigmask = configparse.getValueOfKey<int>("path_trigmask");
    cut2_b14_low = configparse.getValueOfKey<double>("cut2_b14_low");
    cut2_b14_high = configparse.getValueOfKey<double>("cut2_b14_high");
    cut2_itr_low = configparse.getValueOfKey<double>("cut2_itr_low");
    cut1_DCmask = configparse.getValueOfKey<int>("cut1_DCmask");
    cut1_trigmask = configparse.getValueOfKey<int>("cut1_trigmask");
    if (dcs_only){
      cut2_DCmask = configparse.getValueOfKey<int>("cut2_DCmask");
      cut2_trigmask = configparse.getValueOfKey<int>("cut2_trigmask");
    }
    if (fits_only){
      cut1_b14_low = configparse.getValueOfKey<double>("cut1_b14_low");
      cut1_b14_high = configparse.getValueOfKey<double>("cut1_b14_high");
      cut1_itr_low = configparse.getValueOfKey<double>("cut1_itr_low");
    }
    E_low = configparse.getValueOfKey<double>("E_low");
    E_high = configparse.getValueOfKey<double>("E_high");
    r_cut = configparse.getValueOfKey<double>("r_cut");

  }  catch (int e) {
    std::cout << "ERROR READING FROM CONFIG FILE." << std::endl;
    return 1;
  };
  //Integers used in bifurcated analysis equation
  int a = 0;  //Passes cut branches 1 and 2
  int b = 0;  //Passes cut 1 but fails cut 2
  int c = 0;  //Passes cut 2 but fails cut 1
  int d = 0;  //Fails cut branches 1 and 2
  std::string outinit = "Files used: ";
  write_to_out_file(outinit);
  //For each ntuple, we'll check if events pass or fail each bifurcation branch cut
  for (int f=1; f<argc; f++)
  {
    const string& filename = string(argv[f]);
    std::string outwrite = string(argv[f]) + " ";
    write_to_out_file(outwrite);
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


    bool cut1_clean;
    bool cut2_clean;
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

      if(!dcs_only && !fits_only){ //We're using a DC branch and a fit branch
        //Next, see if the event is clean accoding to the defined DC branch 
        if(~(dcFlagged) & cut1_DCmask) //is dirty if dcFlagged has any cut_DCmask bits
          cut1_clean = 0;
        if((trigWord) & cut1_trigmask) //is dirty if it has the OwlEHi bit
          cut1_clean = 0;
        //Label the event dirty if it fails
        if((beta14 < cut2_b14_low) | (cut2_b14_high < beta14))
          cut2_clean = 0;
        if(ITR < cut2_itr_low)
          cut2_clean = 0;
      }
      else if(dcs_only){   //We're only using the DC cuts
        //Next, see if the event is clean accoding to the defined DC branch 
        if(~(dcFlagged) & cut1_DCmask) 
          cut1_clean = 0;
        if((trigWord) & cut1_trigmask) //is dirty if it has the OwlEHi bit
          cut1_clean = 0;
       if(~(dcFlagged) & cut2_DCmask) 
          cut2_clean = 0;
       if((trigWord) & cut2_trigmask) //is dirty if it has the OwlEHi bit
          cut1_clean = 0;
     }

     else if(fits_only)
     {
       if((beta14 < cut1_b14_low) | (cut1_b14_high < beta14))
         cut1_clean = 0;
       if(ITR < cut1_itr_low)
         cut1_clean = 0;

       if((beta14 < cut2_b14_low) | (cut2_b14_high < beta14))
         cut2_clean = 0;
       if(ITR < cut2_itr_low)
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
  configuration::converter convert;
  std::string space = "\n\n\n";
  write_to_out_file(space);
  std::string aresult = "a=" + convert.T_to_string(a);
  std::string bresult = "b=" + convert.T_to_string(b);
  std::string cresult = "c=" + convert.T_to_string(c);
  std::string dresult = "d=" + convert.T_to_string(d);
  write_to_out_file(aresult);
  write_to_out_file(bresult);
  write_to_out_file(cresult);
  write_to_out_file(dresult);
  cout << "a: " << a << endl;
  cout << "b: " << b << endl;
  cout << "c: " << c << endl;
  cout << "d: " << d << endl;
} //End main
