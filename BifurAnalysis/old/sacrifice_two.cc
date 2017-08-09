////////////////////////////////////////////////////////////////////////////////
// Sacrifice studies for data cleaning
////////////////////////////////////////////////////////////////////////////////

#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/DataCleaningBits.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/DataQCFlags.hh>
#include <RAT/DS/BitMask.hh>
#include <RAT/DS/Meta.hh>

#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TChain.h>
#include <TFile.h>

#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <map>

using namespace std;

#define dbug

/// Overview
// Loop over all data cleaning words
// For each event loop over params: nhits, energy, radius, direction, isotropy
// Each param_word combination will have 2 histograms
//     1. [cut] Fill with events that fail a particular word
//     2. [denom] Denominator histogram: All events
// An extra histogram will be filled which is the OR of the cut histogram
//     1. [total] Total histogram: Fill if ANY word exists
// Plots:
//     -

vector<string> dc_vector()
{
  // List of all of the cuts to perform sacrifice studies on
  
  const char *vinit[] = {"ftscut", "flashergeocut", "itctimespreadcut", 
                         "junkcut", "muontag", "neckcut",
                         "owlcut", "qcluster", "qvnshit", "qvt", "ringoffire",
                         "crateisotropy", "zerozerocut", //"prescalecut", 
                         "tpmuonfollowercut", "caencut", "waterBlind" };
  
  //const char *vinit[] = {"ftscut"};
  int length = sizeof(vinit)/sizeof(vinit[0]);
  vector<string> vee(vinit, vinit+length);
  return vee;
}

struct HistParams
{
  int params[3];
};

map<string, HistParams> hist_param_map()
{
  // Map of histogram parameters for those found in param_vector
  // using compound literals of the HistParams struct
  // Parameters are: nhits, position, direction, energy, radius, and isotropy
  map<string, HistParams> hpmap;
  hpmap["nhits"]      = (HistParams){300,  0,  300};
  hpmap["position"]   = (HistParams){100,   0,  100};
  hpmap["direction"]  = (HistParams){100,   -2,  2};
  hpmap["energy"]     = (HistParams){400,   0,  40};
  hpmap["radius"]     = (HistParams){800,   0,  pow(8,3)};
  hpmap["isotropy"]   = (HistParams){100,   0,  100};
  return hpmap;
}

vector<string> param_vector()
{
  // List of all of the cuts to perform sacrifice studies on
  /*
  const char *vinit[] = {"nhits", "position", "direction", 
                         "energy", "radius", "isotropy" };
  */
  const char *vinit[] = {"nhits", "energy", "radius", "direction"};
  int length = sizeof(vinit)/sizeof(vinit[0]);
  vector<string> vee(vinit, vinit+length);
  return vee;
}

map< string, map<string, TH1D*> > cut_mapper(vector<string> dcvector, 
                                           vector<string> pvector)
{
  map<string, HistParams> hpmap = hist_param_map();
  map<string,map<string, TH1D*> > mymap;
  for(vector<string>::iterator dcit=dcvector.begin(); 
      dcit!=dcvector.end(); ++dcit)
  {
    map<string, TH1D*> insidemap;
    for(vector<string>::iterator pit=pvector.begin();
        pit!=pvector.end(); ++pit)
    {
      char hname[50];
      sprintf(hname, "cut_%s_%s",(*dcit).c_str(), (*pit).c_str());
      HistParams hp = hpmap[(*pit)];
      TH1D* nhist = new TH1D(hname, hname, hp.params[0], hp.params[1], hp.params[2]);
      nhist->Sumw2();
      insidemap[*pit] = nhist;
    }
    mymap[*dcit] = insidemap;
  }

  return mymap;
}

map<string, TH1D*> denom_mapper(vector<string> pvector)
{
  map<string, HistParams> hpmap = hist_param_map();
  map<string, TH1D*> histmap;
  for(vector<string>::iterator pit=pvector.begin();
      pit!=pvector.end(); ++pit)
  {
    char hname[50];
    sprintf(hname, "denom_%s", (*pit).c_str());
    HistParams hp = hpmap[(*pit)];
    TH1D* nhist = new TH1D(hname, hname, hp.params[0], hp.params[1], hp.params[2]);
    nhist->Sumw2();
    histmap[*pit] = nhist;
  }
  return histmap;
}

map<string, TH1D*> total_mapper(vector<string> pvector)
{
  map<string, HistParams> hpmap = hist_param_map();
  map<string, TH1D*> histmap;
  for(vector<string>::iterator pit=pvector.begin();
      pit!=pvector.end(); ++pit)
  {
    char hname[50];
    sprintf(hname, "total_%s", (*pit).c_str());
    HistParams hp = hpmap[(*pit)];
    TH1D* nhist = new TH1D(hname, hname, hp.params[0], hp.params[1], hp.params[2]);
    nhist->Sumw2();
    histmap[*pit] = nhist;
  }
  return histmap;
}

void add_unique(vector<string>& tvect, string element)
{
  if(find(tvect.begin(), tvect.end(), element) == tvect.end())
    tvect.push_back(element);
}

int main(int argc, char** argv)
{
  string tname = string(argv[1]) + "_sacrifice";
  TFile* tfile = new TFile(tname.c_str(), "recreate");

#ifdef dbug
  printf("arguments given %i\n", argc);
  for(int i=0; i<argc; ++i)
    printf("arg %i: %s\n", i, argv[i]);
#endif
  // Define all of the data cleaning parameters
  vector<string> dcvector = dc_vector();
#ifdef dbug
  printf("Running on the following data cleaning cuts\n");
  for(vector<string>::iterator it=dcvector.begin(); it!=dcvector.end(); ++it)
    cout << *it << endl;
#endif
  // Define the plotting parameters
  vector<string> pvector = param_vector();
#ifdef dbug
  printf("Running on the following parameters\n");
  for(vector<string>::iterator it=pvector.begin(); it!=pvector.end(); ++it)
    cout << *it << endl;
#endif
  // Double loop, for each dcvector, for each pvector :: Make a cut histogram
  // Structure: map<dccut, map<param, hist>>
  map< string, map<string, TH1D*> > cmap = cut_mapper(dcvector, pvector);
  // Single loop, for each pvector :: Make a denominator and total histogram
  // Structure: map<param, hist>
  map<string, TH1D*> denommap = denom_mapper(pvector);
  map<string, TH1D*> totalmap = total_mapper(pvector);

  // Time to fill the histograms looping over the root file (or files) given
  // Tree is T
  vector<string> filenames;
  for(int i=1; i<argc; i++)
    filenames.push_back(argv[i]);
  // Time to RAT
  RAT::DU::DSReader dsReader( filenames );

  int events = dsReader.GetEntryCount();
#ifdef dbug
  printf("Events: %i\n", events);
#endif
  // events = 5; // For testing purpose, will remove, todo

  RAT::DU::DataCleaningBits rDataCleaningBits = 
    RAT::DU::Utility::Get()->GetDataCleaningBits();
#ifdef dbug
  vector<string> unique_dctests;
#endif

  // Pass information
  const RAT::DS::Meta& rMeta = dsReader.GetMeta();
  const UInt_t nPasses = rMeta.GetPassCount();
  UInt_t pass = rMeta.GetCurrentPass();
  bool setinfo = false;

  // Loop through entries to fill histograms
  for( int ev=0; ev < events; ev++ )
  {
    const RAT::DS::Entry& rDS = dsReader.GetEntry( ev );
    // Only look at triggered events
    if( !(ev%1000) )
      std::cerr << "Processed " << ev << " of " << events << "\r" << std::flush;
    for( int iEV = 0; iEV < rDS.GetEVCount(); iEV++ )
    {
      if( iEV != 0 ) {printf("not iEV 0!\n"); continue;}

      const RAT::DS::EV& rEV = rDS.GetEV(iEV);
      const RAT::DS::DataQCFlags& rDataQCFlags = rEV.GetDataCleaningFlags();
      //const UInt_t pass = rMeta.GetCurrentPass();
      // For the first event only, get Datacleaning pass to use
      if( !setinfo )
      {
        for(int passnum=0; passnum < nPasses; ++passnum)
        {
          const RAT::DS::BitMask& rBitMaskApplied = rDataQCFlags.GetApplied(passnum);
          if(rBitMaskApplied.Get(0))
            pass = passnum;
        }
        printf("Last data cleaning pass: %i\n", pass);
        setinfo = true;
      }

      // Check if there exists a fit result
      if( !rEV.DefaultFitVertexExists() ) continue;
      if( !rEV.GetDefaultFitVertex().ContainsEnergy() ) continue;
      if( !rEV.GetDefaultFitVertex().ValidEnergy() ) continue;
      if( !rEV.GetDefaultFitVertex().ContainsPosition() ) continue;
      if( !rEV.GetDefaultFitVertex().ValidPosition() ) continue;
      if( !rEV.GetDefaultFitVertex().ContainsDirection() ) continue;
      if( !rEV.GetDefaultFitVertex().ValidDirection() ) continue;
      if( !rEV.ClassifierResultExists("isotropy:waterFitter") ) continue;
      if( !rEV.GetClassifierResult("isotropy:waterFitter").GetValid() ) continue;
      if( !rEV.ClassifierResultExists("ITR:waterFitter") ) continue;
      if( !rEV.GetClassifierResult("ITR:waterFitter").GetValid() ) continue;

      // Needed parameters
      int nhits = rEV.GetCalPMTs().GetAllCount();
      double energy = rEV.GetDefaultFitVertex().GetEnergy();
      TVector3 position = rEV.GetDefaultFitVertex().GetPosition();
      double radius = pow(position.Mag()/1000.0,3);
      // This "direction" is weird to me, need to look into this
      double direction = position.Unit().Dot(rEV.GetDefaultFitVertex().GetDirection().Unit());
      double isotropy = rEV.GetClassifierResult("isotropy:waterFitter").GetClassification("snobeta14");
      double itr = rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR");

      // Denom map
      for(vector<string>::iterator pit=pvector.begin();
          pit!=pvector.end(); ++pit)
      {
        if( !(*pit).compare("energy") )
          denommap[*pit]->Fill(energy);
        if( !(*pit).compare("nhits") )
          denommap[*pit]->Fill(double(nhits));
        if( !(*pit).compare("radius") )
          denommap[*pit]->Fill(radius);
        if( !(*pit).compare("direction") )
          denommap[*pit]->Fill(direction);
      }

      // Test if the event passes everything
      bool passes = true;

      // put isotropy cuts here? That's not a good idea, here for later

      // Check the dc flags exist
      if(!rDataQCFlags.ExistFlags(pass)) continue;

      const RAT::DS::BitMask& rBitMaskFlags = rDataQCFlags.GetFlags(pass);
      const RAT::DS::BitMask& rBitMaskApplied = rDataQCFlags.GetApplied(pass);
      for( map<size_t, string>::iterator iDCBits = rDataCleaningBits.GetInverseMapBegin();
          iDCBits != rDataCleaningBits.GetInverseMapEnd(); ++iDCBits )
      {
        // Get the data cleaning words
        string cutName = rDataCleaningBits.GetBitName( iDCBits->first );
        if( !(rBitMaskApplied.Get( iDCBits->first )) ) continue;
        if( rBitMaskFlags.Get( iDCBits->first ) ) continue;
#ifdef debug
        add_unique(unique_dctests, cutName);
#endif
        // Check if the word is in list of words to plot and test for
        bool fill_total=true;
        if( cmap.find(cutName) != cmap.end() )
        {
          for(vector<string>::iterator pit=pvector.begin();
              pit!=pvector.end(); ++pit)
          {
            if( !(*pit).compare("energy") )
            {
              cmap[cutName][*pit]->Fill(energy);
              if(fill_total)
                totalmap[*pit]->Fill(energy);
            }
            if( !(*pit).compare("nhits") )
            {
              cmap[cutName][*pit]->Fill(double(nhits));
              if(fill_total)
                totalmap[*pit]->Fill(double(nhits));
            }
            if( !(*pit).compare("radius") )
            {
              cmap[cutName][*pit]->Fill(radius);
              if(fill_total)
                totalmap[*pit]->Fill(radius);
            }
            if( !(*pit).compare("direction") )
            {
              cmap[cutName][*pit]->Fill(direction);
              if(fill_total)
                totalmap[*pit]->Fill(direction);
            }
          } // Loop over parameter vector
          fill_total = false;
        } // Adding information to cut histograms
      } // Iterating over data cleaning bits
    } // Loop over triggered sub events
  } // Loop over all events in root files
  double sacrifice_count = (totalmap.begin()->second)->GetEntries();
  printf("Analysis done, total sacrifice: %.0f of %i (%.3f\%)\n", 
         sacrifice_count, events, 100*double(sacrifice_count)/events);
  /*
  // Divide all the histograms by the appropriate denominator cmap/denommap
  for( map<string, map<string, TH1D*> >::iterator cmapit=cmap.begin();
       cmapit!=cmap.end(); ++cmapit )
  {
    string cutname = cmapit->first;
    for( map<string, TH1D*>::iterator pit=(cmapit->second).begin();
         pit!=(cmapit->second).end(); ++pit )
    {
      string paramname = pit->first;
      (pit->second)->Divide(denommap[paramname]);
    }
  } 
  // Divide the total map as well
  for( map<string, TH1D*>::iterator tmapit=totalmap.begin();
       tmapit!=totalmap.end(); ++tmapit )
  {
    string paramname = tmapit->first;
    (tmapit->second)->Divide(denommap[paramname]);
  }
  */

  tfile->Write();
#ifdef debug
  for(vector<string>::iterator it=unique_dctests.begin();
      it!=unique_dctests.end(); ++it)
    printf("%s\n", (*it).c_str());
  printf("Total Unique: %i\n", unique_dctests.size() );
#endif

  return 0;
}
