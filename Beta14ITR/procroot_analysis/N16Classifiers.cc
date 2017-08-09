////////////////////////////////////////////////////////////////////
/// \file N16Sac.cc
///
/// \Program to read out data quality flags from a processing ROOT file
////////////////////////////////////////////////////////////////////
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DataCleaningUtility.hh>
#include <RAT/DU/DataCleaningBits.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/DataQCFlags.hh>
#include <RAT/DS/BitMask.hh>
#include <RAT/DS/Meta.hh>

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

//Taken from some code by Eric Marzec
bool is_n16_event(const RAT::DS::EV & ev)
{
  const RAT::DS::CalPMTs& PMTs = ev.GetCalPMTs();
  for( size_t iPMT = 0; iPMT < PMTs.GetFECDCount(); iPMT++ )
  {
    const RAT::DS::PMTCal& pmt = PMTs.GetFECDPMT( iPMT );
    int id = pmt.GetID();
    if(id == 9188)
      return 1;
  }
}

vector<string> pathological_events()
{
    // List of some known cuts that are not the TPMuonFollower.
    // Want to weed out the values so we're not printing so much.
    const char *vinit[] = {"junkcut", "muoncut","itctimespreadcut"};
    int length = sizeof(vinit)/sizeof(vinit[0]);
    vector<string> vec(vinit, vinit+length);
    return vec;
}

bool is_pathological(string& cutname)
{
  vector<string> pev = pathological_events();
  bool isqvcut = false;
  for (int i=0; i<pev.size(); i++)
  {
    if(cutname == pev[i])
      isqvcut = true;
  }
  return isqvcut;
}

bool has_pathologicalcut(const RAT::DS::EV & ev)
{
  bool haspecut = false;
  RAT::DU::DataCleaningBits rDataCleaningBits = RAT::DU::Utility::Get()->GetDataCleaningBits();
  const RAT::DS::DataQCFlags& rDataQCFlags = ev.GetDataCleaningFlags();
  const UInt_t pass = rDataQCFlags.GetLatestPass();
  const RAT::DS::BitMask& rBitMaskFlags = rDataQCFlags.GetFlags(pass);
  const RAT::DS::BitMask& rBitMaskApplied = rDataQCFlags.GetApplied(pass);
  // Loop through bits, search for qv-like events
  for( map<size_t, string>::iterator iDCBits = rDataCleaningBits.GetInverseMapBegin();
          iDCBits != rDataCleaningBits.GetInverseMapEnd(); iDCBits++ )
  {
    Bool_t rApplied = rBitMaskApplied.Get( iDCBits->first );
    if( !rApplied ) // Checks cut has been applied on this run
      continue;
    Bool_t rFlag = rBitMaskFlags.Get( iDCBits->first ); //see if event is clean
    if( rFlag )
      continue;
    //This bit was flagged; see if the bit name is a qv name
    string cutName = rDataCleaningBits.GetBitName( iDCBits->first ).c_str();
    if(is_pathological(cutName))
    {
      cout << "HAS CUT " << cutName << endl;
      haspecut = true;
    }
  }
  return haspecut;
}

int main(int argc, char** argv)
{
  // Define our histograms

  //Define histograms to fill in
  TH2F* h_N16_dirty = new TH2F("h_N16_dirty", "h_N16_dirty", 50,0.,1.,50,-0.5,2.);
  TH2F* h_N16_clean = new TH2F("h_N16_clean", "h_N16_clean", 50,0.,1.,50,-0.5,2.);

  //some cut selections you could apply
  //To build this mask, see snopl.us/docs/rat/user_manual/html/node226.html
  //int prescaleonly
  bool plotClean = true;
  bool plotDirty = true;

  int pathTrigMask = 0b110001100000;  //ESum , PED, and PulseGT bits in trigmask
  //If splitting dirty and clean according to beta14 and ITR values, these
  //Are the values according to cuts from SNO required
  double b14_low = -0.12;
  double b14_high = 0.95;
  double itr_low = 0.55;

  int nofit = 0;
  int pathological = 0;

  for (int f=1; f<argc; f++)
  {
    //load up the next root file
    RAT::DU::DSReader dsReader( argv[f] );

    //Get the analysis mask as defined in ROOT
    ULong64_t rAnalysisMask = RAT::GetDataCleaningWord ( "analysis_mask" );
    for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
    {
      if( iEntry % 100000 == 0 ){ cout << "iEntry: " << iEntry << " / " << dsReader.GetEntryCount() << endl; }
      const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
      // Look at each triggered event in the entry
      for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ )
      {
        //REALLY WANT THIS SPOT?
        if( iEV != 0 )
          continue;
        const RAT::DS::EV& rEV = rDS.GetEV( iEV );
        bool isN16 = is_n16_event(rEV);
        if (!(isN16))
          continue;

        //We're at an N16 event.  Get the beta14 and ITR parameters.
        //FIRST, see if the event is a pathological event.  If so, continue.
        bool ispath = has_pathologicalcut(rEV);
        if (ispath)
        {
          pathological+=1;
          continue;
        }
        Int_t trigWord = rEV.GetTrigType();
        if (trigWord & pathTrigMask)
        {
          pathological+=1;
          continue;
        }

        Bool_t b14exists = rEV.ClassifierResultExists("isotropy:waterFitter");
        Bool_t ITRexists = rEV.ClassifierResultExists("ITR:waterFitter");
        if((!(b14exists)) | (!(ITRexists)))
        {
            nofit+=1;
            continue;
        }
        try{
          const RAT::DS::ClassifierResult& b14Result = rEV.GetClassifierResult("isotropy:waterFitter");
          const RAT::DS::ClassifierResult& itrResult = rEV.GetClassifierResult("ITR:waterFitter");
//      const RAT::DS::DataQCFlags& rDataQCFlags = rEV.GetDataCleaningFlags();
          double b14Value = b14Result.GetClassification("snobeta14");
          double itrValue = itrResult.GetClassification("ITR");
          //check if the event is clean or dirty according to the analysi mask
          //const RAT::DS::Meta& rMeta = dsReader.GetMeta();
          //Bool_t eventIsClean = RAT::EventIsClean( rEV, rAnalysisMask );
          //if (eventIsClean)
          if ((b14Value > b14_high) | (b14Value<b14_low) | (itrValue < itr_low))
            h_N16_dirty->Fill(itrValue,b14Value);
          else
            h_N16_clean->Fill(itrValue,b14Value);
        } catch( RAT::DS::ClassifierResult::NoClassificationError& e)
        {
          cout << "Classifier failed an event.  continuing."  << endl;
        }
      }//end looking at events in each entry
    }//end loop through entries in root file
  }//end loop through all rootfiles

  cout << "NO. N16s WITH NO BETA14 OR ITR: " << nofit << endl;
  cout << "NO. N16s WITH A PATHOLOGICAL CUT: " << pathological << endl;
  h_N16_dirty->SetMarkerStyle(20);
  h_N16_dirty->SetMarkerColor(2);
  h_N16_dirty->SetMarkerSize(0.7);
  h_N16_dirty->GetXaxis()->SetTitle("ITR");
  h_N16_dirty->GetYaxis()->SetTitle("Beta14");

  h_N16_clean->SetMarkerStyle(20);
  h_N16_clean->SetMarkerColor(4);
  h_N16_clean->SetMarkerSize(0.7);
  h_N16_clean->GetXaxis()->SetTitle("ITR");
  h_N16_clean->GetYaxis()->SetTitle("Beta14");

  TFile* thehists = new TFile("booyaw.root","CREATE");
  if(plotDirty)
    thehists->Add(h_N16_dirty);
  if(plotClean)
    thehists->Add(h_N16_clean);
  thehists->Write();
  thehists->Close();
//  delete h_N16_clean;
//  delete h_N16_dirty;
} //End main
