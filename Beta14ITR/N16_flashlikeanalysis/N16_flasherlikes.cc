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

bool is_flashercut(string& cutname)
{
  const char *vinit[] = {"qvnhit","qvt","neckcut","crateisotropy","qcluster",
                         "ftscut","flashergeocut","itctimespreadcut"};
  int length = sizeof(vinit)/sizeof(vinit[0]);
  vector<string> flashercuts(vinit, vinit+length);
  bool isflashcut = false;
  for (int i=0; i<flashercuts.size(); i++)
  {
    if(cutname == flashercuts[i])
      isflashcut = true;
  }
  return isflashcut;
}

bool is_qvcut(string& cutname)
{
  const char *vinit[] = {"qvnhit","qvt"};
  int length = sizeof(vinit)/sizeof(vinit[0]);
  vector<string> qvcuts(vinit, vinit+length);
  bool isqvcut = false;
  for (int i=0; i<qvcuts.size(); i++)
  {
    if(cutname == qvcuts[i])
      isqvcut = true;
  }
  return isqvcut;
}

bool has_qvcut(const RAT::DS::EV & ev)
{
  bool hasqvcut = false;
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
    if(is_qvcut(cutName))
      cout << "HAS CUT " << cutName << endl;
      hasqvcut = true;
  }
  return hasqvcut;
}

bool has_flashercut(const RAT::DS::EV & ev)
{
  bool hasflashercut = false;
  RAT::DU::DataCleaningBits rDataCleaningBits = RAT::DU::Utility::Get()->GetDataCleaningBits();
  const RAT::DS::DataQCFlags& rDataQCFlags = ev.GetDataCleaningFlags();
  const UInt_t pass = rDataQCFlags.GetLatestPass();
  const RAT::DS::BitMask& rBitMaskFlags = rDataQCFlags.GetFlags(pass);
  const RAT::DS::BitMask& rBitMaskApplied = rDataQCFlags.GetApplied(pass);
  // Loop through bits, search for flasher-like events
  for( map<size_t, string>::iterator iDCBits = rDataCleaningBits.GetInverseMapBegin();
          iDCBits != rDataCleaningBits.GetInverseMapEnd(); iDCBits++ )
  {
    Bool_t rApplied = rBitMaskApplied.Get( iDCBits->first );
    if( !rApplied ) // Checks cut has been applied on this run
      continue;
    Bool_t rFlag = rBitMaskFlags.Get( iDCBits->first ); //see if event is clean
    if( rFlag )
      continue;
    //This bit was flagged; see if the bit name is a flasher name
    string cutName = rDataCleaningBits.GetBitName( iDCBits->first ).c_str();
    if(is_flashercut(cutName))
      hasflashercut = true;
  }
  return hasflashercut;
}


//Duplicate of code in /rat/example/root/PlotDataCleaningCuts.cc, with
//The addition of selecting only the FECD events
TH1D* MakeDCFHist_FECD(const string& fileName)
{
  RAT::DU::DSReader dsReader( fileName );

  RAT::DU::DataCleaningBits rDataCleaningBits = RAT::DU::Utility::Get()->GetDataCleaningBits(); // To get the data cleaning bits

  size_t numberOfCuts = rDataCleaningBits.GetInverseMapLast()->first+1;
  TH1D* hCutFlags = new TH1D( "hCutFlags", "Data cleaning cut flags", numberOfCuts, 0, numberOfCuts );

  // Loop through entries in rootfile
  for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ){
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );

    // Look at each triggered event in the entry
    for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ ){
      const RAT::DS::EV& rEV = rDS.GetEV( iEV );
      bool isN16 = is_n16_event(rEV);
      if (!(isN16))
        continue;
      const RAT::DS::DataQCFlags& rDataQCFlags = rEV.GetDataCleaningFlags();
      const UInt_t pass = rDataQCFlags.GetLatestPass();

      // Check flags exist
      if(rDataQCFlags.ExistFlags(pass)){ // Checks flags exist
        const RAT::DS::BitMask& rBitMaskFlags = rDataQCFlags.GetFlags(pass);
        const RAT::DS::BitMask& rBitMaskApplied = rDataQCFlags.GetApplied(pass);

        // Loop through bits
        for( map<size_t, string>::iterator iDCBits = rDataCleaningBits.GetInverseMapBegin();
             iDCBits != rDataCleaningBits.GetInverseMapEnd(); iDCBits++ )
          {
            Bool_t rApplied = rBitMaskApplied.Get( iDCBits->first );
            if( rApplied ){ // Checks cut has been applied
              Bool_t rFlag = rBitMaskFlags.Get( iDCBits->first );
              if( !rFlag ) hCutFlags->Fill( iDCBits->first ); // Fill if event fails the cut
            }
          }

      }
    } // EV
  } // Entry

  hCutFlags->GetYaxis()->SetTitle( "Number of events" );
  for( map<size_t, string>::iterator iDCBits = rDataCleaningBits.GetInverseMapBegin();
       iDCBits != rDataCleaningBits.GetInverseMapEnd(); iDCBits++ )
    {
      hCutFlags->GetXaxis()->SetBinLabel(iDCBits->first+1, rDataCleaningBits.GetBitName(iDCBits->first).c_str());
    }
  return hCutFlags;
}


int main(int argc, char** argv)
{
  // Define our histograms

  //Define histograms to fill in
  TH2F* h_N16_dirty = new TH2F("h_N16_dirty", "h_N16_dirty", 50,0.,1.,50,-0.5,2.);
  TH2F* h_N16_clean = new TH2F("h_N16_clean", "h_N16_clean", 50,0.,1.,50,-0.5,2.);
  TH2F* h_nonFECD_hasflashcut = new TH2F("h_nonFECD_hasflashcut", "h_nonFECD_hasflashcut", 50,0.,1.,50,-0.5,2.);
  TH1D* h_fitr_hasflashcut = new TH1D("h_fitr_hasflashcut", "h_fitr_hasflashcut",20,0.,9000.);
  TH1D* h_nhits_hasflashcut = new TH1D("h_nhits_hasflashcut", "h_nhits_hasflashcut",40,20.,100);
 
  //To build this mask, see snopl.us/docs/rat/user_manual/html/node226.html
  //int prescaleonly
  bool plotClean = true;
  bool plotDirty = true;

  //Create the DCFlags histogram for the first file only.
//  TH1D* h_N16_DCFlags = MakeDCFHist_FECD(argv[1]);
  for (int f=1; f<argc; f++)
  {
    //load up the next root file
    RAT::DU::DSReader dsReader( argv[f] );
    

    //Get the analysis mask as defined in ROOT
    ULong64_t rAnalysisMask = RAT::GetDataCleaningWord ( "analysis_mask" );

    for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
    {
      const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
      // Look at each triggered event in the entry
      for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ )
      {
        //REALLY WANT THIS SPOT?
        if( iEV != 0 )
          continue;
        const RAT::DS::EV& rEV = rDS.GetEV( iEV );
        bool isN16 = is_n16_event(rEV);
        if (isN16)
        {
          //We're at an N16 event.  Get the beta14 and ITR parameters.
          Bool_t b14exists = rEV.ClassifierResultExists("isotropy:waterFitter");
          Bool_t ITRexists = rEV.ClassifierResultExists("ITR:waterFitter");
          if((!(b14exists)) | (!(ITRexists)))
            continue;
          try{
            const RAT::DS::ClassifierResult& b14Result = rEV.GetClassifierResult("isotropy:waterFitter");
            const RAT::DS::ClassifierResult& itrResult = rEV.GetClassifierResult("ITR:waterFitter");
            double b14Value = b14Result.GetClassification("snobeta14");
            double itrValue = itrResult.GetClassification("ITR");
            //check if the event is clean or dirty according to the analysi mask
            Bool_t eventIsClean = RAT::EventIsClean( rEV, rAnalysisMask );
            if (eventIsClean)
              h_N16_clean->Fill(itrValue,b14Value);
            else
            {
              h_N16_dirty->Fill(itrValue,b14Value);
              cout << "DIRTY IN FECD GTID: " << rEV.GetGTID() << endl;
              if (has_qvcut(rEV))
                cout << "STILL GTID: " << rEV.GetGTID() << endl;
            }
          } catch( RAT::DS::ClassifierResult::NoClassificationError& e)
          {
            cout << "Classifier failed this event.  continuing." << endl;
            continue;
          }
        } else
        {
          //only look at nhit > 20
          if( rEV.GetCalPMTs().GetAllCount() < 21)
            continue;
          Bool_t b14exists = rEV.ClassifierResultExists("isotropy:waterFitter");
          Bool_t ITRexists = rEV.ClassifierResultExists("ITR:waterFitter");
          if((!(b14exists)) | (!(ITRexists)))
            continue;
          try{
            const RAT::DS::ClassifierResult& b14Result = rEV.GetClassifierResult("isotropy:waterFitter");
            const RAT::DS::ClassifierResult& itrResult = rEV.GetClassifierResult("ITR:waterFitter");
            double b14Value = b14Result.GetClassification("snobeta14");
            double itrValue = itrResult.GetClassification("ITR");
            //check if the event is clean or dirty according to the analysi mask
            if (has_flashercut(rEV))
            {
              h_nonFECD_hasflashcut->Fill(itrValue,b14Value);
              //Also fill in the fitted radius histogram
              if( !rEV.GetFitResult("waterFitter").GetVertex(0).ContainsPosition() )
                continue;
              if( !rEV.GetFitResult("waterFitter").GetVertex(0).ValidPosition() )
                continue;
              const RAT::DS::FitResult FitRes = rEV.GetFitResult("waterFitter");
              const RAT::DS::FitVertex FitVer = FitRes.GetVertex(0); //get first fit
              TVector3 posn = FitVer.GetPosition();
              h_fitr_hasflashcut->Fill(posn.Mag());
              h_nhits_hasflashcut->Fill(rEV.GetCalPMTs().GetAllCount());
            }
          } catch( RAT::DS::ClassifierResult::NoClassificationError& e)
          {
            cout << "Classifier failed this event.  continuing." << endl;
            continue;
          }
        }
      }//end looking at events in each entry
    }//end loop through entries in root file
  }//end loop through all rootfiles

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

  h_nonFECD_hasflashcut->SetMarkerStyle(20);
  h_nonFECD_hasflashcut->SetMarkerColor(15);
  h_nonFECD_hasflashcut->SetMarkerSize(0.7);
  h_nonFECD_hasflashcut->GetXaxis()->SetTitle("ITR");
  h_nonFECD_hasflashcut->GetYaxis()->SetTitle("Beta14");

  h_fitr_hasflashcut->GetXaxis()->SetTitle("Radius (cm)");
  h_fitr_hasflashcut->GetYaxis()->SetTitle("Events");

  h_nhits_hasflashcut->GetXaxis()->SetTitle("nhit");
  h_nhits_hasflashcut->GetYaxis()->SetTitle("Events");

  TFile* thehists = new TFile("booyaw.root","RECREATE");
  if(plotDirty)
    thehists->Add(h_N16_dirty);
  if(plotClean)
    thehists->Add(h_N16_clean);
//  thehists->Add(h_N16_DCFlags);
  thehists->Add(h_nonFECD_hasflashcut);
  thehists->Add(h_fitr_hasflashcut);
  thehists->Add(h_nhits_hasflashcut);
  thehists->Write();
  thehists->Close();
  return 0;
} //End main
