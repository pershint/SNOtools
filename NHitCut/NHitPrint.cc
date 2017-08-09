////////////////////////////////////////////////////////////////////
/// \file ReadDCFlags_NoPlots.cc
///
/// \Program to read out data quality flags from a processing ROOT file
////////////////////////////////////////////////////////////////////

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

#include <typeinfo>
#include <string>
#include <map>
#include <iterator>
using namespace std;


vector<int> NHitCuts()
{
    // List of some known cuts that are not the TPMuonFollower.
    // Want to weed out the values so we're not printing so much.
    int vinit[] = {40,60};
    int length = sizeof(vinit)/sizeof(vinit[0]);
    vector<int> vec(vinit, vinit+length);
    return vec;
}

int main(int argc, char** argv)
{
  const string& filename = string(argv[1]);
  const string material = string(argv[2]);


  // Read root files
  RAT::DU::DSReader dsReader( filename );
  cout << "PRINTING YOUR NHITS FOR FILE " << filename << endl;
  vector<int> nhitcuts = NHitCuts();
  for (int i = 0; i < nhitcuts.size(); i++){
    int yaknow = nhitcuts[i];
    cout << "THE FOLLOWING HAVE NHIT < " << yaknow << endl;
    // Loop through entries in rootfile
    for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ){
      if( iEntry % 10000 == 0 ){ cout << "iEntry: " << iEntry << " / " << dsReader.GetEntryCount() << endl; }
      const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );

      // Look at each triggered event in the entry
      for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ ){
        const RAT::DS::EV& rEV = rDS.GetEV( iEV );
        int nhits = rEV.GetCalPMTs().GetAllCount();
        if( nhits > yaknow){
          cout << "EVENT: " << iEV << endl;
          cout << "NHITS: " << nhits << endl;
        }
      }
    }
  }
}
