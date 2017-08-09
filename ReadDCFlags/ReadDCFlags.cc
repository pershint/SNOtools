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


vector<string> known_cuts()
{
    // List of some known cuts that are not the TPMuonFollower.
    // Want to weed out the values so we're not printing so much.
    const char *vinit[] = {"prescalecut", "crateisotropy","qcluster",
                           "ringoffire","waterblindlow0","waterblindhigh0",
                           "itctimespreadcut","qvnhit","flashergeocut","owlcut",
                           "zerozerocut","qvt","ftscut","junkcut", 
                           "waterblindhigh4","waterblindlow9"};
    int length = sizeof(vinit)/sizeof(vinit[0]);
    vector<string> vec(vinit, vinit+length);
    return vec;
}

bool isKnownDCFlag(string input)
{
    /* Compares a fed in DCFlag against the known_cuts vector. If in the
     * list, returns true.
    */
    bool isknown = 0;
    vector<string> knowncuts = known_cuts();
    for (int i = 0; i < knowncuts.size(); i++){
      if (input == knowncuts[i])
        isknown = 1;
    }
    return isknown;
}

int main(int argc, char** argv)
{
  const string& filename = string(argv[1]);

  // Read root files
  RAT::DU::DSReader dsReader( filename );
  cout << "CHECKING FOR TPMFC IN FILE " << filename << endl;
  RAT::DU::DataCleaningBits rDataCleaningBits = RAT::DU::Utility::Get()->GetDataCleaningBits(); // To get the data cleaning bits

  // Loop through entries in rootfile
  for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ){
    if( iEntry % 10000 == 0 ){ cout << "iEntry: " << iEntry << " / " << dsReader.GetEntryCount() << endl; }
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );

    // Look at each triggered event in the entry
    for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ ){
      if( iEV != 0 )
        continue;
      const RAT::DS::EV& rEV = rDS.GetEV( iEV );
      const RAT::DS::DataQCFlags& rDataQCFlags = rEV.GetDataCleaningFlags();
      const RAT::DS::Meta& rMeta = dsReader.GetMeta();
      const UInt_t pass = rMeta.GetCurrentPass();


      // Get Nhits for all and OWL tubes
      int calnormnhits = rEV.GetCalPMTs().GetAllCount();
      int owlhits = rEV.GetCalPMTs().GetOWLCount();
      int neckhits = rEV.GetCalPMTs().GetNeckCount();

      // Check flags exist
      if(!(rDataQCFlags.ExistFlags(pass))) // Checks flags exist
        continue;
      // Flags exist! For this event, get the flags
      const RAT::DS::BitMask& rBitMaskFlags = rDataQCFlags.GetFlags(pass);
      const RAT::DS::BitMask& rBitMaskApplied = rDataQCFlags.GetApplied(pass);

      // Loop through bits
      for( map<size_t, string>::iterator iDCBits = rDataCleaningBits.GetInverseMapBegin();
           iDCBits != rDataCleaningBits.GetInverseMapEnd(); iDCBits++ ){
        Bool_t rApplied = rBitMaskApplied.Get( iDCBits->first );
        if( !rApplied ) // Checks cut has been applied
          continue;
        Bool_t rFlag = rBitMaskFlags.Get( iDCBits->first );
        if( rFlag )  //I can't remember what this check is for...
          continue;
        string cutName = rDataCleaningBits.GetBitName( iDCBits->first ).c_str();
        if (calnormnhits > 100) {
          cout << "# total inward pmts hit:" << calnormnhits << endl;
          cout << "# total owl pmts hit:" << owlhits << endl;
          cout << "# total neck pmts hit:" << neckhits << endl;
// Uncomment the below lines to only see cuts NOT defined in the
// isKnownDCFlag function
//        }
//        if(!isKnownDCFlag(cutName)){
          cout << "NEW BIT CUT NAME:" << cutName << endl;
          cout << "IN EVENT " << iEV << endl;
          cout << "IN PASS " << pass << endl;
          cout << "IN ENTRY "<< iEntry << endl;
        }
      } // bit loop

      // Divide cut events by denominator to get sacrifice
    } // EV
  } // Entry
}
