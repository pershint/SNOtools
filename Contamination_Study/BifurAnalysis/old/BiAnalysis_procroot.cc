////////////////////////////////////////////////////////////////////
/// \file BifurcatedAnalysis.cc
///
/// \Estimate the contamination for an input ROOT file with a data set
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

vector<string> DCCuts()
{
    // Returns a list of cuts that will be used as the "Data Cleaning
    // Total Cut".  All assumed perpendicular to Beta14 and ITR.
    const char *vinit[] = {"prescalecut", "crateisotropy","qcluster",
                           "ringoffire","waterblindlow0","waterblindhigh0",
                           "itctimespreadcut","qvnhit","flashergeocut","owlcut",
                           "zerozerocut","qvt","ftscut","junkcut", 
                           "waterblindhigh4","waterblindlow9"};
    int length = sizeof(vinit)/sizeof(vinit[0]);
    vector<string> vec(vinit, vinit+length);
    return vec;
}

vector<string> PrelimCuts()
{
    // Returns a list of the pathological event flags we do not want to
    // Consider when doing the bifurcated analysis. 
    const char *vinit[] = {"prescalecut", "crateisotropy","qcluster",
                           "ringoffire","waterblindlow0","waterblindhigh0",
                           "itctimespreadcut","qvnhit","flashergeocut","owlcut",
                           "zerozerocut","qvt","ftscut","junkcut", 
                           "waterblindhigh4","waterblindlow9"};
    int length = sizeof(vinit)/sizeof(vinit[0]);
    vector<string> vec(vinit, vinit+length);
    return vec;
}

bool hasValidClassifiers(const RAT::DS::EV & ev)
{
    // Checks whether the event has valid Beta14 and ITR parameters.
    // If either are invalid, return false.
    bool areValid = 1;
    Bool_t b14exists = ev.ClassifierResultExists("isotropy:waterFitter");
    Bool_t ITRexists = ev.ClassifierResultExists("ITR:waterFitter");
    if((!(b14exists)) | (!(ITRexists)))
        areValid = 0;
    return areValid;
}

bool isInCutList(string input, vector<string> cuts_list)
{
    /* Compares a fed in DCFlag against the fed in vector of cuts. If in the
     * list, returns true.
    */
    bool isknown = 0;
    for (int i = 0; i < cuts_list.size(); i++){
      if (input == cuts_list[i])
        isknown = 1;
    }
    return isknown;
}

int main(int argc, char** argv)
{
  const string& filename = string(argv[1]);

  //Constants defining range of valid Beta14 and ITR cuts
  //Data Cleaning cut applications should be changed in functions above.
  double b14_low = -0.12;
  double b14_high = 0.95;
  double itr_low = 0.5;

  //Here, we have the variables used to solve the bifurcated analysis
  //Equation and calculate our contamination.
  //a: passes both cuts, b: passes DC, fails Classifiers
  //c: fails DC, passes classifiers, d: fails both cuts
  int a = 0;
  int b = 0;
  int c = 0;
  int d = 0;

  // Read root files
  RAT::DU::DSReader dsReader( filename );

  cout << "RUNNING BIFURCATION ANALYSIS ON " << filename << endl;
  cout << "CHOSEN BETA14 RANGES: " << b14_low << "," << b14_high << endl;
  cout << "CHOSEN ITR RANGES: " << itr_low << ",None" << endl;
  
  RAT::DU::DataCleaningBits rDataCleaningBits = RAT::DU::Utility::Get()->GetDataCleaningBits(); // To get the data cleaning bits

  // Loop through entries in rootfile
  for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ){
    if( iEntry % 10000 == 0 ){ cout << "iEntry: " << iEntry << " / " << dsReader.GetEntryCount() << endl; }
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
    
    //These flags determine which box the event goes into for the
    //Bifurcated analysis.
    bool DC_clean = 1;
    bool Class_clean = 1;

    // Look at each triggered event in the entry
    for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ ){
      if( iEV != 0 )
        continue;
      const RAT::DS::EV& rEV = rDS.GetEV( iEV );


      //Let's first look to make sure a valid B14 and ITR exist.
      if (!hasValidClassifiers(rEV))
        continue;
      //They're valid!  So, now check if the event's B14 and ITR are
      //Within the defined ranges for this analysis.
      try{
        const RAT::DS::ClassifierResult& b14Result = 
            rEV.GetClassifierResult("isotropy:waterFitter");
        const RAT::DS::ClassifierResult& itrResult = 
            rEV.GetClassifierResult("ITR:waterFitter");
        double b14 = b14Result.GetClassification("snobeta14");
        double itr = itrResult.GetClassification("ITR");
      } catch( RAT::DS::ClassifierResult::NoClassificationError& e){
        cout << "Classifier failed this readout.  Going to next event" << endl;
        continue;
      }
       //FIXME: make this a function so you're not copypasta-ing code?
      const RAT::DS::ClassifierResult& b14Result = 
          rEV.GetClassifierResult("isotropy:waterFitter");
      const RAT::DS::ClassifierResult& itrResult = 
          rEV.GetClassifierResult("ITR:waterFitter");
      double b14 = b14Result.GetClassification("snobeta14");
      double itr = itrResult.GetClassification("ITR");
      if (!((b14 > b14_low) & (b14 < b14_high)))
        Class_clean = 0;
      else if (!(itr > itr_low))
        Class_clean = 0;

      //Cool.  Now, we move to the DCflags part of the event characterization.
      const RAT::DS::DataQCFlags& rDataQCFlags = rEV.GetDataCleaningFlags();
      const UInt_t pass = rDataQCFlags.GetLatestPass();
      // Check flags exist
      if(!(rDataQCFlags.ExistFlags(pass))) // Checks flags exist
        continue;

      // Flags exist! For this event, get the flags
      const RAT::DS::BitMask& rBitMaskFlags = rDataQCFlags.GetFlags(pass);
      const RAT::DS::BitMask& rBitMaskApplied = rDataQCFlags.GetApplied(pass);

      //FIXME: Need to have a check that the Flags are all a part of the
      //Applied mask
      
      // Loop through bits; we will determine if the event has a preliminary
      // Cut, or one of the DC cuts that make it "dirty" in the Bif. Analysis
      for( map<size_t, string>::iterator iDCBits = rDataCleaningBits.GetInverseMapBegin();
           iDCBits != rDataCleaningBits.GetInverseMapEnd(); iDCBits++ ){
        Bool_t rApplied = rBitMaskApplied.Get( iDCBits->first );
        if( !rApplied ) // Checks cut was applied in processing
          continue;
        Bool_t rFlag = rBitMaskFlags.Get( iDCBits->first );
        if( rFlag )  //Checks if this particular event was flagged
          continue;
        string cutName = rDataCleaningBits.GetBitName( iDCBits->first ).c_str();
        //Check if this cut is in the preliminary cuts mask
       if (isInCutList(cutName, PrelimCuts()))
         continue;
         //FIXME: This will only break the bit for loop; want to get out
         //And move to the next event here
       if (isInCutList(cutName, DCCuts()))
         DC_clean = 0;
         break;  //it's dirty; no need to look further in bitword
      } // bit loop


      //Next, we look at the DC_clean and Class_clean bools to see how to
      //Classify the event.  Depending on which box it goes to, increment
      //Some integers that correspond to a,b,c, and d in the SNO data
      //cleaning salt document.
      if ((DC_clean) & (Class_clean))
        a+=1;
      else if ((DC_clean) & !(Class_clean))
        b+=1;
      else if (!(DC_clean) & (Class_clean))
        c+=1;
      else if (!(DC_clean) & !(Class_clean))
        d+=1;
      else
        cout << "How did you get here? Better check for bugs" << endl;
    } // EV
  } // Entry
  //Just print a, b, c, and d for debugging here
  cout << "a: " << a << endl;
  cout << "b: " << b << endl;
  cout << "c: " << c << endl;
  cout << "d: " << d << endl;
}
