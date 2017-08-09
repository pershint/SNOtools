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

#include <TApplication.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TFile.h>

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
                           "zerozerocut","qvt","ftscut","junkcut"};
    int length = sizeof(vinit)/sizeof(vinit[0]);
    vector<string> vec(vinit, vinit+length);
    return vec;
}

string Run_Number(string filename)
{
    //Gets the run number out of a DCAProc filename. 
    //Only works for files with format Procr####_DCAPROC.root
    string name = filename;
    int namesize = name.size();
    name.string::erase(namesize-13,13); //Erase _DCAPROC.root at end of filename
    name.string::erase(0,5); //Erase Procr at start of string
    return name;
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
  TApplication* myapp = new TApplication("myapp",0,0);

  // Define histogram ranges
  // Define Nhits histograms
  int mfcl_time_min = 3870;
  int mfcl_time_max = 3880;
  int mfcs_time_min = 4270;
  int mfcs_time_max = 4280;
  int mttime_min = 4470;
  int mttime_max = 4480;

  // Read root files
  TFile* tfile = TFile::Open(filename.c_str(), "UPDATE");
  cout << "CHECKING FOR TPMFC IN FILE " << filename << endl;

  TCanvas *c1 = new TCanvas("c1","c1",800,1200);
  c1->Divide(2,1);
  TCanvas *c2 = new TCanvas("c2","c2",800,1200);
  c2->Divide(1,1);
  TCanvas *c3 = new TCanvas("c3","c3",800,1200);
  c3->Divide(2,1);

  //Get Time difference histograms from DCAProc file
  cout << "TRYING TO GET THE HISTOGRAM NOW" << endl;
  TH1F* h_mfcs_time = (TH1F*)tfile->Get("tpmuonfollowercut-short_time");
  TH1F* h_mfcl_time = (TH1F*)tfile->Get("tpmuonfollowercut-long_time");
  TH1F* h_mt_time = (TH1F*)tfile->Get("muontag_time");
  TH1F* h_mfcs_nhit = (TH1F*)tfile->Get("tpmuonfollowercut-short_nhit");
  TH1F* h_mfcl_nhit = (TH1F*)tfile->Get("tpmuonfollowercut-long_nhit");
  h_mfcs_time->GetXaxis()->SetRangeUser(mfcs_time_min,mfcs_time_max);
  h_mfcl_time->GetXaxis()->SetRangeUser(mfcl_time_min,mfcl_time_max);
  h_mt_time->GetXaxis()->SetRangeUser(mttime_min,mttime_max);

  //Generate titles for time histograms
  string run_num = Run_Number(filename);
  string short_descr = "Time Diff. between short muon flag events, run ";
  string long_descr = "Time Diff. between long muon flag events, run ";
  string mttime_descr = "Time Diff. between muon tagged events, run ";
  string shortnhit_descr = "Nhit distribution for muon flag events, run ";
  string longnhit_descr = "Nhit distribution for long muon flag events, run ";
  short_descr.append(run_num);
  long_descr.append(run_num);
  shortnhit_descr.append(run_num);
  longnhit_descr.append(run_num);
  mttime_descr.append(run_num);
  const char* short_title = short_descr.c_str();
  const char* long_title = long_descr.c_str();
  const char* shortnhit_title = shortnhit_descr.c_str();
  const char* longnhit_title = longnhit_descr.c_str();
  const char* mttime_title = mttime_descr.c_str();
  h_mfcs_time->SetTitle(short_title);
  h_mfcl_time->SetTitle(long_title);
  h_mfcs_nhit->SetTitle(shortnhit_title);
  h_mfcl_nhit->SetTitle(longnhit_title);
  h_mt_time->SetTitle(mttime_title);

  cout << "DRAWING" << endl;
  c1->cd(1);
  h_mfcs_time->Draw();
  c1->cd(2);
  h_mfcl_time->Draw();
  c2->cd(1);
  h_mt_time->Draw();
  c3->cd(1);
  h_mfcs_nhit->Draw();
  c3->cd(2);
  h_mfcl_nhit->Draw();

  myapp->Run();
  //delete histograms
  delete h_mfcl_time;
  delete h_mfcs_time;
  delete h_mfcl_nhit;
  delete h_mfcs_nhit;
  delete h_mt_time;
  return 0;
}
