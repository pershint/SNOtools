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


vector<int> GTIDsOfInterest()
{
    // List of some known cuts that are not the TPMuonFollower.
    // Want to weed out the values so we're not printing so much.
    int vinit[] = {7432879,7432880,7432881,7432882,7432883,7432884,7432885,7432886};
//10112680,10112684,10112685,10112693,10112697,10112723,
//    10112726,10112727,10112728,10112729,10112731,10112738,10112739,10112745};
    int length = sizeof(vinit)/sizeof(vinit[0]);
    vector<int> vec(vinit, vinit+length);
    return vec;
}


void RunFTSCut(const RAT::DS::EV& ev)
{
  RAT::DBLinkPtr fLClean = RAT::DB::Get()->GetLink("DATACLEANING","ftscut");
  int fClusterMode = fLClean->GetI("clustermode");
  int fDeltaTThresh = fLClean->GetI("deltatthresh");
  int fDistThresh = fLClean->GetI("distthresh");
  int fCountThresh = fLClean->GetI("countthresh");
  int fMedianCut = fLClean->GetD("mediancut");
  int fMaxPairs = 10000;

  int cal_count,pair_count,cluster;
  int iccc,jccc;
  float tubedist, deltat, dtmedian;
  vector<const RAT::DS::PMTCal *> pmt_list;
  vector<int> ccc_hits(19*16*32,0);
  vector<float> deltat_list;

  bool fPassFlag = true;

  // get a list of the pmts with good
  // time calibrations
  cal_count = ev.GetCalPMTs().GetCount();
  cout << "#NUM OF CALIBRATED PMTS: " << cal_count << endl;
  for (int ipmt=0;ipmt<cal_count;ipmt++){
    const RAT::DS::PMTCal* pmt = &ev.GetCalPMTs().GetPMT(ipmt);
    // if time is calibrated
    // add to the list
    if (abs(pmt->GetTime()) < 8000){
      pmt_list.push_back(pmt);
      ccc_hits[pmt->GetID()] = 1;
    }
  }

  pair_count = 0;
  cluster = 0;
  // loop over all pairs of hit pmts
  for (size_t ipmt=0;ipmt<pmt_list.size();ipmt++){
    if (pair_count >= fMaxPairs)
      break;
    for (size_t jpmt=ipmt;jpmt<pmt_list.size();jpmt++){
      if (pair_count >= fMaxPairs)
        break;
      if (ipmt == jpmt)
        continue;

      iccc = pmt_list[ipmt]->GetID();
      jccc = pmt_list[jpmt]->GetID();
      TVector3 ipos = RAT::DU::Utility::Get()->GetPMTInfo().GetPosition(iccc);
      TVector3 jpos = RAT::DU::Utility::Get()->GetPMTInfo().GetPosition(jccc);

      // compute the linear distance between the
      // pmts and the absolute value of their time
      // difference
      tubedist = sqrt(pow((ipos.X() - jpos.X()),2) + pow((ipos.Y() - jpos.Y()),2) + pow((ipos.Z() - jpos.Z()),2));
      deltat = abs(pmt_list[ipmt]->GetTime() - pmt_list[jpmt]->GetTime());

      // count this pmt pair as a valid pair and store the time diff if
      // deltat < threshold in ns
      // pmts closer than threshold in cm
      // pmt channels not in a cluster
      if ((deltat < fDeltaTThresh) and (tubedist < fDistThresh) and (!cluster)){


        // if cluster pair removal is enabled, check to see if
        // every tube between these two has been hit
        if (fClusterMode){
          int high_ccc,low_ccc;
          if (iccc > jccc){
            high_ccc = iccc;
            low_ccc = jccc;
          }else{
            high_ccc = jccc;
            low_ccc = iccc;
          }

          cluster = 1;
          for (int iclust=low_ccc+1;iclust<high_ccc;iclust++){
            if (ccc_hits[iclust] == 0){
              cluster = 0;
              break;
            }
          }
        } // end if cluster mode

        if (!cluster){
          pair_count++;
          deltat_list.push_back(deltat);
        }
      }
    } // end loop over j pmts
  } // end loop over i pmts


  // now we have constructed our list of deltats.
  // if we have enough statistics, compute the median.
  // fail the event if median is out of range
  if (pair_count >= fCountThresh){
    int middle = deltat_list.size()/2;
    nth_element(deltat_list.begin(),deltat_list.begin()+middle,deltat_list.end());
    if (deltat_list.size()%2){
      nth_element(deltat_list.begin(),deltat_list.begin()+middle-1,deltat_list.end());
      dtmedian = (deltat_list[middle]+deltat_list[middle-1])/2.0;
    }else{
      dtmedian = deltat_list[middle];
    }
    if (dtmedian > fMedianCut){
      fPassFlag = false;
    }
  }
  cout << "NUMBER OF PAIRS: " << pair_count << endl;
  cout << "DELTAT MEDIAN(NS): " << dtmedian << endl;
  cout << "IS CLEAN?" << fPassFlag << endl;
}

int main(int argc, char** argv)
{
  const string& filename = string(argv[1]);
  const string material = string(argv[2]);
  // Read root files
  RAT::DU::DSReader dsReader( filename );
  cout << "PRINTING YOUR NHITS FOR FILE " << filename << endl;
  RAT::DU::DataCleaningBits rDataCleaningBits = RAT::DU::Utility::Get()->GetDataCleaningBits(); // To get the data cleaning bits

  vector<int> GTIDoIs = GTIDsOfInterest();
  // Loop through entries in rootfile
  for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ){
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );

    // Look at each triggered event in the entry
    for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ ){
      const RAT::DS::EV& rEV = rDS.GetEV( iEV );
      const RAT::DS::Meta& rMeta = dsReader.GetMeta();
      const UInt_t pass = rMeta.GetCurrentPass();
      int isGTIDMatch = 0;
      for(int j = 0; j < GTIDoIs.size(); j++){
        if(rEV.GetGTID() == GTIDoIs[j])
          isGTIDMatch = 1;
      }
      if (isGTIDMatch==0)
        continue;
      cout << "FOUND A GTID OF INTEREST: IS" << rEV.GetGTID() << endl;
      cout << "TOTAL NUM. CALPMTs IN EVENT: " << rEV.GetCalPMTs().GetCount();
      cout << "RUNNING FTS CUT..."  << endl;
      RunFTSCut(rEV);
      const RAT::DS::DataQCFlags& rDataQCFlags = rEV.GetDataCleaningFlags();
      const RAT::DS::BitMask& rBitMaskApplied = rDataQCFlags.GetApplied(pass);
      const RAT::DS::BitMask& rBitMaskFlags = rDataQCFlags.GetFlags(pass);
      for( map<size_t, string>::iterator iDCBits = rDataCleaningBits.GetInverseMapBegin();
        iDCBits != rDataCleaningBits.GetInverseMapEnd(); iDCBits++ ){
        Bool_t rApplied = rBitMaskApplied.Get( iDCBits->first );
        if( !rApplied ) // Checks cut has been applied
          continue;
        Bool_t rFlag = rBitMaskFlags.Get( iDCBits->first );
        if( rFlag )
          continue;
        string cutName = rDataCleaningBits.GetBitName( iDCBits->first ).c_str();
        cout << "FLAG APPLIED:" << cutName << endl;
      }
    }
  }
}
