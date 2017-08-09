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

double avg1(vector<double> const& v)
{
  int n = 0;
  double mean = 0.0;
  for (int n = 0; n<v.size(); n++) 
  {
    double delta = v[n] - mean;
    mean += delta/(n+1);
  }
  return mean;
}

double pointdist(vector<double> const& v1, vector<double> const& v2)
{
  return sqrt(pow(v1[0]-v2[0],2) + pow(v1[1]-v2[1],2) + pow(v1[2]-v2[2],2));
}

vector<double> SCCut(const RAT::DS::EV& ev)
{
//  RAT::DBLinkPtr fLClean = RAT::DB::Get()->GetLink("DATACLEANING","sclustercut");
//  int fNHitMin = fLClean->GetI("nhit_min");
//  int fDeltaTThresh = fLClean->GetI("deltatthresh");
//  int fDistThresh = fLClean->GetI("distthresh");
//  int fCountThresh = fLClean->GetI("countthresh");
//  int fMaxPairs = 10000;

// THE DB IS BEING A PAIN.  JUST HARD CODE FOR NOW

  int fNHitMin = 20;
  int fDeltaTThresh = 200;
  int fDistThresh = 6000;
  int fCountThresh = 15;
  int fMaxPairs = 10000;
  Double_t fRval = 0.6;

  int cal_count,pair_count,cluster;
  int iccc,jccc;
  float tubedist, deltat, dtmedian;
  vector<const RAT::DS::PMTCal *> pmt_list;
  vector<float> deltat_list;
  vector<double> cluster_posn;
  bool fPassFlag = true;


  cal_count = ev.GetCalPMTs().GetCount();
  if (cal_count <= fNHitMin)
    return cluster_posn;

  // Get your list of hit PMT IDs
  for (int ipmt=0;ipmt<cal_count;ipmt++){
    const RAT::DS::PMTCal* pmt = &ev.GetCalPMTs().GetPMT(ipmt);
    // if time is calibrated
    // add to the list
    if (abs(pmt->GetTime()) < 8000)
      pmt_list.push_back(pmt);
  }


  vector<double> clustx;
  vector<double> clusty;
  vector<double> clustz;
  pair_count = 0;
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
      // pmts closer than threshold in mm

      // NOTE; Could just add pmt positions for each valid pair to the position
      // list; the more a PMT is valid in pairs, the more it weighs the avg.pos

      if ((deltat < fDeltaTThresh) and (tubedist < fDistThresh)){
        clustx.push_back(ipos.X());
        clustx.push_back(jpos.X());
        clusty.push_back(ipos.Y());
        clusty.push_back(jpos.Y());
        clustz.push_back(ipos.Z());
        clustz.push_back(jpos.Z());
        pair_count++;
      }
    } // end loop over j pmts
  } // end loop over i pmts


  // now we have constructed our number of clustered pairs.
  // if over 75% of PMT pairs match the criteria, the event is flagged
  if (pair_count >= fCountThresh){
    Double_t n2_full = (pow(cal_count,2) - cal_count)/2.;
    Double_t Valid34 = (fRval) * ((fRval*(double)cal_count) - 1.)/((double)cal_count -1.);
   cout << "VALID35 VALUE/REQUIREMENT FOR " << ev.GetGTID() << endl;
   cout << (pair_count/n2_full) << "/" << Valid34 << endl;
   if ((pair_count/n2_full) > Valid34)
    {
      //Our flag has been triggered!  Return the average position of valid cluster ids
      cout << "AAAAAND TRIGGERED FOR " << ev.GetGTID() << endl;
      cluster_posn.push_back(avg1(clustx));
      cluster_posn.push_back(avg1(clusty));
      cluster_posn.push_back(avg1(clustz));
      fPassFlag = false;
    }
  }
//  cout << "NUMBER OF PAIRS: " << pair_count << endl;
//  cout << "IS CLEAN?" << fPassFlag << endl;
  return cluster_posn;  //fPassFlag;
}





int main(int argc, char** argv)
{
  //Initialize interevent checks
  RAT::DBLinkPtr fLClean = RAT::DB::Get()->GetLink("DATACLEANING","sclustercut");
//  fEventDt = fLClean->GetI("dtevent");
  int fEventDt = 1000000000;   //0.5 sec in time
  int fNumclust = 3; //Number of clusters that must be w/in fEventdt
  int fClustDist = 3000; //Distance the cluster pmt avg posns. must be within, in mm
  cout << "YEEEEP" << endl;
  const string& filename = string(argv[1]);
  const string material = string(argv[2]);
  // Read root files
  cout << "NOT IN FILE" << filename << endl;
  RAT::DU::DSReader dsReader( filename );
  RAT::DU::DataCleaningBits rDataCleaningBits = RAT::DU::Utility::Get()->GetDataCleaningBits(); // To get the data cleaning bits

  vector<int> Event_GTIDs;  //GTIDs meeting sclustercut criteria
  vector<ULong64_t> Event_UnivNSTime; //Universal times of events meeting scluster criteria
  vector< vector<double> > Event_Spot;
  // Loop through entries in rootfile
  for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ ){
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );

    // Look at each triggered event in the entry
    for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ ){
      const RAT::DS::EV& rEV = rDS.GetEV( iEV );
      const RAT::DS::Meta& rMeta = dsReader.GetMeta();
      const UInt_t pass = rMeta.GetCurrentPass();
      bool PassesCut = true;
      vector<double> clusterspot;
      clusterspot = SCCut(rEV);
      if(!clusterspot.empty()) {
        Event_GTIDs.push_back(rEV.GetGTID());
        Event_UnivNSTime.push_back(rEV.GetClockCount10() * 100);
        Event_Spot.push_back(clusterspot);
//        Event_10MHzCLock.push_back(rEV.GetUniversalTime().GetSeconds());
//      Event_Utimes.push_back(rEV.universalTime.GetNanoSeconds());
      }
    }
  }
  // Now, use the pmt_weighted position and interevent times
  // to determine if events are at the same rough place and
  // occur in a burst
  for (int i=0; i<Event_GTIDs.size(); i++)
  {
    cout << "CLOCK: " << Event_UnivNSTime[i] << endl;
    cout << "GTID: " << Event_GTIDs[i] << endl;
    cout << "CLUST_SPOT" << Event_Spot[i][0] << " " << Event_Spot[i][1] <<
        " " << Event_Spot[i][2] << endl;
  }

  cout << "NOW, APPLY SECOND PASS CUTS" << endl;
  vector<int> Event_GTIDs_failSP;
  for (int j=2; j<Event_GTIDs.size(); j++)
  {
    if ((Event_UnivNSTime[j] - Event_UnivNSTime[j-2]) < fEventDt)
    {
      //check the interdistance of all three of the events
      bool allclose = true;
      for (int l=(j-2); l<(j+1); l++)
      {
        for (int k=(l+1); k<(j+1); k++)
        {
          if(l==k)
            continue;
          double distance = pointdist(Event_Spot[l],Event_Spot[k]);
          if (distance > fClustDist)
            allclose = false;
        }
      }

      // if they're all withing the required distance, events are flagged
      if(allclose)
      {
        //Not quite right; getting pointers in the array, not the values themselves
        Event_GTIDs_failSP.push_back(Event_GTIDs[j]);
        Event_GTIDs_failSP.push_back(Event_GTIDs[j-1]);
        Event_GTIDs_failSP.push_back(Event_GTIDs[j-2]);
        cout << "GTIDFAILSP: " << Event_GTIDs[j] << " " <<
            Event_GTIDs[j-1] << " " << Event_GTIDs[j-2] << endl;
      }
    }
  }
}
