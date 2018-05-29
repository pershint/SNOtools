#This code reads out N16 metadata from the DB directory, and will
#produce sacrifice plots that show either:
#  a) How the fit/DC sacrifice varies with position
#  b_ How the fit/DC sacrifice varies over time
# To use:
#  1) Add a database entry for the calibration run you are using
#     Follow the format of "N16_Positions_1.json" to do this
#  2) Make sure you have your root histograms output from running processed
#     data in DC_Sacrifice and Beta14ITR_Sacrifice in the directories
#     Defined below in DCDIR and FITDIR
#  3) Modify the variables below in the VARIABLES section
#  4) Let it rip


import os,sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import json
import ROOT
import glob

class SacrificeHistAnalyzer(object):
    def __init__(self, Sacrifice_Histograms=None, config_dict = {},source='N/A'):
        self.source = Sacrifice_Histograms.source
        self.hist_rootfiles = Sacrifice_Histograms.histogram_files  #All sacrifice root files
        self.hist_runnums = Sacrifice_Histograms.histogram_runnums
        self.source=source
        self.cdict = config_dict
        self.calib_positions = {}
        self.sacrifice_byrun = {}
        self.sacrifice_summary = {}
        self.sacrifice_summary['sourcetype'] = self.source
        self.sacrifice_byrun['sourcetype'] = self.source
        if self.source == "MC":
            comment = ('Monte carlos pathological cuts are only defined by'
            'the energy, radius, and zcuts.  No trigger cuts or pathological'
            'data cleaning cuts included (assumed to not be pathological)')
            self.sacrifice_summary['comment'] = comment
            self.sacrifice_byrun['comment'] = comment
        self._Initialize_SacrificeSummaries()

    def ShowSacrificeResults(self):
        print("SACRIFICE AND UNCERTAINTIES FOR CUT BRANCH 1")
        print(self.sacrifice_summary['cut1'])
        print("SACRIFICE AND UNCERTAINTIES FOR CUT BRANCH 2")
        print(self.sacrifice_summary['cut2'])
    
    def _Initialize_SacrificeSummaries(self):
        for cut in ['cut1','cut2']:
            self.sacrifice_byrun[cut] = {}
            self.sacrifice_summary[cut] = {}
            self.sacrifice_summary[cut]["total_events"] = 0
            self.sacrifice_summary[cut]["nonpath_events"] = 0
            self.sacrifice_summary[cut]["total_evtsflagged"] = 0
            self.sacrifice_summary[cut]["total_fracsac"] = 0.0
            self.sacrifice_summary[cut]["total_fracsac_unc"] = 0.0
            self.sacrifice_byrun[cut]["position"] = []
            self.sacrifice_byrun[cut]["run"] = []
            self.sacrifice_byrun[cut]["total_events"] = []
            self.sacrifice_byrun[cut]["events_flagged"] = []
            self.sacrifice_byrun[cut]["fractional_sac"] = []
            self.sacrifice_byrun[cut]["fractional_sac_unc"] = []
            self.sacrifice_byrun[cut]["nonpath_events"] = []

    def LoadCalibrationPositions(self,db_directory):
        all_cdicts = []
        calibration_filenames = glob.glob(db_directory+"/"+self.source+"*.json")
        for cfile in calibration_filenames:
            with open(cfile,"r") as f:
                all_cdicts.append(json.load(f))
        for cdict in all_cdicts:
            for run in cdict:
                self.calib_positions[run] = cdict[run]

    def CalculateSacrifices(self):
        '''Clears current sacrifice summary dictionaries, and runs the by-run
        and total summary calculator functions'''
        self._Initialize_SacrificeSummaries()
        self.Calculate_SacrificesByRun()
        self.Calculate_SacrificeSummary()

    def Calculate_SacrificeSummary(self):
        '''Summarizes information found after running Calculate_SacrificesByRun
        into one master dictionary'''

        for cut in ['cut1','cut2']:
            #Add fractional sacrifice estimates for all runs investigated
            total = float(self.sacrifice_summary[cut]["nonpath_events"])
            totalflagged = float(self.sacrifice_summary[cut]["total_evtsflagged"])
            if total <= 0:
                print("No nonpathological events in any run analyzed")
                self.sacrifice_summary[cut]["total_fracsac"] = -1
                self.sacrifice_summary[cut]["total_fracsac_unc"] = -1
                continue
            self.sacrifice_summary[cut]["total_fracsac"] = totalflagged/total
            self.sacrifice_summary[cut]["total_fracsac_unc"] = np.sqrt((np.sqrt(totalflagged) /\
                                total)**2+((totalflagged*np.sqrt(total))/(total**2))**2)

    def Calculate_SacrificesByRun(self):
        '''populates the self.sacrifice_byrun dictionary with run-by-run
        sacrifice information for each data cleaning cut branch'''
        for j,filename in enumerate(self.hist_rootfiles):
            thisrun = str(self.hist_runnums[j])
            position = None
            for run in self.calib_positions:
                if run == thisrun:
                    #Found the filename corresponding to run in meta
                    position = self.calib_positions[run]["position"]
                    break
            if position is None:
                position = "N/A"
            #Now we get this run's sacrifice information
            _file0 = ROOT.TFile.Open(filename,"READ")
            h_allnonpath = copy.deepcopy(_file0.Get("h_AllNonpathEvents"))
            h_all = copy.deepcopy(_file0.Get("h_AllEvents"))
            numall = float(h_all.GetEntries())
            numnonpath = float(h_allnonpath.GetEntries())
            for cut in ['cut1','cut2']:
                h_cutflagged = copy.deepcopy(_file0.Get("h_"+cut+"_FlaggedEvents"))
                cut_numflagged= float(h_cutflagged.GetEntries())
                self.sacrifice_summary[cut]["total_events"]+=numall
                self.sacrifice_summary[cut]["nonpath_events"]+=numnonpath
                self.sacrifice_summary[cut]["total_evtsflagged"]+=cut_numflagged
                if thisrun in self.sacrifice_byrun[cut]["run"]:
                    #Have a subrun; add into the previous run's contents
                    #Get the index where self.sacrifice_byrun[cut]["run"] == run
                    currentruns = np.array(self.sacrifice_byrun[cut]["run"])
                    runindex = np.where(currentruns == str(thisrun))[0]
                    runindex = runindex[0]
                    self.sacrifice_byrun[cut]["total_events"][runindex]+=numall
                    self.sacrifice_byrun[cut]["nonpath_events"][runindex]+=numnonpath
                    self.sacrifice_byrun[cut]["events_flagged"][runindex]+=\
                            cut_numflagged
                else:
                    #new run; append all the contents here
                    self.sacrifice_byrun[cut]["run"].append(thisrun)
                    self.sacrifice_byrun[cut]["position"].append(position)
                    self.sacrifice_byrun[cut]["total_events"].append(numall)
                    self.sacrifice_byrun[cut]["nonpath_events"].append(numnonpath)
                    self.sacrifice_byrun[cut]["events_flagged"].append(cut_numflagged)
        for j,run in enumerate(self.sacrifice_byrun[cut]["run"]):
            numflagged = self.sacrifice_byrun[cut]["events_flagged"][j]
            numnonpath = self.sacrifice_byrun[cut]["nonpath_events"][j]
            if numnonpath == 0:
                self.sacrifice_byrun[cut]["fractional_sac_unc"].append(-1)
                self.sacrifice_byrun[cut]["fractional_sac"].append(-1)
                continue
            cut_fractional_sac = ( numflagged / numnonpath)
            cut_frac_sac_unc = np.sqrt((np.sqrt(numflagged) /\
                    numnonpath)**2+((numflagged*np.sqrt(numnonpath))/(numnonpath**2))**2)
            self.sacrifice_byrun[cut]["fractional_sac_unc"].append(cut_frac_sac_unc)
            self.sacrifice_byrun[cut]["fractional_sac"].append(cut_fractional_sac)
    
    def SaveSacrificeSummary(self,savedir,savename):
        savesacloc = savedir+"/"+savename
        with open(savesacloc,"w") as f:
            json.dump(self.sacrifice_summary, f, sort_keys=True,indent=4)

    def SaveSacrificeByRun(self,savedir,savename):
        savesacloc = savedir+"/"+savename
        with open(savesacloc,"w") as f:
            json.dump(self.sacrifice_byrun, f, sort_keys=True,indent=4)

