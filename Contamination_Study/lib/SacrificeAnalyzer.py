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
    def __init__(self, Sacrifice_Histograms=None, config_dict = {},source="N/A"):
        self.hist_rootfiles = Sacrifice_Histograms.histogram_files  #All sacrifice root files
        #FIXME: Eventually store this metadata in the histogram files as a tree?
        self.source = "N/A"
        self.cdict = config_dict
        self.calib_positions = {}
        self.cut_sacrifices = {"sourcetype": self.source}
        self.sacrifice_summary = {"sourcetype": self.source}

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
        #Given a calibration run dictionary and the corresponding roots
        #output from either DC_Sacrifice or Beta14ITR_Sacrifice, return an array of
        #positions, sacrifice associated with each position, the uncertainty of that,
        #and the total number of events flagged and total number of events total
        for cut in ['cut1','cut2']:
            self.cut_sacrifices[cut] = {}
            self.sacrifice_summary[cut] = {}
            self.sacrifice_summary[cut]["total_events"] = 0
            self.sacrifice_summary[cut]["nonpath_events"] = 0
            self.sacrifice_summary[cut]["total_evtsflagged"] = 0
            self.sacrifice_summary[cut]["total_fracsac"] = 0.0
            self.sacrifice_summary[cut]["total_fracsac_unc"] = 0.0
            self.cut_sacrifices[cut]["position"] = []
            self.cut_sacrifices[cut]["run"] = []
            self.cut_sacrifices[cut]["total_events"] = []
            self.cut_sacrifices[cut]["events_flagged"] = []
            self.cut_sacrifices[cut]["fractional_sac"] = []
            self.cut_sacrifices[cut]["fractional_sac_unc"] = []
            self.cut_sacrifices[cut]["nonpath_events"] = []
        for run in self.calib_positions:
            for filename in self.hist_rootfiles:
                if run in filename:
                    #Found the filename corresponding to run in meta
                    position = self.calib_positions[run]["position"]
                    #Now we get this run's sacrifice information
                    _file0 = ROOT.TFile.Open(filename,"READ")
                    h_allnonpath = copy.deepcopy(_file0.Get("h_AllNonpathEvents"))
                    h_all = copy.deepcopy(_file0.Get("h_AllEvents"))
                    numall = float(h_all.GetEntries())
                    numnonpath = float(h_allnonpath.GetEntries())
                    if numnonpath == 0:
                        print("RUN %s HAS NO NON-PATHOLOGICAL EVENTS." % (_file0.GetName()))
                    for cut in ['cut1','cut2']:
                        h_cutflagged = copy.deepcopy(_file0.Get("h_"+cut+"_FlaggedEvents"))
                        cut_numflagged= float(h_cutflagged.GetEntries())
                        self.sacrifice_summary[cut]["total_events"]+=numall
                        self.sacrifice_summary[cut]["nonpath_events"]+=numnonpath
                        self.sacrifice_summary[cut]["total_evtsflagged"]+=cut_numflagged
                        if int(run) in self.cut_sacrifices[cut]["run"]:
                            #Have a subrun; add into the previous run's contents
                            #Get the index where self.cut_sacrifices[cut]["run"] == run
                            currentruns = np.array(self.cut_sacrifices[cut]["run"])
                            runindex = np.where(currentruns == int(run))[0]
                            runindex = runindex[0]
                            self.cut_sacrifices[cut]["total_events"][runindex]+=numall
                            self.cut_sacrifices[cut]["nonpath_events"][runindex]+=numnonpath
                            self.cut_sacrifices[cut]["events_flagged"][runindex]+=\
                                    cut_numflagged
                        else:
                            #new run; append all the contents here
                            self.cut_sacrifices[cut]["run"].append(int(run))
                            self.cut_sacrifices[cut]["position"].append(position)
                            self.cut_sacrifices[cut]["total_events"].append(numall)
                            self.cut_sacrifices[cut]["nonpath_events"].append(numnonpath)
                            self.cut_sacrifices[cut]["events_flagged"].append(cut_numflagged)

        #Calculate the total fractional acceptance and uncertainty, and
        #run-by-run sacrifice now
        print(self.cut_sacrifices)
        for cut in ['cut1','cut2']:
            for j,run in enumerate(self.cut_sacrifices[cut]["run"]):
                numflagged = self.cut_sacrifices[cut]["events_flagged"][j]
                numnonpath = self.cut_sacrifices[cut]["nonpath_events"][j]
                cut_fractional_sac = ( numflagged / numnonpath)
                cut_frac_sac_unc = np.sqrt((np.sqrt(numflagged) /\
                        numnonpath)**2+((numflagged*np.sqrt(numnonpath))/(numnonpath**2))**2)
                self.cut_sacrifices[cut]["fractional_sac_unc"].append(cut_frac_sac_unc)
                self.cut_sacrifices[cut]["fractional_sac"].append(cut_fractional_sac)
            #Add fractional sacrifice estimates for all runs investigated
            total = float(self.sacrifice_summary[cut]["nonpath_events"])
            totalflagged = float(self.sacrifice_summary[cut]["total_evtsflagged"])
            self.sacrifice_summary[cut]["total_fracsac"] = totalflagged/total
            self.sacrifice_summary[cut]["total_fracsac_unc"] = np.sqrt((np.sqrt(totalflagged) /\
                                total)**2+((totalflagged*np.sqrt(total))/(total**2))**2)


    def ShowSacrificeResults(self):
        print("SACRIFICE AND UNCERTAINTIES FOR CUT BRANCH 1")
        print(self.sacrifice_summary['cut1'])
        print("SACRIFICE AND UNCERTAINTIES FOR CUT BRANCH 2")
        print(self.sacrifice_summary['cut2'])

    def SaveSacrificeSummary(self,savedir,savename):
        savesacloc = savedir+"/"+savename
        with open(savesacloc,"w") as f:
            json.dump(self.sacrifice_summary, f, sort_keys=True,indent=4)

    def SaveSacrificeByRun(self,savedir,savename):
        savesacloc = savedir+"/"+savename
        with open(savesacloc,"w") as f:
            json.dump(self.cut_sacrifices, f, sort_keys=True,indent=4)
