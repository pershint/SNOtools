#TODO: Write a class that takes in a list of N16 files, and can return a
#Subset of N16 files based on date (From the calibration database), a 
#single run, or a set of run ranges.  
#Generalize this so it can be used for physics ntuple files too; we will
#Need to do the bifurcation analysis in different time bins as well
import glob
import json

def LoadCalibrationDicts(calibdir,source):
        all_calibdicts = []
        calibration_files = glob.glob(calibdir+"/"+source+"*.json")
        print(calibration_files)
        for cfile in calibration_files:
            with open(cfile,"r") as f:
                all_calibdicts.append(json.load(f))
        return all_calibdicts

def ApplyZCut(calibdir,source,zcutval,filelist):
        all_calibdicts = LoadCalibrationDicts(calibdir,source)
        trimmed_filelist = []
        for f in filelist:
            for cdict in all_calibdicts:
                for run in cdict:
                    if run in f:
                        if cdict[run]['position'][2] < zcutval:
                            trimmed_filelist.append(f)
        return trimmed_filelist
                        


#    def _checkEnergyRange(self):
#        if self.energy_range[0] < self.hist_erange[0]:
#            print("WARNING: Configuration file low energy range is outside"+\
#                    "range defined in sacrifice histograms.")
#        if self.energy_range[1] > self.hist_erange[1]:
#            print("WARNING: Configuration file high energy range higher"+\
#                    "than range defined in sacrifice histograms.")


#    def GetDate(self, date):
#        #load position dict with files only associated with given date
#        datedict = {}
#        for rundict in self.all_calib_positions:
#            for run in rundict:
#                if rundict[run]["date"] == str(date):
#                    datedict[run] = rundict[run]
#        self.positions_to_analyze = GetDate
#    
#    def GetRunRange(self,runrange):
#        if runrange[1] < runrange[0]:
#            print("order your run range correctly, come on...")
#            return None
#        self.runrange= runrange
#        rrdict={}
#        for rundict in self.all_calib_positions:
#            for run in rundict:
#                if runrange[0] <= int(run) <= runrange[1]:
#                    rrdict[run] = rundict[run]
#        print(rrdict)
#        self.positions_to_analyze = rrdict
#
#    def GetRun(runnum):
#        self.runrange=runnum
#        datedict={}
#        for rundict in self.all_calib_positions:
#            for run in rundict:
#                if int(run)==runnum:
#                    datedict[run] = rundict[run]
#        self.positions_to_analyze = datedict
#
