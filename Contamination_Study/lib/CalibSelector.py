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

def ApplyZCut(calibdir,source,zcutrange,filelist):
        all_calibdicts = LoadCalibrationDicts(calibdir,source)
        trimmed_filelist = []
        for f in filelist:
            for cdict in all_calibdicts:
                for run in cdict:
                    if run in f:
                        print(cdict[run]['position'][2])
                        print(zcutrange)
                        if cdict[run]['position'][2] < float(zcutrange[0]) and \
                                cdict[run]['position'][2] > float(zcutrange[1]):
                            trimmed_filelist.append(f)
        return trimmed_filelist 
