#Some functions that can be used to Load the calibration dictionaries
#Defined in ../DB/ and remove runs from a filelist based on the run number
#zposition associated with the run

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
                        if cdict[run]['position'][2] < float(zcutrange[0]) and \
                                cdict[run]['position'][2] > float(zcutrange[1]):
                            trimmed_filelist.append(f)
        return trimmed_filelist 
