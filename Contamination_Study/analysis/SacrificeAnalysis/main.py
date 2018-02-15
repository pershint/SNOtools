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

MAINDIR = os.path.dirname(__file__)
CONFIGDIR = os.path.abspath(os.path.join(MAINDIR,"..","..","config"))
DBDIR = os.path.abspath(os.path.join(MAINDIR, "DB"))
SACDIR = os.path.abspath(os.path.join(MAINDIR, "sacrifice_rootfiles", "old", "may2016_N16"))
#DCDIR = os.path.abspath(os.path.join(MAINDIR, "..", "N16_DCsacs", "4_9_MeV"))
#FITDIR = os.path.abspath(os.path.join(MAINDIR, "..", "N16_Fitsacs","4_9_MeV"))

######### VARIABLES ########
energy_range = [5.5,9.0]
branch = "DC"   #Either DC or Fit
source = "N16"  #For plot labels
pcolor = 'b'
DATE = None  #Input a date to only look at those files. Type as None otherwise
RUN = None      #Calculate the sacrifice for this run only.  Type None if not in use
RUNRANGE = [100934,100951]  #Look at the sacrifice for all runs in this range.  Type
#none if not in use
######## /VARIABLES #######

def GetDate(rundicts, date):
    #Returns the run entries associated with only a specific date
    datedict = {}
    for rundict in rundicts:
        for run in rundict:
            if rundict[run]["date"] == str(date):
                datedict[run] = rundict[run]
    return datedict

def GetRunRange(rundicts,runrange):
    if runrange[1] < runrange[0]:
        print("order your run range correctly, come on...")
        return None
    rrdict={}
    for rundict in rundicts:
        for run in rundict:
            if runrange[0] <= int(run) <= runrange[1]:
                rrdict[run] = rundict[run]
    print(rrdict)
    return rrdict

def GetRun(rundicts,runnum):
    datedict={}
    for rundict in rundicts:
        for run in rundict:
            if int(run)==runnum:
                datedict[run] = rundict[run]
    return datedict


def calculate_sacrifice(calib_run_dict, sacrifice_filenames, branch):
    #Given a calibration run dictionary and the corresponding roots
    #output from either DC_Sacrifice or Beta14ITR_Sacrifice, return an array of
    #positions, sacrifice associated with each position, the uncertainty of that,
    #and the total number of events flagged and total number of events total
    positions = []
    source = None
    total_events = 0
    total_flagged_events = 0
    fractional_sac = []
    fractional_sac_unc = []
    for run in calib_run_dict:
        for filename in sacrifice_filenames:
            if run in filename:
                print(run)
                if source is None:
                    source = calib_run_dict[run]["source"]
                elif source != calib_run_dict[run]["source"]:
                    print("WARNING: not all sources are of the same...")
                #Found the filename corresponding to run in meta
                position = calib_run_dict[run]["position"]
                positions.append(position)

                #Now we get this run's sacrifice information
                _file0 = ROOT.TFile.Open(filename,"READ")
                h_all = copy.deepcopy(_file0.Get("h_AllEvents"))
                if branch == "DC":
                    h_flagged = copy.deepcopy(_file0.Get("h_DC_FlaggedEvents"))
                elif branch == "Fit":
                    h_flagged = copy.deepcopy(_file0.Get("h_BI_FlaggedEvents"))
                #FIXME: For the specified energy range, get only the entries from
                #Those bins
                numall = float(h_all.GetEntries())
                total_events = total_events + numall
                numflagged= float(h_flagged.GetEntries())
                total_flagged_events = total_flagged_events + numflagged
                fractional_sac.append( numflagged / numall)
                fractional_sac_unc.append( np.sqrt((np.sqrt(numflagged) / numall)**2+\
                        ((numflagged*np.sqrt(numall))/(numall**2))**2))
    fractional_sac = np.array(fractional_sac)
    fractional_sac_unc = np.array(fractional_sac_unc)
    return positions, fractional_sac, fractional_sac_unc, numall, numflagged

def plot_sacrificevsZ(positions, fractional_sac, fractional_sac_unc):
    z_positions = []
    for position in positions:
        z_positions.append(position[2])
    #We've got our position and sacrifice information for this calibration set.
    #Now, just plot it.
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(z_positions,fractional_sac, xerr=0, yerr=fractional_sac_unc, \
            marker='o', linestyle='none', color = pcolor, alpha=0.6, \
            label = 'Fractional Sacrifice')
    ax.set_xlabel(str(source) + " Z Position (m)")
    ax.set_ylabel("Fraction of events sacrificed")
    ax.set_title("Fractional sacrifice due to " + branch + " as source" + \
            " Z position varies\n" + \
            str(source) + " Source used")
    ax.grid(True)
    plt.show()

def plot_sacrificevsR(positions, fractional_sac, fractional_sac_unc):
    r_positions = []
    for position in positions:
        radius = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
        r_positions.append(radius)
    #We've got our position and sacrifice information for this calibration set.
    #Now, just plot it.
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(r_positions,fractional_sac, xerr=0, yerr=fractional_sac_unc, \
            marker='o', linestyle='none', color = pcolor, alpha=0.6, \
            label = 'Fractional Sacrifice')
    ax.set_xlabel(str(source) + " Radius (m)")
    ax.set_ylabel("Fraction of events sacrificed")
    ax.set_title("Fractional sacrifice due to " + branch + " as source" + \
            " radial position varies\n" + \
            str(source) + " Source used")
    ax.grid(True)
    plt.show()

if __name__ == '__main__':
    #Choose what branch of cuts you want to look at (DC or Fit)
    #First, get all our filenames.
    sacrifice_filenames = glob.glob(SACDIR + "/*")
    #Now, get our calibration dictionary.
    calibration_filenames = glob.glob(DBDIR+"/"+"N16*.json")
    calibration_dicts = []
    for cfile in calibration_filenames:
        with open(cfile,"r") as f:
            calibration_dicts.append(json.load(f))
    if DATE is not None:
        plot_dict = GetDate(calibration_dicts, "05/25/2017")
    if RUN is not None:
        plot_dict = GetRun(calibration_dicts, RUN)
    if RUNRANGE is not None:
        plot_dict = GetRunRange(calibration_dicts, RUNRANGE)
    #Try the plotting out
    positions, fs, fs_unc, total, flagged= calculate_sacrifice(plot_dict, sacrifice_filenames, branch)
    plot_sacrificevsR(positions,fs,fs_unc)
    #Now, let's get the average and standard deviation
    Frac_average = np.average(fs)
    Frac_stdev = np.std(fs)
    stat_uncertainty = np.sqrt((np.sqrt(flagged) / total)**2+\
                        ((flagged*np.sqrt(total))/(total**2))**2)
    print("TOTAL EVENTS: " + str(total))
    print("AVERAGE FRAC. SACRIFICE OF DATA SET USED: " + str(np.average(fs)))
    print("STDEV OF FRAC. SACRIFICE OF DATA SET USED " + str(np.std(fs)))
    print("STATISTICAL UNCERTAINTY OF TOTAL SACRIFICE OF DATA SET: " + str(stat_uncertainty))
