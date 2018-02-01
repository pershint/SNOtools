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
DBDIR = os.path.abspath(os.path.join(MAINDIR, "DB"))
SACDIR = os.path.abspath(os.path.join(MAINDIR, "sacrifice_rootfiles", "old", "may2016_N16"))
#DCDIR = os.path.abspath(os.path.join(MAINDIR, "..", "N16_DCsacs", "4_9_MeV"))
#FITDIR = os.path.abspath(os.path.join(MAINDIR, "..", "N16_Fitsacs","4_9_MeV"))

######### VARIABLES ########
branch = "Fit"   #Either DC or Fit
db_entry = "N16_Positions_1.json"    #Point to which N16 run info you want
pcolor = 'b'
#FIXME: implement this!
date = "05/25/2017"  #Input a date to only look at those files. Type as None otherwise
######## /VARIABLES #######

def GetDate(rundict, date):
    #Returns the run entries associated with only a specific date
    datedict = {}
    for run in rundict:
        if rundict[run]["date"] == str(date):
            datedict[run] = rundict[run]
    return datedict

def plot_sacrificevsZ(calib_run_dict, sacrifice_filenames, branch):
    #Given a calibration run dictionary and the corresponding roots
    #output from either DC_Sacrifice or Beta14ITR_Sacrifice, plot
    #the physics sacrifice as a function of the z position of the source.
    positions = []
    source = None
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
                numall = float(h_all.GetEntries())
                numflagged= float(h_flagged.GetEntries())
                fractional_sac.append( numflagged / numall)
                fractional_sac_unc.append( np.sqrt((numflagged / (numall**2))+\
                        ((numflagged*np.sqrt(numall))/(numall**2))**2))
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
    return z_positions, fractional_sac, fractional_sac_unc

def weighted_avg_and_std(values, weights):
     """
     Returns the weighted average and standard deviation given values and
     Their weights.  The weights are assumed to be ((1/unc(value)**2))
     """
     average = np.average(values, weights=weights)
     variance = np.sqrt((np.sum(weights))/((np.sum(weights))**2))
     weighted_variance = 1./np.sqrt(np.sum(weights))
     tot_variance = np.sqrt((variance**2) + (weighted_variance**2))
     return (average, tot_variance)

if __name__ == '__main__':
    #Choose what branch of cuts you want to look at (DC or Fit)
    #First, get all our filenames.
    filenames = glob.glob(SACDIR + "/*")
    #Now, get our calibration dictionary.
    with open(DBDIR + "/"+db_entry,"r") as f:
        calib_dict = json.load(f)
    plot_dict = GetDate(calib_dict, "05/25/2017")
    #Try the plotting out
    z_pos, fs, fs_unc = plot_sacrificevsZ(plot_dict, filenames, branch)
    z_pos = np.array(z_pos)
    fs = np.array(fs)
    fs_unc = np.array(fs_unc)
    #Now, let's get the average and standard deviation
    Frac_average = np.average(fs)
    Frac_stdev = np.std(fs)
    print("AVERAGE FRAC. SACRIFICE: " + str(Frac_average))
    print("STDEV OF FRAC. SACRIFICE: " + str(Frac_stdev))
