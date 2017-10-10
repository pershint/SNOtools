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
DCDIR = os.path.abspath(os.path.join(MAINDIR, "..", "N16_DCsacs"))
FITDIR = os.path.abspath(os.path.join(MAINDIR, "..", "N16_Fitsacs"))

######### VARIABLES ########
branch = "Fit"
SACDIR = FITDIR
db_entry = "N16_Positions_1.json"
pcolor = 'b'
######## /VARIABLES #######

def plot_sacrificevsZ(calib_run_dict, sacrifice_filenames, branch):
    #Given a calibration run dictionary and the corresponding roots
    #output from either DC_Sacrifice or Beta14ITR_Sacrifice, plot
    #the physics sacrifice as a function of the z position of the source.
    z_positions = []
    fractional_sac = []
    fractional_sac_unc = []
    for run in calib_run_dict["Runs"]:
        for filename in sacrifice_filenames:
            if run in filename:
                #Found the filename corresponding to run in meta
                zposition = calib_run_dict["Runs"][run][2]
                z_positions.append(zposition)

                #Now we get this run's sacrifice information
                _file0 = ROOT.TFile.Open(filename,"READ")
                if branch == "DC":
                    h_all = copy.deepcopy(_file0.Get("h_AllEvents"))
                    h_flagged = copy.deepcopy(_file0.Get("h_FlaggedEvents"))
                elif branch == "Fit":
                    h_all = copy.deepcopy(_file0.Get("h_fit_E"))
                    h_flagged = copy.deepcopy(_file0.Get("h_fit_fail_E"))
                numall = float(h_all.GetEntries())
                numflagged= float(h_flagged.GetEntries())
                fractional_sac.append( numflagged / numall)
                fractional_sac_unc.append( np.sqrt((numflagged / (numall**2))+\
                        ((numflagged*np.sqrt(numall))/(numall**2))**2))
    #We've got our position and sacrifice information for this calibration set.
    #Now, just plot it.
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(z_positions,fractional_sac, xerr=0, yerr=fractional_sac_unc, \
            marker='o', linestyle='none', color = pcolor, alpha=0.6, \
            label = 'Fractional Sacrifice')
    ax.set_xlabel(calib_run_dict["Source"] + " Z Position (m)")
    ax.set_ylabel("Fraction of events sacrificed")
    ax.set_title("Fractional sacrifice due to " + branch + " as source" + \
            " Z position varies\n" + \
            "Calibration on " + calib_run_dict["Date"] + ", " + \
            calib_run_dict["Source"] + " Source used")
    ax.grid(True)
    plt.show()
    return z_positions, fractional_sac, fractional_sac_unc

#Nabbed shamelessly from stack overvlow
def weighted_avg_and_std(values, weights):
     """
     Return the weighted average and standard deviation.
     values, weights -- Numpy ndarrays with the same shape.
     """
     average = np.average(values, weights=weights)
     variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
     return (average, np.sqrt(variance))

if __name__ == '__main__':
    #Choose what branch of cuts you want to look at (DC or Fit)
    #First, get all our filenames.
    filenames = glob.glob(SACDIR + "/*")
    #Now, get our calibration dictionary.
    with open(DBDIR + "/"+db_entry,"r") as f:
        calib_dict = json.load(f)

    #Try the plotting out
    zp, fs, fs_unc = plot_sacrificevsZ(calib_dict, filenames, branch)
    zp = np.array(zp)
    fs = np.array(fs)
    fs_unc = np.array(fs_unc)
    #Now, let's get the average and standard deviation
    Frac_average = np.average(fs)
    Frac_stdev = np.std(fs)
    print("AVERAGE FRAC. SACRIFICE: " + str(Frac_average))
    print("STDEV OF FRAC. SACRIFICE: " + str(Frac_stdev))

    #Try a weighted average.  Weigh each fs term based on the uncertainties
    weights = 1. / fs_unc
    weight_av, weight_std = weighted_avg_and_std(fs, weights)
    print("WEIGHTED AVG. FRAC. SACRIFICE: " + str(weight_av))
    print("WEIGHTED STD. OF AVG: " + str(weight_std))
