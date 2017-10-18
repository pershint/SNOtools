#Main script looks for files from ../results and shows different results from
#Each bifurcation analysis result run.  Can also take in bifurcation analysis
#Class results and plot a fancy grid showing the relative correlations of each.


import numpy as np
import sys, os
import glob
import pylab as pl
import lib.playDarts as pd
import lib.resultgetter as rg
import lib.leakestimator as le
import lib.plots as plots


MAINDIR = os.path.dirname(__file__)
RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "..", "results"))


bifurcation_boxes = {"a": 9., "b": 1., "c": 1., "d": 83.}
bifurcation_boxes_SNO = {"a":369., "b":447., "c":15., "d":94264.}
acceptance_rates = {"DC": 0.9273, "DC_unc": 0.0007, "Fit":0.9984, \
        "Fit_unc": 0.0005}
acceptance_rates_SNO = {"DC": 0.9998, "DC_unc": 0.00005, "Fit": 0.995, "Fit_unc":0.0005}

if __name__ == '__main__':
    #For fun, let's grab all the results

    #Let's grab the list of files in the results directory
    BA = le.BifurAnalysisRun(bifurcation_boxes, acceptance_rates)
    print("y_dc:\n" + str(BA.y_dc()))
    print("y_fit:\n" + str(BA.y_fit()))
    print("uncertainty on y_dc:\n" + str(BA.y_dc_unc()))
    print("uncertainty on y_fit:\n" + str(BA.y_fit_unc()))
    print("Total events fed into B.A.: " + str(BA.total_events))
    print("Contamination assuming all events are background: ")
    print(BA.event_contamination())
    print("Contamination uncertainty: " + str(BA.event_contamination_unc()))

    #Show the boxes and Num. events in each box
    plots.BoxDistribution(BA)
    #Show the best fit y_dc and y_fit values given the bifurcation results
    plots.Contamination_Minimum(BA)

    #Show the histograms for y_dc and y_fit uncertainties with boostrapping
    plots.LeakageUncertainty_Bootstrap_Plot(BA)
