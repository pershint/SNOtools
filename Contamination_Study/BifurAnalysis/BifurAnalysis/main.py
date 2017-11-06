#Main script looks for files from ../results and shows different results from
#Each bifurcation analysis result run.  Can also take in bifurcation analysis
#Class results and plot a fancy grid showing the relative correlations of each.


import numpy as np
import sys, os
import glob
import json
import lib.playDarts as pd
import lib.resultgetter as rg
import lib.leakestimator as le
import lib.correlationestimator as ce
import lib.plots as plots

DEBUG = True

MAINDIR = os.path.dirname(__file__)
RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "..", "results", "OpenGolden_MayProcessed","4_9_MeV"))
ACCDIR = os.path.abspath(os.path.join(MAINDIR,"DB","Acceptance_Rates"))

acceptance_rates_SNO = {"DC": 0.9998, "DC_unc": 0.00005, "Fit": 0.995, "Fit_unc":0.0005}

if __name__ == '__main__':
    #Get the acceptance rate json entry you want
    with open(ACCDIR + "/N16Acceptance_Ewindow.json","r") as f:
        acceptances = json.load(f)
    #FIXME: Have some code that chooses the proper entry for the energy range
    acceptance_rates = acceptances["EnergyRanges"][0]
    #For fun, let's grab all the results
    allresult_filenames = glob.glob(RESULTDIR + "/*_results.out")
    if DEBUG:
        #Let's grab the list of files in the results directory
        test_result = rg.GetResultDict(allresult_filenames[0])
        BA = le.BifurAnalysisRun(test_result, acceptance_rates)
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
    
    #The following takes bifurcation analysis results for single cuts and
    #plots out the pearson coefficients in a grid fashion
    isSymmetric = True
    column_titles, row_titles, PC_rows = ce.buildTitlesAndPhiMatrix(allresult_filenames, acceptance_rates, isSymmetric)
    plots.CorrelationBoxes(PC_rows, column_titles, row_titles)
