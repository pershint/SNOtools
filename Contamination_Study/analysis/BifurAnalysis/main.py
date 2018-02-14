#Main script looks for files from ../results and shows different results from
#Each bifurcation analysis result run.  Can also take in bifurcation analysis
#Class results and plot a fancy grid showing the relative correlations of each.

import argparse
import numpy as np
import sys, os
import glob
import json
import lib.playDarts as pd
import lib.resultgetter as rg
import lib.leakestimator as le
import lib.correlationestimator as ce
import lib.plots as plots
import matplotlib.pyplot as plt
import numpy as np

MAINDIR = os.path.dirname(__file__)
RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "..", "..","oldresults","results", "OpenGolden_MayProcessed"))
ERANGE = "5p5_9_MeV"
FILENAME = "BA_5p5to9MeV_results.out"
ACCDIR = os.path.abspath(os.path.join(MAINDIR,"DB","Acceptance_Rates"))

parser = argparse.ArgumentParser(description='Parser to decide what analysis to do')
parser.add_argument('--debug', dest='debug',action='store_true',
        help='Run code in debug mode')
parser.add_argument('--correlations', dest='CORRELATIONS',action='store_true',
        help='Run the code that plots the correlations of different cuts/classifiers')
parser.add_argument('--vsenergy', dest='CONTVSEN', action='store_true',
        help='Plot the contamination as a function of energy for the energy windows analyzed')
parser.set_defaults(CORRELATIONS=False,CONTVSEN=False,debug=False)
args = parser.parse_args()

CORRELATIONS = args.CORRELATIONS
CONTVSEN = args.CONTVSEN
DEBUG = args.debug

print("CORRELATIONS BOOL: " + str(CORRELATIONS))
print("CONTVSEN BOOL: " + str(CONTVSEN))

if __name__ == '__main__':
    #Get the acceptance rate json entry you want
    with open(ACCDIR + "/N16Acceptance_Ewindow.json","r") as f:
        acceptances = json.load(f)
    #FIXME: Have some code that chooses the proper entry for the energy range
    acceptance_rates = acceptances["EnergyRanges"][0]
    #For fun, let's grab all the results
    allresult_filenames = glob.glob(RESULTDIR+'/'+ERANGE+'/'+FILENAME)
    if DEBUG:
        #Let's grab the first file in the results directory
        test_result = rg.GetResultDict(allresult_filenames[0])
        print(test_result)
        BA = le.BifurAnalysisRun(test_result, acceptance_rates)
        print("y_dc:\n" + str(BA.y_dc()))
        print("y_fit:\n" + str(BA.y_fit()))
        print("Total events fed into B.A.: " + str(BA.total_events))
        print("Contamination assuming all events are background: ")
        print(BA.event_contamination())
        #Use the bootstrap method to find the 90% CL on y1y2beta
        y1y2beta_distribution, Upper_90CL = BA.y1y2beta_bootstrap_unc(0.90,100000)
        print("90% CL ON UPPER BOUND OF CONTAMINATION: " + str(Upper_90CL))
   
        #Show the boxes and Num. events in each box
        plots.BoxDistribution(BA)
        
    
        #Show the histograms for y_dc and y_fit uncertainties with boostrapping
        plots.LeakageUncertainty_Bootstrap_Plot(BA)
   
    ############ CORRELATIONS FOR BIFURCATION ANALYSIS###########
    #The following takes bifurcation analysis results for single cuts and
    #plots out the pearson coefficients in a grid fashion
    isSymmetric = True
    #Cut matrices come out ordered alphabetically in PC and cov matrices
    #TODO: We now need code that will build up distributions of variables
    #From B14, ITR, and data cleaning cuts.  The Pearson Coefficients and
    #variances will be fed in below.
    column_titles, row_titles, PC_matrix, cov_matrix = ce.buildTitlesAndPhiCorrMatrices(allresult_filenames, acceptance_rates, isSymmetric)
    plots.PearsonCoeffMatrix(PC_matrix, column_titles, row_titles)
    plots.CovarianceMatrix(cov_matrix, column_titles, row_titles)

    variables = []
    varibvariances = []
    CovMatrix = ce.CovarianceMatrix(cov_matrix, variables, varibvariances)
    CovMatrix.choleskydecompose()
    n = 94 #FIXME: have this filled with the number from the bifurcation
           #analysis result
    cutmap = prob_dict.keys()
    BoxMaker = ce.BA_Correlations(cutmap,cutmap[2:len(cutmap)],cutmap[0:2]) 
    all_passfails = []
    i = 0
    while i < 500:
        fired_variable_sets = CovMatrix.shoot_corrcuts(n)
        #TODO: need a function that takes in the variables, the cut
        #Labels, and gets the cut thresholds from RATDB. Checks if the
        #Cut passes or fails.
        #passfail_flags = check_passfail()
        #BoxMaker.fillBABoxes(passfail_flags)
        i+=1
           
    ######### END CORRELATIONS CODE #############

    ###########
    #The following is hacky, but will get the job done.  We
    #Load in each of the energy window results and append the
    #Contaminations into arrays.  Using subtraction, we'll Get
    #The expected contamination you get from each additional
    #Step in energy width.
    ###########
    contaminations = []
    contamination_uncs = []
    Result_dirs = ["4_9_MeV","4p5_9_MeV","5_9_MeV","5p5_9_MeV"]
    bin_centers = [4.25, 4.75, 5.25, 7.25]
    bin_widths = [0.25, 0.25, 0.25, 1.75]
    for direc in Result_dirs:
        erange_results = glob.glob(RESULTDIR+'/'+direc+ "/*_results.out")
        erange_result_dict = rg.GetResultDict(erange_results[0])
        BA = le.BifurAnalysisRun(erange_result_dict, acceptance_rates)
        contaminations.append(BA.event_contamination())
        contamination_uncs.append(BA.event_contamination_unc())
    print(contaminations)
    print(contamination_uncs)
    contam_binned = []
    contam_binned_uncs = []
    for j, entry in enumerate(contaminations):
        if (j > 0):
            contam_binned.append(contaminations[j-1] - contaminations[j])
            contam_binned_uncs.append(np.sqrt((contamination_uncs[j-1]**2) + \
                    (contamination_uncs[j]**2)))
    contam_binned.append(contaminations[len(contaminations)-1])
    contam_binned_uncs.append(contamination_uncs[len(contamination_uncs)-1])
    #Dear god move this to a library please
    import matplotlib.pyplot as plt
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(bin_centers, contam_binned, xerr=bin_widths, \
            yerr=contam_binned_uncs, marker = 'o', linestyle='none',\
            color='r', alpha=0.8)
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("# Instrumentals in region")
    ax.set_title("Estimated # Non-cherenkov events that leak through \n"+\
            "Data cleaning and Fits in 11 days of data taking")
    ax.grid(True)
    plt.show()
