#Main script looks for files from ../results and shows different results from
#Each bifurcation analysis result run.  Can also take in bifurcation analysis
#Class results and plot a fancy grid showing the relative correlations of each.


import numpy as np
import sys, os
import glob
import pylab as pl
import lib.playDarts as pd
import lib.resultgetter as rg
import lib.maskbuilder as mb
import lib.leakestimator as le
import lib.plots as plots

DEBUG = False

MAINDIR = os.path.dirname(__file__)
RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "..", "results", "OpenGolden_MayProcessed"))


bifurcation_boxes = {"a": 9., "b": 1., "c": 1., "d": 83.}
bifurcation_boxes_SNO = {"a":369., "b":447., "c":15., "d":94264.}
#acceptance_rates = {"DC": 0.9273, "DC_unc": 0.0007, "Fit":0.9984, \
#        "Fit_unc": 0.0005, "E_range":[5.5,9.0]}
acceptance_rates = {"DC": 0.93031, "DC_unc": 0.0004, "Fit":0.9948, \
        "Fit_unc": 0.0001, "E_range":[4.0,9.0]}
cceptance_rates_SNO = {"DC": 0.9998, "DC_unc": 0.00005, "Fit": 0.995, "Fit_unc":0.0005}

def getPhiRows(cut1_list, cut2_list, PC_list, isSymmetric):
    '''
    Takes in all the phi coefficients from each bifurcation analysis and
    returns what titles go on the x-axis, the y-axis, and the row structure
    of the phi coefficients that will go with them.
    '''
    column_titles = list(set(cut1_list))
    column_titles.sort()
    row_titles = list(set(cut2_list))
    row_titles.sort()
    if isSymmetric:
        #Column and Row titles will be the same, since if you have one set
        #of column/row pairs, you have the mirrored row/column pair value
        column_titles = list(set(column_titles + row_titles))
        column_titles.sort()
        row_titles = list(set(column_titles + row_titles))
        row_titles.sort()
    PC_rows = []
    for y in row_titles:
        row=[]
        for x in column_titles:
            cut1_indices = [i for i,j in enumerate(cut1_list) if j==x]
            cut2_indices = [i for i,j in enumerate(cut2_list) if j==y]
            #Now, get the common index; that is, the index that matches
            #This title pair in the PC_list
            try:
                PC_index = list(set(cut1_indices) & set(cut2_indices))[0]
                row.append(PC_list[PC_index])
            except IndexError:
                #No match in data set.
                #Phi coefficient is symmetric in row/column. Try to get the
                #value from swapping cut1 and cut2
                if isSymmetric:
                    cut1_indices = [i for i,j in enumerate(cut1_list) if j==y]
                    cut2_indices = [i for i,j in enumerate(cut2_list) if j==x]
                    #Now, get the common index; that is, the index that matches
                    #This title pair in the PC_list
                    try:
                        PC_index = list(set(cut1_indices) & set(cut2_indices))[0]
                        row.append(PC_list[PC_index])
                    except IndexError:
                        #No match in mirrored set. set this value as -2
                        row.append(-2)
                else:
                    #No match in mirrored set. set this value as -2
                    row.append(-2)
        PC_rows.append(row)
    return column_titles, row_titles, PC_rows


def getTitles(result_dict):
    '''
    Takes in the results from a bifurcation analysis run and returns what
    the title for cut1 and cut2 should be.  The order for what cut1 and
    cut2 is is dependent on how the analysis is done in ../BifurAnalysis,
    so if you change the analysis/result output logic, CHANGE IT HERE TOO!
    '''
    print(result_dict)
    DCmask_dict = mb.get_dcwords()
    if result_dict["DCandFit"] == 1:
        #There's both data cleaning and fits used
        if result_dict["ITR_only"]==1:
            #Cut1 is a DC mask and Cut2 is ITR
            cut2_title = "ITR"
            cut1_title = DCmask_dict[mb.int_to_binary_bit(int(result_dict["DC_mask1"]))]
        elif result_dict["B14_only"]==1:
            #Cut1 is a DC mask and Cut2 is B14
            cut2_title = "Beta14"
            cut1_title = DCmask_dict[mb.int_to_binary_bit(int(result_dict["DC_mask1"]))]
        else:
            #Cut1 is all DC and Cut2 is both Fitters
            cut1_title = "All DCs"
            cut2_title = "All Fits"
    else:
        #It's only either two DC cuts or the fitters
        if result_dict["DC_mask2"]==None:
            #It's two fitters
            cut1_title = "Beta14"
            cut2_title = "ITR"
        else:
            #It's two DC cuts
            cut1_title = cut1_title = DCmask_dict[mb.int_to_binary_bit(int(result_dict["DC_mask1"]))]
            cut2_title = cut2_title = DCmask_dict[mb.int_to_binary_bit(int(result_dict["DC_mask2"]))]
    return cut1_title, cut2_title


if __name__ == '__main__':
    #For fun, let's grab all the results
    if DEBUG:
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
    
    #The following takes bifurcation analysis results for single cuts and
    #plots out the pearson coefficients in a grid fashion
    allresult_filenames = glob.glob(RESULTDIR + "/*_results.out")
    cut1_list = []
    cut2_list = []
    PC_list = []
    for fname in allresult_filenames:
        result_dict = rg.GetResultDict(fname) #Has the file's results
        BA = le.BifurAnalysisRun(result_dict, acceptance_rates)
        PC_list.append(BA.pearson_coeff())
        cut1_title, cut2_title = getTitles(result_dict)
        cut1_list.append(cut1_title)
        cut2_list.append(cut2_title)
    #Now, if you have all the proper files to fill out a square grid, the
    #increments you should split the lists above into are the sqrt the
    #list length
    numcolumns = np.sqrt(float(len(allresult_filenames)))
    if numcolumns / int(numcolumns) != 1:
        print("It looks like you don't have all your results needed to make" + \
                "a nice nxn grid.  Just warning you...")
    
    isSymmetric = True
    column_titles, row_titles, PC_rows = getPhiRows(cut1_list, cut2_list, PC_list,isSymmetric)
    plots.CorrelationBoxes(PC_rows, column_titles, row_titles)
