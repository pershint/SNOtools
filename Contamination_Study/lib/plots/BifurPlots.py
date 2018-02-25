#File contains plotting functions for use in main.  
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl

def BoxDistribution(Bifurcation_boxes):
    '''
    Takes in the BifurcationAnalysis class and plots out the heatmap giving the
    best fit solution for leakage rates due to branch 1 and branch 2.
    '''
    try:
        entries = np.array([[Bifurcation_boxes['b'], Bifurcation_boxes['d']], \
                [Bifurcation_boxes['a'],Bifurcation_boxes['c']]])
    except KeyError:
        print("Key error in bifurcation box summary given.  Check you are")
        print("inputting a bifurcation summary dictionary")
    #Plot the heatmap showing where the y_dc and y_fit minima are
    im = plt.imshow(entries,interpolation='none', aspect='auto')
    plt.colorbar()
    plt.tick_params(axis='both', which='both', bottom='off', top='off', \
            labelbottom='off', labelleft='off')
    plt.text(0.0, 0.0, str(int(BifurAnalysis.b)), size = '20', \
            backgroundcolor='white')
    plt.text(1.0, 0.0, str(int(BifurAnalysis.d)), size = '20', \
            backgroundcolor='white')
    plt.text(1.0, 1.0, str(int(BifurAnalysis.c)), size = '20', \
            backgroundcolor='white')
    plt.text(0.0, 1.0, str(int(BifurAnalysis.a)), size = '20', \
            backgroundcolor='white')
    plt.xlabel("Data Cleaning Branch (PASS | FAIL)")
    plt.ylabel("Fit Classifier Branch (PASS | FAIL)")
    plt.title("Result of Bifurcation Analysis")
    plt.show()


def Contamination_Minimum(BifurAnalysis):
    '''
    Takes in the BifurcationAnalysis class and plots out the heatmap giving the
    best fit solution for leakage rates due to branch 1 and branch 2.
    '''
    y1_range = [0., 0.1]
    y2_range = [0., 0.1]
    extend_range = y1_range + y2_range
    y1_arr = np.arange(y1_range[0], y1_range[1], 0.0001) #Leakage rate of DC cuts
    y2_arr = np.arange(y2_range[0],y2_range[1], 0.0001) #Leakage rate of Fit Class. Cuts
    X,Y = pl.meshgrid(y1_arr, y2_arr)

    #Plot the heatmap showing where the y_dc and y_fit minima are
    Z = BifurAnalysis(X,Y)
    print("MINIMUM VALUE FOR Y_DC AND Y_FIT BASED ON LEAST SQUARES OF NUISANCE PARAMETER\n")
    im = plt.imshow(Z, aspect='auto', origin='lower',extent=extend_range)
    plt.colorbar()
    plt.xlabel("Data Cleaning Cut Leakage Fraction")
    plt.ylabel("Fit Classifier Cut Leakage Fraction")
    plt.title("Goodness of solution parameter in Contamination\n" +\
            "Contamination numbers from SNO")
    plt.show()

def LeakageUncertainty_Bootstrap_Plot(BifurAnalysis):
    y_dc_spread = BifurAnalysis.y_dc_bootstrap_unc(1000000)
    y_fit_spread = BifurAnalysis.y_fit_bootstrap_unc(1000000)
    y1y2b_spread = BifurAnalysis.y1y2beta_bootstrap_unc(0.90,1000000)
    plt.hist(y_dc_spread, 100, facecolor='green', alpha=0.75)
    plt.title("Distribution of data cleaning cut branch fractional leakage")
    plt.xlabel("Y_DC")
    plt.ylabel("Probability (arb. units)")
    plt.show()
    plt.hist(y_fit_spread, 200, facecolor='blue', alpha=0.75)
    plt.title("Distribution of fit classifier cut branch fractional leakage")
    plt.xlabel("Y_Fit")
    plt.ylabel("Probability (arb. units)")
    plt.show()
    plt.hist(y1y2b_spread, 200, facecolor='red', alpha=0.75)
    plt.title("Distribution of total estimated contamination")
    plt.xlabel("Y_Fit*Y_DC*beta")
    plt.ylabel("Probability (arb. units)")
    plt.show()
