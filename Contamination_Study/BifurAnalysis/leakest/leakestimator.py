#This program takes the variables output from the bifurcation analysis and
#Finds the best estimate on the leakage rates for both the data cleaning cuts
#and the fit classifiers.

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import lib.playDarts as pd

#RESULTS ARE FROM RUNNING BIFURCATION ANALYSIS ON 11 DAYS OPEN GOLD DATA
#FIXME: Could make a sub-directory holding all of these results w/ descriptions

bifurcation_boxes = {"a": 9., "b": 1., "c": 1., "d": 83.}
bifurcation_boxes_SNO = {"a":369., "b":447., "c":15., "d":94264.}
acceptance_rates = {"DC": 0.9273, "DC_unc": 0.0007, "Fit":0.9984, \
        "Fit_unc": 0.0005}
acceptance_rates_SNO = {"DC": 0.9998, "DC_unc": 0.00005, "Fit": 0.995, "Fit_unc":0.0005}

class BifurAnalysisRun(object):
    def __init__(self, bifur_boxes, acceptances):
                    #SNO #SNO+Wat
        self.a = bifur_boxes["a"] #369. #9.
        self.b = bifur_boxes["b"] #447. #1.
        self.c = bifur_boxes["c"] #15.  #1.
        self.d = bifur_boxes["d"] #94264.  #51.

        self.x_dc = acceptances["DC"] #Acceptance rate of N16 events in data cleaning
        self.x_fc = acceptances["Fit"] #Acceptance rate of N16 events in fit classifiers
        self.x_dc_unc = acceptances["DC_unc"]
        self.x_fc_unc = acceptances["Fit_unc"]

        self.bkg_events = self.b+self.c+self.d
        self.phy_events = self.a
        self.total_events = self.bkg_events + self.phy_events

    #Functions defining how "off" each bifurcation analysis eqn is from true
    def p_1(self,y_1):
        return self.a + self.c - (self.x_dc * self.phy_events) - \
                (y_1 * self.bkg_events)

    def p_2(self,y_2):
        return self.a + self.b - (self.x_fc * self.phy_events) - \
                (y_2 * self.bkg_events)

    def p_3(self, y_1, y_2):
        return self.a - (self.x_dc * self.x_fc * self.phy_events) - \
    (y_1 * y_2 * self.bkg_events)

    def wrongness_val(self, y_1, y_2):
        return np.sqrt((self.p_1(y_1))**2 + (self.p_2(y_2))**2 + \
                (self.p_3(y_1, y_2))**2)

    def y_dc(self):
        return ((self.a + self.c) - self.x_dc*self.a)/self.bkg_events

    def y_fit(self):
        return ((self.b + self.a) - self.x_fc*self.a)/self.bkg_events

    def y_dc_unc(self):
        #Calculate the uncertainty in y_1
        t_1 = ((1-self.x_dc)*np.sqrt(self.a))/(self.bkg_events)
        t_2 = (np.sqrt(self.c)/(self.bkg_events))*(1. + \
            ((self.a + self.c - (self.x_dc*self.a))/ \
            (self.bkg_events)))
        t_3 = (self.x_dc_unc*self.a)/(self.bkg_events)
        t_4 = np.sqrt(self.b + self.d)*(((self.a + self.c) - (self.x_dc*self.a)) / \
            (self.bkg_events**2))
        return np.sqrt((t_1**2)+(t_2**2)+(t_3**2)+(t_4**2))

    def y_dc_bootstrap_unc(self, n):
        '''
        Returns an array of y_dc values calculated using values sampled from a 
        poisson of average a,b,c, and d, and a gaussian of mean x_dc and variance
        x_dc_unc.
        '''
        avals = pd.RandShoot_p(self.a, n)
        bvals = pd.RandShoot_p(self.b, n)
        cvals = pd.RandShoot_p(self.c, n)
        dvals = pd.RandShoot_p(self.d, n)
        avals.astype(float)
        bvals.astype(float)
        cvals.astype(float)
        dvals.astype(float)
        x_dc_vals = pd.RandShoot_g0(self.x_dc, self.x_dc_unc, n)
        betavals = dvals + bvals + cvals
        y_dc_unc = ((cvals + avals) - x_dc_vals*avals)/betavals
        return y_dc_unc

    def y_fit_unc(self):
        #Calculate the uncertainty in y_1
        t_1 = ((1-self.x_fc)*np.sqrt(self.a))/(self.bkg_events)
        t_2 = (np.sqrt(self.b)/(self.bkg_events))*(1. + \
            ((self.a + self.b) - (self.x_fc*self.a))/ \
            (self.bkg_events))
        t_3 = (self.x_fc_unc*self.a)/(self.bkg_events)
        t_4 = np.sqrt(self.c + self.d)*(((self.a + self.b) - (self.x_fc*self.a)) / \
            (self.bkg_events**2))
        return np.sqrt((t_1**2)+(t_2**2)+(t_3**2)+(t_4**2))

    def y_fit_bootstrap_unc(self, n):
        '''
        Returns an array of y_fit values calculated using values sampled from a 
        poisson of average a,b,c, and d, and a gaussian of mean x_fit and variance
        x_fit_unc.
        '''
        avals = pd.RandShoot_p(self.a, n)
        bvals = pd.RandShoot_p(self.b, n)
        cvals = pd.RandShoot_p(self.c, n)
        dvals = pd.RandShoot_p(self.d, n)
        avals.astype(float)
        bvals.astype(float)
        cvals.astype(float)
        dvals.astype(float)
        x_fit_vals = pd.RandShoot_g0(self.x_fc, self.x_fc_unc, n)
        betavals = dvals + bvals + cvals
        y_fit_spread = ((bvals + avals) - x_fit_vals*avals)/(betavals)
        return y_fit_spread

    def event_contamination(self):
        return self.y_dc() * self.y_fit() * self.total_events

    def event_contamination_unc(self):
        t_1 = self.y_dc_unc() * self.y_fit() * self.total_events
        t_2 = self.y_dc() * self.y_fit_unc() * self.total_events
        return np.sqrt((t_1)**2 + (t_2)**2)

    def __call__(self, y_1, y_2):
        return self.wrongness_val(y_1, y_2)

if __name__ == "__main__":
        BA = BifurAnalysisRun(bifurcation_boxes, acceptance_rates)
        print("y_dc:\n" + str(BA.y_dc()))
        print("y_fit:\n" + str(BA.y_fit()))
        print("uncertainty on y_dc:\n" + str(BA.y_dc_unc()))
        print("uncertainty on y_fit:\n" + str(BA.y_fit_unc()))
        print("Total events fed into B.A.: " + str(BA.total_events))
        print("Contamination assuming all events are background: ")
        print(BA.event_contamination())
        print("Contamination uncertainty: " + str(BA.event_contamination_unc()))

        y1_range = [0., 0.1]
        y2_range = [0., 0.1]
        extend_range = y1_range + y2_range
        y1_arr = np.arange(y1_range[0], y1_range[1], 0.0001) #Leakage rate of DC cuts
        y2_arr = np.arange(y2_range[0],y2_range[1], 0.0001) #Leakage rate of Fit Class. Cuts
        X,Y = pl.meshgrid(y1_arr, y2_arr)

        #Plot the heatmap showing where the y_dc and y_fit minima are
        Z = BA(X,Y)
        im = plt.imshow(Z, aspect='auto', origin='lower',extent=extend_range)
        plt.colorbar()
        plt.xlabel("Data Cleaning Cut Leakage Fraction")
        plt.ylabel("Fit Classifier Cut Leakage Fraction")
        plt.title("Goodness of solution parameter in Contamination\n" +\
                "Contamination numbers from SNO")
        plt.show()

        #Show the histograms for y_dc and y_fit uncertainties with boostrapping
        y_dc_spread = BA.y_dc_bootstrap_unc(1000000)
        y_fit_spread = BA.y_fit_bootstrap_unc(1000000)
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
