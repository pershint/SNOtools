#Class takes in results from a bifurcation analysis and calculates the
#Estimated contamination.

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import playDarts as pd

class BifurAnalysisRun(object):
    '''
    bifur_boxes is a dictionary of the form output from the lib.resultgetter
    functions.  One of the acceptance dicts should be fed in from above.
    #ASSUMPTIONS: b is the # events that pass the fitter, but fail DC
                  c is the # events that pass DC, but fail the fitter
    #With the above assumptions, by default, y_1 = y_dc and y_2 = y_fit
    '''
    def __init__(self, bifur_boxes, acceptances):
        #SNO #SNO+Wat
        self.a = bifur_boxes['a'] #369. #9.
        self.b = bifur_boxes["b"] #447. #1.
        self.c = bifur_boxes["c"] #15.  #1.
        self.d = bifur_boxes["d"] #94264.  #51.

        self.x_dc = acceptances["DC"] #Acceptance rate of N16 events in data cleaning
        self.x_fc = acceptances["Fit"] #Acceptance rate of N16 events in fit classifiers
        self.x_dc_unc = acceptances["DC_unc"]
        self.x_fc_unc = acceptances["Fit_unc"]

        #background and signal according to the data cleaning cuts
        self.bkg_events = self.b+self.c+self.d
        self.phy_events = self.a
        self.total_events = self.bkg_events + self.phy_events

    def y1y2(self):
        #Third equation in the bifurcation analysis.  Assume signal is box a
        return (self.a - (self.x_dc * self.x_fc * self.a))/self.bkg_events

    def y_dc(self):
        return ((self.a + self.b) - self.x_dc*self.a)/self.bkg_events

    def y_fit(self):
        return ((self.c + self.a) - self.x_fc*self.a)/self.bkg_events
    
    def leastsq_y1y2(self):
        phi_1 = self.y1y2()
        phi_2 = self.y_dc() * self.y_fit()
        if phi_1 == 0:
            print("THE y1y2 TERM IS ZERO.  ASSUMING THE PRODUCT IS THE OTHER")
            print("TERM BECAUSE THE LEAST SQUARES SOLUTION IS INVALID")
            return phi_2
        if phi_2 == 0:

            print("THE Y_DC OR Y_FIT EQUATION IS ZERO.  ASSUMING THE PRODUCT IS THE OTHER")
            print("TERM BECAUSE THE LEAST SQUARES SOLUTION IS INVALID")
            return phi_1
        #Now, we have two equations: p1 - y1y2 = 0, and p2 - y1y2 = 0
        #The analytical solution for the sum of the squares will give the
        #"least wrong" value for y1y2. Ignore the imaginary if you're solving
        #directly...
        return (phi_2*phi_1**2 + phi1*phi_2**2)/(phi_1**2 + phi_2**2)


    def y1y2beta_bootstrap_unc(self, CL, n):
        '''
        Returns the upper limit of y1y2beta to a 90% CL using re-fires of a,b,c, and
        d as well as the acceptance values within their uncertainties.  n is the
        number of fired experiments used to calculate y1y2beta
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
        x_dc_vals = pd.RandShoot_g0(self.x_dc, self.x_dc_unc, n)
        betavals = dvals + bvals + cvals
        y_fit_distribution = ((cvals + avals) - x_fit_vals*avals)/(betavals)
        y_dc_distribution = ((bvals + avals) - x_dc_vals*avals)/(betavals)
        #Now, for the final contamination estimate, assume all sample is bkg
        total = avals + bvals + cvals + dvals
        y1y2beta = np.sort(y_fit_distribution*y_dc_distribution*total)
        y1y2beta_CL = y1y2beta[int(CL*n)]
        return y1y2beta, y1y2beta_CL

    def y_dc_bootstrap_unc(self, n):
        '''
        Returns an array of y_dc values of length n based on random shots of
        a,b,c,and d (from a poisson of the average evaluated in this experiment).
        A gaussian of mean x_dc and variance x_dc_unc is assumed.
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
        y_dc_distribution = ((bvals + avals) - x_dc_vals*avals)/betavals
        return y_dc_CL

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
        y_fit_distribution = ((cvals + avals) - x_fit_vals*avals)/(betavals)
        return y_fit_distribution

    def pearson_coeff(self):
        '''
        Returns the pearson coefficient of our Bifurcation Analysis results.
        Values near -1 or +1 indicate highly correlated binary variables.
        '''
        n11 = float(self.a)
        n00 = float(self.d)
        n10 = float(self.c)
        n01 = float(self.b)
        nA1 = n01 + n11
        n1A = n10 + n11
        n0A = n01 + n00
        nA0 = n10 + n00
        numerator = (n11*n00) - (n10*n01)
        denominator = np.sqrt(nA1*n1A*n0A*nA0)
        if (numerator == 0) and (denominator==0):
            return 0.0
        else:
            return numerator / denominator

    def event_contamination(self):
        return self.leastsq_y1y2() * self.total_events

    def event_contamination_unc(self):
        #Returns the uncertainties based on only the first two equations
        #Uncertainty for total events not included here
        t_1 = self.y_dc_unc() * self.y_fit() * self.total_events
        t_2 = self.y_dc() * self.y_fit_unc() * self.total_events
        return np.sqrt((t_1)**2 + (t_2)**2)
#####DEPRICATED APPROACH TO GRAPHING UNCERTAINTY######

    def y_2_unc(self):
        #Calculate the uncertainty in y_2
        t_1 = ((1-self.x_fc)*np.sqrt(self.a))/(self.bkg_events)
        t_2 = (np.sqrt(self.b)/(self.bkg_events))*(1. + \
            ((self.a + self.b) - (self.x_fc*self.a))/ \
            (self.bkg_events))
        t_3 = (self.x_fc_unc*self.a)/(self.bkg_events)
        t_4 = np.sqrt(self.c + self.d)*(((self.a + self.b) - (self.x_fc*self.a)) / \
            (self.bkg_events**2))
        return np.sqrt((t_1**2)+(t_2**2)+(t_3**2)+(t_4**2))

    def y_1_unc(self):
        #Calculate the uncertainty in y_1
        t_1 = ((1-self.x_dc)*np.sqrt(self.a))/(self.bkg_events)
        t_2 = (np.sqrt(self.c)/(self.bkg_events))*(1. + \
            ((self.a + self.c - (self.x_dc*self.a))/ \
            (self.bkg_events)))
        t_3 = (self.x_dc_unc*self.a)/(self.bkg_events)
        t_4 = np.sqrt(self.b + self.d)*(((self.a + self.c) - (self.x_dc*self.a)) / \
            (self.bkg_events**2))
        return np.sqrt((t_1**2)+(t_2**2)+(t_3**2)+(t_4**2))

   #Functions defining how "off" each bifurcation analysis eqn is from true
    def p_1(self,y_1):
        return self.a + self.b - (self.x_dc * self.phy_events) - \
                (y_1 * self.bkg_events)

    def p_2(self,y_2):
        return self.a + self.c - (self.x_fc * self.phy_events) - \
                (y_2 * self.bkg_events)

    def p_3(self, y_1, y_2):
        return self.a - (self.x_dc * self.x_fc * self.phy_events) - \
    (y_1 * y_2 * self.bkg_events)

    def wrongness_val(self, y_1, y_2):
        return np.sqrt((self.p_1(y_1))**2 + (self.p_2(y_2))**2 + \
                (self.p_3(y_1, y_2))**2)

    def __call__(self, y_1, y_2):
        return self.wrongness_val(y_1, y_2)
####END DEPRICATED#####

if __name__ == "__main__":
    print("NO IMPLEMENTATION OF THE STUFFS HERE")
