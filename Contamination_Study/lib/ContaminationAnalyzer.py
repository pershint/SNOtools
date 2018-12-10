#This code reads acceptance values and the bifurcation box values output from a
#Sacrifice analysis and bifurcation analysis, respectively.

from cStringIO import StringIO
import os,sys
import copy
import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import playDarts as pd
import scipy as sp
from scipy.special import gammaln
from scipy.optimize import fmin
import iminuit as im
import json
import ROOT
import glob

class ContaminationEstimator(object):
    def __init__(self, Bifurcation_Summary=None, Sacrifice_Summary=None):
        self.ss = Sacrifice_Summary
        self.acceptances_loaded = False
        self.minuit_summary = None

        print("SACRIFICE SUMMARY: " + str(Sacrifice_Summary))
        if Sacrifice_Summary is not None:
            self._CalculateSacrifices()
        
        self.a = Bifurcation_Summary['a']
        self.a_unc = np.sqrt(self.a)
        self.b = Bifurcation_Summary['b']
        self.b_unc = np.sqrt(self.b)
        self.c = Bifurcation_Summary['c']
        self.c_unc = np.sqrt(self.c)
        self.d = Bifurcation_Summary['d']
        self.d_unc = np.sqrt(self.d)
        self.N = self.a + self.b + self.c + self.d

        #Build the dictionary that will save various results from the
        #Contamination Estimator
        self.contamination_summary = {}

    def _CalculateSacrifices(self):
        x1_vals = []
        x1_statuncs = []
        x1_sysuncs = []
        x2_vals = []
        x2_statuncs = []
        x2_sysuncs = []
        try:
            for var in self.ss['cut1']:
                x1_vals.append(self.ss['cut1'][var]['sacrifice'])
                x1_statuncs.append(self.ss['cut1'][var]['stat_unc'])
                x1_sysuncs.append(self.ss['cut1'][var]['sys_unc'])
        except KeyError:
            print("No sacrifice information for cut 1.  Continuing to get cut 2 info.")
            pass
        try:
            for var in self.ss['cut2']:
                x2_vals.append(self.ss['cut2'][var]['sacrifice'])
                x2_statuncs.append(self.ss['cut2'][var]['stat_unc'])
                x2_sysuncs.append(self.ss['cut2'][var]['sys_unc'])
        except KeyError:
            print("No sacrifice information for cut 2.  Continuing to get sacrifice summary")
            pass
        self.x1 = 1.0 - np.average(x1_vals)
        self.x2 = 1.0 - np.average(x2_vals)
        self.x1_unc = lg.norm([lg.norm(x1_statuncs),lg.norm(x1_sysuncs)]) 
        self.x2_unc = lg.norm([lg.norm(x2_statuncs),lg.norm(x2_sysuncs)]) 
        print("CUT 1 ACCEPTANCE: " + str(self.x1))
        print("CUT 1 ACC. UNC: " + str(self.x1_unc))
        print("CUT 2 ACCEPTANCE: " + str(self.x2))
        print("CUT 2 ACC. UNC: " + str(self.x2_unc))
        self.acceptances_loaded = True

    def SaveContaminationSummary(self,savedir,savename):
        savesacloc = savedir+"/"+savename
        with open(savesacloc,"w") as f:
            json.dump(self.contamination_summary, f, sort_keys=True,indent=4)
    #Equation 1 of bifur. analysis

class LogLikelihoodContam(ContaminationEstimator):
    '''Estimates found for y1, y2, and B (and their uncertainties) associated with
    minimizing the negative log likelhood function calculating the probability the
    a,b,c, and d box values come from the estimates fitted.'''
    
    def lnfact(self,x):
        '''natural log of the factorial function used in negative
        log likelihood function fit'''
        if x == 0:
            return 1
        else:
            return gammaln(x)
    
    def nll(self,x,a,b,c,d,x1,x2):
        '''Negative log likelihood function that, when maximized,
        finds the y1, y2, and B values to produce the a,b,c, and d
        box entry counts given acceptances x1 and x2'''
        y1, y2, B = x
        N = a + b + c + d
        S = N - B
    
        expected_a = x1*x2*S + y1*y2*B
        expected_b = (1-x2)*x1*S + (1-y2)*y1*B
        expected_c = x2*(1-x1)*S + y2*(1-y1)*B
        expected_d = (1-x1)*(1-x2)*S + (1-y1)*(1-y2)*B
    
        return +expected_a - a*np.log(expected_a) + self.lnfact(a) + \
               +expected_b - b*np.log(expected_b) + self.lnfact(b) + \
               +expected_c - c*np.log(expected_c) + self.lnfact(c) + \
               +expected_d - d*np.log(expected_d) + self.lnfact(d)

    def nll_minuit(self,y1,y2,B):
        '''Negative log likelihood function that, when maximized,
        finds the y1, y2, and B values to produce the a,b,c, and d
        box entry counts given acceptances x1 and x2'''
        N = self.a + self.b + self.c + self.d
        S = N - B
    
        expected_a = self.x1*self.x2*S + y1*y2*B
        expected_b = (1-self.x2)*self.x1*S + (1-y2)*y1*B
        expected_c = self.x2*(1-self.x1)*S + y2*(1-y1)*B
        expected_d = (1-self.x1)*(1-self.x2)*S + (1-y1)*(1-y2)*B
    
        return +expected_a - self.a*np.log(expected_a) + self.lnfact(self.a) + \
               +expected_b - self.b*np.log(expected_b) + self.lnfact(self.b) + \
               +expected_c - self.c*np.log(expected_c) + self.lnfact(self.c) + \
               +expected_d - self.d*np.log(expected_d) + self.lnfact(self.d)

    def FitContaminationParameters(self,useMinuit=True):
        '''Uses scipy optimize to fit the negative log likelihood
        function defined'''
        if not self.acceptances_loaded:
            print("You must first calculate the acceptances associated "+\
                    "with your expected signal ROI using self.CalculateAcceptances")
            return
        if useMinuit:
            self.contamination_summary["Use_minuit"] = True
            im.describe(self.nll_minuit)
            #Save the migrad output to a string
            sys.stdout = stringbuf = StringIO()
            m = im.Minuit(self.nll_minuit, limit_y1=(0.0,1.0), limit_y2=(0.0,1.0), limit_B=(0.0,self.N),
                    y1 = 0.1, y2 = 0.1, B=(self.N/2.0))
            m.migrad()
            m.minos()
            print("MINIMIZATION OUTPUT: " + str(m.values))
            print("MINIMUM VALUE: " + str(m.fval))
            sys.stdout = sys.__stdout__
            self.contamination_summary['y1_bestfit'] = m.values['y1']
            self.contamination_summary['y2_bestfit'] = m.values['y2']
            self.contamination_summary['B_bestfit'] = m.values['B']
            self.minuit_summary = stringbuf
        else:
            self.contamination_summary["Use_minuit"] = False
            xopt = fmin(self.nll,[self.x1,self.x2,self.N/2.],args=(self.a,self.b,self.c,
                self.d,self.x1,self.x2))
            print("fitted y1 = %.2g" % (xopt[0]))
            print("fitted y2 = %.2g" % (xopt[1]))
            print("fitted B  = %.2g" % (xopt[2]))
            self.contamination_summary["y1_bestfit"] = xopt[0]
            self.contamination_summary["y2_bestfit"] = xopt[1]
            self.contamination_summary["B_bestfit"] = xopt[2]
    def SaveMinimizationSummary(self,savedir,savename):
        if self.minuit_summary is None:
            print("Looks like you didn't use minuit.  Not saving Minuit summary")
            return
        savesacloc = savedir+"/"+savename
        with open(savesacloc,"w") as f:
            f.write(str(self.minuit_summary.getvalue()))
            f.close()

class LowEContamination(ContaminationEstimator):
    '''Contamination estimate that approximates that the number of
    signal events in the a-box is much larger than the contamination'''
    def __init__(self, Bifurcation_Summary=None, Sacrifice_Summary=None):
        super(LowEContamination,self).__init__(Bifurcation_Summary=Bifurcation_Summary,
                Sacrifice_Summary=Sacrifice_Summary)
        self.setting_upperlimit = False
        if self.a == 0:
            self.a = 1.15
            self.a_unc = 0
            self.setting_upperlimit = True
        if self.b == 0:
            self.b = 1.15
            self.b_unc = 0
            self.setting_upperlimit = True
        if self.c == 0:
            self.c = 1.15
            self.c_unc = 0
            self.setting_upperlimit = True
        if self.d == 0:
            self.d = 1.15
            self.d_unc = 0
            self.setting_upperlimit = True
        self.N = self.a + self.b + self.c + self.d
        self.signal_estimate = self.a/(self.x1*self.x2)
        
        self.signal_estimate_unc = self.signal_estimate*np.sqrt((self.x1_unc/self.x1)**2 +\
                (self.x2_unc/self.x2)**2 + (self.a_unc/self.a)**2)
        self.contamination_summary['type'] = 'LowE'

    def _N1(self):
        return (self.b - (self.x1*(1.0-self.x2)*self.signal_estimate))
    def _N1_unc(self):

        t2_unc = abs(self.x1*(1.0-self.x2)*self.signal_estimate)*np.sqrt(\
                (self.x1_unc/self.x1)**2 + (self.x2_unc/(1-self.x2))**2 + \
                (self.signal_estimate_unc/self.signal_estimate)**2)
        n1u = np.sqrt((self.b_unc)**2 + (t2_unc)**2)
        return n1u

    def _N2_unc(self):
        t2_unc = abs((self.x2*(1.0-self.x1)*self.signal_estimate))*np.sqrt(\
                (self.x1_unc/(1-self.x1))**2 + (self.x2_unc/(self.x2))**2 + \
                (self.signal_estimate_unc/self.signal_estimate)**2)
        n2u = np.sqrt((self.c_unc)**2 + (t2_unc)**2)
        return n2u

    def _N2(self):
        return (self.c - (self.x2*(1.0-self.x1)*self.signal_estimate))
    
    def _D(self):
        return (self.d - ((1.0-self.x2)*(1.0-self.x1)*self.signal_estimate))
    
    def _D_unc(self):
        t3_unc = (((1.0-self.x2)*(1.0-self.x1)*self.signal_estimate))*np.sqrt(\
                (self.x1_unc/(1-self.x1))**2 + (self.x2_unc/(1-self.x2))**2 + \
                (self.signal_estimate_unc/self.signal_estimate)**2)
        Du = np.sqrt((self.d_unc)**2 + (t3_unc)**2)
        return Du

    def _y1y2B(self):
        return self._N1() * self._N2() / self._D()

    def CalculateContamination(self):
        '''returns estimated y1y2B value based on bifurcation box values & signal estimate'''
        print("N1: " + str(self._N1()))
        print("N2: " + str(self._N2()))
        print("DENOM: " + str(self._D()))
        y1y2B = self._N1() * self._N2() / self._D() 
        self.contamination_summary['y1y2B'] = y1y2B
        self.contamination_summary['setting_68CLbound'] = self.setting_upperlimit
        return y1y2B

    def CalculateContaminationUnc(self):
        print("N1_UNC: " + str(self._N1_unc()))
        print("N2_UNC: " + str(self._N2_unc()))
        print("DENOM_UNC: " + str(self._D_unc()))
        y1y2B_unc = abs(self._y1y2B()) * np.sqrt((self._N1_unc()/self._N1())**2 + \
                (self._N2_unc()/self._N2())**2 + (self._D_unc()/self._D())**2)
        self.contamination_summary['y1y2B_unc'] = y1y2B_unc
        return y1y2B_unc


class NDContamination(ContaminationEstimator):
    '''Contamination estimate that does not require an estimate on the
    number of signal events in the a-box'''

    def __init__(self, Bifurcation_Summary=None, Sacrifice_Summary=None):
        super(NDContamination,self).__init__(Bifurcation_Summary=Bifurcation_Summary,
                Sacrifice_Summary=Sacrifice_Summary)
        self.contamination_summary['type'] = 'NucleonDecayROI'
        #Assume that all events in b, c, and d are background
        self.bkg_events = self.b + self.c + self.d

    def _y_1(self,a,b,x1,bkg):
        return ((a + b) - x1*a)/bkg

    def _y_1_unc(self,a,b,x1,x1_unc,bkg):
        #Calculate the uncertainty in y_1
        t_1 = ((1-x1)*np.sqrt(a))/(bkg)
        t_2 = (np.sqrt(b)/(bkg))*(1. + ((a + b - (x1*a))/ \
            (bkg)))
        t_3 = (x1_unc*a)/(bkg)
        t_4 = np.sqrt(c + d)*(((a + b) - (x1*a)) / (bkg**2))
        return np.sqrt((t_1**2)+(t_2**2)+(t_3**2)+(t_4**2))

    #Equation 2 of bifur. analysis
    def _y_2(self,a,c,x2,bkg):
        return ((a + c) - x2*a)/bkg

    def _y_2_unc(self,a,c,x2,x2_unc,bkg):
        #Calculate the uncertainty in y_2
        t_1 = ((1-x2)*np.sqrt(a))/(bkg)
        t_2 = (np.sqrt(c)/(bkg))*(1. + ((a + c) - (x2*a))/(bkg))
        t_3 = (x2_unc*a)/(bkg)
        t_4 = np.sqrt(b + d)*(((a + c) - (x2*a)) /(bkg**2))
        return np.sqrt((t_1**2)+(t_2**2)+(t_3**2)+(t_4**2))

    #Equation 3 of bifur. analysis
    def _y1y2(self, a, x1, x2, bkg):
        return (a - (x1 * x2 * a))/bkg

    def _avg_y1y2(self,a,b,c,x1,x2,bkg):
        return (self._y1y2(a,x1,x2,bkg) + self._y_1(a,b,x1,bkg)*self._y_2(a,c,x2,bkg))/2.0

    def _highest_y1y2(self,a,b,c,x1,x2,bkg):
        eqn3 = self._y1y2(a, x1, x2, bkg)
        eqn12 = self._y_1(a,b,x1,bkg)*self._y_2(a,c,x2,bkg)
        if type(eqn3) is np.ndarray:
            result = []
            for j, val in enumerate(eqn3):
                if val > eqn12[j]:
                    result.append(val)
                else:
                    result.append(eqn12[j])
            return np.array(result)
        else:
            if eqn3 > eqn12:
                return eqn3
            else:
                return eqn12 

    def _leastsq_y1y2(self,a,b,c,x1,x2,bkg):
        phi_1 = self._y1y2(a,x1,x2,bkg)
        phi_2 = self._y_1(a,b,x1,bkg) * self._y_2(a,c,x2,bkg)
        if type(phi_1) is np.ndarray:
            result =(phi_2*phi_1**2 + phi_1*phi_2**2)/(phi_1**2 + phi_2**2)
            for j, val in enumerate(result):
                if val == 0:
                    #Take the higher of the two values
                    if phi_1[j] > phi_2[j]:
                        result[j] = phi_1[j]
                    else:
                        result[j] = phi_2[j]
            return result
        else:
            if phi_1 == 0:
                print("THE y1y2 TERM IS ZERO.  ASSUMING THE PRODUCT IS THE OTHER")
                print("TERM BECAUSE THE LEAST SQUARES SOLUTION IS INVALID")
                return phi_2
            if phi_2 == 0:
                print("EQNS. 1 OR 2 ARE ZERO.  ASSUMING THE PRODUCT IS THE OTHER")
                print("TERM BECAUSE THE LEAST SQUARES SOLUTION IS INVALID")
                return phi_1
            #Now, we have two equations: p1 - y1y2 = 0, and p2 - y1y2 = 0
            #The analytical solution for the sum of the squares will give the
            #"least wrong" value for y1y2. Ignore the imaginary if you're solving
            #directly...
            return (phi_2*phi_1**2 + phi_1*phi_2**2)/(phi_1**2 + phi_2**2)

    def BootstrapCL(self,CL,n):
        #Using the bifurcation boxes as averages, shoot values for a,b,c, and d
        #assuming poisson fluctuations.  Also shoots from a gaussian for the 
        #acceptances.
        ashot = pd.RandShoot_p(self.a, n)
        bshot = pd.RandShoot_p(self.b, n)
        cshot = pd.RandShoot_p(self.c, n)
        dshot = pd.RandShoot_p(self.d, n)
        bkgevts = bshot + cshot + dshot
        #see how many results have no background (beta=0). These will have an
        #estimated zero contamination. Delete the elements and add zeros to
        #y1y2 later.
        divzeros = np.where(bkgevts==0)[0]
        print("DIVZEROS: " + str(divzeros))
        numdivzeros = len(divzeros)
        ashot = np.delete(ashot,divzeros)
        bshot = np.delete(bshot,divzeros)
        cshot = np.delete(cshot,divzeros)
        dshot = np.delete(dshot,divzeros)
        bkgevts = np.delete(bkgevts,divzeros)
        n_nonz = n - numdivzeros
        x1_shot = pd.RandShoot(self.x1, self.x1_unc, n_nonz)
        x2_shot = pd.RandShoot(self.x2, self.x2_unc, n_nonz)
        _avg_y1y2s = self._avg_y1y2(ashot, bshot, cshot, x1_shot, x2_shot, bkgevts)
        _avg_y1y2Bs = _avg_y1y2s*bkgevts
        #Add estimated background of zero for events with no background back
        zero_results = np.zeros(numdivzeros)
        _avg_y1y2Bs = np.append(_avg_y1y2Bs, zero_results)
        _avg_y1y2s = np.append(_avg_y1y2s, zero_results)
        _avg_y1y2Bs = np.sort(_avg_y1y2Bs)
        _avg_y1y2s = np.sort(_avg_y1y2s)
        y1y2B_CL = _avg_y1y2Bs[int(float(CL)*float(len(_avg_y1y2Bs)))]
        y1y2_CL = _avg_y1y2s[int(float(CL)*float(len(_avg_y1y2s)))]
        self.contamination_summary["CL"] = CL
        self.contamination_summary["y1y2_to_CL"] = y1y2_CL
        self.contamination_summary["y1y2B_to_CL"] = y1y2B_CL
        return _avg_y1y2Bs

    def CalculateContaminationValues(self):
        self.contamination_summary["y1"] = self._y_1(self.a, self.b,\
                self.x1,self.bkg_events)
        self.contamination_summary["y2"] = self._y_2(self.a, self.c,\
                self.x2,self.bkg_events)
        self.contamination_summary["y1*y2"] = self._y1y2(self.a,self.x1,\
                self.x2,self.bkg_events)
        self.contamination_summary["y1*y2"] = self._y1y2(self.a,self.x1,\
                self.x2,self.bkg_events)*self.bkg_events
        self.contamination_summary["highest_y1y2"] = self._highest_y1y2(self.a,\
                self.b,self.c,self.x1,self.x2,self.bkg_events)
        self.contamination_summary["leastsq_y1y2"] = self._leastsq_y1y2(self.a,\
                self.b,self.c,self.x1,self.x2,self.bkg_events)
        self.contamination_summary["avg_y1y2"] = self._avg_y1y2(self.a,
                self.b,self.c,self.x1,self.x2,self.bkg_events)
        self.contamination_summary["est_bkg_evts"] = self.bkg_events
