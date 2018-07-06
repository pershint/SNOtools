#This code reads acceptance values and the bifurcation box values output from a
#Sacrifice analysis and bifurcation analysis, respectively.


import os,sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import playDarts as pd
import scipy as sp
import json
import ROOT
import glob

class ContaminationEstimator(object):
    def __init__(self, Bifurcation_Summary=None, Sacrifice_Summary=None):
        self.ss = Sacrifice_Summary

        self.a = Bifurcation_Summary['a']
        self.b = Bifurcation_Summary['b']
        self.c = Bifurcation_Summary['c']
        self.d = Bifurcation_Summary['d']


        self.x1 = 1.0 - Sacrifice_Summary['cut1']['total_fracsac']
        self.x2 = 1.0 - Sacrifice_Summary['cut2']['total_fracsac']
        self.x1_unc = Sacrifice_Summary['cut1']['total_fracsac_unc']
        self.x2_unc = Sacrifice_Summary['cut2']['total_fracsac_unc']

        #Build the dictionary that will save various results from the
        #Contamination Estimator
        self.contamination_summary = {}
    
    def SaveContaminationSummary(self,savedir,savename):
        savesacloc = savedir+"/"+savename
        with open(savesacloc,"w") as f:
            json.dump(self.contamination_summary, f, sort_keys=True,indent=4)
    #Equation 1 of bifur. analysis

class LowEContamination(ContaminationEstimator):
    '''Contamination estimate that requires an estimate on the
    number of signal events in the a-box'''
    def __init__(self, Bifurcation_Summary=None, Sacrifice_Summary=None,sigdict=None):
        super(LowEContamination,self).__init__(Bifurcation_Summary=Bifurcation_Summary,
                Sacrifice_Summary=Sacrifice_Summary)
        self.signal_estimate = sigdict['cut1']['nonpath_events']
        self.signal_estimate_unc = np.sqrt(self.signal_estimate)
        self.contamination_summary['type'] = 'LowE'

    def _N1(self):
        return (self.b - (self.x1*(1.0-self.x2)*self.signal_estimate))

    def _N1_unc(self):
        n1u = np.sqrt((np.sqrt(self.b) -(self.x1*(1.0-self.x2)*self.signal_estimate))**2 +\
                (self.b -(self.x1_unc*(1.0-self.x2)*self.signal_estimate))**2 +\
                (self.b + (self.x1*self.x2_unc*self.signal_estimate))**2 + \
                (self.b - (self.x1*(1.0-self.x2)*self.signal_estimate_unc))**2)
        return n1u

    def _N2_unc(self):
        n2u = np.sqrt((np.sqrt(self.c) -(self.x2*(1.0-self.x1)*self.signal_estimate))**2 +\
                (self.c -(self.x2_unc*(1.0-self.x1)*self.signal_estimate))**2 +\
                (self.c + (self.x2*self.x1_unc*self.signal_estimate))**2 + \
                (self.c - (self.x2*(1.0-self.x1)*self.signal_estimate_unc))**2)
        return n2u

    def _N2(self):
        return (self.c - (self.x2*(1.0-self.x1)*self.signal_estimate))
    
    def _D(self):
        return (self.d - ((1.0-self.x2)*(1.0-self.x1)*self.signal_estimate))
    def _D_unc(self):
        Du = np.sqrt((np.sqrt(self.d) - ((1.0-self.x2)*(1.0-self.x1)*self.signal_estimate))**2 +\
                (self.d +(self.x2_unc*(1-self.x1)*self.signal_estimate))**2 +\
                (self.d + ((1.0-self.x2)*self.x1_unc*self.signal_estimate))**2 + \
                (self.d - ((1.0-self.x1)*(1.0-self.x2)*self.signal_estimate_unc))**2)
        return Du

    def _y1y2B(self):
        return self._N1() * self._N2() / self._D()

    def CalculateContamination(self):
        '''returns estimated y1y2B value based on bifurcation box values & signal estimate'''
        y1y2B = self._N1() * self._N2() / self._D() 
        self.contamination_summary['y1y2B'] = y1y2B
        return y1y2B

    def CalculateContaminationUnc(self):
        y1y2B_unc = self._y1y2B() * np.sqrt((self._N1_unc()/self._N1())**2 + \
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
        #Add estimated background of zero for events with no background back
        zero_results = np.zeros(numdivzeros)
        _avg_y1y2s = np.append(_avg_y1y2s, zero_results)
        np.sort(_avg_y1y2s)
        y1y2_CL = _avg_y1y2s[int(float(CL)*float(len(_avg_y1y2s)))]
        self.contamination_summary["CL"] = CL
        self.contamination_summary["y1y2_to_CL"] = y1y2_CL
        return _avg_y1y2s

    def CalculateContaminationValues(self):
        self.contamination_summary["y1"] = self._y_1(self.a, self.b,\
                self.x1,self.bkg_events)
        self.contamination_summary["y2"] = self._y_2(self.a, self.c,\
                self.x2,self.bkg_events)
        self.contamination_summary["y1*y2"] = self._y1y2(self.a,self.x1,\
                self.x2,self.bkg_events)
        self.contamination_summary["highest_y1y2"] = self._highest_y1y2(self.a,\
                self.b,self.c,self.x1,self.x2,self.bkg_events)
        self.contamination_summary["leastsq_y1y2"] = self._leastsq_y1y2(self.a,\
                self.b,self.c,self.x1,self.x2,self.bkg_events)
        self.contamination_summary["avg_y1y2"] = self._avg_y1y2(self.a,
                self.b,self.c,self.x1,self.x2,self.bkg_events)
        self.contamination_summary["est_bkg_evts"] = self.bkg_events
