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

        #Assume that all events in b, c, and d are background
        self.bkg_events = self.b + self.c + self.d

        self.x1 = 1.0 - Sacrifice_Summary['cut1']['total_fracsac']
        self.x2 = 1.0 - Sacrifice_Summary['cut2']['total_fracsac']
        self.x1_unc = Sacrifice_Summary['cut1']['total_fracsac_unc']
        self.x2_unc = Sacrifice_Summary['cut2']['total_fracsac_unc']

        #Build the dictionary that will save various results from the
        #Contamination Estimator
        self.contamination_summary = {}

    #Equation 1 of bifur. analysis
    def __y_1(self,a,b,x1,bkg):
        return ((a + b) - x1*a)/bkg
    #Equation 2 of bifur. analysis
    def __y_2(self,a,c,x2,bkg):
        return ((a + c) - x2*a)/bkg
    #Equation 3 of bifur. analysis
    def __y1y2(self, a, x1, x2, bkg):
        return (a - (x1 * x2 * a))/bkg

    def avg_y1y2(self,a,b,c,x1,x2,bkg):
        return (self.__y1y2(a,x1,x2,bkg) + self.__y_1(a,b,x1,bkg)*self.__y_2(a,c,x2,bkg))/2.0

    #FIXME: Need to have these be functional fills
    def highest_y1y2(self,a,b,c,x1,x2,bkg):
        eqn3 = self.__y1y2(a, x1, x2, bkg)
        eqn12 = self.__y_1(a,b,x1,bkg)*self.__y_2(a,c,x2,bkg)
        if eqn3 > eqn12:
            return eqn3
        else:
            return eqn12 

    def leastsq_y1y2(self,a,b,c,x1,x2,bkg):
        phi_1 = self.__y1y2(a,x1,x2,bkg)
        phi_2 = self.__y_1(a,b,x1,bkg) * self.__y_2(a,c,x2,bkg)
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
        avg_y1y2s = self.avg_y1y2(ashot, bshot, cshot, x1_shot, x2_shot, bkgevts)
        #Add estimated background of zero for events with no background back
        zero_results = np.zeros(numdivzeros)
        avg_y1y2s = np.append(avg_y1y2s, zero_results)
        np.sort(avg_y1y2s)
        y1y2_CL = avg_y1y2s[int(float(CL)*float(len(avg_y1y2s)))]
        self.contamination_summary["CL"] = CL
        self.contamination_summary["y1y2_to_CL"] = y1y2_CL
        return avg_y1y2s

    def CalculateContaminationValues(self):
        self.contamination_summary["y1"] = self.__y_1(self.a, self.b,\
                self.x1,self.bkg_events)
        self.contamination_summary["y2"] = self.__y_2(self.a, self.c,\
                self.x2,self.bkg_events)
        self.contamination_summary["y1*y2"] = self.__y1y2(self.a,self.x1,\
                self.x2,self.bkg_events)
        self.contamination_summary["highest_y1y2"] = self.highest_y1y2(self.a,\
                self.b,self.c,self.x1,self.x2,self.bkg_events)
        self.contamination_summary["leastsq_y1y2"] = self.leastsq_y1y2(self.a,\
                self.b,self.c,self.x1,self.x2,self.bkg_events)
        self.contamination_summary["avg_y1y2"] = self.avg_y1y2(self.a,
                self.b,self.c,self.x1,self.x2,self.bkg_events)
        self.contamination_summary["est_bkg_evts"] = self.bkg_events

    def SaveContaminationSummary(self,savedir,savename):
        savesacloc = savedir+"/"+savename
        with open(savesacloc,"w") as f:
            json.dump(self.contamination_summary, f, sort_keys=True,indent=4)
