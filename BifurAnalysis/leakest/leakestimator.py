#This program takes the variables output from the bifurcation analysis and
#Finds the best estimate on the leakage rates for both the data cleaning cuts
#and the fit classifiers.

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl


class BifurAnalysisRun(object):
    def __init__(self):
                    #SNO #SNO+Wat
        self.a = 9. #369. #9.
        self.b = 1. #447. #1.
        self.c = 1. #15.  #1.
        self.d = 51. #94264.  #51.

        self.x_dc = 0.997 #Acceptance rate of N16 events in data cleaning
        self.x_fc = 0.998 #Acceptance rate of N16 events in fit classifiers

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

    def __call__(self, y_1, y_2):
        return self.wrongness_val(y_1, y_2)

if __name__ == "__main__":
        y1_range = [0., 0.1]
        y2_range = [0., 0.1]
        extend_range = y1_range + y2_range
        y1_arr = np.arange(y1_range[0], y1_range[1], 0.0001) #Leakage rate of DC cuts
        y2_arr = np.arange(y2_range[0],y2_range[1], 0.0001) #Leakage rate of Fit Class. Cuts
        X,Y = pl.meshgrid(y1_arr, y2_arr)
        BA = BifurAnalysisRun()
        Z = BA(X,Y)
        im = plt.imshow(Z, aspect='auto', origin='lower',extent=extend_range)

        plt.colorbar()
        plt.xlabel("Data Cleaning Cut Leakage Fraction")
        plt.ylabel("Fit Classifier Cut Leakage Fraction")
        plt.title("Goodness of solution parameter in Contamination\n" +\
                "Contamination numbers from SNO")
        plt.show()
