#This program is a one-shot program that outputs the results of the SNO+
#Conamination study.  First, the data cleanning and fit sacrifice are estimated
#Using the given calibration files.  Then, the bifurcation analysis is run over
#the given data physics files.  Finally, an estimate on the contamination for the
#chosen configuration of cuts/ROI choice is output.
import matplotlib.pyplot as plt

import ROOT
import numpy as np
import lib.Bifurcator as bi
import lib.SacrificeHists as sh
import lib.plots.SacrificePlots as sp
import lib.plots.BifurPlots as bp
import lib.SacrificeAnalyzer as sa
import lib.ContaminationAnalyzer as ca
import lib.ConfigParser as cp
import lib.CalibSelector as cs
import lib.ResultUtils as ru
import os,sys
import glob
import json

#FIXME: set up a simpler argparser here.  you can choose a config file, override
#some of the configuration values with flags (changes dictionary values), and
#Could have flags to only run sacrifice, bifurcation, bifurcation analysis, or all.
CONFIGFILE='cuts_def_oldschool.json'
JOBNUM=1
ZCUT=600.0
SOURCE="N16"
PLOTS=False
SACANALYSIS=False
BIFURCATE=False
ESTIMATECONTAMINATION=True
ERANGE=None
DEBUG=True

MAINDIR = os.path.dirname(__file__)
RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "output","results_j"+str(JOBNUM)))
if not os.path.exists(RESULTDIR):
    os.makedirs(RESULTDIR)
DBDIR = os.path.abspath(os.path.join(MAINDIR, "DB"))
CALIBDIR = '/home/onetrueteal/share/May2016_N16'#os.path.abspath(os.path.join(MAINDIR, "ntuples", "N16"))
PHYSDIR = os.path.abspath(os.path.join(MAINDIR, "ntuples", "physics_data"))
CONFIGDIR = os.path.abspath(os.path.join(MAINDIR, 'config'))

if __name__ == '__main__':
    #Get your run files to load
    #Load the configuration file to use
    if SACANALYSIS or BIFURCATE is True:
        ConfigParser = cp.ConfigParser(CONFIGDIR+"/"+CONFIGFILE)
        config_dict = ConfigParser.Load_JsonConfig()
        if ERANGE is not None:
            config_dict["E_low"] = ERANGE[0]
            config_dict["E_high"] = ERANGE[1]
        ConfigParser.SaveConfiguration(config_dict,RESULTDIR,"used_configuration.json")

    if SACANALYSIS is True:
        N16_roots = glob.glob(CALIBDIR+"/*.ntuple.root")
        print("LEN OF N16: " + str(len(N16_roots)))
        print("N16_ROOTS: " + str(N16_roots))
        if ZCUT is not None:
            N16_roots = cs.ApplyZCut(DBDIR,SOURCE,ZCUT,N16_roots)
        print("LEN AFTER ZCUT: " + str(len(N16_roots)))
        ru.save_calib_list(RESULTDIR, N16_roots)       
        SacHists = sh.SacrificeHistGen(rootfiles=N16_roots,config_dict=config_dict,
                sourcetype=SOURCE)
        SacHists.GenerateHistograms()
        SacHists.SaveHistograms(RESULTDIR)
    
        SacSysUnc = sa.SacrificeSystematics(Sacrifice_Histograms=SacHists,\
                config_dict=config_dict)
        SacSysUnc.LoadCalibrationPositions(DBDIR)
        SacSysUnc.CalculateSacrifices()
        if DEBUG is True:
            SacSysUnc.ShowSacrificeResults()
        SacSysUnc.SaveSacrificeByRun(RESULTDIR,"cut_sacrifices_byrun.json")
        SacSysUnc.SaveSacrificeSummary(RESULTDIR,"cut_sacrifices_total.json")

    if PLOTS is True:
       cut_sacrifices = ru.LoadJson(RESULTDIR,"cut_sacrifices_byrun.json")
       cut_sac_summary = ru.LoadJson(RESULTDIR,"cut_sacrifices_total.json")
       for axis in ['x','y','z']:
           sp.plot_sacrificevsCart(cut_sacrifices, 'cut1', axis)
           sp.plot_sacrificevsCart(cut_sacrifices, 'cut2', axis)
       for cut in ['cut1','cut2']:
           sp.plot_sacrificevsR(cut_sacrifices, cut)
           #sp.plot_XYSacrifice(cut_sacrifices, cut)

    #Run bifurcation analysis on Physics files
    if BIFURCATE is True:
        physics_roots = glob.glob(PHYSDIR+"/*.ntuple.root")
        ru.save_physics_list(RESULTDIR,physics_roots)
        print("PHYS_ROOTS: " + str(physics_roots))
        Bifurcator = bi.Bifurcator(rootfiles=physics_roots,config_dict=config_dict)
        Bifurcator.Bifurcate()
        Bifurcator.SaveBifurcationSummary(RESULTDIR,"bifurcation_boxes.json")
    if PLOTS is True:
        bifurcation_summary = ru.LoadJson(RESULTDIR,"bifurcation_boxes.json")
        ru.BoxDistribution(bifurcation_summary)
    if ESTIMATECONTAMINATION is True:
        bifurcation_summary = ru.LoadJson(RESULTDIR,"bifurcation_boxes.json")
        cut_sac_summary = ru.LoadJson(RESULTDIR,"cut_sacrifices_total.json")
        CE = ca.ContaminationEstimator(bifurcation_summary,cut_sac_summary)
        CE.CalculateContaminationValues() #Calculate contamination eqns.
        if PLOTS is True:
            values = CE.BootstrapCL(0.90,100000) #Estimate upper end of y1y2
            values=values*CE.contamination_summary['est_bkg_evts']
            plt.hist(values,100,range=(min(values),max(values)))
            plt.show()
        CE.SaveContaminationSummary(RESULTDIR,"contamination_summary.json")
