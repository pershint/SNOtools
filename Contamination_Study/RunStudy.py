#This program is a one-shot program that outputs the results of the SNO+
#Conamination study.  First, the data cleanning and fit sacrifice are estimated
#Using the given calibration files.  Then, the bifurcation analysis is run over
#the given data physics files.  Finally, an estimate on the contamination for the
#chosen configuration of cuts/ROI choice is output.
import matplotlib.pyplot as plt

import numpy as np
import argparse
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
parser = argparse.ArgumentParser(description='Parser to decide what analysis to do')
parser.add_argument('--debug', dest='debug',action='store_true',
        help='Run code in debug mode')
parser.add_argument('--sacrifice', dest='SACANALYSIS',action='store_true',
        help='Run the code that plots the correlations of different cuts/classifiers')
parser.add_argument('--bifurcate', dest='BIFURCATE',action='store_true',
        help='Run the bifurcation analysis on physics data.  Saves a bifurcation summary')
parser.add_argument('--contamination', dest='ESTIMATECONTAMINATION',
        action='store_true', help='Run the contamination estimation using'+\
                'bifurcation and sacrifice results.  Save a summary.')
parser.add_argument('--plots', dest='PLOTS', action='store_true',
        help='Show plots resulting from sacrifice and contamination studies.  If no sacrifice or contamination study, loads results from the result directory and plots what is available.')
parser.add_argument('--jobnum', dest='JOBNUM', action='store',
        help='Specify this jobs number among others.  Will save results'+\
                'to ./output/results_jN')
parser.add_argument('--source', dest='SOURCE', action='store',
        help='Specify the source to use for sacrifice estimation.  Currently'+\
                'supported: N16 or AmBe')
parser.add_argument('--erange', dest='ERANGE', action='store',nargs='+',
        help='Specify an energy range to run all analyses over.  If running'+\
                '--plots only, will check range matches that in results'+\
                'directory (usage: --erange 2.0 5.0)')
parser.set_defaults(SACANALYSIS=False,BIFURCATE=False,debug=False,
        ESTIMATECONTAMINATION=False,JOBNUM=0,PLOTS=False,erange=None,SOURCE='N16')
args = parser.parse_args()

DEBUG = args.debug
PLOTS=args.PLOTS
SACANALYSIS=args.SACANALYSIS
BIFURCATE=args.BIFURCATE
ESTIMATECONTAMINATION=args.ESTIMATECONTAMINATION
ERANGE=args.ERANGE
JOBNUM=args.JOBNUM
SOURCE=args.SOURCE

import ROOT

print(SOURCE)
print(ERANGE)

CONFIGFILE='cuts_def_oldschool.json'
ZCUT=600.0

MAINDIR = os.path.dirname(__file__)
RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "output","results_j"+str(JOBNUM)))
if not os.path.exists(RESULTDIR):
    os.makedirs(RESULTDIR)
DBDIR = os.path.abspath(os.path.join(MAINDIR, "DB"))
CALIBDIR = os.path.abspath(os.path.join(MAINDIR, "ntuples", "N16"))
PHYSDIR = os.path.abspath(os.path.join(MAINDIR, "ntuples", "physics_data"))
CONFIGDIR = os.path.abspath(os.path.join(MAINDIR, 'config'))

if __name__ == '__main__':
    #Get your run files to load
    #Load the configuration file to use
    if SACANALYSIS or BIFURCATE is True:
        ConfigParser = cp.ConfigParser(CONFIGDIR+"/"+CONFIGFILE)
        config_dict = ConfigParser.Load_JsonConfig()
        if ERANGE is not None:
            config_dict["E_low"] = float(ERANGE[0])
            config_dict["E_high"] = float(ERANGE[1])
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
