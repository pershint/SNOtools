#This program is a one-shot program that outputs the results of the SNO+
#Conamination study.  First, the data cleanning and fit sacrifice are estimated
#Using the given calibration files.  Then, the bifurcation analysis is run over
#the given data physics files.  Finally, an estimate on the contamination for the
#chosen configuration of cuts/ROI choice is output.
import matplotlib.pyplot as plt

import numpy as np
import argparse
import os,sys
import glob
import json

#FIXME: set up a simpler argparser here.  you can choose a config file, override
#some of the configuration values with flags (changes dictionary values), and
parser = argparse.ArgumentParser(description='Parser to decide what analysis to do')
parser.add_argument('--nosave', dest='NOSAVE',action='store_true',
        help='Do not save any outputs; ensures no writing is done')
parser.add_argument('--debug', dest='debug',action='store_true',
        help='Run code in debug mode')
parser.add_argument('--jobnum', dest='JOBNUM', action='store',
        help='Specify this jobs number among others.  Will save results'+\
                'to ./output/results_jN')
parser.add_argument('--configfile', dest='CONFIGDIR',action='store',
        type=str,help='specify the config file that will be used (JSON format)')
parser.add_argument('--resultdir', dest='RESULTDIR',action='store',
        type=str,help='specify the location and filename for results to be read'+\
                'or written to.  No job number support with this flag called.')
parser.add_argument('--sacrifice', dest='SACANALYSIS',action='store_true',
        help='Run the code that plots the correlations of different cuts/classifiers')
parser.add_argument('--bifurcate', dest='BIFURCATE',action='store_true',
        help='Run the bifurcation analysis on physics data.  Saves a bifurcation summary')
parser.add_argument('--contamination', dest='ESTIMATECONTAMINATION',
        action='store_true', help='Run the contamination estimation using'+\
                'bifurcation and sacrifice results.  Save a summary.')
parser.add_argument('--plots', dest='PLOTS', action='store_true',
        help='Show plots resulting from sacrifice and contamination studies.  If no sacrifice or contamination study, loads results from the result directory and plots what is available.')
parser.add_argument('--source', dest='SOURCE', action='store',
        help='Specify the source to use for sacrifice estimation.  Currently'+\
                'supported: N16 or AmBe')
parser.add_argument('--erange', dest='ERANGE', action='store',nargs='+',
        help='Specify an energy range to run all analyses over.  If running'+\
                '--plots only, will check range matches that in results'+\
                'directory (usage: --erange 2.0 5.0)')
parser.add_argument('--zrange', dest='ZRANGE', action='store',nargs='+',
        help='Specify an upper and lower zcut in cm range to perform the analysis'+\
                'over.  Applied in sacrifice and comtanimation studies.'+\
                '(usage: --zrange 600 -500)')

MAINDIR = os.path.dirname(__file__)
CALIBDIR = '/home/onetrueteal/share/May2016_N16_2'#os.path.abspath(os.path.join(MAINDIR, "ntuples", "N16"))
cfg_default = os.path.abspath(os.path.join(MAINDIR, 'config','cuts_default.json'))
parser.set_defaults(NOSAVE=False,SACANALYSIS=False,BIFURCATE=False,debug=False,
        ESTIMATECONTAMINATION=False,JOBNUM=0,PLOTS=False,erange=None,SOURCE='N16',
        RESULTDIR=None,CONFIGDIR=cfg_default,ZRANGE=None)
args = parser.parse_args()

DEBUG = args.debug
NOSAVE=args.NOSAVE
PLOTS=args.PLOTS
SACANALYSIS=args.SACANALYSIS
BIFURCATE=args.BIFURCATE
ESTIMATECONTAMINATION=args.ESTIMATECONTAMINATION
ERANGE=args.ERANGE
ZRANGE=args.ZRANGE
JOBNUM=args.JOBNUM
SOURCE=args.SOURCE
RESULTDIR=args.RESULTDIR
CONFIGDIR=args.CONFIGDIR

print("ZRANGE: " + str(ZRANGE))

import ROOT
import lib.Bifurcator as bi
import lib.SacrificeHists as sh
import plots.SacrificePlots as sp
import plots.BifurPlots as bp
import lib.SacrificeAnalyzer as sa
import lib.ContaminationAnalyzer as ca
import lib.ConfigParser as cp
import lib.CalibSelector as cs
import lib.ResultUtils as ru

if RESULTDIR is None:
    RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "output","results_j"+str(JOBNUM)))
if not os.path.exists(RESULTDIR):
    os.makedirs(RESULTDIR)
DBDIR = os.path.abspath(os.path.join(MAINDIR, "DB"))
PHYSDIR = os.path.abspath(os.path.join(MAINDIR, "ntuples", "physics_data"))

if __name__ == '__main__':
    #Get your run files to load
    #Load the configuration file to use
    if SACANALYSIS or BIFURCATE is True:
        ConfigParser = cp.ConfigParser(CONFIGDIR)
        config_dict = ConfigParser.Load_JsonConfig()
        if ZRANGE is not None:
            config_dict["Z_high"] = float(ZRANGE[0])
            config_dict["Z_low"] = float(ZRANGE[1])
        if ERANGE is not None:
            config_dict["E_low"] = float(ERANGE[0])
            config_dict["E_high"] = float(ERANGE[1])
        if NOSAVE is False:
            ConfigParser.SaveConfiguration(config_dict,RESULTDIR,"used_configuration.json")

    if SACANALYSIS is True:
        N16_roots = glob.glob(CALIBDIR+"/*.ntuple.root")
        print("LEN OF N16: " + str(len(N16_roots)))
        print("N16_ROOTS: " + str(N16_roots))
        if ZRANGE is not None:
            N16_roots = cs.ApplyZCut(DBDIR,SOURCE,ZRANGE,N16_roots)
        print("LEN AFTER ZCUT: " + str(len(N16_roots)))
        ru.save_calib_list(RESULTDIR, N16_roots)       
        SacHists = sh.SacrificeHistGen(rootfiles=N16_roots,config_dict=config_dict,
                sourcetype=SOURCE)
        SacHists.GenerateHistograms()
        if NOSAVE is False:
            SacHists.SaveHistograms(RESULTDIR)
    
        SacSysUnc = sa.SacrificeSystematics(Sacrifice_Histograms=SacHists,\
                config_dict=config_dict)
        SacSysUnc.LoadCalibrationPositions(DBDIR)
        SacSysUnc.CalculateSacrifices()
        if DEBUG is True:
            SacSysUnc.ShowSacrificeResults()
        if NOSAVE is False:
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
        if NOSAVE is False:
            Bifurcator.SaveBifurcationSummary(RESULTDIR,"bifurcation_boxes.json")
    if PLOTS is True:
        bifurcation_summary = ru.LoadJson(RESULTDIR,"bifurcation_boxes.json")
        bp.BoxDistribution(bifurcation_summary)
    if ESTIMATECONTAMINATION is True:
        bifurcation_summary = ru.LoadJson(RESULTDIR,"bifurcation_boxes.json")
        cut_sac_summary = ru.LoadJson(RESULTDIR,"cut_sacrifices_total.json")
        CE = ca.ContaminationEstimator(bifurcation_summary,cut_sac_summary)
        CE.CalculateContaminationValues() #Calculate contamination eqns.
        values = CE.BootstrapCL(0.90,100000) #Estimate upper end of y1y2
        if PLOTS is True:
            values=values*CE.contamination_summary['est_bkg_evts']
            plt.hist(values,100,range=(min(values),max(values)))
            plt.xlabel(r"Total estimated contamination (y1y2$\beta$)")
            plt.ylabel(r"Relative probability (unitless)")
            plt.title("Distribution of estimated contamination after\n"+\
                    "re-firing variables with statistical uncertainties")
            plt.show()
        if NOSAVE is False:
            CE.SaveContaminationSummary(RESULTDIR,"contamination_summary.json")
