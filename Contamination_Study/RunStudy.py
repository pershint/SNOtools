#This program is a one-shot program that outputs the results of the SNO+
#Conamination study.  First, the data cleanning and fit sacrifice are estimated
#Using the given calibration files.  Then, the bifurcation analysis is run over
#the given data physics files.  Finally, an estimate on the contamination for the
#chosen configuration of cuts/ROI choice is output.

import ROOT
import numpy as np
import lib.Bifurcator as bi
import lib.SacrificeHists as sh
import lib.SacrificePlots as sp
import lib.SacrificeAnalyzer as sa
import lib.ConfigParser as cp
import lib.ResultUtils as ru
import os,sys
import glob
import json

#FIXME: set up a simpler argparser here.  you can choose a config file, override
#some of the configuration values with flags (changes dictionary values), and
#Could have flags to only run sacrifice, bifurcation, bifurcation analysis, or all.
CONFIGFILE='cuts_default.json'
JOBNUM=0
SOURCE="N16"
PLOTS=True
SACANALYSIS=True
BIFURCATE=False
ERANGE=None
DEBUG=True

MAINDIR = os.path.dirname(__file__)
RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "output","results_j"+str(JOBNUM)))
if not os.path.exists(RESULTDIR):
    os.makedirs(RESULTDIR)
DBDIR = os.path.abspath(os.path.join(MAINDIR, "DB"))
CALIBDIR = '/home/onetrueteal/share/May2016_N16_2'#os.path.abspath(os.path.join(MAINDIR, "ntuples", "N16"))
PHYSDIR = os.path.abspath(os.path.join(MAINDIR, "ntuples", "OpenGolden"))
CONFIGDIR = os.path.abspath(os.path.join(MAINDIR, 'config'))

def save_physics_list(directory,fullpath):
    #Strips off the directory, places the root names in a list, and saves it
    #in the RESULTDIR
    physics_roots = []
    for fullname in fullpath:
        rootname = fullname.lstrip(fullpath+"/")
        physics_roots.append(rootname)
    analysis_list = {}
    analysis_list["runs_used_in_bifurcation"] = physics_roots
    with open(RESULTDIR+"/physics_list.json","w") as f:
        json.dump(analysis_list,f,sort_keys=True,indent=4)

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
        print("N16_ROOTS: " + str(N16_roots))
        #Generate Histograms of fractional sacrifice for Calibration data
        #These histograms are used to output a sacrifice and uncertainty below
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
           sp.plot_XYSacrifice(cut_sacrifices, cut)
    #Run bifurcation analysis on Physics files
    if BIFURCATE is True:
        physics_roots = glob.glob(PHYSDIR+"/*.ntuple.root")
        save_physics_list(PHYSDIR,physics_roots)
        print("PHYS_ROOTS: " + str(physics_roots))
       #Bifurcator = bi.Bifurcator(rootfiles=physics_roots,config_dict=config_dict,
        #        save_directory=RESULTDIR)
        #Bifurcator.bifurcate()
        #Bifurcator.SaveBifurcationRoot()

