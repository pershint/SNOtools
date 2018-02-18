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
import lib.configparser as cp
import os,sys
import glob
import json

#FIXME: set up a simpler argparser here.  you can choose a config file, override
#some of the configuration values with flags (changes dictionary values), and
#Could have flags to only run sacrifice, bifurcation, bifurcation analysis, or all.
CONFIGFILE='cuts_default.json'
JOBNUM=0
SOURCE="N16"
ERANGE=None
DEBUG=True

MAINDIR = os.path.dirname(__file__)
SAVEDIR = os.path.abspath(os.path.join(MAINDIR, "output","results_j"+str(JOBNUM)))
if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)
DBDIR = os.path.abspath(os.path.join(MAINDIR, "DB"))
CALIBDIR = os.path.abspath(os.path.join(MAINDIR, "ntuples", "N16"))
PHYSDIR = os.path.abspath(os.path.join(MAINDIR, "ntuples", "OpenGolden"))
CONFIGDIR = os.path.abspath(os.path.join(MAINDIR, 'config'))

def save_physics_list(directory,fullpath):
    #Strips off the directory, places the root names in a list, and saves it
    #in the SAVEDIR
    physics_roots = []
    for fullname in fullpath:
        rootname = fullname.lstrip(fullpath+"/")
        physics_roots.append(rootname)
    analysis_list = {}
    analysis_list["runs_used_in_bifurcation"] = physics_roots
    with open(SAVEDIR+"/physics_list.json","w") as f:
        json.dump(analysis_list,f,sort_keys=True,indent=4)

if __name__ == '__main__':
    #Get your run files to load
    N16_roots = glob.glob(CALIBDIR+"/*.ntuple.root")
    physics_roots = glob.glob(PHYSDIR+"/*.ntuple.root")
    save_physics_list(PHYSDIR,physics_roots)
    print("N16_ROOTS: " + str(N16_roots))
    print("PHYS_ROOTS: " + str(physics_roots))
    #Load the configuration file to use
    ConfigParser = cp.ConfigParser(CONFIGDIR+"/"+CONFIGFILE)
    config_dict = ConfigParser.Load_JsonConfig()
    if ERANGE is not None:
        config_dict["E_low"] = ERANGE[0]
        config_dict["E_high"] = ERANGE[1]
    ConfigParser.SaveConfiguration(config_dict,SAVEDIR,"used_configuration.json")

    #Generate Histograms of fractional sacrifice for Calibration data
    #These histograms are used to output a sacrifice and uncertainty below
    SacHists = sh.SacrificeHistGen(rootfiles=N16_roots,config_dict=config_dict,
            sourcetype=SOURCE)
    SacHists.GenerateHistograms()
    SacHists.SaveHistograms(SAVEDIR)

    SacSysUnc = sa.SacrificeSystematics(Sacrifice_Histograms=SacHists,\
            config_dict=config_dict)
    SacSysUnc.LoadCalibrationPositions(DBDIR)
    SacSysUnc.CalculateAcceptances()
    if DEBUG is True:
        SacSysUnc.ShowAcceptanceResults()
    SacSysUnc.SaveAcceptanceByRun(SAVEDIR,"cut_acceptances_byrun.json")
    SacSysUnc.SaveAcceptanceSummary(SAVEDIR,"cut_acceptances_total.json")

    #Run bifurcation analysis on Physics files
    #Bifurcator = bi.Bifurcator(rootfiles=physics_roots,config_dict=config_dict,
    #        save_directory=SAVEDIR)
    #Bifurcator.bifurcate()
    #Bifurcator.SaveBifurcationRoot()

