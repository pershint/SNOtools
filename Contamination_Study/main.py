#This program is a one-shot program that outputs the results of the SNO+
#Conamination study.  First, the data cleanning and fit sacrifice are estimated
#Using the given calibration files.  Then, the bifurcation analysis is run over
#the given data physics files.  Finally, an estimate on the contamination for the
#chosen configuration of cuts/ROI choice is output.
import matplotlib.pyplot as plt
import lib.ArgParser as ap
import numpy as np
import os,sys
import glob
import json

args = ap.args

DEBUG = args.debug
NOSAVE=args.NOSAVE
CALIBSACANALYSIS=args.CALIBSACANALYSIS
MCSACANALYSIS=args.MCSACANALYSIS
BIFURCATE=args.BIFURCATE
LOWECONTAM=args.LOWECONTAM
ESTIMATECONTAMINATION=args.ESTIMATECONTAMINATION
ERANGE=args.ERANGE
ZRANGE=args.ZRANGE
JOBNUM=args.JOBNUM
CALIBDIR=args.CALIBDIR
PHYSDIR=args.PHYSDIR
MCSIGNALDIR=args.MCSIGNALDIR

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

basepath = os.path.dirname(__file__)
MAINDIR = os.path.dirname(__file__)
CONFIGDIR = os.path.abspath(os.path.join(MAINDIR,"config"))

RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "output","results_j"+str(JOBNUM)))
if not os.path.exists(RESULTDIR):
    os.makedirs(RESULTDIR)
DBDIR = os.path.abspath(os.path.join(MAINDIR, "DB"))

if __name__ == '__main__':
    #Get your run files to load
    #Load the configuration file to use
    if MCSACANALYSIS or CALIBSACANALYSIS or BIFURCATE is True:
        ConfigParser = cp.ConfigParser(CONFIGDIR)
        print("CONFIGDIR: " + str(CONFIGDIR))
        config_dict = ConfigParser.Load_JsonConfig('cuts_default.json')
        setup_dict = ConfigParser.Load_JsonConfig('setup.json')
        if ZRANGE is not None:
            config_dict["Z_high"] = float(ZRANGE[0])
            config_dict["Z_low"] = float(ZRANGE[1])
        if ERANGE is not None:
            config_dict["E_low"] = float(ERANGE[0])
            config_dict["E_high"] = float(ERANGE[1])
        if NOSAVE is False:
            ConfigParser.SaveConfiguration(config_dict,RESULTDIR,"used_cutsconfig.json")

    if MCSACANALYSIS is True:
        mc_data_all = glob.glob("%s/*.ntuple.root"%(MCSIGNALDIR))
        print("NUMBER OF SIGNAL MC FILES: " + str(len(mc_data_all)))
        if DEBUG is True:
            print("MONTECARLO_ROOTS: " + str(mc_data_all))
        ru.save_sacrifice_list(RESULTDIR, mc_data_all,"sacestimate_mcfiles.json")
        SacHists = sh.SacrificeHistGen(rootfiles=mc_data_all,config_dict=config_dict)
        SacHists.GenerateHistograms()
        if NOSAVE is False:
            SacHists.SaveHistograms(savedir="%s/MC_sachists"%(RESULTDIR))
    
        SacSysUnc = sa.SacrificeHistAnalyzer(Sacrifice_Histograms=SacHists,\
                config_dict=config_dict,source="Monte Carlo")
        SacSysUnc.LoadCalibrationPositions(DBDIR)
        SacSysUnc.CalculateSacrifices()
        if DEBUG is True:
            SacSysUnc.ShowSacrificeResults()
        if NOSAVE is False:
            SacSysUnc.SaveSacrificeByRun(RESULTDIR,"MC_cut_sacrifices_byrun.json")
            SacSysUnc.SaveSacrificeSummary(RESULTDIR,"MC_cut_sacrifices_total.json")

    if CALIBSACANALYSIS is True:
        CALIBSOURCE=setup_dict['CALIBSOURCE']
        calib_data_all = glob.glob("%s/%s/*.ntuple.root"%(CALIBDIR,SOURCE))
        print("NUMBER OF N16 FILES: " + str(len(calib_data_all)))
        if DEBUG is True:
            print("N16_ROOTS: " + str(calib_data_all))
        if ZRANGE is not None:
            calib_data = cs.ApplyZCut(DBDIR,SOURCE,ZRANGE,calib_data_all)
        print("NUMBER OF FILES AFTER ZCUT: " + str(len(calib_data)))
        ru.save_sacrifice_list(RESULTDIR, calib_data,'sacestimate_calibfiles.json')
        SacHists = sh.SacrificeHistGen(rootfiles=calib_data,config_dict=config_dict,
                source=CALIBSOURCE)
        SacHists.GenerateHistograms()
        if NOSAVE is False:
            SacHists.SaveHistograms(RESULTDIR)
    
        SacSysUnc = sa.SacrificeHistAnalyzer(Sacrifice_Histograms=SacHists,\
                config_dict=config_dict)
        SacSysUnc.LoadCalibrationPositions(DBDIR)
        SacSysUnc.CalculateSacrifices()
        if DEBUG is True:
            SacSysUnc.ShowSacrificeResults()
        if NOSAVE is False:
            SacSysUnc.SaveSacrificeByRun(RESULTDIR,"calib_cut_sacrifices_byrun.json")
            SacSysUnc.SaveSacrificeSummary(RESULTDIR,"calib_cut_sacrifices_total.json")

#    if PLOTS is True:
#       cut_sacrifices = ru.LoadJson(RESULTDIR,"cut_sacrifices_byrun.json")
#       cut_sac_summary = ru.LoadJson(RESULTDIR,"cut_sacrifices_total.json")
#       for axis in ['x','y','z']:
#           sp.plot_sacrificevsCart(cut_sacrifices, 'cut1', axis)
#           sp.plot_sacrificevsCart(cut_sacrifices, 'cut2', axis)
#       for cut in ['cut1','cut2']:
#           sp.plot_sacrificevsR(cut_sacrifices, cut)
#           #sp.plot_XYSacrifice(cut_sacrifices, cut)

    #Run bifurcation analysis on Physics files
    if BIFURCATE is True:
        physics_roots = glob.glob(PHYSDIR+"/*.ntuple.root")
        ru.save_bifurcation_list(RESULTDIR,physics_roots,'bifurcation_analysisfiles.json')
        if DEBUG is True:
            print("PHYS_ROOTS: " + str(physics_roots))
        Bifurcator = bi.Bifurcator(rootfiles=physics_roots,config_dict=config_dict)
        Bifurcator.Bifurcate()
        if NOSAVE is False:
            Bifurcator.SaveBifurcationSummary(RESULTDIR,"bifurcation_boxes.json")
     
    if ESTIMATECONTAMINATION is True:
        bifurcation_summary = None
        cut_sac_summary = None
        signal_estimate_summary = None
        try:
            bifurcation_summary = ru.LoadJson(RESULTDIR,"bifurcation_boxes.json")
            cut_sac_summary = ru.LoadJson(RESULTDIR,"cut_sacrifices_total.json")
        except IOError:
            print("Bifurcation or Cut Sacrifice Summary loading error.  Were these"+\
                    "analyses run?")
            raise
        if LOWECONTAM is True:
            print("CALCULATING CONTAMINATION ASSUMING LOW ENERGY REGION")
            try:
                signal_estimate_summary = ru.LoadJson(RESULTDIR,"signal_estimate.json")
            except IOError:
                print("Error loading signal estimate file... does it exist?")
                raise
            CE = ca.LowEContamination(bifurcation_summary,cut_sac_summary,\
                    signal_estimate_summary)
            contam = CE.CalculateContamination()
            contam_unc = CE.CalculateContaminationUnc()
            if DEBUG is True:
                print("Estimated contamination: %f"%(contam))
                print("Estimate's uncertainty: %f"%(contam_unc))
        else:
            CE = ca.NDContamination(bifurcation_summary,cut_sac_summary)
            CE.CalculateContaminationValues() #Calculate contamination eqns.
            values = CE.BootstrapCL(0.90,100000) #Estimate upper end of y1y2
            #if PLOTS is True:
            #    values=values*CE.contamination_summary['est_bkg_evts']
            #    plt.hist(values,100,range=(min(values),max(values)))
            #    plt.xlabel(r"Total estimated contamination (y1y2$\beta$)")
            #    plt.ylabel(r"Relative probability (unitless)")
            #    plt.title("Distribution of estimated contamination after\n"+\
            #        "re-firing variables with statistical uncertainties")
            #    plt.show()
        if NOSAVE is False:
            CE.SaveContaminationSummary(RESULTDIR,"contamination_summary.json")
