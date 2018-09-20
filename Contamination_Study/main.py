#This program is a one-shot program that outputs the results of the SNO+
#Conamination study.  First, the data cleanning and fit sacrifice are estimated
#Using the given calibration files.  Then, the bifurcation analysis is run over
#the given data physics files.  Finally, an estimate on the contamination for the
#chosen configuration of cuts/ROI choice is output.
import lib.ArgParser as ap

args = ap.args

DEBUG = args.debug
NOSAVE=args.NOSAVE
CONFIGFILE=args.CONFIGFILE
CALIBSACANALYSIS=args.CALIBSACANALYSIS
BIFURCATE=args.BIFURCATE
LETACONTAM=args.LETACONTAM
ESTIMATECONTAMINATION=args.ESTIMATECONTAMINATION
PLOTS=args.PLOTS
ERANGE=args.ERANGE
ZRANGE=args.ZRANGE
JOBNAME=args.JOBNAME
CALIBDIR=args.CALIBDIR
CALIBMCDIR=args.CALIBMCDIR
ANALYSISDIR=args.ANALYSISDIR

import numpy as np
import os,sys,time
import glob
import json
import matplotlib.pyplot as plt
import ROOT
import lib.Bifurcator as bi
import plots.BifurPlots as bp
import lib.DCClassSacrifice as sa
import lib.ContaminationAnalyzer as ca
import lib.ConfigParser as cp
import lib.ResultUtils as ru

basepath = os.path.dirname(__file__)
MAINDIR = os.path.dirname(__file__)
CONFIGDIR = os.path.abspath(os.path.join(MAINDIR,"config"))

results_exist = True 
RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "output","results_j"+str(JOBNAME)))
if not os.path.exists(RESULTDIR):
    results_exist = False
    os.makedirs(RESULTDIR)
if not os.path.exists(RESULTDIR+"/plots"):
    os.makedirs(RESULTDIR+"/plots")
DBDIR = os.path.abspath(os.path.join(MAINDIR, "DB"))

if __name__ == '__main__':
    #Get your run files to load
    #Load the setup JSON defining what plots/variables to use
    #for sacrifice estimate
    with open("./config/setup.json","r") as f:
        setup_dict = json.load(f)
    if CALIBSACANALYSIS or BIFURCATE is True:
        if results_exist is False:
            ConfigParser = cp.ConfigParser(CONFIGDIR)
            print("CONFIGDIR: " + str(CONFIGDIR))
            config_dict = ConfigParser.Load_JsonConfig(CONFIGFILE)
            if ZRANGE is not None:
                config_dict["Z_high"] = float(ZRANGE[0])
                config_dict["Z_low"] = float(ZRANGE[1])
            if ERANGE is not None:
                config_dict["E_low"] = float(ERANGE[0])
                config_dict["E_high"] = float(ERANGE[1])
            if NOSAVE is False:
                ConfigParser.SaveConfiguration(config_dict,RESULTDIR,"used_cutsconfig.json")
        else:
            config_dict = ru.LoadJson(RESULTDIR,"used_cutsconfig.json")

    if CALIBSACANALYSIS is True:
        calib_data = glob.glob("%s/*.ntuple.root"%(CALIBDIR))
        calib_mc = glob.glob("%s/*.ntuple.root"%(CALIBMCDIR))
        print("NUMBER OF N16 FILES: " + str(len(calib_data)))
        print("NUMBER OF N16 MC FILES: " + str(len(calib_mc)))
        if DEBUG is True:
            print("N16_ROOTS: " + str(calib_data))
        ru.save_sacrifice_list(RESULTDIR,calib_data,"calibration_data_used.json")
        ru.save_sacrifice_list(RESULTDIR,calib_mc,"calibration_MCdata_used.json") 
        CUTS_TODO = setup_dict['CUTS_TODO']
        SACRIFICE_VARIABLES = setup_dict['SACRIFICE_VARIABLES']
        SacComp_results = {} 
        try:
            SacComp_results = ru.LoadJson(RESULTDIR,"calib_DCSacClassComp_totals.json")
        except IOError:
            print("No Sacrifice estimates done previously. Initializing empty dict")
            pass 
        for cut in CUTS_TODO:
            if cut not in SacComp_results.keys():
                SacComp_results[cut] = {} 
        for variable in SACRIFICE_VARIABLES:
            for cut in CUTS_TODO:
                if cut == 'cut1': 
                    DCSacs = sa.DCSacrificeAnalyzer(rootfiles_data=calib_data, cuts_dict=config_dict)
                    DCSacs.SetBinNumber(setup_dict["BINNUM_CUT1"])
                    DCSacs.AnalyzeData(var=variable)
                    SacComp_results['cut1'][variable] = DCSacs.GetFitTotalAndUncertainties() 
                    if PLOTS is True:
                        DCSacs.ShowPlottedSacrifice(title=setup_dict["PLOT_TITLE_CUT1"],
                                savedir="%s/%s"%(RESULTDIR,"plots"))
                elif cut == 'cut2':
                    ClassSacs = sa.ClassSacrificeAnalyzer(rootfiles_data=calib_data,
                             cuts_dict=config_dict)
                    ClassSacs.SetBinNumber(setup_dict["BINNUM_CUT2"])
                    ClassSacs.AnalyzeData(var=variable)
                    SacComp_results['cut2'][variable] = ClassSacs.GetFitTotalAndUncertainties()
                elif cut == 'cut2_DataMCComp':
                    ClassComps = sa.DataMCClassAnalyzer(rootfiles_data=calib_data,
                            rootfiles_mc=calib_mc, cuts_dict=config_dict)
                    ClassComps.SetBinNumber(setup_dict["BINNUM_CUT2"])
                    ClassComps.AnalyzeData(var=variable)
                    SacComp_results['cut2_DataMCComp'][variable] = ClassComps.GetFitTotalAndUncertainties() 
                    if variable == "energy":
                        ClassComps.SaveRatioToCSV(RESULTDIR,"RatioVsEnergy.csv")
                    if PLOTS is True: 
                        ClassComps.PlotRatio(title=setup_dict["PLOT_TITLE_CUT2"],
                                savedir="%s/%s"%(RESULTDIR,"plots"))
                else:
                    print("Cut type not supported.  Please use only cut1 and/or cut2"+\
                            " in your setup file.")
        if NOSAVE is False:
            with open("%s/calib_DCSacClassComp_totals.json"%(RESULTDIR),"w") as f:
                json.dump(SacComp_results,f,sort_keys=True,indent=4)


    #Run bifurcation analysis on Physics files
    if BIFURCATE is True:
        physics_roots = glob.glob(ANALYSISDIR+"*tuple.root")
        ru.save_bifurcation_list(RESULTDIR,physics_roots,'bifurcation_analysisfiles.json')
        if DEBUG is True:
            print("ANALYSIS_ROOTS: " + str(physics_roots))
        Bifurcator = bi.Bifurcator(rootfiles=physics_roots,config_dict=config_dict)
        Bifurcator.Bifurcate()
        if NOSAVE is False:
            Bifurcator.SaveBifurcationSummary(RESULTDIR,"bifurcation_boxes.json")
     
    if ESTIMATECONTAMINATION is True:
        bifurcation_summary = None
        cut_sac_summary = None
        try:
            bifurcation_summary = ru.LoadJson(RESULTDIR,"bifurcation_boxes.json")
        except IOError:
            print("Bifurcation Summary loading error.  Was this "+\
                    " analysis run?")
            pass
        try:
            cut_sac_summary = ru.LoadJson(RESULTDIR,"calib_DCSacClassComp_totals.json")
        except IOError:
            print("Cut sacrifice Summary loading error.  Was this "+\
                    " analysis run?")
            pass
        if LETACONTAM is True:
            print("CALCULATING CONTAMINATION ASSUMING LOW ENERGY REGION")
            CE = ca.LowEContamination(bifurcation_summary,cut_sac_summary)
            contam = CE.CalculateContamination()
            contam_unc = CE.CalculateContaminationUnc()
            if DEBUG is True:
                print("Estimated contamination: %f"%(contam))
                print("Estimate's uncertainty: %f"%(contam_unc))
            CE.SaveContaminationSummary(RESULTDIR,"contamination_summary_LETA.json")
        else:
            CE = ca.NDContamination(bifurcation_summary,cut_sac_summary)
            CE.CalculateContaminationValues() #Calculate contamination eqns.
            values = CE.BootstrapCL(0.683,100000) #Estimate upper end of y1y2
            if PLOTS is True:
                values=values
                plt.hist(values,100,range=(min(values),max(values)))
                plt.xlabel(r"Total estimated contamination (y1y2$\beta$)")
                plt.ylabel(r"Relative probability (unitless)")
                plt.title("Distribution of estimated contamination after\n"+\
                    "re-firing variables with statistical uncertainties")
                plt.show()
        if NOSAVE is False:
            CE.SaveContaminationSummary(RESULTDIR,"contamination_summary.json")
