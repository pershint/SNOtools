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
LOWECONTAM=args.LOWECONTAM
ESTIMATECONTAMINATION=args.ESTIMATECONTAMINATION
PLOTS=args.PLOTS
ERANGE=args.ERANGE
ZRANGE=args.ZRANGE
JOBNUM=args.JOBNUM
CALIBDIR=args.CALIBDIR
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

RESULTDIR = os.path.abspath(os.path.join(MAINDIR, "output","results_j"+str(JOBNUM)))
if not os.path.exists(RESULTDIR):
    os.makedirs(RESULTDIR)
DBDIR = os.path.abspath(os.path.join(MAINDIR, "DB"))

if __name__ == '__main__':
    #Get your run files to load
    #Load the setup JSON defining what plots/variables to use
    #for sacrifice estimate
    with open("./config/setup.json","r") as f:
        setup_dict = json.load(f)
    if CALIBSACANALYSIS or BIFURCATE is True:
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

    if CALIBSACANALYSIS is True:
        calib_data = glob.glob("%s/*.ntuple.root"%(CALIBDIR))
        print("NUMBER OF N16 FILES: " + str(len(calib_data)))
        if DEBUG is True:
            print("N16_ROOTS: " + str(calib_data))
        ru.save_sacrifice_list(RESULTDIR,calib_data,"calibration_data_used.json") 
        CUTS_TODO = setup_dict['CUTS_TODO']
        SACRIFICE_VARIABLES = setup_dict['SACRIFICE_VARIABLES']
        sacrifice_results = {} 
        for cut in CUTS_TODO:
            sacrifice_results[cut] = {} 
        for variable in SACRIFICE_VARIABLES:
            for cut in sacrifice_results:
                if cut == 'cut1': 
                    DCSacs = sa.DCSacrificeAnalyzer(rootfiles=calib_data, cuts_dict=config_dict)
                    DCSacs.AnalyzeData(var=variable,nbins=9)
                    sacrifice_results['cut1'][variable] = DCSacs.GetFitTotalAndUncertainties() 
                    if PLOTS is True:
                        tit="Fractional sacrifice due to data cleaning\n "+\
                                "internal ROI, Nov 2017 N16 scan"
                        DCSacs.ShowPlottedSacrifice(title=tit)
                elif cut == 'cut2': 
                    ClassSacs = sa.ClassSacrificeAnalyzer(rootfiles=calib_data, cuts_dict=config_dict)
                    ClassSacs.AnalyzeData(var=variable,nbins=9)
                    sacrifice_results['cut2'][variable] = ClassSacs.GetFitTotalAndUncertainties() 
                    if PLOTS is True: 
                        tit="Fractional sacrifice due to classifiers\n "+\
                                "internal ROI, March 2018 N16 external"
                        ClassSacs.ShowPlottedSacrifice(title=tit)
                else:
                    print("Cut type not supported.  Please use only cut1 and/or cut2"+\
                            " in your setup file.")
        if NOSAVE is False:
            with open("%s/calib_cut_sacrifices_total.json"%(RESULTDIR),"w") as f:
                json.dump(sacrifice_results,f,sort_keys=True,indent=4)


    #Run bifurcation analysis on Physics files
    if BIFURCATE is True:
        physics_roots = glob.glob(ANALYSISDIR+"*.ntuple.root")
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
            cut_sac_summary = ru.LoadJson(RESULTDIR,"calib_cut_sacrifices_total.json")
        except IOError:
            print("Bifurcation or Cut Sacrifice Summary loading error.  Were these"+\
                    "analyses run?")
            raise
        if LOWECONTAM is True:
            print("CALCULATING CONTAMINATION ASSUMING LOW ENERGY REGION")
            CE = ca.LowEContamination(bifurcation_summary,cut_sac_summary)
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
