from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import maskbuilder as mb
import seaborn as sns
import matplotlib.pyplot as plt
import pandas
import time

def GetAnalysisMask():
    util = RAT.DU.Utility.Get()
    util.LoadDBAndBeginRun()
    dcmask = RAT.GetDataCleaningWord('analysis_mask')
    return dcmask

def SacVSClass_Data(ambefiles=[],n16files=[],var="nhits",nbins=10,xmin=10.0,xmax=100.0,dcmask=None):
    '''
    Takes in a list of N16 and AmBe root files and plots the sacrifice versus
    variable chosen for N16, prompts, and delayeds.  '''
    n16data = ROOT.TChain("output")
    for rf in rootfiles:
        n16data.Add(rf)
    ambedata = ROOT.TChain("CombinedOutput")
    for af in ambefiles:
        ambedata.Add(af)
    #Now, make a new dictionary object: key is cutname, value is histogram
    n16data.Draw("%s>>h_alln16events(%i,%f,%f)"% (var,nbins,xmin,xmax),
            "fitValid==1&&posr<5000&&isCal==1","goff")
    ambedata.Draw("%s>>h_allambe_p(%i,%f,%f)"% (var,nbins,xmin,xmax),
            "fitValid_p==1&&posr_p<5000","goff")
    ambedata.Draw("%s>>h_allambe_d(%i,%f,%f)"% (var,nbins,xmin,xmax),
            "fitValid_d==1&&posr_d<5000","goff")
    h_alln16events = gDirectory.Get("h_alln16events")
    h_alln16events.Sumw2()
    h_allambe_p = gDirectory.Get("h_allambe_p")
    h_allambe_p.Sumw2()
    h_allambe_d = gDirectory.Get("h_allambe_d")
    h_allambe_d.Sumw2()
    
    allclasssacs = {}
    labeldict = ["AmBe_Prompt","AmBe_Delayed","N16"]
    for dat in labeldict:
        graphdict={}
        vardat, fs,fs_unc =(), (), () #pandas wants ntuples
        h_cut_FracFlagged = ROOT.TH1D("h_cut_FracFlagged", "h_cut_FracFlagged", nbins,xmin,xmax)
        h_cut_FracFlagged.Sumw2()
        if dat == "AmBe_Prompt":
            n16data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "fitValid_p==1&&posr_p<5000&&((dcFlagged_p&%s)!=%s)"%(dcmask,dcmask),"goff")
            h_all = h_allambe_p
        if dat == "AmBe_Delayed":
            ambedata.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "fitValid_d==1&&posr_d<5000&&((dcFlagged_d&%s)!=%s)"%(dcmask,dcmask),"goff")
            h_all = h_allambe_d
        if dat == "N16":
            ambedata.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "fitValid==1&&posr<5000&&isCal==1&&((dcFlagged&%s)!=%s)"%(dcmask,dcmask) ,"goff")
            h_all = h_alln16events
        print("CUT: " + str(dat))
        h_flagged = gDirectory.Get("h_flagged")
        h_flagged.Sumw2()
        h_cut_FracFlagged.Divide(h_flagged,h_all,1.,1.,"b")
        for i in xrange(int(h_cut_FracFlagged.GetNbinsX())):
            vardat =  vardat + (float(h_cut_FracFlagged.GetBinWidth(i)) + float(h_cut_FracFlagged.GetBinLowEdge(i)),)
            fs = fs + (h_cut_FracFlagged.GetBinContent(i),)
            fs_unc = fs_unc + (h_cut_FracFlagged.GetBinError(i),)
        graphdict["vardat"] = vardat
        graphdict["fractional_sacrifice"] = fs
        graphdict["fs_uncertainty"] = fs_unc
        allclasssacs[dat]= pandas.DataFrame(data=graphdict)
        del h_cut_FracFlagged
        del h_flagged
    #graphdict has the sacrifice information for each cut. Now, let's plot it.
    return allclasssacs,var 



def SacVSVar_Plot(allclasssacs,variable="nhits"):
    sns.set_style("whitegrid")
    sns.axes_style("whitegrid")
    xkcd_colors = ['slate blue', 'fluro green', 'brown', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','twilight']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(allclasssacs)))
    for cut in allclasssacs:
        plt.errorbar(x=allclasssacs[cut].vardat, 
                y=allclasssacs[cut].fractional_sacrifice,
                yerr=allclasssacs[cut].fs_uncertainty,
                linestyle='none', marker='o', label=cut, markersize=8,
                elinewidth=3, capsize=0)
    plt.legend(loc=3)
    plt.yscale("log")
    plt.ylabel("Fractional sacrifice",fontsize=34)
    plt.xlabel(variable,fontsize=34)
    plt.tick_params(labelsize=32)
    plt.title("Fractional sacrifice of central tagged N16 and prompt/delayed"+\
            "AmBe events with Valid Fit and radius<5 m",fontsize=36)
    plt.ion()
    plt.show()
