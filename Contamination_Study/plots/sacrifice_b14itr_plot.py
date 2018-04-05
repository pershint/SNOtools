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

def SacVSClass_Data(rootfiles=[],var="nhits",nbins=10,xmin=10.0,xmax=100.0,dcmask=None):
    '''
    Takes in a rootfile and returns a PandaFrame object that can be used
    for plotting in matplotlib.  Returns var vs. fractional sacrifice for
    events w/ b14<b14_low, b14>b14_high, and itr < itr_low.
    '''
    data = ROOT.TChain("output")
    for rf in rootfiles:
        data.Add(rf)
    #Now, make a new dictionary object: key is cutname, value is histogram
    if dcmask is not None:
        data.Draw("%s>>h_allevents(%i,%f,%f)"% (var,nbins,xmin,xmax),
                "fitValid==1&&isCal==1&&((dcFlagged&%i)!=%i)" % (dcmask,dcmask),"goff")
    else:
        data.Draw("%s>>h_allevents(%i,%f,%f)"% (var,nbins,xmin,xmax),
                "fitValid==1&&isCal==1","goff")
    h_allevents = gDirectory.Get("h_allevents")
    h_allevents.Sumw2()
    cutnames = ()
    allclasssacs = {}
    labeldict = ["b14_low", "b14_high","itr_low","total"]
    for cut in labeldict:
        graphdict={}
        vardat, fs,fs_unc =(), (), () #pandas wants ntuples
        h_cut_FracFlagged = ROOT.TH1D("h_cut_FracFlagged", "h_cut_FracFlagged", nbins,xmin,xmax)
        h_cut_FracFlagged.Sumw2()
        if cut == "b14_low":
            data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "energy>5.5&&fitValid==1&&isCal==1&&beta14<-0.12","goff")
        if cut == "b14_high":
            data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "energy>5.5&&fitValid==1&&isCal==1&&beta14>0.95" ,"goff")
        if cut == "itr_low":
            data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "energy>5.5&&fitValid==1&&isCal==1&&itr<0.55" ,"goff")
        if cut == "total":
            data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "energy>5.5&&fitValid==1&&isCal==1&&(beta14<-0.12 || beta14>0.95 || itr<0.55)" ,"goff")
        print("CUT: " + str(cut))
        h_flagged = gDirectory.Get("h_flagged")
        h_flagged.Sumw2()
        h_cut_FracFlagged.Divide(h_flagged,h_allevents,1.,1.,"b")
        for i in xrange(int(h_cut_FracFlagged.GetNbinsX())):
            vardat =  vardat + (float(h_cut_FracFlagged.GetBinWidth(i)) + float(h_cut_FracFlagged.GetBinLowEdge(i)),)
            fs = fs + (h_cut_FracFlagged.GetBinContent(i),)
            fs_unc = fs_unc + (h_cut_FracFlagged.GetBinError(i),)
        graphdict["vardat"] = vardat
        graphdict["fractional_sacrifice"] = fs
        graphdict["fs_uncertainty"] = fs_unc
        allclasssacs[cut]= pandas.DataFrame(data=graphdict)
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
    plt.title("Fractional sacrifice of Tagged N16 Events with Valid Fit\n"+\
            "due to classifiers, November 2017 scan",fontsize=36)
    plt.ion()
    plt.show()
