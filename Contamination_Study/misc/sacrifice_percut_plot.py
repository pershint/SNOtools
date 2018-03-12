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

def SacVSNhit_Data(rootfile,nbins=10,nmin=10.0,nmax=100.0,dcmask=None):
    '''
    Takes in a rootfile and returns a PandaFrame object that can be used
    for plotting in matplotlib.  Returns nhits vs. fractional sacrifice for
    Each cut in the given dcmask.
    '''
    if dcmask == None:
        print("No DC mask given.  Please give a DC mask of cuts to plot")
        return

    rf = ROOT.TFile(rootfile,"READ")
    data = rf.Get("output")
    fullmask = mb.get_dcwords()
    plotmask = {}
    for cut in fullmask:
        if 'prescale' in fullmask[cut]:
            continue
        if 'waterblind' in fullmask[cut]:
            continue
        if cut&dcmask==cut:
            plotmask[cut] = fullmask[cut]
    #Now, make a new dictionary object: key is cutname, value is histogram
    data.Draw("nhits>>h_allevents(%i,%f,%f)"% (nbins,nmin,nmax),
            "","goff")
    h_allevents = gDirectory.Get("h_allevents")
    h_allevents.Sumw2()
    cutnames = ()
    allcutsacs = {}
    for cut in plotmask:
        graphdict={}
        nhit, fs,fs_unc =(), (), () #pandas wants ntuples
        h_cut_FracFlagged = ROOT.TH1D("h_cut_FracFlagged", "h_cut_FracFlagged", nbins,nmin,nmax)
        h_cut_FracFlagged.Sumw2()
        data.Draw("nhits>>h_flagged(%i,%f,%f)" % (nbins,nmin,nmax),
                "((dcFlagged&%i)!=%i)" % (cut,cut),"goff")
        h_flagged = gDirectory.Get("h_flagged")
        h_flagged.Sumw2()
        h_cut_FracFlagged.Divide(h_flagged,h_allevents,1.,1.,"b")
        for i in xrange(h_cut_FracFlagged.GetNbinsX()):
            nhit =  nhit + (int(h_cut_FracFlagged.GetBinWidth(i) + h_cut_FracFlagged.GetBinLowEdge(i)),)
            fs = fs + (h_cut_FracFlagged.GetBinContent(i),)
            fs_unc = fs_unc + (h_cut_FracFlagged.GetBinError(i),)
        graphdict["nhit"] = nhit
        graphdict["fractional_sacrifice"] = fs
        graphdict["fs_uncertainty"] = fs_unc
        allcutsacs[plotmask[cut]]= pandas.DataFrame(data=graphdict)
        del h_cut_FracFlagged
        del h_flagged
    #graphdict has the sacrifice information for each cut. Now, let's plot it.
   
    return allcutsacs 

def GetTopSacs(allcutsacs,topnumber=5):
    if len(allcutsacs) <= topnumber:
        print("You already only have five cuts.  Not combining any")
        return
    names = []
    sacrifices = []
    for cut in allcutsacs:
        sacrifices.append(np.average(allcutsacs[cut]["fractional_sacrifice"]))
        names.append(cut)
    sortednames = [x for _,x in sorted(zip(sacrifices,names))]
    topnames = sortednames[len(sortednames)-topnumber:len(sortednames)]
    return topnames
     

def SacVSNHit_Plot(allcutsacs,topcuts=None):
    sns.set_style("whitegrid")
    sns.axes_style("whitegrid")
    xkcd_colors = ['slate blue', 'fluro green', 'twilight', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','yellowgreen']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(allcutsacs)))
    if topcuts is None:
        for cut in allcutsacs:
            plt.errorbar(x=allcutsacs[cut].nhit, 
                    y=allcutsacs[cut].fEractional_sacrifice,
                    yerr=allcutsacs[cut].fs_uncertainty,
                    linestyle='none', marker='+', label=cut, markersize=8,
                    elinewidth=2, capsize=0)
    else:
        fracsum,fracuncsum=[],[]
        for cut in allcutsacs:
            if cut in topcuts:
                plt.errorbar(x=allcutsacs[cut].nhit, y=allcutsacs[cut].fractional_sacrifice,
                        yerr=allcutsacs[cut].fs_uncertainty, linestyle='none',
                        marker='o', label=cut, capsize=0, elinewidth=2, markersize=5)
            else:
                fracsum.append(allcutsacs[cut].fractional_sacrifice)
                fracuncsum.append(allcutsacs[cut].fs_uncertainty)
        plt.errorbar(x=allcutsacs[cut].nhit,y=sum(fracsum),yerr=sum(fracuncsum),
                linestyle='none',marker='o', capsize=0, elinewidth=2, label='All other cuts', markersize=5)
    plt.legend(loc=3)
    plt.yscale("log")
    plt.ylabel("Fractional sacrifice",fontsize=22)
    plt.xlabel("NHit",fontsize=22)
    plt.tick_params(labelsize=20)
    plt.title("Fractional sacrifice of dataset by data cleaning cuts",fontsize=24)
    plt.ion()
    plt.show()
