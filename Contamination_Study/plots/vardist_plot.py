from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import maskbuilder as mb
import seaborn as sns
import matplotlib.pyplot as plt
import pandas
import time

def EVDist_Data(rootfiles=[],var="nhits",nbins=10,xmin=10.0,xmax=100.0,dcmask=None):
    '''
    Takes in a rootfile and returns a PandaFrame object that can be used
    for plotting in matplotlib.  Returns var's event distribution for
    events in a dataset.
    '''
    data = ROOT.TChain("output")
    for rf in rootfiles:
        data.Add(rf)
    #Now, make a new dictionary object: key is cutname, value is histogram
    if dcmask is not None:
        data.Draw("%s>>h_vardist(%i,%f,%f)"% (var,nbins,xmin,xmax),
                "fitValid==1&&isCal==1&&((dcFlagged&%i)!=%i)" % (dcmask,dcmask),"goff")
    else:
        data.Draw("%s>>h_vardist(%i,%f,%f)"% (var,nbins,xmin,xmax),
                "fitValid==1&&isCal==1","goff")
    h_vardist = gDirectory.Get("h_vardist")
    h_vardist.Sumw2()
    print("GETENTRIESRESULT: " + str(h_vardist.GetEntries()))
    h_vardist.Scale(1.0/float(h_vardist.GetEntries()))
    graphdict = {}
    vardat, evdist, evdist_err = (), (), ()
    for i in xrange(h_vardist.GetNbinsX()):
        vardat =  vardat + (float(h_vardist.GetBinWidth(i) + h_vardist.GetBinLowEdge(i)),)
        evdist = evdist + (h_vardist.GetBinContent(i),)
        evdist_err = evdist_err + (h_vardist.GetBinError(i),)
    graphdict[var] = vardat
    graphdict["evdistribution"] = evdist
    graphdict["evdist_unc"] = evdist_err
    distdict= pandas.DataFrame(data=graphdict)
    distdict.var = var
    distdict.xwidth = (xmax-xmin)/(2.*float(nbins))
    del h_vardist
   
    return distdict 

def EVDist_Plot(distdict):
    sns.set_style("whitegrid")
    sns.axes_style("whitegrid")
    xkcd_colors = ['purple', 'fluro green', 'brown', 'slate blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','twilight']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(distdict)))
    plt.errorbar(x=getattr(distdict,distdict.var), 
        y=distdict.evdistribution,
        yerr=distdict.evdist_unc,xerr=distdict.xwidth,
        linestyle='none', marker='o', label="N16 tagged events", markersize=2,
        elinewidth=2, capsize=0)
    plt.legend(loc=3,fontsize=30)
    plt.yscale("log")
    #plt.ylabel("Fractional sacrifice",fontsize=22)
    plt.xlabel(distdict.var,fontsize=30)
    plt.tick_params(labelsize=28)
    plt.title("Event distribution of Tagged N16 Events with Valid Fit\n"+\
            "November 2017 scan",fontsize=32)
    plt.ion()
    plt.show()
