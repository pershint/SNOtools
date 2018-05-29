from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import maskbuilder as mb
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.optimize as spc
import pandas
import time

def GetAnalysisMask():
    util = RAT.DU.Utility.Get()
    util.LoadDBAndBeginRun()
    dcmask = RAT.GetDataCleaningWord('analysis_mask')
    return dcmask

def flatline(x,b):
    return b

def weighted_stdev(vals,valavg, uncs):
    weights = 1/(uncs**2)
    return np.sqrt((len(weights)*np.sum(weights*(vals-valavg)**2))/((len(weights)-1)*np.sum(weights)))

def PrepData_ClassSac(rootfiles=[],var="nhits",nbins=10,xmin=10.0,xmax=100.0,dcmask=None):
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
                "fitValid==1&&isCal==1&&energy>5.5&&energy<9.0&&posr<5500&&((dcFlagged&%i)!=%i)" % (dcmask,dcmask),"goff")
    else:
        data.Draw("%s>>h_allevents(%i,%f,%f)"% (var,nbins,xmin,xmax),
                "fitValid==1&&isCal==1&&energy>5.5&&energy<9.0&&posr<5500","goff")
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
                    "energy>5.5&&energy<9.0&&posr<5500&&fitValid==1&&isCal==1&&beta14<-0.12","goff")
        if cut == "b14_high":
            data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "energy>5.5&&energy<9.0&&posr<5500&&fitValid==1&&isCal==1&&beta14>0.95" ,"goff")
        if cut == "itr_low":
            data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "energy>5.5&&energy<9.0&&posr<5500&&fitValid==1&&isCal==1&&itr<0.55" ,"goff")
        if cut == "total":
            data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "energy>5.5&&energy<9.0&&posr<5500&&fitValid==1&&isCal==1&&(beta14<-0.12 || beta14>0.95 || itr<0.55)" ,"goff")
        print("CUT: " + str(cut))
        h_flagged = gDirectory.Get("h_flagged")
        h_flagged.Sumw2()
        h_cut_FracFlagged.Divide(h_flagged,h_allevents,1.,1.,"b")
        for i in xrange(int(h_cut_FracFlagged.GetNbinsX())):
            if i==0:
                continue
            vardat =  vardat + (float(h_cut_FracFlagged.GetBinWidth(i)) + float(h_cut_FracFlagged.GetBinLowEdge(i)),)
            fs = fs + (h_cut_FracFlagged.GetBinContent(i),)
            fs_unc = fs_unc + (h_cut_FracFlagged.GetBinError(i),)
        graphdict["vardat"] = vardat
        graphdict["fractional_sacrifice"] = fs
        graphdict["fs_uncertainty"] = fs_unc
        allclasssacs[cut]= pandas.DataFrame(data=graphdict)
        del h_cut_FracFlagged
        del h_flagged

    meta = {"binwidth":((xmax-xmin)/float(nbins)),"variable":var}
    #graphdict has the sacrifice information for each cut. Now, let's plot it.
    return allclasssacs,meta 



def Plot_ClassSac(allclasssacs,meta={"binwidth":10.0, "variable":"nhits"},
        fittotal=True):
    sns.set_style("whitegrid")
    xkcd_colors = ['slate blue', 'black', 'brown', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','twilight']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(allclasssacs)))
    for cut in allclasssacs:
        plt.errorbar(x=allclasssacs[cut].vardat, 
                y=allclasssacs[cut].fractional_sacrifice,
                yerr=allclasssacs[cut].fs_uncertainty,
                linestyle='none', marker='o', label=cut, markersize=6,
                elinewidth=3, capsize=0)
    if fittotal is True:
        popt, pcov = spc.curve_fit(flatline, allclasssacs["total"].vardat, allclasssacs["total"].fractional_sacrifice,
                p0=[0.02], sigma=allclasssacs["total"].fs_uncertainty)
        #one standard deviation
        print("BEST FIT: " + str(popt))
        print("PCOVARIANCE: " + str(pcov))
        stdev = np.sqrt(np.diag(pcov))
        plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
        print("WEIGHTED STANDARD DEVIATION FROM FLAT: " + str(weighted_stdev(allclasssacs["total"].fractional_sacrifice,
            float(popt[0]),allclasssacs["total"].fs_uncertainty)))
    plt.legend(loc=3)
    plt.yscale("log")
    plt.ylabel("Fractional sacrifice",fontsize=34)
    plt.xlabel(meta["variable"],fontsize=34)
    plt.tick_params(labelsize=32)
    plt.title("Fractional sacrifice of Tagged N16 Events \n"+\
            "due to classifiers, November 2017 scan",fontsize=36)
    plt.ion()
    plt.show()
