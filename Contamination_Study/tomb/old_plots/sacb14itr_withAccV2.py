from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import maskbuilder as mb
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
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


def DeleteEmptyBins(data):
    '''Deletes any bins with zero events and zero uncertainty (i.e. bins
    with zero events in them) from the self.sacrifices dataset'''
    total_uncertainties = data.fa_uncertainty
    zero_unc_entries = [] 
    zero_unc_entries = [j for j in xrange(len(total_uncertainties)) if
            float(total_uncertainties[j])==0]
    print("INDICES: " + str(zero_unc_entries))
    if len(zero_unc_entries)>0:
        print("WARNING: Graphing region is being trimmed; it seems you"+\
                " have some regions in your ROI with no events accepted in it.")
    data = data.drop(data.index[\
            zero_unc_entries])
    data = data.reset_index()
    return data

def DeleteEmptyBins_all(datdict):
    '''Deletes any bins from all data sets in array where the uncertainty is zero (this
    happens when the bin has no data in it at all)'''
    allzeros = [] 
    for data in datdict: 
        total_uncertainties = datdict[data].fa_uncertainty
        zero_unc_entries = []
        zero_unc_entries = [j for j in xrange(len(total_uncertainties)) if
                float(total_uncertainties[j])==0]
        allzeros = allzeros + zero_unc_entries
    allzeros = list(set(allzeros)) 
    print("INDICES: " + str(allzeros))
    if len(allzeros)>0:
        print("WARNING: Graphing region is being trimmed; it seems "+\
                "one of the acceptance fractions in the region graphed "+\
                "has no data in it.  Trimming those data points")
    for data in datdict: 
        datdict[data] = datdict[data].drop(datdict[data].index[\
                allzeros])
        datdict[data] = datdict[data].reset_index()
    return datdict

def PrepData_ClassSac(datafiles=[],MCfiles=[],var="nhits",precuts=None,nbins=10,xmin=10.0,xmax=100.0):
    '''
    Takes in a rootfile and returns a PandaFrame object that can be used
    for plotting in matplotlib.  Returns var vs. fractional sacrifice for
    events w/ b14<b14_low, b14>b14_high, and itr < itr_low.
    '''
    data = ROOT.TChain("output")
    MC = ROOT.TChain("output") 
    for rf in datafiles:
        data.Add(rf)
    for mf in MCfiles:
        MC.Add(mf)
    #Now, make a new dictionary object: key is cutname, value is histogram
    if precuts is not None:
        data.Draw("%s>>h_dallevents(%i,%f,%f)"% (var,nbins,xmin,xmax),
                "fitValid==1&&isCal==1&&%s" % (precuts),"goff")
        numd = data.GetEntries("fitValid==1&&%s" % (precuts))
        print("NUMBER OF DATA EVENTS W PRECUTS: " + str(numd))
        MC.Draw("%s>>h_mallevents(%i,%f,%f)"% (var,nbins,xmin,xmax),
                "fitValid==1&&%s" % (precuts),"goff")
        numMC = MC.GetEntries("fitValid==1&&%s" % (precuts))
        print("NUMBER OF MC EVENTS W PRECUTS: " + str(numMC))
    else:
        data.Draw("%s>>h_dallevents(%i,%f,%f)"% (var,nbins,xmin,xmax),
                "fitValid==1&&isCal==1","goff")
        MC.Draw("%s>>h_mallevents(%i,%f,%f)"% (var,nbins,xmin,xmax),
                "fitValid==1","goff")
    h_dallevents = gDirectory.Get("h_dallevents")
    h_dallevents.Sumw2()
    h_mallevents = gDirectory.Get("h_mallevents")
    h_mallevents.Sumw2()
    cutnames = ()
    allclasssacs = {}
    labeldict = ["Data total","MC total"]
    if precuts is None:
        precuts = "&&"
    else:
        precuts ="&&%s&&"%(precuts)
    for label in labeldict:
        graphdict={}
        vardat, fs,fs_unc =(), (), () #pandas wants ntuples
        fa,fa_unc = (),()
        h_cut_FracFlagged = ROOT.TH1D("h_cut_FracFlagged", "h_cut_FracFlagged", nbins,xmin,xmax)
        h_cut_FracFlagged.Sumw2()
        h_fracClean = ROOT.TH1D("h_fracClean", "h_fracClean", nbins, xmin, xmax)
        h_fracClean.Sumw2() 
        if label == "Data total":
            data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "fitValid==1%sisCal==1&&(beta14<-0.12 || beta14>0.95 || itr<0.55)"%(precuts) ,"goff")
            data.Draw("%s>>h_clean(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "fitValid==1%sisCal==1&&(beta14>-0.12 && beta14<0.95 && itr>0.55)"%(precuts) ,"goff")
            numclean = data.GetEntries("fitValid==1%sisCal==1&&(beta14>-0.12 && beta14<0.95 && itr>0.55)"%(precuts))
            print("NUMBER OF CLEAN EVENTS IN DATA: " + str(numclean))
        if label == "MC total":
            MC.Draw("%s>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "fitValid==1%s(beta14<-0.12 || beta14>0.95 || itr<0.55)"%(precuts) ,"goff")
            MC.Draw("%s>>h_clean(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "fitValid==1%s(beta14>-0.12 && beta14<0.95 && itr>0.55)"%(precuts) ,"goff")
            numclean = MC.GetEntries("fitValid==1%s(beta14>-0.12 && beta14<0.95 && itr>0.55)"%(precuts))
            print("NUMBER OF CLEAN EVENTS IN MC: " + str(numclean))
        print("LABEL: " + str(label))
        h_flagged = gDirectory.Get("h_flagged")
        h_flagged.Sumw2()
        h_clean = gDirectory.Get("h_clean")
        h_clean.Sumw2()
        if label == "Data total": 
            h_clean.Divide(h_dallevents)
            h_flagged.Divide(h_dallevents)
        if label == "MC total": 
            h_clean.Divide(h_mallevents)
            h_flagged.Divide(h_mallevents)
        for i in xrange(int(h_flagged.GetNbinsX())):
            vardat =  vardat + (float(h_flagged.GetBinWidth(i)) + float(h_flagged.GetBinLowEdge(i)),)
            fs = fs + (h_flagged.GetBinContent(i),)
            fs_unc = fs_unc + (h_flagged.GetBinError(i),)
        for i in xrange(int(h_clean.GetNbinsX())):
            fa = fa + (h_clean.GetBinContent(i),)
            fa_unc = fa_unc + (h_clean.GetBinError(i),)
        #for i in xrange(int(h_flagged.GetNbinsX())):
        #    vardat =  vardat + (float(h_flagged.GetBinWidth(i)) + float(h_flagged.GetBinLowEdge(i)),)
        #    fs = fs + (h_flagged.GetBinContent(i),)
        #    fs_unc = fs_unc + (h_flagged.GetBinError(i),)
        #for i in xrange(int(h_clean.GetNbinsX())):
        #    fa = fa + (h_clean.GetBinContent(i),)
        #    fa_unc = fa_unc + (h_clean.GetBinError(i),)
        graphdict["vardat"] = vardat
        graphdict["fractional_sacrifice"] = fs
        graphdict["fractional_acceptance"] = fa
        graphdict["fs_uncertainty"] = fs_unc
        graphdict["fa_uncertainty"] = fa_unc
        allclasssacs[label]= pandas.DataFrame(data=graphdict)
        del h_cut_FracFlagged
        del h_flagged
        del h_clean
        del h_fracClean
    meta = {"binwidth":((xmax-xmin)/float(nbins)),"variable":var}
    #graphdict has the sacrifice information for each cut. Now, let's plot it.
    print("THEDATA:\n")
    print(allclasssacs)
    return allclasssacs,meta 



def Plot_SacComparison(allclasssacs,meta={"binwidth":10.0, "variable":"nhits"},
        title="Heres a plot title",xlabel=None):
    sns.set_style("whitegrid")
    xkcd_colors = ['black', 'red', 'brown', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','twilight']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(allclasssacs)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for cut in allclasssacs:
        plt.errorbar(x=allclasssacs[cut].vardat, 
                y=allclasssacs[cut].fractional_sacrifice,
                yerr=allclasssacs[cut].fs_uncertainty,
                linestyle='none', marker='o', label=cut, markersize=6,
                elinewidth=3, capsize=0)
    legend = plt.legend(loc=3,frameon=1)
    frame = legend.get_frame()
    frame.set_facecolor("white")
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    plt.yscale("log")
    plt.ylabel("Fractional sacrifice",fontsize=34)
    if xlabel is None:
        plt.xlabel(meta["variable"],fontsize=34)
    else:
        plt.xlabel(xlabel,fontsize=34)
    plt.tick_params(labelsize=32)
    plt.title(title, fontsize=36)
    #plt.ion()
    plt.show()

def Plot_AccComparison(allclasssacs,meta={"binwidth":10.0, "variable":"nhits"},
        title="Heres a plot title",xlabel=None):
    sns.set_style("whitegrid")
    xkcd_colors = ['black', 'red', 'brown', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','twilight']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(allclasssacs)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for cut in allclasssacs:
        plotterdat = DeleteEmptyBins(allclasssacs[cut]) 
        print(plotterdat) 
        plt.errorbar(x=plotterdat.vardat, 
                y=plotterdat.fractional_acceptance,
                yerr=plotterdat.fa_uncertainty,
                linestyle='none', marker='o', label=cut, markersize=6,
                elinewidth=3, capsize=0)
        #plt.errorbar(x=allclasssacs[cut].vardat, 
        #        y=allclasssacs[cut].fractional_acceptance,
        #        yerr=allclasssacs[cut].fs_uncertainty,
        #        linestyle='none', marker='o', label=cut, markersize=6,
        #        elinewidth=3, capsize=0)
    legend = plt.legend(loc=3,frameon=1)
    frame = legend.get_frame()
    frame.set_facecolor("white")
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    plt.ylabel("Fractional acceptance",fontsize=34)
    if xlabel is None:
        plt.xlabel(meta["variable"],fontsize=34)
    else:
        plt.xlabel(xlabel,fontsize=34)
    plt.tick_params(labelsize=32)
    plt.title(title, fontsize=36)
    #plt.ion()
    plt.show()

def Plot_Ratio(allclasssacs,meta={"binwidth":10.0, "variable":"nhits"},
        fittotal=True,title="Title for graph",xlabel="nhits"):
    sns.set_style("whitegrid")
    xkcd_colors = ['blue','black',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','twilight']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(allclasssacs)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    allclasssacs = DeleteEmptyBins_all(allclasssacs) 
    ratio = pandas.Series(allclasssacs["Data total"].fractional_acceptance /
            allclasssacs["MC total"].fractional_acceptance) 
    ratio_unc = np.sqrt(allclasssacs["MC total"].fa_uncertainty**2 + \
            allclasssacs["Data total"].fa_uncertainty**2)
    ratio_unc = pandas.Series(ratio_unc)
    allclasssacs["ratio"] = ratio
    allclasssacs["ratio_unc"] = ratio_unc
    plt.errorbar(x=allclasssacs["MC total"].vardat, 
            y=allclasssacs["ratio"],
            yerr=allclasssacs["ratio_unc"],
            linestyle='none', marker='o', label="Data/MC Ratio", markersize=6,
            elinewidth=3, capsize=0)
    if fittotal is True:
        popt, pcov = spc.curve_fit(flatline, allclasssacs["MC total"].vardat, allclasssacs["ratio"],
                p0=[0.95], sigma=allclasssacs["ratio_unc"])
        #one standard deviation
        print("BEST FIT: " + str(popt))
        print("PCOVARIANCE: " + str(pcov))
        stdev = np.sqrt(np.diag(pcov))
        plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
        print("WEIGHTED STANDARD DEVIATION FROM FLAT: " + str(weighted_stdev(allclasssacs["ratio"],
            float(popt[0]),allclasssacs["ratio_unc"])))
    legend = plt.legend(loc=3,frameon=1)
    frame = legend.get_frame()
    frame.set_facecolor("white")
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
    #plt.yscale("log")
    plt.ylabel("Ratio",fontsize=34)
    if xlabel is None:
        plt.xlabel(meta["variable"],fontsize=34)
    else:
        plt.xlabel(xlabel,fontsize=34)
    plt.tick_params(labelsize=32)
    plt.title(title, fontsize=36)
    plt.show()
