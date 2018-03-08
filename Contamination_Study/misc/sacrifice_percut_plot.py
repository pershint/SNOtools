from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import maskbuilder as mb
import seaborn
import pandas
import time

def GetAnalysisMask():
    util = RAT.DU.Utility.Get()
    util.LoadDBAndBeginRun()
    dcmask = RAT.GetDataCleaningWord('analysis_mask')
    return dcmask

def SacVSNhit(rootfile,nbins=10,nmin=10.0,nmax=100.0,dcmask=None):
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
    x,unc,y = (), (),() #pandas wants ntuples
    for cut in plotmask:
        h_cut_FracPassed = ROOT.TH1D("h_cut_FracPassed", "h_cut_FracPassed", nbins,nmin,nmax)
        h_cut_FracPassed.Sumw2()
        data.Draw("nhits>>h_passed(%i,%f,%f)" % (nbins,nmin,nmax),
                "((dcFlagged&%i)==%i)" % (cut,cut),"goff")
        h_passed = gDirectory.Get("h_passed")
        h_passed.Sumw2()
        h_cut_FracPassed.Divide(h_passed,h_allevents,1.,1.,"b")
        for i in xrange(h_cut_FracPassed.GetNbinsX()):
            cutnames = cutnames + (plotmask[cut],) 
            x =  x + (int(h_cut_FracPassed.GetBinWidth(i) + h_cut_FracPassed.GetBinLowEdge(i)),)
            y = y + (h_cut_FracPassed.GetBinContent(i),)
            unc = unc + (h_cut_FracPassed.GetBinError(i),)
        del h_cut_FracPassed
        del h_passed
    #graphdict has the sacrifice information for each cut. Now, let's plot it.
    graphdict = {"cutnames": cutnames, "nhits": x, "frac_passed": y, "frac_passed_unc": unc}
   
    return pandas.DataFrame(data=graphdict)
