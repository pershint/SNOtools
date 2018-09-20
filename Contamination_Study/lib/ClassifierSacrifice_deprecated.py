#Methods that were used when the statistics of the AVBox
#N16 MC were low.  Got all the runs, and statistics are fine
#Now.  
from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import maskbuilder as mb
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.optimize as spc
import pandas
import sys,time
    
    def AnalyzeData_FullAcc(self,var="nhits",xmin=None,xmax=None):
        '''
        Takes in a rootfile and returns a PandaFrame object that can be used
        for plotting in matplotlib.  Calculates the acceptance directly, not
        using (1-sacrifice).  Necessary when sacrifice statistics are low.
        '''
        print("PRECUTS: " + str(self.precuts_data))
        havebounds = False 
        for v in self.analyze_range:
            if (xmin is None or xmax is None) and v == var: 
                xmin,xmax = self.analyze_range[var][0],self.analyze_range[var][1]
                havebounds = True
        if not havebounds:
            if xmin is None or xmax is None:
                print("You do not have bounds defined for analyzing this variable"+\
                        "defined.  Check your config file.")
                sys.exit(0)
        varname = var
        if var=='posr3':
            xmin = (float(xmin)/6000.0)**3
            xmax = (float(xmax)/6000.0)**3
        print("OURVAR: " + str(var))

        print("######ANALYZING DATA/MC CLASSIFIER COMPARISON#######")
        print("XMIN, XMAX: %s,%s"%(str(xmin),str(xmax))) 
        data = ROOT.TChain("output")
        MC = ROOT.TChain("output") 
        for rf in self.rootfiles_data:
            data.Add(rf)
        for mf in self.rootfiles_mc:
            MC.Add(mf)
        #Now, make a new dictionary object: key is cutname, value is histogram
        if self.precuts_data is not None:
            data.Draw("%s>>h_dallevents(%i,%f,%f)"% (var,self.nbins,xmin,xmax),
                    "%s" % (self.precuts_data),"goff")
            MC.Draw("%s>>h_mallevents(%i,%f,%f)"% (var,self.nbins,xmin,xmax),
                    "%s" % (self.precuts_mc),"goff")
        else:
            data.Draw("%s>>h_dallevents(%i,%f,%f)"% (var,self.nbins,xmin,xmax),
                    "fitValid==1&&isCal==1","goff")
            MC.Draw("%s>>h_mallevents(%i,%f,%f)"% (var,self.nbins,xmin,xmax),
                    "fitValid==1","goff")
        h_dallevents = gDirectory.Get("h_dallevents")
        #h_dallevents.Sumw2()
        h_mallevents = gDirectory.Get("h_mallevents")
        #h_mallevents.Sumw2()
        cutnames = ()
        allclasssacs = {}
        labeldict = ["Data total","MC total"]
        for label in labeldict:
            graphdict={}
            vardat, fa,fa_unc =(), (), () #pandas wants ntuples
            fs = ()
            h_cut_FracAcc = ROOT.TH1D("h_cut_FracAcc", "h_cut_FracAcc", self.nbins,xmin,xmax)
            h_cut_FracAcc.Sumw2()
            if label == "Data total":
                data.Draw("%s>>h_accepted(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                        "%s&&(beta14>%f && beta14<%f && itr>%f)"%(self.precuts_data,
                            self.cdict["cut2_b14_low"], self.cdict["cut2_b14_high"],
                            self.cdict["cut2_itr_low"]) ,"goff")
            if label == "MC total":
                MC.Draw("%s>>h_accepted(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                        "%s&&(beta14>%f && beta14<%f && itr>%f)"%(self.precuts_mc,
                            self.cdict["cut2_b14_low"], self.cdict["cut2_b14_high"],
                            self.cdict["cut2_itr_low"]) ,"goff")
            h_accepted = gDirectory.Get("h_accepted")
            if label == "Data total": 
                h_cut_FracAcc.Divide(h_accepted,h_dallevents,1.,1.,"b")
            if label == "MC total": 
                h_cut_FracAcc.Divide(h_accepted,h_mallevents,1.,1.,"b")
            for i in xrange(int(h_accepted.GetNbinsX())+1):
                if i==0:
                    continue
                num = float(h_accepted.GetBinContent(i))
                if label == "Data total":
                    denom = float(h_dallevents.GetBinContent(i))
                if label == "MC total":
                    denom = float(h_mallevents.GetBinContent(i))
                #Assume poisson errors.
                binerror=None
                if denom==0:
                    print("no data in this bin. setting error to zero")
                    binerror = 0
                elif num==0:
                    binerror = 0
                else:
                    binerror = (num/denom)*np.sqrt((1./num)+(1./denom))
                vardat =  vardat + ((float(h_accepted.GetBinWidth(i))/2.0) + float(h_accepted.GetBinLowEdge(i)),)
                fa = fa + (h_cut_FracAcc.GetBinContent(i),)
                fa_unc = fa_unc + (binerror,)
            h_fracDirty = ROOT.TH1D("h_fracDirty", "h_fracDirty", self.nbins, xmin, xmax)
            h_fracDirty.Sumw2()
            for j in xrange(int(h_cut_FracAcc.GetNbinsX())+1):
                if j==0:
                    continue
                h_fracDirty.SetBinContent(j,1.0)
            h_fracDirty.Add(h_cut_FracAcc,-1.0)
            for i in xrange(int(h_fracDirty.GetNbinsX())+1):
                if i==0:
                    continue
                fs = fs + (h_fracDirty.GetBinContent(i),)
            graphdict["vardat"] = vardat
            graphdict["fractional_sacrifice"] = fs
            graphdict["fractional_acceptance"] = fa
            graphdict["fs_uncertainty"] = fa_unc
            print("THEGRAPHDICT: " + str(graphdict))
            allclasssacs[label]= pandas.DataFrame(data=graphdict)
            del h_cut_FracAcc
            del h_accepted
            del h_fracDirty
    
        meta = {"binwidth":((xmax-xmin)/float(self.nbins)),"variable":var}
        print("SACPERCUT DATA: " + str(self.sac_percut)) 
        self.sac_percut = allclasssacs
        self.sac_percut_metadata = meta
        #If any bins had no statistics, remove them
        self.sac_percut = self._DeleteEmptyBins_all_FullAcc(self.sac_percut) 


    def _DeleteEmptyBins_all_FullAcc(self,datdict):
        '''Deletes any bins from all data sets in array where the uncertainty is zero (this
        happens when the bin has no data in it at all)'''
        allzeros = []
        print(datdict)
        for data in datdict:
            if data=='ratio' or data=='ratio_unc':
                continue
            print("WE HAVE KEY: " + str(data))
            fractional_acceptance = datdict[data]["fractional_acceptance"]
            print(fractional_acceptance)
            zero_unc_entries = []
            zero_unc_entries = [j for j in xrange(len(fractional_acceptance)) if
                    float(fractional_acceptance[j])==0]
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
            print("%s AFTER TRIM: %s"%(str(data),str(datdict[data])))
        return datdict

