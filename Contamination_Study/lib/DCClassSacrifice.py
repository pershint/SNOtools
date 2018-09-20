from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import maskbuilder as mb
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as spc
import pandas
import sys,time

def GetAnalysisMask():
    util = RAT.DU.Utility.Get()
    util.LoadDBAndBeginRun()
    dcmask = RAT.GetDataCleaningWord('analysis_mask')
    return dcmask


class SacrificeAnalyzer(object):
    def __init__(self, rootfiles_data=None,rootfiles_mc=None,cuts_dict=None):
        self.rootfiles_data = rootfiles_data
        self.rootfiles_mc = rootfiles_mc
        self.cdict = cuts_dict
        self.sac_percut = {} #Contains binned sacrifice for each of the highest cuts
        self.sac_percut_metadata = {} #Metadata, like variable fitted to and bin width
        self.nbins = 11
        self.analyze_range = {}
        self._SetAnalyzeRanges() 
        self._DefinePrecuts()
        self.xlabel_dict = {'energy': 'Energy (MeV)', 'udotr': r'$U . R$',
                'posr3': r'$(R/R_{AV})^{3}$','posr':'Radius (mm)'}
    
    def SetBinNumber(self, nbins):
        self.nbins = nbins + 1
    
    def _flatline(self,x,b):
        return b

    def _weighted_stdev(self,vals,valavg, uncs):
        weights = 1/(uncs**2)
        return np.sqrt((len(weights)*np.sum(weights*(vals-valavg)**2))/((len(weights)-1)*np.sum(weights)))

    def _SetAnalyzeRanges(self):
        '''Uses variable ranges defined in cut dictionary and sets the analysis range
        to these windows'''
        self.analyze_range['energy'] = [self.cdict['E_low'],self.cdict['E_high']]
        self.analyze_range['udotr'] = [self.cdict['udotr_low'],self.cdict['udotr_high']]
        self.analyze_range['beta14'] = [self.cdict['cut2_b14_low'],self.cdict['cut2_b14_high']]
        self.analyze_range['itr'] = [self.cdict['cut2_itr_low'],1.0]
        self.analyze_range['posr'] = [self.cdict['r_low'],self.cdict['r_high']] 
        self.analyze_range['posr3'] = [self.cdict['r_low'],self.cdict['r_high']] 
        self.analyze_range['nhits'] = [self.cdict['nhits_low'],self.cdict['nhits_high']] 
    
    def _DefinePrecuts(self):
        self.precuts_data=[]
        self.precuts_mc=[]
        if self.cdict['r_high'] is not None:
            self.precuts_data.append("posr<"+str(self.cdict['r_high']))
            self.precuts_mc.append("posr<"+str(self.cdict['r_high']))
        if self.cdict['r_low'] is not None:
            self.precuts_data.append("posr>"+str(self.cdict['r_low']))
            self.precuts_mc.append("posr>"+str(self.cdict['r_low']))
        if self.cdict['nhits_high'] is not None:
            self.precuts_data.append('nhits'+"<"+str(self.cdict['nhits_high']))
            self.precuts_mc.append('nhits'+"<"+str(self.cdict['nhits_high']))
        if self.cdict['nhits_low'] is not None:
            self.precuts_data.append('nhits'+">"+str(self.cdict['nhits_low']))
            self.precuts_mc.append('nhits'+">"+str(self.cdict['nhits_low']))
        if self.cdict['udotr_high'] is not None:
            self.precuts_data.append("udotr<"+str(self.cdict['udotr_high']))
            self.precuts_mc.append("udotr<"+str(self.cdict['udotr_high']))
        if self.cdict['udotr_low'] is not None:
            self.precuts_data.append("udotr>"+str(self.cdict['udotr_low']))
            self.precuts_mc.append("udotr>"+str(self.cdict['udotr_low']))
        if self.cdict['E_high'] is not None: 
            self.precuts_data.append("energy<"+str(self.cdict['E_high']))
            self.precuts_mc.append("energy<"+str(self.cdict['E_high']))
        if self.cdict['E_low'] is not None: 
            self.precuts_data.append("energy>"+str(self.cdict['E_low']))
            self.precuts_mc.append("energy>"+str(self.cdict['E_low']))
        if self.cdict['Z_low'] is not None:
            self.precuts_data.append("posz>"+str(self.cdict['Z_low']*10.0))
            self.precuts_mc.append("posz>"+str(self.cdict['Z_low']*10.0))
        if self.cdict['Z_high'] is not None:
            self.precuts_data.append("posz<"+str(self.cdict['Z_high']*10.0))
            self.precuts_mc.append("posz<"+str(self.cdict['Z_high']*10.0))
        if self.cdict['fitValid'] is not None: 
            if self.cdict['fitValid'] is True: 
                self.precuts_data.append("fitValid==1")
                self.precuts_mc.append("fitValid==1")
            else:
                self.precuts_data.append("fitValid==0")
                self.precuts_mc.append("fitValid==0")
        if self.cdict['AVudotrCut'] is not None:
            if self.cdict['AVudotrCut'] is True:
                rcut = self.cdict["AVudotrCut_rcut"]
                cutparam = (rcut/6000.0)**3
                self.precuts_data.append("(udotr > (1.0 - 12 *((posr3-%f)**2)))"%(cutparam))
                self.precuts_mc.append("(udotr > (1.0 - 12 *((posr3-%f)**2)))"%(cutparam))
        if self.cdict['isCal'] is not None:
            if self.cdict['isCal'] is True: self.precuts_data.append("isCal==1")
            else: self.precuts_data.append("isCal==0")
        self.precuts_data.append("((dcFlagged&%s)==%s)" % (self.cdict["sacpath_DCmask"],\
            self.cdict["sacpath_DCmask"]))
        self.precuts_data.append("((triggerWord&%s)==0)" % (self.cdict["path_trigmask"]))
        self.precuts_data = self._cutfuse(self.precuts_data, "&&")
        self.precuts_mc = self._cutfuse(self.precuts_mc, "&&")

    def _cutfuse(self,stringlist,delim):
        outstring = ""
        for j,s in enumerate(stringlist):
            if j == 0:
                outstring = s
            else:
                outstring=outstring+delim+s
        return outstring


    def _DeleteEmptyBins(self):
        '''Deletes any bins with zero events and zero uncertainty (i.e. bins
        with zero events in them) from the self.sacrifices dataset'''
        total_uncertainties = self.sac_percut["total"].fs_uncertainty
        zero_unc_entries = [] 
        zero_unc_entries = [j for j in xrange(len(total_uncertainties)) if
                float(total_uncertainties[j])==0]
        if len(zero_unc_entries)>0:
            print("WARNING: Graphing region is being trimmed; it seems you"+\
                    " have some regions in your ROI with no data in it.")
        for cut in self.sac_percut:
            self.sac_percut[cut] = self.sac_percut[cut].drop(self.sac_percut[cut].index[\
                    zero_unc_entries])
            self.sac_percut[cut] = self.sac_percut[cut].reset_index()
            #index(\
            #        index=range(0,len(self.sac_percut[cut])-1))

    def GetFitTotalAndUncertainties(self):
        '''Finds the best fit total sacrifice and uncertainties (statistical and
        systematic is weighted stdev. from flat).'''
        popt, pcov = spc.curve_fit(self._flatline, self.sac_percut["total"].vardat, self.sac_percut["total"].fractional_sacrifice,
                p0=[0.02], sigma=self.sac_percut["total"].fs_uncertainty)
        #one standard deviation
        print("BEST FIT: " + str(popt))
        print("PCOVARIANCE: " + str(pcov))
        statistical_unc = None 
        try:
            statistical_unc = np.sqrt(np.diag(pcov))
        except ValueError:
            print("Likely that your fit failed to converge.  You may have had"+\
                    " bins with no events in them in the dataset. Searching for "+\
                    "empty bins and deleting them before calculation...")
            return None
        systematic_unc = self._weighted_stdev(self.sac_percut["total"].fractional_sacrifice,
            float(popt[0]),self.sac_percut["total"].fs_uncertainty)
        print("WEIGHTED STANDARD DEVIATION FROM FLAT: " + str(systematic_unc)) 
        self.total_sacrifice['sacrifice'] = float(popt[0])
        self.total_sacrifice['stat_unc'] = float(statistical_unc)
        self.total_sacrifice['sys_unc'] = float(systematic_unc)
        print(self.total_sacrifice) 
        return self.total_sacrifice
   

class DCSacrificeAnalyzer(SacrificeAnalyzer):
    def __init__(self, rootfiles_data=None, cuts_dict=None):
        super(DCSacrificeAnalyzer, self).__init__(rootfiles_data=rootfiles_data, rootfiles_mc=None,cuts_dict=cuts_dict)
        self.cut1_mask = {'dcmask':self.cdict['cut1_sacDCmask'],'dcmask_cutnames': None}
        self.total_sacrifice = {} #Total sacrifice, statistical unc, and sys unc.
    
    def AnalyzeData(self,var="nhits",xmin=None,xmax=None):
        '''
        Takes in a rootfile and returns a PandaFrame object that can be used
        for plotting in matplotlib.  Returns nhits vs. fractional sacrifice for
        Each cut in the given dcmask. if precuts==None, isCal&&fitValid enforced as
        precuts.  Otherwise, make sure to put them in the precuts defined for N16!
        '''
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
        
        print("######ANALYZING DATA CLEANING SACRIFICE#######")
        print("XMIN, XMAX: %s,%s"%(str(xmin),str(xmax))) 
        dcmask = self.cdict['cut1_sacDCmask']
        if dcmask == None:
            print("No DC mask in cuts config.  Please give a DC mask of cuts to plot")
            sys.exit(0)
        data = ROOT.TChain("output")
        for rf in self.rootfiles_data:
            data.Add(rf)
        fullmask = mb.get_dcwords()
        print(fullmask)
        plotmask = {}
        for cut in fullmask:
            if 'prescale' in fullmask[cut]:
                continue
            if 'waterblind' in fullmask[cut]:
                continue
            if mb.binary_bit_to_int(cut)&dcmask==mb.binary_bit_to_int(cut):
                print("WE ARE ADDING IN CUT: " + str(cut))
                print("THAT IS... CUT: " + str(fullmask[cut]))
                plotmask[cut] = fullmask[cut]
        plotmask[dcmask] = "total"
        #Now, make a new dictionary object: key is cutname, value is histogram
        if self.precuts_data is not None:
            data.Draw("%s>>h_allevents(%i,%f,%f)"% (var,self.nbins,xmin,xmax),
                    "%s"%(self.precuts_data),"goff")

        else:
            data.Draw("%s>>h_allevents(%i,%f,%f)"% (var,self.nbins,xmin,xmax),
                    "","goff")
        h_allevents = gDirectory.Get("h_allevents")
        h_allevents.Sumw2()
        cutnames = ()
        self.sac_percut = {}
        for cut in plotmask:
            if cut!= dcmask:
                cutint = mb.binary_bit_to_int(cut)
            else:
                cutint = cut
            graphdict={}
            vardat, fs,fs_unc =(), (), () #pandas wants ntuples
            h_cut_FracFlagged = ROOT.TH1D("h_cut_FracFlagged", "h_cut_FracFlagged", self.nbins,xmin,xmax)
            h_cut_FracFlagged.Sumw2()
            if self.precuts_data is not None:
                data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                        "%s"%(self.precuts_data)+\
                             "&&((dcFlagged&%i)!=%i)" % (cutint,cutint),"goff")
            else:
                data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                        "((dcFlagged&%i)!=%i)" % (cutint,cutint),"goff")
            h_flagged = gDirectory.Get("h_flagged")
            h_flagged.Sumw2()
            h_cut_FracFlagged.Divide(h_flagged,h_allevents,1.,1.,"b")
            for i in xrange(int(h_cut_FracFlagged.GetNbinsX()+1)):
                if i==0:
                    continue
                vardat =  vardat + ((float(h_cut_FracFlagged.GetBinWidth(i))/2.0) + float(h_cut_FracFlagged.GetBinLowEdge(i)),)
                fs = fs + (h_cut_FracFlagged.GetBinContent(i),)
                fs_unc = fs_unc + (h_cut_FracFlagged.GetBinError(i),)
            graphdict["vardat"] = vardat
            graphdict["fractional_sacrifice"] = fs
            graphdict["fs_uncertainty"] = fs_unc
            self.sac_percut[plotmask[cut]]= pandas.DataFrame(data=graphdict)
            del h_cut_FracFlagged
            del h_flagged
        #graphdict has the sacrifice information for each cut. Now, let's plot it.
        self.sac_percut_metadata = {"binwidth":((xmax-xmin)/float(self.nbins)),"variable":varname}
        self._GetTopSacs()
        self._DeleteEmptyBins()
        print self.sac_percut

    def _GetTopSacs(self,topnumber=7):
        dcmask = self.cdict['cut1_sacDCmask']
        if len(self.sac_percut) <= topnumber:
            print("You already only have five cuts.  Not combining any")
            return
        names = []
        sacrifices = []
        for cut in self.sac_percut:
            sacrifices.append(np.average(self.sac_percut[cut]["fractional_sacrifice"]))
            names.append(cut)
        sortednames = [x for _,x in sorted(zip(sacrifices,names))]
        topnames = sortednames[(len(sortednames)-(topnumber+1)):len(sortednames)]
        self.cut1_mask['dcmask_cutnames'] = topnames

    def ShowPlottedSacrifice(self,fittotal=True,title=None,savedir="."):
        sns.set_style("whitegrid")
        xkcd_colors = ['black','slate blue', 'fluro green', 'brown', 'blue',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'leaf',
                'aqua blue','vomit', 'red','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        fig = plt.figure()
        fig.set_size_inches(18.5,11.5)
        ax = fig.add_subplot(1,1,1)
        if self.cut1_mask['dcmask_cutnames'] is None:
            for cut in self.sac_percut:
                plt.errorbar(x=self.sac_percut[cut].vardat, 
                        y=self.sac_percut[cut].fractional_sacrifice,
                        yerr=self.sac_percut[cut].fs_uncertainty,
                        linestyle='none', marker='+', label=cut, markersize=7,
                        elinewidth=2, capsize=0)
        else:
            fracsum,fracuncsum=[],[]
            for cut in self.sac_percut:
                if cut in self.cut1_mask['dcmask_cutnames']:
                    print(cut)
                    plt.errorbar(x=self.sac_percut[cut].vardat, y=self.sac_percut[cut].fractional_sacrifice,
                            yerr=self.sac_percut[cut].fs_uncertainty, linestyle='none',
                            marker='o', label=cut, capsize=0, elinewidth=2, markersize=8)
                else:
                    print("CUT GOING INTO ALL OTHER CUTS: " + str(cut))
                    fracsum.append(self.sac_percut[cut].fractional_sacrifice)
                    print("FRACSUM ARRAY: " + str(fracsum))
                    fracuncsum.append(self.sac_percut[cut].fs_uncertainty)
            plt.errorbar(x=self.sac_percut[cut].vardat,y=sum(fracsum),yerr=sum(fracuncsum),
                    linestyle='none',marker='o', capsize=0, elinewidth=2, label='All other cuts', markersize=8)
    
        if fittotal is True:
            popt, pcov = spc.curve_fit(self._flatline, self.sac_percut["total"].vardat, self.sac_percut["total"].fractional_sacrifice,
                    p0=[0.02], sigma=self.sac_percut["total"].fs_uncertainty)
            #one standard deviation
            try:
                stdev = np.sqrt(np.diag(pcov))
            except ValueError:
                print("Fit likely failed.  You may not have enough data"+\
                        "in this region to fit a straight line.")
                raise
            plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
        legend = plt.legend(loc=3,frameon=1)
        frame = legend.get_frame()
        frame.set_facecolor("white")
        plt.yscale("log")
        plt.ylabel("Fractional sacrifice",fontsize=32)
      
        varindict = False
        for var in self.xlabel_dict:
            if self.sac_percut_metadata["variable"] == var:
                plt.xlabel(self.xlabel_dict[var],fontsize=32)
                varindict = True
            else:
                if varindict is False: 
                    plt.xlabel(self.sac_percut_metadata["variable"],fontsize=32)
        plt.tick_params(labelsize=30)
        if title is None:
            plt.title("Fractional sacrifice due to data cleaning cuts",fontsize=36)
        else: 
            plt.title(title,fontsize=36)
        variable = self.sac_percut_metadata["variable"]
        plt.savefig(savedir+"/DCSac_%s.pdf"%(variable))
        #plt.show()
        plt.close()

class DataMCClassAnalyzer(SacrificeAnalyzer):
    def __init__(self, rootfiles_data=None, rootfiles_mc=None,cuts_dict=None):
        super(DataMCClassAnalyzer, self).__init__(rootfiles_data=rootfiles_data, 
                rootfiles_mc=rootfiles_mc, cuts_dict=cuts_dict)
        self.datamc_ratio = {} 
    
    def GetFitTotalAndUncertainties(self):
        '''Finds the best fit total sacrifice and uncertainties (statistical and
        systematic is weighted stdev. from flat).'''
        if self.sac_percut is None:
            print("You must analyze your data before you can do this fit!")
            return
        ratio = pandas.Series(self.sac_percut["Data total"].fractional_acceptance /
                self.sac_percut["MC total"].fractional_acceptance) 
        ratio_unc = np.sqrt(self.sac_percut["MC total"].fs_uncertainty**2 + \
                self.sac_percut["Data total"].fs_uncertainty**2)
        ratio_unc = pandas.Series(ratio_unc)
        self.sac_percut["ratio"] = ratio
        self.sac_percut["ratio_unc"] = ratio_unc
        popt, pcov = spc.curve_fit(self._flatline, self.sac_percut["MC total"].vardat, self.sac_percut["ratio"],
                p0=[0.98], sigma=self.sac_percut["ratio_unc"])
        #one standard deviation
        print("BEST FIT: " + str(popt))
        print("PCOVARIANCE: " + str(pcov))
        statistical_unc = None 
        try:
            statistical_unc = np.sqrt(np.diag(pcov))
        except ValueError:
            print("Likely that your fit failed to converge.  You may have had"+\
                    " bins with no events in them in the dataset. Searching for "+\
                    "empty bins and deleting them before calculation...")
            return None
        systematic_unc = self._weighted_stdev(self.sac_percut["ratio"],
            float(popt[0]),self.sac_percut["ratio_unc"])
        print("WEIGHTED STANDARD DEVIATION FROM FLAT: " + str(systematic_unc)) 
        self.datamc_ratio['Data/MC ratio'] = float(popt[0])
        self.datamc_ratio['stat_unc'] = float(statistical_unc)
        self.datamc_ratio['sys_unc'] = float(systematic_unc)
        print(self.datamc_ratio) 
        return self.datamc_ratio

    def SaveRatioToCSV(self,loc,name):
        '''Saves the ratio, ratio_unc, and bin centers for the current dataset to
        a CSV file at the specified location. Values are the right-hand edge of
        the bins for the variable'''
        ratio = pandas.Series(self.sac_percut["ratio"]) 
        variable = pandas.Series(self.sac_percut["Data total"].vardat)
        variable = variable + (self.sac_percut_metadata["binwidth"]/2.0)
        ratio_unc = pandas.Series(self.sac_percut["ratio_unc"])
        xlabel = self.sac_percut_metadata["variable"]
        thegoods = {"ratio": ratio.round(4), "ratio_unc": ratio_unc.round(4), xlabel: variable.round(2)}
        thegoods_pd = pandas.DataFrame(thegoods)
        thegoods_pd.to_csv("%s/%s"%(loc,name))

    def _DeleteEmptyBins(self,data):
        '''Deletes any bins with zero events and zero uncertainty (i.e. bins
        with zero events in them) from the self.sacrifices dataset'''
        total_uncertainties = data.fs_uncertainty
        zero_unc_entries = [] 
        zero_unc_entries = [j for j in xrange(len(total_uncertainties)) if
                float(total_uncertainties[j])==0]
        print("INDICES: " + str(zero_unc_entries))
        if len(zero_unc_entries)>0:
            print("WARNING: Graphing region is being trimmed; it seems you"+\
                    " have some regions in your ROI with no data in it.")
        data = data.drop(data.index[\
                zero_unc_entries])
        data = data.reset_index()
        return data
    
    def _DeleteEmptyBins_all(self,datdict):
        '''Deletes any bins from all data sets in array where the uncertainty is zero (this
        happens when the bin has no data in it at all)'''
        allzeros = []
        print(datdict)
        for data in datdict:
            if data=='ratio' or data=='ratio_unc':
                continue
            print("WE HAVE KEY: " + str(data))
            total_uncertainties = datdict[data]["fs_uncertainty"]
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
    def AnalyzeData(self,var="nhits",xmin=None,xmax=None):
        '''
        Takes in a rootfile and returns a PandaFrame object that can be used
        for plotting in matplotlib.  Returns var vs. fractional sacrifice for
        events w/ b14<b14_low, b14>b14_high, and itr < itr_low.
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
        h_dallevents.Sumw2()
        h_mallevents = gDirectory.Get("h_mallevents")
        h_mallevents.Sumw2()
        cutnames = ()
        allclasssacs = {}
        labeldict = ["Data total","MC total"]
        for label in labeldict:
            graphdict={}
            vardat, fs,fs_unc =(), (), () #pandas wants ntuples
            fa = () 
            h_cut_FracFlagged = ROOT.TH1D("h_cut_FracFlagged", "h_cut_FracFlagged", self.nbins,xmin,xmax)
            h_cut_FracFlagged.Sumw2()
            if label == "Data total":
                data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                        "%s&&(beta14<%f || beta14>%f || itr<%f)"%(self.precuts_data,
                            self.cdict["cut2_b14_low"], self.cdict["cut2_b14_high"],
                            self.cdict["cut2_itr_low"]) ,"goff")
            if label == "MC total":
                MC.Draw("%s>>h_flagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                        "%s&&(beta14<%f || beta14>%f || itr<%f)"%(self.precuts_mc,
                            self.cdict["cut2_b14_low"], self.cdict["cut2_b14_high"],
                            self.cdict["cut2_itr_low"]) ,"goff")
            print("LABEL: " + str(label))
            h_flagged = gDirectory.Get("h_flagged")
            h_flagged.Sumw2()
            if label == "Data total": 
                h_cut_FracFlagged.Divide(h_flagged,h_dallevents,1.,1.,"b")
            if label == "MC total": 
                h_cut_FracFlagged.Divide(h_flagged,h_mallevents,1.,1.,"b")
            for i in xrange(int(h_cut_FracFlagged.GetNbinsX())+1):
                if i==0:
                    continue
                vardat =  vardat + ((float(h_cut_FracFlagged.GetBinWidth(i))/2.0) + float(h_cut_FracFlagged.GetBinLowEdge(i)),)
                fs = fs + (h_cut_FracFlagged.GetBinContent(i),)
                fs_unc = fs_unc + (h_cut_FracFlagged.GetBinError(i),)
            h_fracClean = ROOT.TH1D("h_fracClean", "h_fracClean", self.nbins, xmin, xmax)
            for j in xrange(int(h_cut_FracFlagged.GetNbinsX())+1):
                if j==0:
                    continue
                h_fracClean.SetBinContent(j,1.0)
            h_fracClean.Add(h_cut_FracFlagged,-1.0)
            for i in xrange(int(h_fracClean.GetNbinsX())+1):
                if i==0:
                    continue
                fa = fa + (h_fracClean.GetBinContent(i),)
            graphdict["vardat"] = vardat
            graphdict["fractional_sacrifice"] = fs
            graphdict["fractional_acceptance"] = fa
            graphdict["fs_uncertainty"] = fs_unc
            print("THEGRAPHDICT: " + str(graphdict))
            allclasssacs[label]= pandas.DataFrame(data=graphdict)
            del h_cut_FracFlagged
            del h_flagged
    
        meta = {"binwidth":((xmax-xmin)/float(self.nbins)),"variable":var}
        print("SACPERCUT DATA: " + str(self.sac_percut))
        self.sac_percut = allclasssacs
        self.sac_percut_metadata = meta
        #If any bins did not have any data, remove them
        self.sac_percut = self._DeleteEmptyBins_all(self.sac_percut) 

    def Plot_SacComparison(self,title="Heres a plot title",xlabel=None):
        if self.sac_percut is None:
            print("You must first analyze the input root data!")
            return
        sns.set_style("whitegrid")
        xkcd_colors = ['black', 'red', 'brown', 'blue',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
                'aqua blue','vomit', 'black','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        fig = plt.figure()
        fig.set_size_inches(18.5,11.5)
        ax = fig.add_subplot(1,1,1)
        for cut in self.sac_percut:
            plt.errorbar(x=self.sac_percut[cut].vardat, 
                    y=self.sac_percut[cut].fractional_sacrifice,
                    yerr=self.sac_percut[cut].fs_uncertainty,
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
            plt.xlabel(self.sac_percut_metadata["variable"],fontsize=34)
        else:
            plt.xlabel(xlabel,fontsize=34)
        plt.tick_params(labelsize=32)
        plt.title(title, fontsize=36)
        #plt.show()
        plt.close()
    
    def Plot_AccComparison(self,title="Heres a plot title",xlabel=None):
        if self.sac_percut is None:
            print("You must first analyze the input root data!")
            return
        sns.set_style("whitegrid")
        xkcd_colors = ['black', 'red', 'brown', 'blue',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
                'aqua blue','vomit', 'black','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for cut in self.sac_percut:
            plotterdat = self._DeleteEmptyBins(self.sac_percut[cut]) 
            print(plotterdat) 
            plt.errorbar(x=plotterdat.vardat, 
                    y=plotterdat.fractional_acceptance,
                    yerr=plotterdat.fs_uncertainty,
                    linestyle='none', marker='o', label=cut, markersize=6,
                    elinewidth=3, capsize=0)
        legend = plt.legend(loc=3,frameon=1)
        frame = legend.get_frame()
        frame.set_facecolor("white")
        y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
        plt.ylabel("Fractional acceptance",fontsize=34)
        if xlabel is None:
            plt.xlabel(self.sac_percut_metadata["variable"],fontsize=34)
        else:
            plt.xlabel(xlabel,fontsize=34)
        plt.tick_params(labelsize=32)
        plt.title(title, fontsize=36)
        #plt.show()
        plt.close()
    
    def PlotRatio(self,fittotal=True,title="Title for graph",xlabel="nhits",savedir="."):
        if self.sac_percut is None:
            print("You must first analyze the input root data!")
            return
        sns.set_style("whitegrid")
        xkcd_colors = ['blue','black',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
                'aqua blue','vomit', 'black','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        fig = plt.figure()
        fig.set_size_inches(18.5,11.5)
        ax = fig.add_subplot(1,1,1)
        self.sac_percut = self._DeleteEmptyBins_all(self.sac_percut) 
        ratio = pandas.Series(self.sac_percut["Data total"].fractional_acceptance /
                self.sac_percut["MC total"].fractional_acceptance) 
        ratio_unc = np.sqrt(self.sac_percut["MC total"].fs_uncertainty**2 + \
                self.sac_percut["Data total"].fs_uncertainty**2)
        ratio_unc = pandas.Series(ratio_unc)
        self.sac_percut["ratio"] = ratio
        self.sac_percut["ratio_unc"] = ratio_unc
        plt.errorbar(x=self.sac_percut["MC total"].vardat, 
                y=self.sac_percut["ratio"],
                yerr=self.sac_percut["ratio_unc"],
                linestyle='none', marker='o', label="Data/MC Ratio", markersize=6,
                elinewidth=3, capsize=0)
        if fittotal is True:
            popt, pcov = spc.curve_fit(self._flatline, self.sac_percut["MC total"].vardat, self.sac_percut["ratio"],
                    p0=[0.98], sigma=self.sac_percut["ratio_unc"])
            #one standard deviation
            print("BEST FIT: " + str(popt))
            print("PCOVARIANCE: " + str(pcov))
            stdev = np.sqrt(np.diag(pcov))
            plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
            print("WEIGHTED STANDARD DEVIATION FROM FLAT: " + str(self._weighted_stdev(self.sac_percut["ratio"],
                float(popt[0]),self.sac_percut["ratio_unc"])))
        legend = plt.legend(loc=3,frameon=1)
        frame = legend.get_frame()
        frame.set_facecolor("white")
        y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
        #plt.yscale("log")
        plt.ylabel("Ratio",fontsize=34)
        varindict = False
        for var in self.xlabel_dict:
            if self.sac_percut_metadata["variable"] == var:
                plt.xlabel(self.xlabel_dict[var],fontsize=32)
                varindict = True
            else:
                if varindict is False: 
                    plt.xlabel(self.sac_percut_metadata["variable"],fontsize=32)
        plt.tick_params(labelsize=32)
        plt.title(title, fontsize=36)
        variable = self.sac_percut_metadata["variable"]
        plt.savefig(savedir+"/DataMCComp_%s.pdf"%(variable))
        #plt.show()
        plt.close()

class ClassSacrificeAnalyzer(SacrificeAnalyzer):
    def __init__(self, rootfiles_data=None, cuts_dict=None):
        super(ClassSacrificeAnalyzer, self).__init__(rootfiles_data=rootfiles_data, cuts_dict=cuts_dict)
        self.total_sacrifice = {} #Total sacrifice, statistical unc, and sys unc.

    def AnalyzeData(self,var="nhits",xmin=None,xmax=None):
        '''
        Takes in a rootfile and returns a PandaFrame object that can be used
        for plotting in matplotlib.  Returns var vs. fractional sacrifice for
        events w/ b14<b14_low, b14>b14_high, and itr < itr_low.
        '''
         
        havebounds = False 
        for v in self.analyze_range:
            if xmin is None or xmax is None and v == var: 
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
        print("######ESTIMATING CLASSIFIER SACRIFICE#######") 
        data = ROOT.TChain("output")
        for rf in self.rootfiles_data:
            data.Add(rf)
        #Now, make a new dictionary object: key is cutname, value is histogram
        print("PRECUTS IN CLASSIFIER SAC ESTIMATOR: " + str(self.precuts_data))
        if self.precuts_data is not None:
            data.Draw("%s>>h_allevents(%i,%f,%f)"% (var,self.nbins,xmin,xmax),
                    "%s" % (self.precuts_data),"goff")
        else:
            data.Draw("%s>>h_allevents(%i,%f,%f)"% (var,self.nbins,xmin,xmax),
                    "","goff")
        h_allevents = gDirectory.Get("h_allevents")
        h_allevents.Sumw2()
        cutnames = ()
        labeldict = ["b14_low", "b14_high","itr_low","total"]
        for cut in labeldict:
            graphdict={}
            vardat, fs,fs_unc =(), (), () #pandas wants ntuples
            h_cut_FracFlagged = ROOT.TH1D("h_cut_FracFlagged", "h_cut_FracFlagged", self.nbins,xmin,xmax)
            h_cut_FracFlagged.Sumw2()
            if cut == "b14_low":
                data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                        "%s&&beta14<%s"%(self.precuts_data,str(self.cdict['cut2_b14_low'])),"goff")
            if cut == "b14_high":
                data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                        "%s&&beta14>%s"%(self.precuts_data,str(self.cdict['cut2_b14_high'])) ,"goff")
            if cut == "itr_low":
                data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                        "%s&&itr<%s"%(self.precuts_data,str(self.cdict['cut2_itr_low'])),"goff")
            if cut == "total":
                data.Draw("%s>>h_flagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                        "%s&&(beta14<%f || beta14>%f || itr<%f)"%(self.precuts_data,
                            self.cdict['cut2_b14_low'],self.cdict['cut2_b14_high'],
                            self.cdict['cut2_itr_low']) ,"goff")
            print("CUT: " + str(cut))
            h_flagged = gDirectory.Get("h_flagged")
            h_flagged.Sumw2()
            h_cut_FracFlagged.Divide(h_flagged,h_allevents,1.,1.,"b")
            for i in xrange(int(h_cut_FracFlagged.GetNbinsX())+1):
                if i==0:
                    continue
                vardat =  vardat + ((float(h_cut_FracFlagged.GetBinWidth(i))/2.0) + float(h_cut_FracFlagged.GetBinLowEdge(i)),)
                fs = fs + (h_cut_FracFlagged.GetBinContent(i),)
                fs_unc = fs_unc + (h_cut_FracFlagged.GetBinError(i),)
            graphdict["vardat"] = vardat
            graphdict["fractional_sacrifice"] = fs
            graphdict["fs_uncertainty"] = fs_unc
            self.sac_percut[cut]= pandas.DataFrame(data=graphdict)
            del h_cut_FracFlagged
            del h_flagged
        print(self.sac_percut) 
        self._DeleteEmptyBins()
        self.sac_percut_metadata = {"binwidth":((xmax-xmin)/float(self.nbins)),"variable":varname}
        #graphdict has the sacrifice information for each cut. Now, let's plot it.

    def ShowPlottedSacrifice(self,fittotal=True,title=None,savedir="."):
        sns.set_style("whitegrid")
        xkcd_colors = ['slate blue', 'black', 'brown', 'blue',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
                'aqua blue','vomit', 'black','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        
        for cut in self.sac_percut:
            plt.errorbar(x=self.sac_percut[cut].vardat, 
                    y=self.sac_percut[cut].fractional_sacrifice,
                    yerr=self.sac_percut[cut].fs_uncertainty,
                    linestyle='none', marker='o', label=cut, markersize=6,
                    elinewidth=3, capsize=0)
        if fittotal is True:
            print(self.sac_percut)
            popt, pcov = spc.curve_fit(self._flatline, self.sac_percut["total"].vardat, self.sac_percut["total"].fractional_sacrifice,
                    p0=[0.02], sigma=self.sac_percut["total"].fs_uncertainty)
            #one standard deviation
            print("THE PCOV: %s"%(str(pcov))) 
            try:
                stdev = np.sqrt(np.diag(pcov))
            except ValueError:
                print("Fit likely failed.  You may not have enough data"+\
                        "in this region to fit a straight line.")
                raise
            plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
        plt.legend(loc=3)
        plt.yscale("log")
        plt.ylabel("Fractional sacrifice",fontsize=34)
        varindict = False
        for var in self.xlabel_dict:
            if self.sac_percut_metadata["variable"] == var:
                plt.xlabel(self.xlabel_dict[var],fontsize=32)
                varindict = True
            else:
                if varindict is False: 
                    plt.xlabel(self.sac_percut_metadata["variable"],fontsize=32)
        plt.tick_params(labelsize=32)
        if title is None:
            plt.title("Fractional sacrifice due to classifiers",fontsize=36)
        else: 
            plt.title(title,fontsize=36)
        #plt.savefig("sacrifice_b14itr_%s.pdf"%(self.sac_percut_metadata["variable"]))
        variable = self.sac_percut_metadata["variable"]
        plt.savefig(savedir+"/ClassSac_%s.pdf"%(variable))
        #plt.show()
        plt.close()

