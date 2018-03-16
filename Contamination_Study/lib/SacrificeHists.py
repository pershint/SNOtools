#This class takes in a list of N16 root files and outputs sacrifice histograms
#Associated with each one.  We're writing this to replace Sacrifice.cc, and have
#Everything united under one master Python code.  

import ROOT
from ROOT import gDirectory
import copy
import os,sys

class SacrificeHistGen(object):
    def __init__(self, rootfiles=[], config_dict={},sourcetype=None,zcut=None):
        self.sourcetype = sourcetype
        self.rootfile_list = rootfiles
        self.zcut = zcut
        self.nbins = 14
        self.elow = config_dict["E_low"]
        self.ehigh = config_dict["E_high"]
        self.cdict = config_dict
        self.sacrifice_histograms = []
        self.histogram_files = []

    def add_rootfile(self, rootfile):
        #Add a rootfile to the list of rootfiles to perform analysis on
        self.rootfile_list.append(rootfile)

    def clear_rootfiles(self):
        self.rootfile_list = []

    def load_configdict(self,config_file):
        #Lets you run the sacrifice with a new configuration dictionary
        ConParser = cp.ConfigParser(config_file)
        self.cdict = ConParser.parse_file()

    def GenerateHistograms(self):
        #We will generate sacrifice histograms for each rootfile loaded in
        #Clears currenet sacrifice_histogram list first
        self.sacrifice_histograms = []
        basecuts = []
        if self.cdict['r_cut'] is not None:
            basecuts.append("posr<"+str(self.cdict['r_cut']))
        if self.cdict['E_high'] is not None: 
            basecuts.append("energy<"+str(self.cdict['E_high']))
        if self.cdict['E_low'] is not None: 
            basecuts.append("energy>"+str(self.cdict['E_low']))
        if self.cdict['Z_low'] is not None:
            basecuts.append("posz>"+str(self.cdict['Z_low']*10.0))
        if self.cdict['Z_high'] is not None:
            basecuts.append("posz<"+str(self.cdict['Z_high']*10.0))
        basecuts.append("fitValid==1")
        basecuts.append("isCal==1")
        basecuts.append("((dcFlagged&%s)==%s)" % (self.cdict["sacpath_DCmask"],\
                self.cdict["sacpath_DCmask"]))
        basecuts.append("((triggerWord&%s)==0)" % (self.cdict["path_trigmask"]))
        cut1list = []
        cut1list.append("((dcFlagged&%s)!=%s)" % (self.cdict["cut1_DCmask"],\
                self.cdict["cut1_DCmask"]))
        cut1list.append("((triggerWord&%s)!=0)" % (self.cdict["cut1_trigmask"]))
        cut2list = []
        cut2list.append("beta14>%s" % (self.cdict["cut2_b14_high"]))
        cut2list.append("beta14<%s" % (self.cdict["cut2_b14_low"]))
        cut2list.append("itr<%s" % (self.cdict["cut2_itr_low"]))
        base = self.cutfuse(basecuts, "&&")
        b1 = base + "&&(%s)" % self.cutfuse(cut1list,"||")
        b2 = base + "&&(%s)" % self.cutfuse(cut2list,"||")
        for rf in self.rootfile_list:
            rootfile = ROOT.TFile(rf,"READ")
            #Need the tree that has the data, then xrange it
            h_cut1_FracFlagged = ROOT.TH1D("h_cut1_FracFlagged", "h_cut1_FracFlagged", self.nbins,self.elow,self.ehigh)
            h_cut2_FracFlagged = ROOT.TH1D("h_cut2_FracFlagged", "h_cut2_FracFlagged", self.nbins,self.elow,self.ehigh)
            rootfile.cd()
            datatree=rootfile.Get("output")
            hist_dim = "%d,%f,%f" % (self.nbins,self.elow,self.ehigh)
            datatree.Draw("energy>>h_AllEvents("+hist_dim+")","isCal==1","goff")
            h_AllEvents = gDirectory.Get("h_AllEvents")
            datatree.Draw("energy>>h_AllNonpathEvents("+hist_dim+")",base,"goff")
            h_AllNonpathEvents = gDirectory.Get("h_AllNonpathEvents")
            datatree.Draw("energy>>h_cut1_FlaggedEvents("+hist_dim+")",b1,"goff")
            h_cut1_FlaggedEvents = gDirectory.Get("h_cut1_FlaggedEvents")
            datatree.Draw("energy>>h_cut2_FlaggedEvents("+hist_dim+")",b2,"goff")
            h_cut2_FlaggedEvents = gDirectory.Get("h_cut2_FlaggedEvents")
            h_cut1_FlaggedEvents.Sumw2()
            h_AllEvents.Sumw2()
            h_AllNonpathEvents.Sumw2()
            h_cut2_FlaggedEvents.Sumw2()
            h_cut1_FracFlagged.Divide(h_cut1_FlaggedEvents,h_AllNonpathEvents,1.,1.,"b")
            h_cut2_FracFlagged.Divide(h_cut2_FlaggedEvents,h_AllNonpathEvents,1.,1.,"b")

            h_cut1_FlaggedEvents.GetXaxis().SetTitle("Energy(MeV)")
            h_cut1_FlaggedEvents.GetYaxis().SetTitle("Events")
            h_cut1_FracFlagged.GetXaxis().SetTitle("Energy(MeV)")
            h_cut1_FracFlagged.GetYaxis().SetTitle("Fractional Sacrifice")

            h_cut2_FlaggedEvents.GetXaxis().SetTitle("Energy(MeV)")
            h_cut2_FlaggedEvents.GetYaxis().SetTitle("Events")
            h_cut2_FracFlagged.GetXaxis().SetTitle("Energy(MeV)")
            h_cut2_FracFlagged.GetYaxis().SetTitle("Fractional Sacrifice")

            sachists_thisfile = [copy.deepcopy(h_AllEvents),copy.deepcopy(h_AllNonpathEvents),
                    copy.deepcopy(h_cut2_FracFlagged),
                    copy.deepcopy(h_cut2_FlaggedEvents),copy.deepcopy(h_cut1_FracFlagged),
                    copy.deepcopy(h_cut1_FlaggedEvents)]
            self.sacrifice_histograms.append(sachists_thisfile)
            del h_AllEvents,h_AllNonpathEvents,h_cut1_FracFlagged,h_cut1_FlaggedEvents,h_cut2_FlaggedEvents,\
                    h_cut2_FracFlagged
        
    def cutfuse(self,stringlist,delim):
        outstring = ""
        for j,s in enumerate(stringlist):
            if j == 0:
                outstring = s
            else:
                outstring=outstring+delim+s
        return outstring

    def SaveHistograms(self,savedir):
        self.histogram_files = []
        for j,rf in enumerate(self.rootfile_list):
            outfilename=rf.split("/")
            print("ROOTFILE,SPLIT IN ARRAY" + str(outfilename))
            outfilename=outfilename[len(outfilename)-1].rstrip(".root")+"_sachists.root"
            print("OUTFILENAME: " + outfilename)
            outfiledir=savedir+"/sachists"
            if not os.path.exists(outfiledir):
                os.makedirs(outfiledir)
            outfile = ROOT.TFile(outfiledir+"/"+outfilename,"CREATE")
            for histogram in self.sacrifice_histograms[j]:
                outfile.Add(histogram)
            outfile.Write()
            outfile.Close()
            self.histogram_files.append(outfiledir+"/"+outfilename)

