#This class takes in a list of N16 root files and outputs sacrifice histograms
#Associated with each one.  We're writing this to replace Sacrifice.cc, and have
#Everything united under one master Python code.  

import ROOT
import copy
import os,sys

class SacrificeHistGen(object):
    def __init__(self, rootfiles=[], config_dict={},sourcetype=None):
        self.sourcetype = sourcetype
        self.rootfile_list = rootfiles
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
        for rf in self.rootfile_list:
            rootfile = ROOT.TFile(rf,"READ")
            h_AllEvents = ROOT.TH1D("h_AllEvents", "h_AllEvents", self.nbins,
                    self.elow, self.ehigh)
            h_cut1_FracFlagged = ROOT.TH1D("h_cut1_FracFlagged", "h_cut1_FracFlagged", self.nbins,self.elow,self.ehigh)
            h_cut1_FlaggedEvents = ROOT.TH1D("h_cut1_FlaggedEvents", "h_cut1_FlaggedEvents", self.nbins,self.elow,self.ehigh)
    
            h_cut2_FlaggedEvents = ROOT.TH1D("h_cut2_FlaggedEvents", "h_cut2_FlaggedEvents", self.nbins,self.elow,self.ehigh)
    
            h_cut2_FracFlagged = ROOT.TH1D("h_cut2_FracFlagged", "h_cut2_FracFlagged", self.nbins,self.elow,self.ehigh)
            h_AllEvents.Sumw2()
            h_cut1_FracFlagged.Sumw2()
            h_cut1_FlaggedEvents.Sumw2()
            h_cut2_FracFlagged.Sumw2()
            h_cut2_FlaggedEvents.Sumw2()
            #Need the tree that has the data, then xrange it
            rootfile.cd()
            datatree=rootfile.Get("output")
            for i in xrange(datatree.GetEntries()):
                datatree.GetEntry(i)
                if datatree.posr > self.cdict['r_cut']:
                    continue
                if datatree.energy > self.cdict["E_high"] or datatree.energy < \
                        self.cdict["E_low"]:
                    continue
                if datatree.fitValid is False:
                    continue
                if datatree.isCal is False:
                    continue
                if ((~datatree.dcFlagged) & self.cdict["path_DCmask"]) > 0:
                    continue
                if ((datatree.triggerWord) & self.cdict["path_trigmask"]) > 0:
                    continue
                h_AllEvents.Fill(datatree.energy);
                if (((~datatree.dcFlagged) & self.cdict["cut1_DCmask"]) or \
                        (datatree.triggerWord & self.cdict["cut1_trigmask"])):
                    h_cut1_FlaggedEvents.Fill(datatree.energy);
                if ((datatree.beta14 > self.cdict["cut2_b14_high"]) or \
                        (datatree.beta14 < self.cdict["cut2_b14_low"]) or \
                        (datatree.itr < self.cdict["cut2_itr_low"])):
                    h_cut2_FlaggedEvents.Fill(datatree.energy);
            h_cut1_FracFlagged.Divide(h_cut1_FlaggedEvents,h_AllEvents,1.,1.,"b")
            h_cut2_FracFlagged.Divide(h_cut2_FlaggedEvents,h_AllEvents,1.,1.,"b")

            h_cut1_FlaggedEvents.GetXaxis().SetTitle("Energy(MeV)")
            h_cut1_FlaggedEvents.GetYaxis().SetTitle("Events")
            h_cut1_FracFlagged.GetXaxis().SetTitle("Energy(MeV)")
            h_cut1_FracFlagged.GetYaxis().SetTitle("Fractional Sacrifice")

            h_cut2_FlaggedEvents.GetXaxis().SetTitle("Energy(MeV)")
            h_cut2_FlaggedEvents.GetYaxis().SetTitle("Events")
            h_cut2_FracFlagged.GetXaxis().SetTitle("Energy(MeV)")
            h_cut2_FracFlagged.GetYaxis().SetTitle("Fractional Sacrifice")

            sachists_thisfile = [copy.deepcopy(h_AllEvents),copy.deepcopy(h_cut2_FracFlagged),
                    copy.deepcopy(h_cut2_FlaggedEvents),copy.deepcopy(h_cut1_FracFlagged),
                    copy.deepcopy(h_cut1_FlaggedEvents)]
            self.sacrifice_histograms.append(sachists_thisfile)
            del h_AllEvents,h_cut1_FracFlagged,h_cut1_FlaggedEvents,h_cut2_FlaggedEvents,\
                    h_cut2_FracFlagged
        

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

