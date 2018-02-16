#This class takes in a list of N16 root files and outputs sacrifice histograms
#Associated with each one.  We're writing this to replace Sacrifice.cc, and have
#Everything united under one master Python code.  

import ROOT
import copy
import os,sys

class SacrificeEstimator(object):
    def __init__(self, rootfiles=[], config_dict={}, save_directory=None):
        self.rootfile_list = rootfiles
        self.cdict = config_dict
        self.sacrifice_histograms = []
        self.save_directory = save_directory
        if not os.path.exists(self.save_directory):
            os.makedirs(self.save_directory)

    def set_savedirectory(self,directory):
        if not os.path.exists(directory):
            os.makedirs(directory)
        self.save_directory = directory

    def add_rootfile(self, rootfile):
        #Add a rootfile to the list of rootfiles to perform analysis on
        self.rootfile_list.append(rootfile)

    def clear_rootfiles(self):
        self.rootfile_list = []

    def load_configdict(self,config_file):
        #Lets you run the sacrifice with a new configuration dictionary
        ConParser = cp.ConfigParser(config_file)
        self.cdict = ConParser.parse_file()

    def EstimateSacrifice(self):
        #We will generate sacrifice histograms for each rootfile loaded in
        #Clears currenet sacrifice_histogram list first
        self.sacrifice_histograms = []
        for rf in self.rootfile_list:
            rootfile = ROOT.TFile(rf,"READ")
            h_AllEvents = ROOT.TH1D("h_AllEvents", "h_AllEvents", 58, 0.0, 11.6)
            h_DC_FracFlagged = ROOT.TH1D("h_DC_FracFlagged", "h_DC_FracFlagged", 58, 0.0, 11.6)
            h_DC_FlaggedEvents = ROOT.TH1D("h_DC_FlaggedEvents", "h_DC_FlaggedEvents", 58, 0.0, 11.6)
    
            h_BI_FlaggedEvents = ROOT.TH1D("h_BI_FlaggedEvents", "h_BI_FlaggedEvents", 58, 0.0, 11.6)
    
            h_BI_FracFlagged = ROOT.TH1D("h_BI_FracFlagged", "h_BI_FracFlagged", 58, 0.0, 11.6)
            h_AllEvents.Sumw2()
            h_DC_FracFlagged.Sumw2()
            h_DC_FlaggedEvents.Sumw2()
            h_BI_FracFlagged.Sumw2()
            h_BI_FlaggedEvents.Sumw2()
            #Need the tree that has the data, then xrange it
            rootfile.cd()
            datatree=rootfile.Get("output")
            for i in xrange(datatree.GetEntries()):
                datatree.GetEntry(i)
                if datatree.posr > self.cdict['r_cut']:
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
                    h_DC_FlaggedEvents.Fill(datatree.energy);
                if ((datatree.beta14 > self.cdict["cut2_b14_high"]) or \
                        (datatree.beta14 < self.cdict["cut2_b14_low"]) or \
                        (datatree.itr < self.cdict["cut2_itr_low"])):
                    h_BI_FlaggedEvents.Fill(datatree.energy);
            h_DC_FracFlagged.Divide(h_DC_FlaggedEvents,h_AllEvents,1.,1.,"b")
            h_BI_FracFlagged.Divide(h_BI_FlaggedEvents,h_AllEvents,1.,1.,"b")

            h_DC_FlaggedEvents.GetXaxis().SetTitle("Energy(MeV)")
            h_DC_FlaggedEvents.GetYaxis().SetTitle("Events")
            h_DC_FracFlagged.GetXaxis().SetTitle("Energy(MeV)")
            h_DC_FracFlagged.GetYaxis().SetTitle("Fractional Sacrifice")

            h_BI_FlaggedEvents.GetXaxis().SetTitle("Energy(MeV)")
            h_BI_FlaggedEvents.GetYaxis().SetTitle("Events")
            h_BI_FracFlagged.GetXaxis().SetTitle("Energy(MeV)")
            h_BI_FracFlagged.GetYaxis().SetTitle("Fractional Sacrifice")

            sachists_thisfile = [copy.deepcopy(h_AllEvents),copy.deepcopy(h_BI_FracFlagged),
                    copy.deepcopy(h_BI_FlaggedEvents),copy.deepcopy(h_DC_FracFlagged),
                    copy.deepcopy(h_DC_FlaggedEvents)]
            self.sacrifice_histograms.append(sachists_thisfile)
            del h_AllEvents,h_DC_FracFlagged,h_DC_FlaggedEvents,h_BI_FlaggedEvents,\
                    h_BI_FracFlagged
        

    def SaveHistograms(self):
        for j,rf in enumerate(self.rootfile_list):
            outfilename=rf.split("/")
            print("ROOTFILE,SPLIT IN ARRAY" + str(outfilename))
            outfilename=outfilename[len(outfilename)-1].rstrip(".root")+"_sachists.root"
            print("OUTFILENAME: " + outfilename)
            outfiledir=self.save_directory+"/sachists"
            if not os.path.exists(outfiledir):
                os.makedirs(outfiledir)
            outfile = ROOT.TFile(outfiledir+"/"+outfilename,"CREATE")
            for histogram in self.sacrifice_histograms[j]:
                outfile.Add(histogram)
            outfile.Write()
            outfile.Close()

