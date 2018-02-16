#This class takes in a list of N16 root files and outputs sacrifice histograms
#Associated with each one.  We're writing this to replace Sacrifice.cc, and have
#Everything united under one master Python code.  

class SacrificeEstimator(object):
    def __init__(self, rootfiles=[], config_dict={}, save_directory=None):
        self.rootfile_list = rootfiles
        self.cdict = config_dict
        self.sacrifice_histograms = []
        self.save_directory = save_directory

    def set_savedirectory(self,directory):
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
    
            h_DC_FracFlagged = ROOT.TH1D("h_DC_FracFlagged", "h_DC_FracFlagged", 58, 0.0, 11.6)
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
                if datatree.posr > self.cdict['posr']:
                    continue
                if datatree.fitValid is False:
                    continue
                if datatree.isN16 is False:
                    continue
                if ((~datatree.dcFlagged) & self.cdict["path_DC_DCmask"]) > 0:
                    continue
                if ((datatree.triggerWord) & self.cdict["path_DC_trigmask"]) > 0:
                    continue
                h_AllEvents.Fill(datatree.energy);
                if (((~datatree.dcFlagged) & self.cdict["cut_DCmask"]) or \
                        (datatree.triggerWord & self.cdict["cut_trigmask"])):
                    h_DC_FlaggedEvents.Fill(energy);
                if ((datatree.beta14 > self.cdict["b14_high"]) or \
                        (datatree.beta14 < self.cdict["b14_low"]) or \
                        (datatree.ITR < self.cdict["itr_low"])):
                    h_BI_FlaggedEvents.Fill(energy);
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

    def SaveHistograms(self):
        for j,rf in enumerate(self.rootfile_list):
            outfilename=rf.split("/")
            outfilename=outfilename[len(outfilename)-1].rstrip(".root")+"_sachists.root"
            outfiledir=self.save_directory+"/sachists/"+outfilename
            outfile = ROOT.TFile(outfiledir+"/"+outfilename,"CREATE")
            for histogram in self.sacrifice_histograms[j]:
                outfile.Add(histogram)
            outfile.Write()
            outfile.Close()

    def SaveConfiguration(self,config_outname):
        saveconfigloc = self.save_directory+"/"+config_outname
        json.dump(self.cdict, saveconfigloc, sort_keys=True,indent=4)
