import sys
import ROOT

#Python Class that takes in a configuration file and a list of root files.
#The bifurcation analysis is run over the files, and a root file is output with
#A "bifurcated" tree.  The tree has all the ntuple information for events that
#Pass the preliminary cuts defined in the config file.

class Bifurcator(object):
    def __init__(self,rootfiles=[], config_dict={}, save_directory=None):
        self.rootfile_list = rootfiles
        self.cdict = config_dict
        self.sacrifice_histograms = []
        self.save_directory = save_directory
        if not os.path.exists(self.save_directory):
            os.makedirs(self.save_directory)
        self.bifurcation_rootfile = None

    def set_savedirectory(self,directory):
        if not os.path.exists(directory):
            os.makedirs(self.save_directory)
        self.save_directory = directory

    def add_rootfile(self, rootfile):
        #Add a rootfile to the list of rootfiles to perform analysis on
        self.rootfile_list.append(rootfile)

    def clear_rootfiles(self):
        self.rootfile_list = []

    def load_configdict(self,config_file):
        #Lets you run the bifurcation analysis with a new configuration dictionary
        ConParser = cp.ConfigParser(config_file)
        self.cdict = ConParser.parse_file()

    def Bifurcate(self):
        self.bifurcation_rootfile = None
        bfile = ROOT.TFile(self.save_directory+"/bifurcation_result.ntuple.root","CREATE")
        b_root = ROOT.TTree("output","Events involved with bifurcation analysis")
        result = ROOT.TTree("boxes","values for a,b,c, and d boxes")
        #Initialize variables to fill with events passing preliminary cuts
        runID = np.zeros(1,dtype=int)
        ratversion = np.zeros(1,dtype=str)
        nhits = np.zeros(1,dtype=int)
        nhitsCleaned = np.zeros(1,dtype=int)
        triggerWord = np.zeros(1,dtype=int)
        tubiiWord = np.zeros(1,dtype=int)
        dcApplied = np.zeros(1,dtype=float64)
        dcFlagged = np.zeros(1,dtype=float64)
        uTDays = np.zeros(1,dtype=int)
        uTSecs = np.zeros(1,dtype=int)
        uTNSecs = np.zeros(1,dtype=int)
        fitValid = np.zeros(1,dtype=bool)
        energy = np.zeros(1,dtype=float64)
        itr = np.zeros(1,dtype=float64)
        beta14 = np.zeros(1,dtype=float64)
        cut1BranchFail = np.zeros(1,dtype=bool)
        cut2BranchFail = np.zeros(1,dtype=bool)

        a = np.zeros(1,dtype=int)
        b = np.zeros(1,dtype=int)
        c = np.zeros(1,dtype=int)
        d = np.zeros(1,dtype=int)


        result.Branch('a', a, 'a/I')
        result.Branch('b', b, 'b/I')
        result.Branch('c', c, 'c/I')
        result.Branch('d', d, 'd/I')

        #Now, associate initialized variables with the TTree
        b_root.Branch('runID',runID, 'runID/I')
        b_root.Branch('ratversion',ratversion, 'ratversion/C')
        b_root.Branch('nhits',nhits, 'nhits/I')
        b_root.Branch('nhitsCleaned',nhitsCleaned, 'nhitsCleaned/I')
        b_root.Branch('triggerWord', triggerWord, 'triggerWord/I')
        b_root.Branch('tubiiWord', tubiiWord, 'tubiiWord/I')
        b_root.Branch('dcApplied', dcApplied, 'dcApplied/D')
        b_root.Branch('dcFlagged', dcFlagged, 'dcFlagged/D')
        b_root.Branch('uTDays',uTDays, 'uTDays/I')
        b_root.Branch('uTSecs',uTSecs, 'uTSecs/I')
        b_root.Branch('uTNSecs',uTNsecs, 'uTNsecs/I')
        b_root.Branch('fitValid',fitValid, 'fitValid/O')
        b_root.Branch('energy',energy, 'energy/D')
        b_root.Branch('itr', itr, 'itr/D')
        b_root.Branch('beta14', beta14, 'beta14/D')
        b_root.Branch('cut1BranchFail',cut1BranchFail, 'cut1BranchFail/O')
        b_root.Branch('cut2BranchFail',cut2BranchFail, 'cut2BranchFail/O')

        a[0],b[0],c[0],d[0] = 0,0,0,0
        for rf in self.rootfile_list:
            rootfile=ROOT.TFile(rf,"READ")
            rootfile.cd()
            datatree=rootfile.Get("output")
            for i in xrange(datatree.GetEntries()):
                datatree.GetEntry(i)
                if datatree.posr > self.cdict['r_cut']:
                    continue
                if datatree.fitValid is False:
                    continue
                if datatree.energy > self.cdict['E_high'] or \
                        datatree.energy < self.cdict['E_low']:
                    continue
                if ((~datatree.dcFlagged) & self.cdict["path_DCmask"]) > 0:
                    continue
                if ((datatree.triggerWord) & self.cdict["path_trigmask"]) > 0:
                    continue
                cut1_clean = True
                cut2_clean = True
                if self.cdict["dcs_only"] is false and self.cdict["fits_only"] is false:
                    #We're using a dc branch and a fit branch
                    if (((~datatree.dcFlagged) & self.cdict["cut1_DCmask"]) > 0):
                        cut1_clean = False
                    if ((datatree.trigWord & self.cdict["cut1_trigmask"]) > 0):
                        cut1_clean = False
                    if (datatree.beta14 < self.cdict["cut2_b14_low"] or \
                            self.cdict["cut2_b14_high"] < datatree.beta14):
                        cut2_clean = False
                    if (datatree.itr < self.cdict["cut2_itr_low"]):
                        cut2_clean = False
                else:
                    print("OPTIONS FOR DC AND FIT ONLY NOT YET IMPLEMENTED.")
                    print("JUST LOOK AT BIFURCATOR.CC AND ADD THE LOGIC.")
                    sys.exit()
                #Increment the box values and save the ntuple info
                if cut1_clean and cut2_clean:
                    a[0]+=1
                elif cut1_clean and not cut2_clean:
                    b[0]+=1
                elif cut2_clean and not cut1_clean:
                    c[0]+=1
                elif not cut2_clean and not cut1_clean:
                    d[0]+=1
                runID[0] = datatree.runID
                ratversion[0] = datatree.ratversion
                nhits[0] = datatree.nhits
                nhitsCleaned[0] = datatree.nhitsCleaned
                triggerWord[0] = datatree.triggerWord
                tubiiWord[0] = datatree.tubiiWord
                dcApplied[0] = datatree.dcApplied
                dcFlagged[0] = datatree.dcFlagged
                uTDays[0] = datatree.uTDays
                uTSecs[0] = datatree.uTSecs
                uTNSecs[0] = datatree.uTNSecs
                fitValid[0] = datatree.fitValid
                energy[0] = datatree.energy
                itr[0] = datatree.itr
                beta14[0] = datatree.beta14
                cut1BranchClean[0] = cut1_clean
                cut2BranchClean[0] = cut2_clean
                b_root.Fill()
            result.Fill() #saves a, b, c, and d
            self.bifurcation_rootfile = bfile

    def SaveBifurcationRoot(self):
        self.bifurcation_rootfile.Write()
        self.bifurcation_rootfile.Close()
