#When run in a directory containing all timebin results (outputs from 
#main.py), produces a table showing the data cleaning sacrifice and
#Data/MC ratio for each timebin.  

import numpy as np
import os,sys,time
import glob
import json

basepath = os.path.dirname(__file__)
MAINDIR = os.path.abspath(basepath)

class LaTeXTable(object):
    def __init__(self):
        self.tlines = []

    def ClearTable(self):
        self.tlines = []

    def LoadTable(self,table):
        '''Takes in a filename and loads the file into the table lines'''
        self.tlines = table

    def AddTableLine(self,line):
        self.tlines.append(line)

    def AddCommentLine(self,string):
        self.tlines.append("%% %s \n"%(string))

    def SaveTable(self,name):
        f = open(name,'w')
        for line in self.tlines:
            f.write(line)
        f.close()

class SacClassTable(LaTeXTable):
    def __init__(self,datadicts=None):
        super(SacClassTable,self).__init__()
        self.datadicts = datadicts
    
    def EditSysValues(self):
        '''Takes values found in the datadicts object and
        updates the systematic uncertianties for lines by using
        self.datadict's keys identified with first object in line'''
        for j,line in enumerate(self.tlines):
            for key in self.datadicts:
                if line.find(key) != -1:
                    #Found our string. Break line into entries
                    sline = line.split(" & ")
                    #Dc/class values are 3rd and 4th index elements
                    thislinesdata=self.datadicts[key]
                    Classinfo = thedata['cut2_DataMCComp']
                    classcomps = []
                    comp_stat_uncs = []
                    comp_sys_uncs = []
                    for var in Classinfo:
                        classcomps.append(Classinfo[var]['Data/MC ratio'])
                        comp_stat_uncs.append(Classinfo[var]['stat_unc'])
                        comp_sys_uncs.append(Classinfo[var]['sys_unc'])
                    classcomps = np.array(classcomps)
                    comp_stat_uncs = np.array(comp_stat_uncs)
                    comp_sys_uncs = np.array(comp_sys_uncs)
                    avg_comp = np.round(np.average(classcomps),4)
                    stat_unc = np.sqrt(np.sum(sac_stat_uncs**2))
                    sys_unc = np.sqrt(np.sum(sac_sys_uncs**2))
                    tot_unc = np.round(np.sqrt(stat_unc**2 + sys_unc**2),4)
                    sline[3] = avg_comp
                    sline[4] = tot_unc
                    #Then, rejoin
                    newtline = " & ".join(sline)
                    print(newtline)
                    #when it works...
                    #self.tline[j] = newtline

    def MakeTitle(self,sepstatsys=False):
        '''Add the title line that I want for the tables given to'''
        '''the SNO+ collaborators'''
        if sepstatsys:
            thetitle="TimeBin & DC Sacrifice & DC Sac. Unc. & "+\
                    "Data/MC Class Comp. & Data/MC Class Comp. unc \\\\ \n"
        else:
            thetitle="TimeBin & DC Sacrifice & DC Sac stat Unc. & "+\
                    "DC Sac sys unc & Data/MC Class Comp. & Data/MC stat unc & "+\
                    "Data/MC sys unc\\\\ \n"
        self.tlines.append(thetitle)

    def WriteTableLine(self,key,sepstatsys=False):
        '''For a given key, searches the datadict and writes the
        data line with the average DC/Classifier info and the 
        uncertainties added in quadrature'''
        thedata = None
        theline = "%s & "%(key)
        try:
            thedata = self.datadicts[key]
        except KeyError:
            print("The key given does not exist in the data "+\
                    "dictionary given to the class.")
            return
        DCinfo = thedata['cut1']
        sacrifices = []
        sac_stat_uncs = []
        sac_sys_uncs = []
        for var in DCinfo:
            sacrifices.append(DCinfo[var]['sacrifice'])
            sac_stat_uncs.append(DCinfo[var]['stat_unc'])
            sac_sys_uncs.append(DCinfo[var]['sys_unc'])
        sacrifices = np.array(sacrifices)
        sac_stat_uncs = np.array(sac_stat_uncs)
        sac_sys_uncs = np.array(sac_sys_uncs)
        avg_sac = np.round(np.average(sacrifices),4)
        stat_unc = np.sqrt(np.sum(sac_stat_uncs**2))
        sys_unc = np.sqrt(np.sum(sac_sys_uncs**2))
        tot_unc = np.round(np.sqrt(stat_unc**2 + sys_unc**2),4)
        if sepstatsys:
            theline += "%s & %s & %s & "%(str(avg_sac),str(stat_unc),
                    str(sys_unc))
        else:
            theline += "%s & %s & "%(str(avg_sac),str(tot_unc))
        Classinfo = thedata['cut2_DataMCComp']
        classcomps = []
        comp_stat_uncs = []
        comp_sys_uncs = []
        for var in Classinfo:
            classcomps.append(Classinfo[var]['Data/MC ratio'])
            comp_stat_uncs.append(Classinfo[var]['stat_unc'])
            comp_sys_uncs.append(Classinfo[var]['sys_unc'])
        classcomps = np.array(classcomps)
        comp_stat_uncs = np.array(comp_stat_uncs)
        comp_sys_uncs = np.array(comp_sys_uncs)
        avg_comp = np.round(np.average(classcomps),4)
        stat_unc = np.sqrt(np.sum(sac_stat_uncs**2))
        sys_unc = np.sqrt(np.sum(sac_sys_uncs**2))
        tot_unc = np.round(np.sqrt(stat_unc**2 + sys_unc**2),4)
        if sepstatsys:
            theline += "%s & %s & %s & "%(str(avg_comp),str(stat_unc),
                    str(sys_unc))
        else:
            theline += "%s & %s & "%(str(avg_comp),str(tot_unc))
        self.tlines.append(theline)

def getTBRowLabel(tbstr):
    '''Given what number is at the end of the results, return
    the TB label to be put in the latex table row label'''
    thelabel = None
    if tbstr=='1':
        thelabel = 'TB1'
    elif tbstr=='2lo':
        thelabel = 'TB2 (z<0)'
    elif tbstr=='2up':
        thelabel = 'TB2 (z>0)'
    elif tbstr=='3':
        thelabel = 'TB3'
    elif tbstr=='4':
        thelabel = 'TB4'
    elif tbstr=='5':
        thelabel = 'TB5'
    elif tbstr=='6':
        thelabel = 'TB6'
    else:
        print("job suffix given not recognized...")
    return thelabel

def getAllDirJsons(jsonname):
    '''Returns a dictionary where the keys are TB#, and
    the values are the JSON loaded as a dict for each entry'''
    alldirjsons = {}
    dirnames = glob.glob("%s/*"%(MAINDIR))
    indicestodelete = []
    for j,entry in enumerate(dirnames):
        print("ENTRY IS: " + str(entry))
        istablemaker = entry.find('results_j')
        print(istablemaker)
        if istablemaker==-1:
           indicestodelete.append(j)
    for i in sorted(indicestodelete)[::-1]:
        del dirnames[i]

    print("DIRNAMES AFTER DELETE: " + str(dirnames))
    TBstrs = ['1','2lo','2up','3','4','5','6']
    for name in dirnames:
        theTBstr = None
        name = name.replace(MAINDIR+"/","")
        for j,TB in enumerate(TBstrs):
            istimebin = name.find(TB)
            if istimebin!=-1:
                theTBstr = TB
                break
        theRowLabel = getTBRowLabel(theTBstr)
        with open("%s/%s"%(dirnames[j],jsonname),"r") as f:
            thisTBjson = json.load(f)
            alldirjsons[theRowLabel] = thisTBjson
    return alldirjsons

if __name__=='__main__':
    EDITSYS = True
    SYSFILETOEDIT = "testing.tex"
    print("HERE'S OUR MAIN CODE")
    alldirjsons = getAllDirJsons("calib_DCSacClassComp_totals.json")
    print(alldirjsons)
    if EDITSYS:
        #use alldirjsons, but load in the table we want to modify
        ETable = SacClassTable(datadicts=alldirjsons)
        ETable.LoadLines(SYSFILETOEDIT)
        ETable.EditSysValues()
    else:
        #That will be each column in the table
        TheTable = SacClassTable(datadicts=alldirjsons)
        TheTable.AddCommentLine("DC/Class corrections for AV Box Timebins")
        TheTable.MakeTitle(sepstatsys=True)
        for key in alldirjsons:
            TheTable.WriteTableLine(key,sepstatsys=True)
    TheTable.SaveTable(SYSFILETOEDIT)
