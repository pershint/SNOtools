#This program is a one-shot program that outputs the results of the SNO+
#Conamination study.  First, the data cleanning and fit sacrifice are estimated
#Using the given calibration files.  Then, the bifurcation analysis is run over
#the given data physics files.  Finally, an estimate on the contamination for the
#chosen configuration of cuts/ROI choice is output.

import ROOT
import numpy as np
import lib.Bifurcator as bi
import lib.SacHists as s
import lib.configparser as cp
import os,sys

MAINDIR = os.path.dirname(__file__)
SAVEDIR = os.path.abspath(os.path.join(MAINDIR, "output"))
CALIBDIR = os.path.abspath(os.path.join(MAINDIR, "data", "N16"))
PHYSDIR = os.path.abspath(os.path.join(MAINDIR, "data", "OpenGolden"))
CONFIGDIR = os.path.abspath(os.path.join(MAINDIR, 'config'))

#FIXME: set up a simpler argparser here.  you can choose a config file, override
#some of the configuration values with flags (changes dictionary values), and
#Could have flags to only run sacrifice, bifurcation, bifurcation analysis, or all.
CONFIGFILE='cuts_default.ini'

if __name__ == '__main__':
    #Load the configuration file
    ConfigParser = cp.ConfigParser(CONFIGDIR+"/"+CONFIGFILE)
    config_dict = ConfigParser.Parse()
    print(config_dict)
