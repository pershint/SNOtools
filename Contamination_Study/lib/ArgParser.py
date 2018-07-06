import argparse
import os

basepath = os.path.dirname(__file__)

#FIXME: set up a simpler argparser here.  you can choose a config file, override
#some of the configuration values with flags (changes dictionary values), and
parser = argparse.ArgumentParser(description='Sacrifice and contamination analysis code '+\
        'for SNO+ ntuple data.  Results of all analyses are loaded from/saved to '+\
        './output/results_j{JOBNUM}, where JOBNUM is specified with --jobnum flag. ')
parser.add_argument('--nosave', dest='NOSAVE',action='store_true',
        help='Do not save any outputs; ensures no writing is done')
parser.add_argument('--debug', dest='debug',action='store_true',
        help='Run code in debug mode')
parser.add_argument('--jobnum', dest='JOBNUM', action='store',
        help='Specify this jobs number among others.  Will save results'+\
                'to ./output/results_jN')
parser.add_argument('--lowE', dest='LOWECONTAM',action='store_true',
        help='Tell contamination study whether or not to use analysis for'+\
                'low energy contamination.  For ND, lowE best for E < 5.5 MeV')
parser.add_argument('--signalmcdir', dest='MCSIGNALDIR',action='store',
        type=str,help='Specify the directory where the signal monte carlo '+\
                'root files are stored.  For this analysis, signal="cherenkov events"'+\
                'Will read all files ending with .ntuple.root'+\
                'from the directory. Default: ./ntuples/signal_mc/')
parser.add_argument('--analysisdir', dest='ANALYSISDIR',action='store',
        type=str,help='Specify the directory where the analysis files'+\
                'are stored.  Will read all files ending with .ntuple.root'+\
                'from the directory. Default: ./ntuples/physics_data/')
parser.add_argument('--calibdir', dest='CALIBDIR',action='store',
        type=str,help='Specify the directory where the calibration files'+\
                'are stored.  Will read all files ending with .ntuple.root'+\
                'from the directory. Default: ./ntuples/N16/')
parser.add_argument('--calibsacrifice', dest='CALIBSACANALYSIS',action='store_true',
        help='Runs the sacrifice estimate processor on calibration data selected')
parser.add_argument('--signalmcsacrifice', dest='MCSACANALYSIS',action='store_true',
        help='Runs the sacrifice estimate processor on signal MC files')
parser.add_argument('--bifurcate', dest='BIFURCATE',action='store_true',
        help='Run the bifurcation analysis on files in --analysisdir.')
parser.add_argument('--contamination', dest='ESTIMATECONTAMINATION',
        action='store_true', help='Run the contamination estimation using'+\
                'bifurcation and sacrifice results.')
parser.add_argument('--erange', dest='ERANGE', action='store',nargs='+',
        help='Specify an energy range to run sacrifice/bifurcation analyses over.  If running'+\
                '--plots only, will check range matches that in results'+\
                'directory (usage: --erange 2.0 5.0)')
parser.add_argument('--zrange', dest='ZRANGE', action='store',nargs='+',
        help='Specify an upper and lower zcut in cm range to perform the '+\
                'sacrifice/bifurcation analyses over.'+\
                '(usage: --zrange 600 -500)')

THISDIR = os.path.dirname(__file__)
pd_default = os.path.abspath(os.path.join(THISDIR,"..", "ntuples", "physics_data"))
cal_default = os.path.abspath(os.path.join(THISDIR,"..", "ntuples","calibration"))
signalmc_def = os.path.abspath(os.path.join(THISDIR,"..","ntuples","signal_mc"))
parser.set_defaults(NOSAVE=False,debug=False,CALIBSACANALYSIS=False,
        MCSACANALYSIS=False,BIFURCATE=False,ESTIMATECONTAMINATION=False,
        JOBNUM=0,erange=None,
        ANALYSISDIR=pd_default, LOWECONTAM=False,MCSIGNALDIR=signalmc_def,
        CALIBDIR=cal_default,ZRANGE=None)
args = parser.parse_args()
