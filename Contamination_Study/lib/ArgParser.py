import argparse
import os

basepath = os.path.dirname(__file__)

#FIXME: set up a simpler argparser here.  you can choose a config file, override
#some of the configuration values with flags (changes dictionary values), and
parser = argparse.ArgumentParser(description='Sacrifice and contamination analysis code '+\
        'for SNO+ ntuple data.  Results of all analyses are loaded from/saved to '+\
        './output/results_j{JOBNAME}, where JOBNAME is specified with --jobname flag. ')
parser.add_argument('--nosave', dest='NOSAVE',action='store_true',
        help='Do not save any outputs; ensures no writing is done')
parser.add_argument('--debug', dest='debug',action='store_true',
        help='Run code in debug mode')
parser.add_argument('--showplots', dest='SHOWPLOTS',action='store_true',
        help='Show sacrifice distribution and fit plots')
parser.add_argument('--saveplots', dest='SAVEPLOTS',action='store_true',
        help='Save sacrifice distribution and fit plots')
parser.add_argument('--jobname', dest='JOBNAME', action='store',
        help='Specify this jobs number among others.  Will save results'+\
                'to ./output/results_jN')
parser.add_argument('--configfile', dest='CONFIGFILE',action='store',
        type=str,help='Specify the JSON file in ./config/ that has all cut and ROI"+\
                "selection desired for analysis.  Default is "cuts_default.json"')
parser.add_argument('--LETA', dest='LETACONTAM',action='store_true',
        help='Tell contamination study whether or not to use analysis used '+\
                'in the SNO LETA approach.')
parser.add_argument('--analysisdir', dest='ANALYSISDIR',action='store',
        type=str,help='Specify the directory where the analysis files'+\
                'are stored.  Will read all files ending with .ntuple.root'+\
                'from the directory. Default: ./ntuples/physics_data/')
parser.add_argument('--sacdir_data', dest='SACRIFICEDIR',action='store',
        type=str,help='Specify the directory where the files that will be analyzed'+\
                ' for sacrifice are stored.  Will read all files ending with .ntuple.root'+\
                'from the directory. Default: ./ntuples/sacrifice/data/')
parser.add_argument('--sacdir_mc', dest='SACRIFICEMCDIR',action='store',
        type=str,help='Specify the directory where the simulated sacrifice files '+\
                'are stored.  Will read all files ending with .ntuple.root'+\
                'from the directory. Default: ./ntuples/sacrifice/MC/')
parser.add_argument('--sacrifice', dest='SACANALYSIS',action='store_true',
        help='Runs the sacrifice estimate processor on the data directory selected')
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
calmc_default = os.path.abspath(os.path.join(THISDIR,"..", "ntuples", "calibration", "mc"))
parser.set_defaults(NOSAVE=False,debug=False,CALIBSACANALYSIS=False,
        MCSACANALYSIS=False,BIFURCATE=False,ESTIMATECONTAMINATION=False,
        JOBNAME="0",erange=None,ANALYSISDIR=pd_default,CONFIGFILE='cuts_default.json', 
        LETACONTAM=False,CALIBDIR=cal_default,CALIBMCDIR=calmc_default,ZRANGE=None,PLOTS=False)
args = parser.parse_args()
