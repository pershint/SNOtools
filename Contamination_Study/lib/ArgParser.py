import argparse
import os

basepath = os.path.dirname(__file__)

#FIXME: set up a simpler argparser here.  you can choose a config file, override
#some of the configuration values with flags (changes dictionary values), and
parser = argparse.ArgumentParser(description='Parser to decide what analysis to do')
parser.add_argument('--nosave', dest='NOSAVE',action='store_true',
        help='Do not save any outputs; ensures no writing is done')
parser.add_argument('--debug', dest='debug',action='store_true',
        help='Run code in debug mode')
parser.add_argument('--jobnum', dest='JOBNUM', action='store',
        help='Specify this jobs number among others.  Will save results'+\
                'to ./output/results_jN')
parser.add_argument('--analysisdir', dest='ANALYSISDIR',action='store',
        type=str,help='Specify the directory where the analysis files'+\
                'are stored.  Will read all files ending with .ntuple.root'+\
                'from the directory. Default: ./ntuples/physics_data/')
parser.add_argument('--calibdir', dest='CALIBDIR',action='store',
        type=str,help='Specify the directory where the calibration files'+\
                'are stored.  Will read all files ending with .ntuple.root'+\
                'from the directory. Default: ./ntuples/N16/')
parser.add_argument('--resultdir', dest='RESULTDIR',action='store',
        type=str,help='specify the location and filename for results to be read'+\
                'or written to.  No job number support with this flag called.')
parser.add_argument('--sacrifice', dest='SACANALYSIS',action='store_true',
        help='Run the code that plots the correlations of different cuts/classifiers')
parser.add_argument('--bifurcate', dest='BIFURCATE',action='store_true',
        help='Run the bifurcation analysis on physics data.  Saves a bifurcation summary')
parser.add_argument('--contamination', dest='ESTIMATECONTAMINATION',
        action='store_true', help='Run the contamination estimation using'+\
                'bifurcation and sacrifice results.  Save a summary.')
parser.add_argument('--plots', dest='PLOTS', action='store_true',
        help='Show plots resulting from sacrifice and contamination studies.  If no sacrifice or contamination study, loads results from the result directory and plots what is available.')
parser.add_argument('--erange', dest='ERANGE', action='store',nargs='+',
        help='Specify an energy range to run all analyses over.  If running'+\
                '--plots only, will check range matches that in results'+\
                'directory (usage: --erange 2.0 5.0)')
parser.add_argument('--zrange', dest='ZRANGE', action='store',nargs='+',
        help='Specify an upper and lower zcut in cm range to perform the analysis'+\
                'over.  Applied in sacrifice and comtanimation studies.'+\
                '(usage: --zrange 600 -500)')

THISDIR = os.path.dirname(__file__)
pd_default = os.path.abspath(os.path.join(THISDIR,"..", "ntuples", "physics_data"))
cal_default = os.path.abspath(os.path.join(THISDIR,"..", "ntuples","calibration"))
parser.set_defaults(NOSAVE=False,SACANALYSIS=False,BIFURCATE=False,debug=False,
        ESTIMATECONTAMINATION=False,JOBNUM=0,PLOTS=False,erange=None,
        RESULTDIR=None,PHYSDIR=pd_default,
        CALIBDIR=cal_default,ZRANGE=None)
args = parser.parse_args()
