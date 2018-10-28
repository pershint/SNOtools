#!/bin/bash -l

#This script gives an example of how to load in a config file, select a directory
#of data, and run the sacrifice analyzer

RUNDIR=/home/onetrueteal/Programs/SNOrat/rat/tools/SNOtools/Contamination_Study
RATSRC=/home/onetrueteal/Programs/SNOrat/rat/rat-dev.sh
SACRIFICEDIR=${RUNDIR}/ntuples/FirstPartialFill
CALIBMCDIR=/home/onetrueteal/share_wd/ND_AnalysisData/BenN16MC/Nov2017_IntScan/subtuple
source $RATSRC
cd $RUNDIR

python main.py --configfile cuts_PartialFill.json --sacrifice --sacdir_data ${SACRIFICEDIR} --jobname 200032_partialsac --showplots --saveplots --debug
