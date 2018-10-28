#!/bin/bash -l

#This bash script feeds in an array of energy ranges and starts Contamination Study jobs for each one.

TB=6

ROITYPE=NDROI/contamination
FILEEXTENS=ND
RUNDIR=/home/onetrueteal/Programs/SNOrat/rat/tools/SNOtools/Contamination_Study
#CALIBDIR=/home/onetrueteal/share/SNOPlus/N16/Data/Nov2017_N16/subtuple
#CALIBMCDIR=/home/onetrueteal/share_wd/ND_AnalysisData/BenN16MC/Nov2017_IntScan/subtuple
ANALYSISDIR=/home/onetrueteal/share_wd/ND_AnalysisData/TB${TB}/subtuple/
cd $RUNDIR


if ((${TB}==2));
  then
    python main.py --jobname ${TB}lo --analysisdir ${ANALYSISDIR} --bifurcate --debug
    python main.py --jobname ${TB}up --analysisdir ${ANALYSISDIR} --bifurcate --debug
  fi
if ((${TB} !=2 ));
  then
    python main.py --jobname ${TB} --analysisdir ${ANALYSISDIR} --bifurcate --debug 
  fi
done;
