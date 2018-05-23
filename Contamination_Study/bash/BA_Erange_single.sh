#/bin/bash -l

#This bash script feeds in an array of energy ranges and starts Contamination Study jobs for each one.


RUNDIR=/home/onetrueteal/Programs/SNOrat/rat/tools/SNOtools/Contamination_Study/
ROOTSRC=/home/onetrueteal/Programs/SNOrat/rat/rat-dev.sh
CALIBDIR=/home/onetrueteal/share/May2016_N16_2
DATADIR=/home/onetrueteal/share_sg/lowEContamination_afterThresh
ELOW="4.0"
EHIGH="5.5"
source $ROOTSRC
cd $RUNDIR
python main.py --sacrifice --bifurcate --contamination --calibdir ${CALIBDIR} --analysisdir ${DATADIR} --erange $ELOW $EHIGH --source N16 --jobnum 6
