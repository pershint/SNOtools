#!/bin/bash -l

#This bash script feeds in an array of energy ranges and starts Contamination Study jobs for each one.


RUNDIR=/home/onetrueteal/Programs/SNOrat/rat/tools/SNOtools/Contamination_Study/
ROOTSRC=/home/onetrueteal/Programs/SNOrat/rat/rat-dev.sh
CALIBDIR=/home/onetrueteal/share/May2016_N16_2
DATADIR=/home/onetrueteal/share_sg/lowEContamination_afterThresh
source $ROOTSRC
cd $RUNDIR

ENERGY_ARRAY=("1.5" "2.0" "2.5" "3.0" "3.5" "4.0" "4.5" "5.0" "5.5") #This format to feed values to BA_Single.sh
for maskindex in ${!ENERGY_ARRAY[@]};
do
  echo $maskindex
  if (($maskindex>0));
    then
      python main.py --sacrifice --bifurcate --contamination --calibdir ${CALIBDIR} --analysisdir ${DATADIR} --erange ${ENERGY_ARRAY[(($maskindex-1))]} ${ENERGY_ARRAY[$maskindex]} --source N16 --jobnum $maskindex
      ((currenti++))
    fi
done;
