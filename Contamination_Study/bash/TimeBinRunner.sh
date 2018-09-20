#!/bin/bash -l

#This bash script feeds in an array of energy ranges and starts Contamination Study jobs for each one.

ROITYPE=ExtNonPMT/sacrifice/inAV
FILEEXTENS=extnonPMT
RUNDIR=/home/onetrueteal/Programs/SNOrat/rat/tools/SNOtools/Contamination_Study
ROOTSRC=/home/onetrueteal/Programs/SNOrat/rat/rat-dev.sh
CALIBDIR=/home/onetrueteal/share/SNOPlus/N16/Data/Nov2017_N16/subtuple
CALIBMCDIR=/home/onetrueteal/share_wd/ND_AnalysisData/BenN16MC/Nov2017_IntScan/subtuple
source $ROOTSRC
cd $RUNDIR

TBARRAY=("2" "3" "4" "5" "6") #This format to feed values to BA_Single.sh
#TBARRAY=("1") #This format to feed values to BA_Single.sh
for maskindex in ${!TBARRAY[@]};
do
  echo $maskindex
  if ((${TBARRAY[maskindex]}==2));
    then
      python main.py --configfile ND_configs/${ROITYPE}/cuts_N16sac_${FILEEXTENS}TB${TBARRAY[maskindex]}_lowhem.json --calibsacrifice --calibdir_data ${CALIBDIR} --calibdir_mc ${CALIBMCDIR} --jobname ${TBARRAY[maskindex]}lo --plots
      python main.py --configfile ND_configs/${ROITYPE}/cuts_N16sac_${FILEEXTENS}TB${TBARRAY[maskindex]}_uphem.json --calibsacrifice --calibdir_data ${CALIBDIR} --calibdir_mc ${CALIBMCDIR} --jobname ${TBARRAY[maskindex]}up --plots
  fi
  if ((${TBARRAY[maskindex]} !=2 ));
    then
      python main.py --configfile ND_configs/${ROITYPE}/cuts_N16sac_${FILEEXTENS}TB${TBARRAY[maskindex]}.json --calibsacrifice --calibdir_data ${CALIBDIR} --calibdir_mc ${CALIBMCDIR} --jobname ${TBARRAY[maskindex]} --plots
  fi
done;
