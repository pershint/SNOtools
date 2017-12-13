#Runs current .plotsac over all files in processing_roots
#TO USE THIS SCRIPT:
#  - Make sure only the root files you want to read through
#    are in ./processing_roots/
#  - Change the name you pipe to to whatever runs you are reading

FILEDIR=/home/pershint/recentrun_data/n16scan/ntuples_newqvt/
FILELIST=$(find $FILEDIR -name \*.nt.root | sort)
FILEARR=(${FILEDIR}*.nt.root)
for i in ${!FILEARR[@]};
do
    OUTFILE=${FILEARR[${i}]}_fracflags.root
    echo $OUTFILE
    echo ${FILEARR[${i}]}
    #Pass the list of files as argu{ments to be read from for the analysis
    ./CutVSNHit $OUTFILE ${FILEARR[${i}]} >> tester.out
done
