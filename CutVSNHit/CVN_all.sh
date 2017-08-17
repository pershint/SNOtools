#Runs current .plotsac over all files in processing_roots
#TO USE THIS SCRIPT:
#  - Make sure only the root files you want to read through
#    are in ./processing_roots/
#  - Change the name you pipe to to whatever runs you are reading

FILEDIR=/home/maskins/golden_reproc/ntuple/
#FILEDIR=/home/pershint/snoing/install/rat-dev/tools/SNOtools/BifurAnalysis/procntuples/
FILELIST=$(find $FILEDIR -name \*.ntuple.root | sort)

#Pass the list of files as arguments to be read from for the analysis
./CutVSNHit $FILELIST >> tester.out
