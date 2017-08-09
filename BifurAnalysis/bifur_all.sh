#Runs current .plotsac over all files in processing_roots
#TO USE THIS SCRIPT:
#  - Make sure only the root files you want to read through
#    are in ./processing_roots/
#  - Change the name you pipe to to whatever runs you are reading

FILEDIR=/home/pershint/snoing/install/rat-dev/tools/BifurAnalysis/procroots/ntuples/
FILELIST=$(find $FILEDIR -name \*.ntuple.root | sort)

#Pass the list of files as arguments to be read from for the analysis
./BifurcatedAnalysis $FILELIST >> tester.out
