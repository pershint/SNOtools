#Runs current .plotsac over all files in processing_roots
#TO USE THIS SCRIPT:
#  - Make sure only the root files you want to read through
#    are in ./processing_roots/
#  - Change the name you pipe to to whatever runs you are reading

for file in ./processing_roots/*.root
do
    ./SCCutRunner $file water >> SCCUTFPLOOK.out
done
