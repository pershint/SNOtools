#Runs current .plotsac over all files in processing_roots
#TO USE THIS SCRIPT:
#  - Make sure only the root files you want to read through
#    are in ./processing_roots/
#  - Change the name you pipe to to whatever runs you are reading

for file in ./procroots/*.root
do
    ./ReadDC $file water >> newresults_r10545.out
done
