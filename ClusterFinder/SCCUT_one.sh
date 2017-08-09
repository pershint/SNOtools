#Runs current .plotsac over all files in procroots
#TO USE THIS SCRIPT:
#  - Make sure only the root files you want to read through
#    are in ./procroots/
#  - Change the name you pipe to to whatever runs you are reading

./SCCutRunner ./procroots/$1 water >> SCCUT_${1}.out

