The following directory contains tools for processing a SNO+ ntuple file
with several possible applications in mind:
  - Shifting variables according to necessary corrections in analysis
  - Correcting variables with RAT processors
  - Defining new entries in the output ntuple file

There are several .cpp files here that can provide a guideline for 
writing your own subtupler.  To compile, modify the Makefile in this directory
to compile your subtupler, then in the command line in this directory do:

make

The subtupleMakerScript.sh script provides a simple script for running your
compiled subtupler on all files in a directory.  Alternatively, you can also
give a list of files directly to the subtupler program by doing

./subtupler --data file1 file2 file3 --out outputfile.subtuple.root
