Code for analyzing the sacrifice of data cleaning cuts and
Fit classifiers for use in the bifurcation analysis.

First, define your cuts and regions of interest in the config file. The
different configurables have a description in cuts_default_notes.txt.  You can
also write your own configuration json, and use it using the --config flag.

Then, a sacrifice analysis, bifurcation analysis, and contamination analysis
 can be run using RunStudy.py.

Some configuration paramters can be adjusted via command line (energy is only
supported for now).  To see available parameters and usage, run

python RunStudy.py --help

By default, analyses run will be output to /output/results_j0 .  The job
number can be changed from the command line, or the output can be sent to
a different file all together.

Note that the python file has configurables inside it for looking at DC only,
B14 only, specific run numbers, etc.

