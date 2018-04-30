Code for analyzing the sacrifice of data cleaning cuts and
Fit classifiers for use in the bifurcation analysis.


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

GUIDE TO A QUICK START OF ESTIMATING SACRIFICE/CONTAMINATION

Preparation steps

i) Define your cuts and regions of interest for estimating sacrifice/contamination
in a config file. The different configurables have a description in 
cuts_default_notes.txt.  A good example structure for the json format is in
cuts_default.json.  Tell the program to use your config file by giving the file 
location to the --config flag.

ii) Place your N16 calibration data and physics_data to estimate the contamination of
into the /ntuples/N16 and /ntuples/physics_data/ directories, respectively.

iii) For each N16 run used, add an entry to one of the N16 run info JSON files in
/DB/.  You can also write a new file with the same format as the others and it
will automatically be read in if the file name has the format /DB/N16_{whatever}.json

ANALYSES

Sacrifice: To calculate the sacrifice of physics events in N16 data, do:

python RunStudy.py --sacrifice

The program will take all files in the calibration directory either specified or
in the default location and estimate the sacrifice due to data cleaning cuts in
the specified energy region, z-region, and radial region.  

Results of the sacrifice for all runs and on a run-by-run detail are saved to the
specified results directory.  

Bifurcation: To generate the "pass:pass", "pass:fail", "fail:pass", and "fail:fail"
values for the data cleaning cuts and fit classifiers, do:

python RunStudy.py --bifurcation

The values for the total number of events analyzed, the total number of events that
passed the defined pathological cut masks, and the number of events that go into
each bifurcation box are output to "bifurcation_boxes.json" in the results directory.

Contamination estimate: First you must either:

1) Also give the --sacrifice and --contamination flags when performing this analysis
2) Give a --jobnum that corresponds to a directory in /output/results_j{JOBNUM} that
already has information from the sacrifice and contamination estimates present.

To run the contamination study, do:

python RunStudy.py --contamination

The contamination due to each cut branch is then calculated using the estimated
sacrifice and bifurcation box analyses.  A description of the algorithm for calculating
the y1 and y2 values ("fraction of instrumentals expected to leak into the ROI for
the classifiers and data cleaning cuts") is found in the SNO+ water phase data
cleaning document's data cleaning section.

Several approaches to calculating y1y2 are presented in the output of this study
"contamination_summary.json".  This is due to the fact that the solution to the
value for y1y2 is over-determined.  The different output values in the contam. JSON
are:

"highest_y1y2": Highest y1y2 value of the two over-determined equations
"lowest_y1y2": Lowest y1y2 value of the two over_determined equations
"avg_y1y2": Average value of the y1y2 values from the two over-determined equations
"leastsq_y1y2": Value found using the least-squares approach to finding the value
of y1y2 from the two over-constrained equations
"y1y2_to_CL": Gives an upper-limit on the value for y1y2 to the specified "CL" in 
the JSON file.  Found by re-shooting values for the a,b,c,d and acceptance values
based on their uncertainties and finding the spread of y1y2 due to these.
