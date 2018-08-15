Code for analyzing the sacrifice of data cleaning cuts and
Fit classifiers for use in the bifurcation analysis.


Then, a sacrifice analysis, bifurcation analysis, and contamination analysis
 can be run using main.py.

Some configuration paramters can be adjusted via command line (energy is only
supported for now).  To see available parameters and usage, run

python main.py --help

By default, analyses run will be output to /output/results_j0 .  The job
number can be changed from the command line, or the output can be sent to
a different file all together.

GUIDE TO A QUICK START OF ESTIMATING DC SACRIFICE

Preparation steps

i) Define your cuts and regions of interest for estimating sacrifice/contamination
in a config file. The different configurables have a description in 
cuts_default_notes.txt.  A good example structure for the json format is in
cuts_default.json.  Tell the program to use your config file by giving the file 
name to the --configfile flag.

ii) Place all of your calibration data ntuples into a single directory.  When
running main.py, you will want to give the directory path using the --calibdir flag. 

iii) For the instrumental contamination estimate, place all of your calibration
data ntuples into another single directory.  When running main.py, you will want to
give this directory path using the --analysis dir flag

iv) Set up your plotting preferences and variables you want to calculate the
DC sacrifice distribution for using ./config/setup.json.  For the DC sacrifice,
you want "CUTSTODO" to have ["cut1"] only, and common sacrifice variables used
for the water analysis are "posr3", "udotr", and "energy".  Note that if you
want to use "posr3", your data must have the "posr3" entry added to the ntuples
using the subtupler code (see ./misc/subtupleMaker\_N16.cpp for an example
of the subtuple processor you could use).


ANALYSES

Sacrifice: To calculate the sacrifice of physics events in N16 data, do:

python main.py --calibsacrifice --calibdir /path/to/calibntuples/ --configfile
theconfigfile.json --plots


Results of the sacrifice for all runs and on a run-by-run detail are saved to the
specified results directory. If --plots is given, the distribution of the 
DC sacrifice for each variable defined in setup.json will be shown. 

Bifurcation: To generate the "pass:pass", "pass:fail", "fail:pass", and "fail:fail"
values for the data cleaning cuts and fit classifiers, do:

python main.py --analysisdir /path/to/physics_data/ --configfile
theconfigfile.json --bifurcate 

The values for the total number of events analyzed, the total number of events that
passed the defined pathological cut masks, and the number of events that go into
each bifurcation box are output to "bifurcation_boxes.json" in the results directory output.

Contamination estimate: First you must either:

1) Also give the --calibsacrifice and --bifurcate flags when performing this analysis
2) Give a --jobnum that corresponds to a directory in /output/results_j{JOBNUM} that
already has information from the sacrifice and contamination estimates present.

To run the contamination study, do:

python main.py --contamination

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

A second approache to calculating the instrumental contamination in the presence
of a large amount of signal (as is common for background ROIs at lower energies
with high amounts of Cherenkov events) can be done with

python main.py --contamination --lowE


