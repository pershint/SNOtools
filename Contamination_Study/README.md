Code for analyzing the sacrifice of data cleaning cuts and
Fit classifiers for use in the bifurcation analysis.

First, define your cuts and regions of interest in the config file. Then,
Run:

./bin/Sacrifice -i /path/to/yourntuple.root -o ./outroots/yourntuple_sacrifice.root

Then, python tools are available in ./python for making nice graphs of the
sacrifice as a function of position.  Implementation of sacrifice vs. time is still
ongoing now.

Note that the python file has configurables inside it for looking at DC only,
B14 only, specific run numbers, etc.

NOTE: There are preliminary cuts made before checking
that Beta14 and ITR are either good or bad.  These preliminary
cuts are as you defined in the config file.  Make sure you use the same config
file when running the bifurcation analysis on any physics runs. 
