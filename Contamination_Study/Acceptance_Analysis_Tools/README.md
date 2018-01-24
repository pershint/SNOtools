Beta14ITR_Sacrifice takes in an N16 calibration file and
returns a root file that contains histograms.
These histograms are the clean and dirty events in the
Event as a function of NHit.

The output root file should be moved to ../../N16_Fitsacs
If you want to use it with the Physics_Acceptance
Calculating python script.

DCSacrifice takes in an N16 calibration file and
returns a root file that contains histograms.
These histograms are the clean and dirty events in the
Event as a function of NHit.

The output here should be moved to ../N16_DCSacs

NOTE: There are preliminary cuts made before checking
that Beta14 and ITR are either good or bad.  These preliminary
cuts should be the same as what you find in the Bifurcated
Analysis program and the DCSacrifice program.
