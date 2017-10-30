The directories contained within are specifically for the contamination study used
in the pure water phase of SNO+ for nucleon decay.

PART 0: PREPPING YOUR N16 CALIBRATION ROOTS FOR USE IN THE ACCEPTANCE ANALYSIS
        TOOLS
  0a) For N16 calibration data specifically:
     Take your processed roots and run them through reducify_p1.cpp.  This will
     output a text file with all events in the roots that are FECD-tagged. Then,
     pass the root through reducify.cpp.  The output will be an ntuple of the
     run which also has an integer indicating if an event is an N16 event (1 it is,
     0 it is not).
  0b)In Acceptance_Analysis_Tools/Physics_Acceptance/DB, add a json that has
     the N16 run information for the date you took your data on.  Follow the
     format of N16_Positions_1.json.
  0b)The filenames output from step 0a) should have names like 
     "RUNNUM_N16.root".  Make sure that the run number in the filename matches 
     the run number you put in the database entries in step 0b).

BEFORE GOING ON: Ensure that the pathological cuts are defined the same in all
three codes used in the following three steps.  Also make sure that your branches
of cuts (DC and Fit classifiers) are defined as you want in the DCSacrifice and
B14ITR_Sacrifice programs, respectively.

  PART 1: USING THE ACCEPTANCE ANALYSIS TOOLS TO CALCULATE YOUR ACCEPTANCE RATES

  1) Run calibration data ntuples through the code in ./DC_sacrifice
     to calculate your acceptance of good physics according
     to the DC cuts used to produce the file. Output should be a ROOT file with
     histograms detailing the sacrifice vs. nhit.
  2) Run calibration data ntuples through the code in ./Beta14ITR_sacrifice
     to calculate your acceptance of good physics according to the fit 
     classifiers. Output should be a ROOT file with histograms detailing the
     sacrifice vs. nhit.
  3) Move the output ROOTS from above to the ./Acceptance_Analysis_Tools/N16_DCsac
     and ""/N16_Fitsacs directories.  Now, run /Acceptance_Analysis_Tools/
     Physics_Acceptance/main.py.  Output will be plots and your acceptance with
     weighted mean uncertainties with variation in z.

  4) Add your acceptance rates to a dictionary as is already done at the top of
     ./BifurAnalysis/BifurcatedAnalysis/main.py.  You will want the code to use
     your acceptance values, so change as needed.

  PART 2: RUN THE BIFURCATED ANALYSIS 
  3) Pass all of the ntuples that are associated with your physics data 
     to the code in ./BifurAnalysis/Bifurcator to populate your
     Pass-fail boxes used in the bifurcated analysis contamination estimation.
     Pipe the output from this to "your_filename.out".
     NOTE: You can also select your cut mask and whether to only use B14 or ITR.
     These should come immediately following the program call.  That is, place
     all filenames to pass in at the end of the program call.

  4) Move the .out files into the "./results" directory; you may want to
     make a new subdirectory that will hold the .out files.

  PART 3: CALCULATE THE LEAKAGE RATES IN THE SAMPLE ANALYZED IN PART 2 USING
  THE ACCEPTANCES FROM PART 1

  5) Run ./BifurAnalysis/BifurcatedAnalysis/main.py.  
     The script will estimate the errors on your
     acceptance rates output from steps 1 and 2, and then apply those along with
     poissonian errors for a,b,c, and d to estimate your contamination rate and
     uncertainty for each cut.
     
     #ADD MORE DETAIL ABOUT GENERATING THE .out FILES USING THE -b, -i, or
     -d OPTIONS.  ADDING THE OUTPUTS OF THESE INTO ./results WOULD BE USED IN
     FINDING THE CORRELATIONS OF EACH FILE.
