#Functions are used to get the pass-pass,pass-fail,fail-pass, and fail-fail
#Values that are output from ./BifurcatedAnalysis.  When BifurcatedAnalysis is
#run, send the results to filename.out.  This code loads in files of that output
#format and grabs the values.


import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import lib.playDarts as pd

def GetResultDict(filename):
    boxes = ["a","b","c","d"]
    result_dict = {}
    f = open(filename,"r")
    for line in f:
        #Get the pass-fail values
        for box in boxes:
            if box+":" in line:
                value = line.split(": ")
                result_dict[box] = value[1]
        if "Data cleaning mask" in line:
            value = line.split(": ")
            result_dict["DC_mask"] = value[1]
        if "Used ITR only?" in line:
            value = line.split("?")
            result_dict["ITR_only"] = value[1]
        if "Used B14 only?" in line:
            value = line.split("?")
            result_dict["B14_only"] = value[1]
    return result_dict

if __name__ == "__main__":
    print("NO MAIN CALL IMPLEMENTED IN THIS LIBRARY")
