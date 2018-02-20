#Utilities for loading results from directories

import json

def LoadJson(result_directory,filename):
    loadfileloc = result_directory+"/"+filename
    with open(loadfileloc,"r") as f:
        return json.load(f)

