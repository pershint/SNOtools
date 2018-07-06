#Utilities for saving and loading results from directories

import json

def save_sacrifice_list(directory,fullpath_list,filename):
    #Strips off the directory, places the root names in a list, and saves it
    #in the RESULTDIR
    physics_roots = []
    for fullname in fullpath_list:
        rootname = fullname.split("/")
        rootname = rootname[len(rootname)-1]
        physics_roots.append(rootname)
    analysis_list = {}
    analysis_list["runs_used_in_sac_estimate"] = physics_roots
    with open(directory+"/"+filename,"w") as f:
        json.dump(analysis_list,f,sort_keys=True,indent=4)

def save_bifurcation_list(directory,fullpath_list,filename):
    #Strips off the directory, places the root names in a list, and saves it
    #in the RESULTDIR
    physics_roots = []
    for fullname in fullpath_list:
        rootname = fullname.split("/")
        rootname = rootname[len(rootname)-1]
        physics_roots.append(rootname)
    analysis_list = {}
    analysis_list["runs_used_in_bifurcation"] = physics_roots
    with open(directory+"/"+filename,"w") as f:
        json.dump(analysis_list,f,sort_keys=True,indent=4)

def LoadJson(result_directory,filename):
    loadfileloc = result_directory+"/"+filename
    with open(loadfileloc,"r") as f:
        return json.load(f)

