import matplotlib.pyplot as plt
import numpy as np
import glob
import json

def PlotContamVsEnergy(DIR):
    #Gets all contamination results from sub-directories inside parent directory
    #DIR.  Plots them as a function of energy range as found in used_configuration.json
    E_lows = []
    E_highs = []
    E_mids = []
    Contamination = []
    result_dirs = glob.glob(DIR+"*")
    for d in result_dirs:
        with open(d+'/used_configuration.json','r') as c:
            config = json.load(c)
        with open(d+'/contamination_summary.json','r') as cs:
            contamsum = json.load(cs)
        E_lows.append(config['E_low'])
        E_highs.append(config['E_high'])
        E_mids.append((config['E_low']+config['E_high'])/2.0)
        Contamination.append(contamsum['y1y2_to_CL']*contamsum['est_bkg_evts'])
    E_lows=np.array(E_lows)
    E_highs=np.array(E_highs)
    E_mids=np.array(E_mids)
    Contamination=np.array(Contamination)
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(E_mids, Contamination, xerr=(E_lows-E_highs)/2.0, \
            yerr=0, marker = 'o', markersize=7,linestyle='none',\
            color='r', alpha=0.8,linewidth=4)
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("# Instrumentals in region")
    ax.set_title("Estimated # Non-cherenkov events that leak through \n"+\
            "Data cleaning and Fits in 11 days of data taking to 90% CL")
    ax.grid(True)
    plt.show()

def PlotContamApproaches(DIR):
    #Gets all contamination results from sub-directories inside parent directory
    #DIR.  Plots them as a function of energy range as found in used_configuration.json
    E_lows = []
    E_highs = []
    E_mids = []
    avg_y1y2 = []
    CL_y1y2 = []
    highest_y1y2 = []
    leastsq_y1y2 = []
    result_dirs = glob.glob(DIR+"*")
    for d in result_dirs:
        with open(d+'/used_configuration.json','r') as c:
            config = json.load(c)
        with open(d+'/contamination_summary.json','r') as cs:
            contamsum = json.load(cs)
        E_lows.append(config['E_low'])
        E_highs.append(config['E_high'])
        E_mids.append((config['E_low']+config['E_high'])/2.0)
        CL_y1y2.append(contamsum['y1y2_to_CL'])
        avg_y1y2.append(contamsum['avg_y1y2'])
        highest_y1y2.append(contamsum['highest_y1y2'])
        leastsq_y1y2.append(contamsum['leastsq_y1y2'])
    E_lows=np.array(E_lows)
    E_highs=np.array(E_highs)
    E_mids=np.array(E_mids)
    CL_y1y2=np.array(CL_y1y2)
    avg_y1y2=np.array(avg_y1y2)
    highest_y1y2=np.array(highest_y1y2)
    leastsq_y1y2=np.array(leastsq_y1y2)
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    colors = ['b','r','g','k']
    labels = ['90% CL','leastsq','highest','average']
    y1y2s = [CL_y1y2,leastsq_y1y2,highest_y1y2,avg_y1y2]
    for j,y1y2 in enumerate(y1y2s):
        ax.errorbar(E_mids, y1y2, xerr=(E_lows-E_highs)/2.0, \
                yerr=0, marker = 'o', markersize=7,linestyle='none',\
                color=colors[j], alpha=0.8,linewidth=4,label=labels[j])
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("Value for y1*y2")
    ax.set_title("Values for the product of the leakage fraction for each "+\
            "bifurcation branch calculated in different ways")
    ax.grid(True)
    ax.legend()
    plt.show()

def PlotRadius(cut_sacrifice_byrun,cut):
    d = cut_sacrifice_byrun[cut]
    if cut == 'cut1':
        col = 'b'
    if cut == 'cut2':
        col = 'r'
    fs, radius, fs_u = [], [], []
    for j,posn in enumerate(d['position']):
        if posn[2] >600.0:
            continue
        fs.append(d['fractional_sac'][j])
        radius.append(np.sqrt(posn[0]**2 + posn[1]**2 + posn[2]**2))
        fs_u.append(d['fractional_sac_unc'][j])
    fs = np.array(fs)
    radius = np.array(radius)
    fs_u = np.array(fs_u)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(radius,fs,yerr=fs_u,xerr=0,marker='o',markersize=8,
            color=col,linestyle='none')
    ax.set_xlabel("Radial position of source (cm)",fontsize=16)
    ax.set_ylabel("Fractional sacrifice",fontsize=16)
    ax.set_title("Fractional sacrifice for "+cut+" branch as radial position of"+\
            " source varies",fontsize=18)
    plt.grid(True)
    print("AVERAGE VALUE: " + str(np.average(fs)))
    print("STD DEV.: " + str(np.std(fs)))
    plt.show()

def PlotAxis(cut_sacrifice_byrun,cut,axis):
    d = cut_sacrifice_byrun[cut]
    if cut == 'cut1':
        col = 'b'
    if cut == 'cut2':
        col = 'r'
    fs, posns, fs_u = [], [], []
    for j,posn in enumerate(d['position']):
        if posn[2] > 600.0:
            continue
        if axis == "X":
            fs.append(d['fractional_sac'][j])
            posns.append(posn[0])
            fs_u.append(d['fractional_sac_unc'][j])
        if axis == "Y":
            fs.append(d['fractional_sac'][j])
            posns.append(posn[1])
            fs_u.append(d['fractional_sac_unc'][j])
        if axis == "Z":
            fs.append(d['fractional_sac'][j])
            posns.append(posn[2])
            fs_u.append(d['fractional_sac_unc'][j])
    fs = np.array(fs)
    posns = np.array(posns)
    fs_u = np.array(fs_u)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(posns,fs,yerr=fs_u,xerr=0,marker='o',markersize=8,
            color=col,linestyle='none')
    ax.set_xlabel("Position of source on "+axis+"-axis (cm)",fontsize=16)
    ax.set_ylabel("Fractional sacrifice",fontsize=16)
    ax.set_title("Fractional sacrifice for "+cut+" branch as "+axis+"-position of"+\
            " source varies",fontsize=18)
    plt.grid(True)
    print("AVERAGE VALUE: " + str(np.average(fs)))
    print("STD DEV.: " + str(np.std(fs)))
    plt.show()

def PlotSys(cut_sacrifice_byrun,cut,axis):
    d = cut_sacrifice_byrun[cut]
    if cut == 'cut1':
        col = 'b'
    if cut == 'cut2':
        col = 'r'
    fs, posns, fs_u = [], [], []
    for j,posn in enumerate(d['position']):
        if posn[2] > 600.0:
            continue
        if axis == "X":
            if abs(posn[0]) > 50.0:
                if np.sqrt((posn[1]**2 + posn[2]**2)) < 100.0:
                    fs.append(d['fractional_sac'][j])
                    posns.append(posn[0])
                    fs_u.append(d['fractional_sac_unc'][j])
        if axis == "Y":
            if abs(posn[1]) > 50.0:
                if np.sqrt((posn[0]**2 + posn[2]**2)) < 100.0:
                    fs.append(d['fractional_sac'][j])
                    posns.append(posn[1])
                    fs_u.append(d['fractional_sac_unc'][j])
        if axis == "Z":
            if abs(posn[2]) > 50.0:
                if np.sqrt((posn[0]**2 + posn[1]**2)) < 100.0:
                    fs.append(d['fractional_sac'][j])
                    posns.append(posn[2])
                    fs_u.append(d['fractional_sac_unc'][j])
    fs = np.array(fs)
    posns = np.array(posns)
    fs_u = np.array(fs_u)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(posns,fs,yerr=fs_u,xerr=0,marker='o',markersize=8,
            color=col,linestyle='none')
    ax.set_xlabel("Position of source on "+axis+"-axis (cm)",fontsize=16)
    ax.set_ylabel("Fractional sacrifice",fontsize=16)
    ax.set_title("Fractional sacrifice for "+cut+" branch as "+axis+"-position of"+\
            " source varies",fontsize=18)
    plt.grid(True)
    print("AVERAGE VALUE: " + str(np.average(fs)))
    print("STD DEV.: " + str(np.std(fs)))
    plt.show()

