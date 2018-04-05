import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.optimize as spc
import seaborn as sns
import numpy as np
import glob
import json

def weighted_stdev(vals,valavg, uncs):
    weights = 1/(uncs**2)
    return np.sqrt((len(weights)*np.sum(weights*(vals-valavg)**2))/((len(weights)-1)*np.sum(weights)))


def flatline(x,b):
    return b

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
    sns.set_style("whitegrid")
    sns.axes_style("whitegrid")
    xkcd_colors = ['slate blue', 'fluro green', 'twilight', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','yellowgreen']
    d = cut_sacrifice_byrun[cut]
    if cut == 'cut1':
        colors = ['twilight']
    if cut == 'cut2':
        colors = ['leaf']
    sns.set_palette(sns.xkcd_palette(colors))#,len(allcutsacs)))
    fs, radius, fs_u = [], [], []
    for j,posn in enumerate(d['position']):
        if posn[2] >600.0:
            continue
        fs.append(d['fractional_sac'][j])
        radius.append(np.sqrt(posn[0]**2 + posn[1]**2 + posn[2]**2))
        fs_u.append(d['fractional_sac_unc'][j])
    fs = np.array(fs)
    r3R3 = np.array(radius)**3/np.full(len(radius),600)**3
    fs_u = np.array(fs_u)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(radius,fs,yerr=fs_u,xerr=0,marker='o',markersize=8,
            linestyle='none',elinewidth=3, capsize=0)
    plt.tick_params(labelsize=20)
    plt.yscale('log')
    ax.set_xlabel(r'$R_{s}^{3}/R_{AV}^{3}$',fontsize=22)
    ax.set_ylabel("Fractional sacrifice",fontsize=22)
    if cut == 'cut1':
        clabel = 'Data Cleaning'
    if cut == 'cut2':
        clabel = 'Fit Classifier'
    ax.set_title("Fractional sacrifice for "+clabel+" as radial position of"+\
            " source varies",fontsize=24)
    plt.grid(True)
    print("AVERAGE VALUE: " + str(np.average(fs)))
    print("STD DEV.: " + str(np.std(fs)))
    plt.show()

def PlotAxis(cut_sacrifice_byrun,cut,axis):
    sns.set_style("whitegrid")
    sns.axes_style("whitegrid")
    xkcd_colors = ['slate blue', 'fluro green', 'twilight', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','yellowgreen']
    d = cut_sacrifice_byrun[cut]
    if cut == 'cut1':
        colors = ['slate blue']
    if cut == 'cut2':
        colors = ['leaf']
    sns.set_palette(sns.xkcd_palette(colors))#,len(allcutsacs)))
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
    if cut == 'cut1':
        clabel = 'Data Cleaning'
    if cut == 'cut2':
        clabel = 'Fit Classifier'
    plt.tick_params(labelsize=18)
    plt.yscale('log')
    ax.errorbar(posns,fs,yerr=fs_u,xerr=0,marker='o',markersize=5,
            linestyle='none',elinewidth=3,capsize=0)
    ax.set_xlabel("Position of source on "+axis+"-axis (cm)",fontsize=22)
    ax.set_ylabel("Fractional sacrifice",fontsize=22)
    ax.set_title("Fractional sacrifice for "+clabel+" as "+axis+"-position of"+\
            " source varies",fontsize=24)
    plt.grid(True)
    print("AVERAGE VALUE: " + str(np.average(fs)))
    print("STD DEV.: " + str(np.std(fs)))
    plt.show()

def PlotSys(cut_sacrifice_byrun,cut,axis):
    sns.set_style("whitegrid")
    sns.axes_style("whitegrid")
    sns.set_context("poster",font_scale = 2)
    xkcd_colors = ['slate blue', 'fluro green', 'twilight', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','yellowgreen']
    d = cut_sacrifice_byrun[cut]
    if cut == 'cut1':
        colors = ['twilight']
    if cut == 'cut2':
        colors = ['light eggplant']
    sns.set_palette(sns.xkcd_palette(colors))#,len(allcutsacs)))
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
    #Fit a flat line to the fractional sacrifices
    popt, pcov = spc.curve_fit(flatline, posns, fs,p0=[np.average(fs)], sigma=fs_u)
    #one standard deviation
    print("PCOVARIANCE: " + str(pcov))
    stdev = np.sqrt(np.diag(pcov))
    #plot it
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    if cut == 'cut1':
        clabel = 'Data Cleaning'
    if cut == 'cut2':
        clabel = 'Fit Classifier'
    #plt.tick_params(labelsize=18)
    plt.yscale('log')
    ax.errorbar(posns,fs,yerr=fs_u,xerr=0,marker='o',markersize=10,
            linestyle='none',elinewidth=6,capsize=0, label='Sacrifice by run')
    ax.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,\sigma = %f$' % (float(popt[0]), float(stdev[0])))
    ax.set_xlabel("Position of source on "+axis+"-axis (cm)")
    ax.set_ylabel("Fractional sacrifice")
    ax.set_title("Fractional sacrifice for "+clabel+" as "+axis+"-position of"+\
            " source varies")
    plt.legend()
    plt.grid(True)
    print("AVERAGE SACRIFICE OF RUNS, NO WEIGHT: " + str(np.average(fs)))
    print("STD DEV. OF SACRIFICE, NO WEIGHTS: " + str(np.std(fs)))
    print("BEST FIT SACRIFICE: " + str(popt))
    print("UNC. OF BEST FIT AVERAGE: " + str(stdev))
    print("WEIGHTED STDEV: " + str(weighted_stdev(fs,popt[0],fs_u)))
    plt.show()

