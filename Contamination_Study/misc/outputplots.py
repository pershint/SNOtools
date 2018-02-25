import matplotlib.pyplot as plt
import numpy as np

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

