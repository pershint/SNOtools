#Functions for plotting results from the SacrificeAnalyzer class.
import matplotlib.pyplot as plt
import numpy as np

def plot_XYSacrifice(cut_acceptances_byrun, cut):
    cdict = cut_acceptances_byrun[cut]
    x_points = []
    y_points = []
    fracsac = []
    for j,positions in enumerate(cdict["position"]):
        if abs(positions[2]) < 100.0:
            x_points.append(positions[0])
            y_points.append(positions[1])
            fracsac.append(cdict["fractional_sac"][j])
    x_points = np.array(x_points)
    y_points = np.array(y_points)
    fracsac = np.array(z_points)

    plt.contourf(x_points,y_points,fracsac)
    plt.colorbar()
    plt.xlabel("x-position of source (cm)")
    plt.ylabel("y-position of source (cm)")
    plt.title("Heatmap of fractional sacrifice for "+cut+" branch at"+\
            "abs(z) < 100 cm")
    plt.show()


def plot_sacrificevsCart(cut_acceptances_byrun, cut,axis):
    cdict = cut_acceptances_byrun
    if axis=="x":
        ax=0
    if axis=="y":
        ax=1
    if axis=="z":
        ax=2
    positions = []
    for position in cdict[cut]["position"]:
        positions.append(position[ax])
    fractional_sac = cdict[cut]["fractional_sac"]
    fractional_sac_unc = cdict[cut]["fractional_sac_unc"]
    #We've got our position and sacrifice information for this calibration set.
    #Now, just plot it.
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(positions,fractional_sac, xerr=0, yerr=fractional_sac_unc, \
            marker='o', linestyle='none', color = 'b', alpha=0.6, \
            label = 'Fractional Sacrifice')
    ax.set_xlabel(str(cdict["sourcetype"])+" "+axis+" Position (cm)")
    ax.set_ylabel("Fraction of events sacrificed")
    ax.set_title("Fractional sacrifice due to " + cut + " branch as source " + \
            axis+" position varies\n" + \
            cdict["sourcetype"] + " Source used")
    ax.grid(True)
    plt.show()

def plot_sacrificevsR(cut_acceptances_byrun, cut):
    cdict = cut_acceptances_byrun
    r_positions = []
    for position in cdict[cut]["position"]:
        radius = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
        r_positions.append(radius)
    fractional_sac = cdict[cut]["fractional_sac"]
    fractional_sac_unc = cdict[cut]["fractional_sac_unc"]
    #We've got our position and sacrifice information for this calibration set.
    #Now, just plot it.
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(r_positions,fractional_sac, xerr=0, yerr=fractional_sac_unc, \
            marker='o', linestyle='none', color = 'b', alpha=0.6, \
            label = 'Fractional Sacrifice')
    ax.set_xlabel(str(cdict["sourcetype"]) + " Radius (cm)")
    ax.set_ylabel("Fraction of events sacrificed")
    ax.set_title("Fractional sacrifice due to " + cut + "branch as source" + \
            " radial position varies\n" + cdict["sourcetype"] + \
            " Source used")
    ax.grid(True)
    plt.show()

