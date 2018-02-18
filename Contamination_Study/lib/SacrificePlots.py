#Functions for plotting results from the SacrificeAnalyzer class.

def plot_sacrificevsZ(positions, fractional_sac, fractional_sac_unc):
    z_positions = []
    for position in positions:
        z_positions.append(position[2])
    #We've got our position and sacrifice information for this calibration set.
    #Now, just plot it.
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(z_positions,fractional_sac, xerr=0, yerr=fractional_sac_unc, \
            marker='o', linestyle='none', color = 'b', alpha=0.6, \
            label = 'Fractional Sacrifice')
    ax.set_xlabel(str(source) + " Z Position (m)")
    ax.set_ylabel("Fraction of events sacrificed")
    ax.set_title("Fractional sacrifice due to " + branch + " as source" + \
            " Z position varies\n" + \
            str(source) + " Source used")
    ax.grid(True)
    plt.show()

def plot_sacrificevsR(positions, fractional_sac, fractional_sac_unc):
    r_positions = []
    for position in positions:
        radius = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
        r_positions.append(radius)
    #We've got our position and sacrifice information for this calibration set.
    #Now, just plot it.
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(r_positions,fractional_sac, xerr=0, yerr=fractional_sac_unc, \
            marker='o', linestyle='none', color = 'b', alpha=0.6, \
            label = 'Fractional Sacrifice')
    ax.set_xlabel(str(source) + " Radius (m)")
    ax.set_ylabel("Fraction of events sacrificed")
    ax.set_title("Fractional sacrifice due to " + branch + " as source" + \
            " radial position varies\n" + \
            str(source) + " Source used")
    ax.grid(True)
    plt.show()

