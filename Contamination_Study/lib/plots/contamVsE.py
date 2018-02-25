import matplotlib.pyplot as plt

def ContaminationVsEnergy():
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.errorbar(bin_centers, contam_binned, xerr=bin_widths, \
            yerr=contam_binned_uncs, marker = 'o', linestyle='none',\
            color='r', alpha=0.8)
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("# Instrumentals in region")
    ax.set_title("Estimated # Non-cherenkov events that leak through \n"+\
            "Data cleaning and Fits in 11 days of data taking")
    ax.grid(True)
    plt.show()
