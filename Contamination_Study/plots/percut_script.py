import sacrifice_percut_plot as spp
import glob
import matplotlib.pyplot as plt
import time

if __name__=='__main__':
    dcm = 0x11E7E
    udotr_av = '(posx*dirx + posy*diry + (posz-108.0)*dirz)/sqrt(posx**2 + ' + \
            "(posz-108.0)**2)"
    udotr_av_cut = udotr_av + "> (1.0 - 12 *(((sqrt((posx**2 + posy**2 + (posz-108.0)**2))/"\
            "6000.0)**3)-0.69)**2)"
    r_av3 = "(sqrt((posx**2 + posy**2 + (posz-108.0)**2))/6000.0)**3"
    n16 = glob.glob("/home/onetrueteal/share/N16/Nov2017_N16/N16/*.ntuple.root") 
    pre = "energy>3.0&&energy<5.6&&0.69<%s&&0.90>%s&&%s"%(r_av3, r_av3, udotr_av_cut)
    
    dat, meta = spp.SacVSVar_Data(rootfiles=n16, var = udotr_av, 
            precuts = pre, nbins=11, xmin = 0.5, xmax = 1.0, dcmask = dcm)
    top = spp.GetTopSacs(dat, topnumber=7, dcmask=dcm)
    plt.ion()
    spp.SacVSVar_Plot(dat, topcuts=top, metadata=meta)
    plt.xlabel(r"$U.R_{AV}$")
    #plt.xlabel(r"$(R_{AV}/R_{Vessel})^{3}$")
    #plt.xlabel("Energy (MeV)")
    plt.title("Fractional sacrifice of internal N16 events by data cleaning\n"+\
            "AV ROI, November 2017 scan", fontsize=26)
    plt.show() 
    plt.savefig("output_File.png")
    time.sleep(360)
