import sacrifice_b14itr_plot as sbp
import glob
import matplotlib.pyplot as plt
import time

if __name__=='__main__':
    udotr_av = '(posx*dirx + posy*diry + (posz-108.0)*dirz)/sqrt(posx**2 + ' + \
            "(posz-108.0)**2)"
    udotr_av_cut = udotr_av + "> (1.0 - 12 *(((sqrt((posx**2 + posy**2 + (posz-108.0)**2))/"\
            "6000.0)**3)-0.69)**2)"
    r_av3 = "(sqrt((posx**2 + posy**2 + (posz-108.0)**2))/6000.0)**3"
    n16 = glob.glob("/home/onetrueteal/share/N16/Nov2017_N16/N16/*.ntuple.root") 
    pre = "energy>3.0&&energy<5.6&&0.69<%s&&0.90>%s&&%s"%(r_av3, r_av3, udotr_av_cut)
    
    dat, meta = sbp.PrepData_ClassSac(rootfiles=n16, var = 'energy', 
            precuts = pre, nbins=14, xmin = 3.0, xmax = 5.6)
    sbp.Plot_ClassSac(dat, meta=meta)
    #plt.xlabel(r"$U.R_{AV}$")
    #plt.xlabel(r"$(R_{AV}/R_{Vessel})^{3}$")
    plt.xlabel("Energy (MeV)")
    plt.title("Fractional sacrifice of internal N16 events by classifiers\n"+\
            "AV ROI, November 2017 scan", fontsize=28)
    plt.show() 
    plt.savefig("output_File_b14itr.png")
    time.sleep(360)
