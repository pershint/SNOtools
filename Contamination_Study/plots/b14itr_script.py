import sacrifice_b14itr_DataMCComp as sbp
import glob
import matplotlib.pyplot as plt
import time

if __name__=='__main__':
    #udotr_av = '(posx*dirx + posy*diry + (posz-108.0)*dirz)/sqrt(posx**2 + ' + \
    #        "(posz-108.0)**2)"
    #udotr_av_cut = "(udotr > (1.0 - 12 *((posr3-0.69)**2)))"
    #r_av3 = "(sqrt((posx**2 + posy**2 + (posz-108.0)**2))/6000.0)**3"
    #n16_data = glob.glob("/home/onetrueteal/share/N16/Data/Nov2017_N16/subtuple/*") 
    #n16_MC = glob.glob("/home/onetrueteal/share/N16/MC/Nov2017_IntScan/subtuple/*") 
    n16_data = glob.glob("/home/onetrueteal/share/N16/Data/Mar2018_N16EXT/subtuple/*") 
    n16_MC = glob.glob("/home/onetrueteal/share/N16/MC/Mar2017_ExtScan/subtuple/*") 
    pre = "energy>3.0&&energy<5.6&&posr>6194&&posr<7161&&udotr>0.5"
    #pre=None
    dat, meta = sbp.PrepData_ClassSac(datafiles=n16_data, MCfiles=n16_MC, 
            var = 'energy', 
            precuts = pre, nbins=16, xmin = 3.0, xmax = 5.6)
    sac_title="Fractional sacrifice due to classifiers \n"+\
            " Ext. Non-PMT ROI, Data-MC comparison, Mar. 2018 N16 scan"
    acc_title="Fractional acceptance of events by classifiers \n"+\
            " Ext. Non-PMT ROI, Data-MC comparison, Mar. 2018 N16 scan"
    ratio_title="Data/MC Ratio of acceptance for classifiers \n "+\
            " Ext. Non-PMT ROI, Mar. 2018 N16 scan"
    #xlabel=r"$U.R$" 
    #xlabel="Nhit"
    xlabel="Energy (MeV)" 
    #xlabel=r"$(R/R_{AV})^{3}$"
    sbp.Plot_SacComparison(dat, meta=meta,title=sac_title, xlabel=xlabel)
    sbp.Plot_AccComparison(dat, meta=meta,title=acc_title, xlabel=xlabel)
    sbp.Plot_Ratio(dat, meta=meta,title=ratio_title, xlabel=xlabel)
    #plt.savefig("output_File_b14itr.png")
    #time.sleep(360)
