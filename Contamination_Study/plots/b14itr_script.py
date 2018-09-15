import sacrifice_b14itr_DataMCComp as sbp
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas
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
    pre = "energy>3.0&&energy<5.6&&posr>5569&&posr<5793&&udotr>0.5"
    #pre=None
    dat, meta = sbp.PrepData_ClassSac(datafiles=n16_data, MCfiles=n16_MC, 
            var = 'energy', 
            precuts = pre, nbins=13, xmin = 3.0, xmax = 5.6)
    sac_title="Fractional sacrifice due to classifiers \n"+\
            " Ext. Bkg. Non-PMT ROI (in AV), Data-MC comparison, Mar. 2018 N16 scan"
    acc_title="Fractional acceptance of events by classifiers \n"+\
            " Ext. Bkg. Non-PMT ROI (in AV), Data-MC comparison, Mar. 2018 N16 scan"
    ratio_title="Data/MC Ratio of acceptance for classifiers \n "+\
            " Ext. Bkg. Non-PMT ROI (in AV), Mar. 2018 N16 scan"
    #xlabel=r"$U.R$" 
    #xlabel="Nhit"
    xlabel="Energy (MeV)" 
    #xlabel=r"$(R/R_{AV})^{3}$"
    print(dat)
    ratio = pandas.Series(dat["Data total"].fractional_acceptance /
            dat["MC total"].fractional_acceptance) 
    variable = pandas.Series(dat["Data total"].vardat)
    ratio_unc = np.sqrt(dat["MC total"].fs_uncertainty**2 + \
            dat["Data total"].fs_uncertainty**2)
    ratio_unc = pandas.Series(ratio_unc)
    thegoods = {"ratio": ratio.round(4), "ratio_unc": ratio_unc.round(4), xlabel: variable.round(1)}
    thegoods_pd = pandas.DataFrame(thegoods)
    thegoods_pd.to_csv("RatioDataOutput.csv")
    sbp.Plot_SacComparison(dat, meta=meta,title=sac_title, xlabel=xlabel)
    sbp.Plot_AccComparison(dat, meta=meta,title=acc_title, xlabel=xlabel)
    sbp.Plot_Ratio(dat, meta=meta,title=ratio_title, xlabel=xlabel)
    #plt.savefig("output_File_b14itr.png")
    #time.sleep(360)
