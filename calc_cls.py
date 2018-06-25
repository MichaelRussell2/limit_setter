#!/usr/bin/env python2.7

import math
import numpy as np
from histo import histo
import sys

#load histograms [ units = pb/ bin ]
xsec_sig=histo()
xsec_bkg=histo()
xsec_sig.from_data("ETmiss_sig.dat")
xsec_bkg.from_data("ETmiss_bkg.dat")

#convert to events/bin
lumi=300e3 #pb^-1
xsec_sig.scale(lumi)
xsec_bkg.scale(lumi)

#calculate cut-and-count significance
def poisson_sig(s,b,alpha,beta):
    return s/math.sqrt(b+(alpha*s)**2+(beta*b)**2)

#convert one-sided Gaussian sig to p-value
def p_value(x):
    return (1-math.erf(x/math.sqrt(2)))/2

def plot(data_H0, data_H1, llr_all):

    H0_med = np.median(data_H0)
    H1_med = np.median(data_H1)
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
    
    plt.clf()

    # set histogram bins from data
    nbins = 100.
    data = np.concatenate((data_H0,data_H1))
    data_min, data_max = data.min(),  data.max()
    width = (data_max - data_min) / nbins
    binrange  = np.arange(data_min,data_max+1,width)

    #before plotting in mpl, save histograms to ascii    
    counts_H0, bins, dump = plt.hist(data_H0, bins=binrange, alpha=0.4, normed=True,label='b smear')
    counts_H1, bins, dump = plt.hist(data_H1, bins=binrange, alpha=0.4, normed=True, label='s+b smear')

    # gaussianH0 = np.linspace(H0_med-5*math.sqrt(abs(H0_med)),H0_med+5*math.sqrt(abs(H0_med)),100)
    # gaussianH1 = np.linspace(H1_med-5*math.sqrt(abs(H1_med)),H1_med+5*math.sqrt(abs(H1_med)),100)    
    # plt.plot(gaussianH0,mlab.normpdf(gaussianH0,H0_med,math.sqrt(abs(H0_med))))
    # plt.plot(gaussianH1,mlab.normpdf(gaussianH1,H1_med,math.sqrt(abs(H1_med))))
    
    xlo, xhi = bins[:-1], bins[1:]
    np.savetxt('test.dat', np.c_[xlo,xhi,counts_H0,counts_H1])
                
    plt.xlabel('-2*(log-likelihood ratio)')
    plt.ylabel('Probability')
    plt.axvline(x=llr_all,color='red')
#    plt.annotate('observed',xy=[-15,0.12],color='red')
    plt.legend()
#    plt.show()
    plt.savefig('test.png')


#@TODO: Expected 1 and 2sigma up/down ->brazil bands
#@TODO: Compare Poisson s/sqrt(b) p-value to CLs on total rate
    
#cls of a full histogram
def cls_binned(sighist,bkghist):

    clb_all, clsb_all, clb_all_smear,  clsb_all_smear = 1., 1., 1., 1.

    #run pseudoexperiments, count those with llr < llr_mc
    n_mc = int(5e4)
    data_llr_H0, data_llr_H1 = np.zeros(n_mc), np.zeros(n_mc)
    
    #get likelihood for each bin, then multiply clsb together
    for ind, (sig, bkg) in enumerate(zip(sighist.bins() , bkghist.bins() )):

        # log-likelihood from data
        obs = bkg

        llr_bin = -2*(-sig + obs*math.log(1+sig/bkg))
        #print "Actual observed LLR in data", llr_bin
        
        n_mc_H0, n_mc_H1, n_mc_H0_smear , n_mc_H1_smear= 0., 0., 0., 0.
        
        syst = 0.05
        for i in xrange(n_mc):
            obs_mc_H0 = np.random.poisson(bkg)
            obs_mc_H1 = np.random.poisson(sig+bkg)

            #@TODO: correlated smearing across neighbouring bins
            bkg_smear = bkg+np.random.normal(0,syst*bkg)
                
            obs_mc_H0_smear = np.random.poisson(bkg_smear)
            obs_mc_H1_smear = np.random.poisson(sig+bkg_smear)

            llr_mc_H0 = -2*(-sig + obs_mc_H0*math.log(1+sig/bkg) )
            llr_mc_H1 = -2*(-sig + obs_mc_H1*math.log(1+sig/bkg) )

            llr_mc_H0_smear = -2*(-sig + obs_mc_H0_smear*math.log(1+sig/bkg))
            llr_mc_H1_smear = -2*(-sig + obs_mc_H1_smear*math.log(1+sig/bkg))

            if llr_bin <= llr_mc_H0:
                n_mc_H0 += 1
            if llr_bin <= llr_mc_H1:
                n_mc_H1 += 1
            if llr_bin <= llr_mc_H0_smear:
                n_mc_H0_smear += 1
            if llr_bin <= llr_mc_H1_smear:
                n_mc_H1_smear += 1

            data_llr_H0[i] = llr_mc_H0_smear
            data_llr_H1[i] = llr_mc_H1_smear

        clb = float(n_mc_H0)/float(n_mc)
        clsb = float(n_mc_H1)/float(n_mc)
        clb_smear = float(n_mc_H0_smear)/float(n_mc)
        clsb_smear = float(n_mc_H1_smear)/float(n_mc)
#        print "Bin ", ind, "CLb ", clb, "CLsb ", clsb,  "CLb_syst", clb_smear, "CLsb_syst ", clsb_smear

        clb_all *= clb
        clsb_all *= clsb 
        clb_all_smear *= clb_smear        
        clsb_all_smear *= clsb_smear

    print "Without systematics:"
    print "Expected CLsb, CLb, CLs: ", clsb_all, clb_all, clsb_all/clb_all
    print
    print "With Gaussian systematics:"
    print "Expected CLsb, CLb, CLs: ", clsb_all_smear, clb_all_smear, clsb_all_smear/clb_all_smear
        

#full shape
sigfile = np.loadtxt("ETmiss_sig.dat",unpack=True)
bkgfile = np.loadtxt("ETmiss_bkg.dat",unpack=True)

cls_binned(xsec_sig, xsec_bkg)

#cut and count
s = (sigfile[2]).sum()*lumi
b = (bkgfile[2]).sum()*lumi
#print (sigfile[0])[skip], p_value(poisson_sig(s,b,0.05,0.05))
