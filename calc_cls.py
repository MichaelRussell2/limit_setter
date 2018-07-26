#!/usr/bin/env python2.7

#@TODO: Expected 1 and 2sigma up/down ->brazil bands
#@TODO: correlated smearing across neighbouring bins

import math
import numpy as np
from histo import histo

n_mc = int(1e5)
syst = 0.02 #systematic uncertainty

def main():

    #load histograms [ units = pb/ bin ]
    xsec_sig=histo.from_dat("ptj1_sig.dat")
    xsec_bkg=histo.from_dat("ptj1_bkg.dat")

    #convert to events/bin
    lumi=300e3 #pb^-1
    xsec_sig.scale(lumi)
    xsec_bkg.scale(lumi)

    cls_binned(xsec_sig, xsec_bkg)

    #compare to cut and count
    s = xsec_sig.integral()*lumi
    b = xsec_bkg.integral()*lumi
    print
    print "Cut and count p-value", p_value(poisson_sig(s,b,0,syst))

#get CLsb and CLb for a single bin
def cls(sig, bkg, obs):

    #distributions of LLR, under each hypothesis 
    llr_H0, llr_H1 = np.zeros(n_mc), np.zeros(n_mc)

    # log-likelihood from data
    obs = bkg
    llr_bin = -2*(-sig + obs*math.log(1+sig/bkg))
    print "LLR in data bin", llr_bin
        
    n_mc_H0, n_mc_H1, n_mc_H0_smear , n_mc_H1_smear= 0, 0, 0, 0

    #run pseudoexperiments, count those with llr < llr_mc
    for i in xrange(n_mc):

        #no systematics case
        obs_mc_H0 = np.random.poisson(bkg)
        obs_mc_H1 = np.random.poisson(sig+bkg)
        llr_mc_H0 = -2*(-sig + obs_mc_H0*math.log(1+sig/bkg) )
        llr_mc_H1 = -2*(-sig + obs_mc_H1*math.log(1+sig/bkg) )
        
        #if data less extreme than PE, add this PE to CLsb, CLb
        if llr_bin <= llr_mc_H0:
            n_mc_H0 += 1
        if llr_bin <= llr_mc_H1:
            n_mc_H1 += 1
    
        llr_H0[i] = llr_mc_H0
        llr_H1[i] = llr_mc_H1

        # Gaussian systematic on background
        bkg_smear = bkg+np.random.normal(0,syst*bkg)

        obs_mc_H0_smear = np.random.poisson(bkg_smear)
        obs_mc_H1_smear = np.random.poisson(sig+bkg_smear)
        llr_mc_H0_smear = -2*(-sig + obs_mc_H0_smear*math.log(1+sig/bkg))
        llr_mc_H1_smear = -2*(-sig + obs_mc_H1_smear*math.log(1+sig/bkg))

        if llr_bin <= llr_mc_H0_smear:
            n_mc_H0_smear += 1
        if llr_bin <= llr_mc_H1_smear:
            n_mc_H1_smear += 1

    clb = float(n_mc_H0)/float(n_mc)
    clsb = float(n_mc_H1)/float(n_mc)
    #plot_llr(llr_H0, llr_H1, llr_bin)
        
    clb_smear = float(n_mc_H0_smear)/float(n_mc)
    clsb_smear = float(n_mc_H1_smear)/float(n_mc)

    return clb, clsb, clb_smear, clsb_smear
            
#cls of a full histogram
def cls_binned(sighist,bkghist):

    clb_all, clsb_all, clb_all_smear,  clsb_all_smear = 1., 1., 1., 1.
    
    #get likelihood for each bin, then multiply clsb together
    for ibin, (sig, bkg) in enumerate(zip(sighist.bins() , bkghist.bins() )):

        clb, clsb, clb_smear, clsb_smear = cls(sig, bkg, bkg)

        if ibin > 0: break
        print "Bin ", ibin, "CLb ", clb, "CLsb ", clsb,  "CLs ", clsb/clb, "CLb_syst", clb_smear, "CLsb_syst ", clsb_smear
        print

        clb_all *= clb
        clsb_all *= clsb

        clb_all_smear *= clb_smear        
        clsb_all_smear *= clsb_smear

    print "--------------------------"
    print "SUMMARY"
    print "Without systematics:"
    print "Expected CLsb, CLb, CLs: ", clsb_all, clb_all, clsb_all/clb_all
    print
    print "With Gaussian systematics:"
    print "Expected CLsb, CLb, CLs: ", clsb_all_smear, clb_all_smear, clsb_all_smear/clb_all_smear
    print "--------------------------"
        

#calculate cut-and-count significance
def poisson_sig(s,b,alpha,beta):
    return s/math.sqrt(b+(alpha*s)**2+(beta*b)**2)

#convert one-sided Gaussian sig to p-value
def p_value(x):
    return (1-math.erf(x/math.sqrt(2)))/2

#plot log-likelihoods for each hypothesis
def plot_llr(H0, H1, LLR):

    import matplotlib.pyplot as plt
    plt.clf()

    # set histogram bins from data
    nbins = 100.
    data = np.concatenate((data_H0,data_H1))
    data_min, data_max = data.min(),  data.max()
    binrange  = np.linspace(data_min,data_max,nbins+1)

    #before plotting in mpl, save histograms to ascii    
    counts_H0, bins, _ = plt.hist(data_H0, bins=binrange, alpha=0.4, normed=True,label='b')
    counts_H1, bins, _ = plt.hist(data_H1, bins=binrange, alpha=0.4, normed=True, label='s+b')

    #import matplotlib.mlab as mlab    
    # H0_med = np.median(data_H0)
    # H1_med = np.median(data_H1)
    # gaussianH0 = np.linspace(H0_med-5*math.sqrt(abs(H0_med)),H0_med+5*math.sqrt(abs(H0_med)),100)
    # gaussianH1 = np.linspace(H1_med-5*math.sqrt(abs(H1_med)),H1_med+5*math.sqrt(abs(H1_med)),100)    
    # plt.plot(gaussianH0,mlab.normpdf(gaussianH0,H0_med,math.sqrt(abs(H0_med))))
    # plt.plot(gaussianH1,mlab.normpdf(gaussianH1,H1_med,math.sqrt(abs(H1_med))))
    
    xlo, xhi = bins[:-1], bins[1:]
    np.savetxt('LLR.dat', np.c_[xlo,xhi,counts_H0,counts_H1])
                
    plt.xlabel('-2*(log-likelihood ratio)')
    plt.ylabel('Probability')
    plt.axvline(x=llr_all,color='red')
#    plt.annotate('observed',xy=[-15,0.12],color='red')
    plt.legend()
#    plt.show()
    plt.savefig('LLR.png')

if __name__ == "__main__": main()
