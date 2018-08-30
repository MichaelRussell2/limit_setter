#!/usr/bin/env python2.7

#@TODO: Expected 1 and 2sigma up/down ->brazil bands
#@TODO: correlated smearing across neighbouring bins

from math import sqrt, log
import numpy as np
from histo import histo
import sys

n_mc = int(1e5)
#syst = 0.02 #systematic uncertainty
verbose = False
#lumi=300e3 #pb^-1

if len(sys.argv) != 3:
    print "Enter lumi [fb^-1], bkg_syst as arguments"
    sys.exit()

lumi = 1000*float(sys.argv[1])
syst = float(sys.argv[2])

def main():

    #load histograms [ units = pb/ bin ]
    xsec_sig=histo.from_dat("ptj1_sig.dat")
    xsec_bkg=histo.from_dat("ptj1_bkg.dat")

    #convert to events/bin
    xsec_sig.scale(lumi)
    xsec_bkg.scale(lumi)

    brinv=0.5
    xsec_sig.scale(brinv)    
    
    cls_binned(xsec_sig, xsec_bkg)

    #compare to cut and count
    s = xsec_sig.integral()*lumi
    b = xsec_bkg.integral()*lumi
#    print
#    print "Cut and count p-value", p_value(poisson_sig(s,b,0,syst))

#get CLsb and CLb for a single bin
def cls(sig, bkg, obs):

    #distributions of LLR, under each hypothesis 
    q_H0, q_H1 = np.zeros(n_mc), np.zeros(n_mc)

    # log-likelihood from data
    obs = bkg
    q_bin = -2*(-sig + obs*log(1+sig/bkg))

    q_bin_p1s = -2*(-sig + (obs+sqrt(obs))*log(1+sig/bkg))
    q_bin_m1s = -2*(-sig + (obs-sqrt(obs))*log(1+sig/bkg))
    q_bin_p2s = -2*(-sig + (obs+2*sqrt(obs))*log(1+sig/bkg))
    q_bin_m2s = -2*(-sig + (obs-2*sqrt(obs))*log(1+sig/bkg))

    if verbose:
        print "LLR in data:",
        print lumi/1e3, q_bin, q_bin_p1s, q_bin_m1s, q_bin_p2s, q_bin_m2s

    #initialize pseudo-experiment counters
    n_mc_H0, n_mc_H1 = 0, 0
    n_mc_H0_p1s, n_mc_H0_p2s, n_mc_H0_m1s, n_mc_H0_m2s = 0, 0, 0, 0
    n_mc_H1_p1s, n_mc_H1_p2s, n_mc_H1_m1s, n_mc_H1_m2s = 0, 0, 0, 0

    #run pseudoexperiments, count those with llr < llr_mc
    for i in xrange(n_mc):

        #smear background likelihood distribution with a Gaussian
        bkg_smear = bkg + np.random.normal(0,syst*bkg)

        obs_mc_H0 = np.random.poisson(bkg_smear)
        obs_mc_H1 = np.random.poisson(sig+bkg_smear)
        q_mc_H0 = -2*(-sig + obs_mc_H0*log(1+sig/bkg) )
        q_mc_H1 = -2*(-sig + obs_mc_H1*log(1+sig/bkg) )

        #if data less extreme than PE, add this PE to CLsb, CLb
        if q_bin <= q_mc_H0:
            n_mc_H0 += 1
        if q_bin <= q_mc_H1:
            n_mc_H1 += 1

        #1-sigma bands
        if q_bin_p1s <= q_mc_H0:
            n_mc_H0_p1s += 1
        if q_bin_m1s <= q_mc_H0:
            n_mc_H0_m1s += 1
        if q_bin_p1s <= q_mc_H1:
            n_mc_H1_p1s += 1
        if q_bin_m1s <= q_mc_H1:
            n_mc_H1_m1s += 1

        #2-sigma bands
        if q_bin_p2s <= q_mc_H0:
            n_mc_H0_p2s += 1
        if q_bin_m2s <= q_mc_H0:
            n_mc_H0_m2s += 1
        if q_bin_p2s <= q_mc_H1:
            n_mc_H1_p2s += 1
        if q_bin_m2s <= q_mc_H1:
            n_mc_H1_m2s += 1

        #for plotting log-likelihood distributions
        q_H0[i] = q_mc_H0
        q_H1[i] = q_mc_H1

    #central values
    clb = float(n_mc_H0)/float(n_mc)
    clsb = float(n_mc_H1)/float(n_mc)

    #1-sigma bands
    clb_p1s = float(n_mc_H0_p1s)/float(n_mc)
    clsb_p1s = float(n_mc_H1_p1s)/float(n_mc)
    
    clb_m1s = float(n_mc_H0_m1s)/float(n_mc)
    clsb_m1s = float(n_mc_H1_m1s)/float(n_mc)

    #2-sigma bands
    clb_p2s = float(n_mc_H0_p2s)/float(n_mc)
    clsb_p2s = float(n_mc_H1_p2s)/float(n_mc)
    
    clb_m2s = float(n_mc_H0_m2s)/float(n_mc)
    clsb_m2s = float(n_mc_H1_m2s)/float(n_mc)

    plot_llr_sigmas(q_H0, q_H1, q_bin, q_bin_p1s, q_bin_m1s)
    return clsb, clb, clsb_p1s, clb_p1s, clsb_m1s, clb_m1s, clsb_p2s, clb_p2s, clsb_m2s, clb_m2s
            
#cls of a full histogram
def cls_binned(sighist,bkghist):

    clb_all, clsb_all, cls_all = 1., 1., 1.
    cls_all_p1s, cls_all_p2s, cls_all_m1s, cls_all_m2s = 1., 1., 1., 1.,
    
    #get likelihood for each bin, then multiply clsb together
    for ibin, (sig, bkg) in enumerate(zip(sighist.bins() , bkghist.bins() )):

        if ibin > 10: break
        clsb_cent, clb_cent, clsb_p1s, clb_p1s, clsb_m1s, clb_m1s, clsb_p2s, clb_p2s, clsb_m2s, clb_m2s = cls(sig, bkg, bkg)

        if verbose:
            print "Bin", ibin
            print "Expected CLb ", clb_cent, "Expected CLsb ", clsb_cent,  "Expected CLs ", clsb_cent/clb_cent
            print
            
        cls_all_p1s *= clsb_p1s/clb_p1s
        cls_all_m1s *= clsb_m1s/clb_m1s
        cls_all_p2s *= clsb_p2s/clb_p2s
        cls_all_m2s *= clsb_m2s/clb_m2s
        cls_all *= clsb_cent/clb_cent
        
    print "--------------------------"
    print "SUMMARY"
    print "With Gaussian systematics of ", 100*syst, "%:", " at a luminosity of ", int(lumi/1e3), "fb^-1"
#    print lumi/1e3, cls_all, cls_all_p1s, cls_all_m1s, cls_all_p2s, cls_all_m2s,
    print "CLs cent +1s +2s -1s -2s: ", cls_all, cls_all_p1s, cls_all_p2s, cls_all_m1s, cls_all_m2s
    print "--------------------------"
        

#calculate cut-and-count significance
def poisson_sig(s,b,alpha,beta):
    return s/math.sqrt(b+(alpha*s)**2+(beta*b)**2)

#convert one-sided Gaussian sig to p-value
def p_value(x):
    return (1-math.erf(x/math.sqrt(2)))/2

#plot log-likelihoods for each hypothesis

def plot_llr_sigmas(H0, H1, llr, llr_p1s, llr_m1s):
    print "Plotting log-likelihood ratio to file LLR.png"
    print "Raw data will be saved in LLR.dat"

    import matplotlib.pyplot as plt

    plt.clf()

    # set histogram bins from data
    nbins = 100.
    data = np.concatenate((H0,H1))
    data_min, data_max = data.min(),  data.max()
    binrange  = np.linspace(data_min,data_max,nbins+1)

    #before plotting in mpl, save histograms to ascii    
    counts_H0, bins, _ = plt.hist(H0, bins=binrange, alpha=0.4, normed=True,label='b')
    counts_H1, bins, _ = plt.hist(H1, bins=binrange, alpha=0.4, normed=True, label='s+b')

    xlo, xhi = bins[:-1], bins[1:]
    np.savetxt('LLR.dat', np.c_[xlo,xhi,counts_H0,counts_H1],header="xmin, xmax, H0, H1",fmt="%.10f", delimiter="  ")
    
    plt.xlabel('-2*(log-likelihood ratio)')
    plt.ylabel('Probability')
    plt.axvline(x=llr,color='red')
    plt.axvline(x=llr_p1s,color='black')
    plt.axvline(x=llr_m1s,color='black')
    plt.legend()
    plt.savefig('LLR.png')
                
if __name__ == "__main__": main()
