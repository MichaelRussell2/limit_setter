#!/usr/bin/env python

import math, sys
import numpy as np

if len(sys.argv) != 3:
    print 'provide signal and background histograms'
    sys.exit()
sigfile=sys.argv[1]
bkgfile=sys.argv[2]

# be sure to get lumi and xsec units right!!! typically pb

#calculate gaussian significance
def gaussianSig(sigmaS, sigmaB, lumi, systS, systB):
    if sigmaB == 0: return 0
    S=sigmaS*lumi
    B=sigmaB*lumi
    return S/math.sqrt(B+(systB*B)**2+(systS*S)**2)

#convert one-sided Gaussian sig to p-value
def ZtoPval(x):
    return (1-math.erf(x/math.sqrt(2)))/2

#input files have format xlo, xhi, y (no errors)
_, _, xS = np.loadtxt(sigfile,unpack=True)
_, _, xB = np.loadtxt(bkgfile,unpack=True)

#10% systematics
alpha=0.1
beta=0.1
lumi=300e3 #pb^-1

#per-bin significances
for ind, val in enumerate(xS):
    Z = gaussianSig(xS[ind],xB[ind],lumi,alpha,beta)
    print ind,'\t', Z
#    print cut[ind],'\t', ZtoPval(Z)

#total significance
print gaussianSig(xS.sum(),xB.sum(),lumi,alpha,beta)
