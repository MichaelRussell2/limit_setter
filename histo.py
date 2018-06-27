#!/usr/bin/env/ python2.7


class histo(object):
    "Simple histogram class"

    #constructor
    def __init__(self, xmin=None, xmax=None, bins=None):
        self._xmin = xmin
        self._xmax = xmax
        self._bins = bins
    
    #load from 3 column file    
    @classmethod
    def from_dat(cls, path=None):
        import numpy as np
        columns = np.loadtxt(path,unpack=True)
        try:
            assert len(columns) == 3
        except:
            print "Error. Histogram has wrong input format"
            print "Should have three columns - x_min x_max y "
            print "Exiting"
            import sys
            sys.exit(1)
        xmin = columns[0]
        xmax = columns[1]
        bins = columns[2]
        return cls(xmin, xmax, bins)

    #load from ROOT file
    @classmethod()
    def from_root(cls,path=None,th1f=None):
        from ROOT import TFile
        import numpy as np

        fin = TFile(path)
        hist = fin.Get(th1f)
        xmin = np.asarray([ hist.GetBinLowEdge(i) for i in xrange(hist.GetSize())])
        xmax = np.asarray([ hist.GetBinLowEdge(i+1) for i in xrange(hist.GetSize())])
        bins = np.asarray([ hist.GetBinContent(i+1) for i in xrange(hist.GetSize())]) #bin0=underflow
        return cls(xmin, xmax, bins)
    
    def xmin(self,binnum=None):
         return self._xmin[binnum] if binnum is not None else self._xmin

    def xmax(self,binnum=None):
         return self._xmax[binnum] if binnum is not None else self._xmax

    def bin(self,binnum=None):
        return self._bins[binnum] if binnum is not None else self._bins

    def bins(self):
        return self._bins
    
    def nbins(self):
        return len(self._bins)

    def integral(self):
        return self._bins.sum()

    def scale(self,scalefactor):
        self._bins *= scalefactor

    #cut out first n bins @TODO: bug here?
    def cut_first(self,n):
        self._xmin = self._xmin[n:]
        self._xmax = self._xmax[n:]
        self._bins = self._bins[n:]

    #cut out last n bins
    def cut_last(self,n):
        self._xmin = self._xmin[:n]
        self._xmax = self._xmax[:n]
        self._bins = self._bins[:n]
