#!/usr/bin/env/ python2.7


class histo(object):
    "Simple histogram class"
    
#    def __init__(self, path=None):
        
        
    def from_data(self, path=None):
        import numpy as np
        self.path = path
        self.columns = np.loadtxt(path,unpack=True)
        try:
            assert len(self.columns) == 3
        except:
            print "Error. Histogram has wrong input format"
            print "Should have three columns - x_min x_max y "
            print "Exiting"
            import sys
            sys.exit(1)
        self._xmin = self.columns[0]
        self._xmax = self.columns[1]
        self._bins = self.columns[2]
     
    def from_root(self,path=None,th1f=None):
        from ROOT import TFile
        fin = TFile(path)
        hist = fin.Get(th1f)
        self._xmin = [ hist.GetBinLowEdge(i) for i in xrange(hist.GetSize() ) ]
        self._xmax = [ hist.GetBinLowEdge(i+1) for i in xrange(hist.GetSize() ) ]
        self._bins = [ hist.GetBinContent(i+1) for i in xrange(hist.GetSize() ) ]                
            
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
