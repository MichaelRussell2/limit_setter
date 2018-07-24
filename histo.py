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
    @classmethod
    def from_root(cls,path=None,th1f=None):
        from ROOT import TFile
        import numpy as np

        fin = TFile(path)
        hist = fin.Get(th1f)
        xmin = np.asarray([ hist.GetBinLowEdge(i) for i in xrange(hist.GetSize())])
        xmax = np.asarray([ hist.GetBinLowEdge(i+1) for i in xrange(hist.GetSize())])
        bins = np.asarray([ hist.GetBinContent(i+1) for i in xrange(hist.GetSize())]) #bin0=underflow
        return cls(xmin, xmax, bins)

    #bin raw data from ROOT tree
    @classmethod
    def from_tree(cls,path=None,tree="tree",branch=None,binwidth=None,xmin=None,xmax=None,weight=None,normed=False,dtype=None):
        from ROOT import TFile, TTree
        import numpy as np
        import array
        
        fin = TFile(path)
        tin = fin.Get(tree)

        #store branch value as pointer
        var = array.array(dtype,[0])
        tin.SetBranchAddress(branch,var)

        entries = tin.GetEntriesFast()
        data = np.zeros(entries)

        print "Looping over entries in tree..."
        for i in range(tin.GetEntries()):
	    tin.GetEntry(i)
            data[i] = var[0]
        print "Done."
            
        #optionally, fill histogram with weights instead of event counts
        weights = np.full(len(data),weight) if weight is not None else np.full(len(data),1)

        #if not specified, set histogram range and binning from data (default 10 bins)
        xmax = xmax if xmax is not None else data.max()  
        xmin = xmin if xmin is not None else data.min()  
        nbins =  abs(xmax-xmin)/float(binwidth) if binwidth is not None else 10
        nbins = int(nbins)

        if normed:
            yvals, binedges = np.histogram(data,nbins, range=(xmin,xmax),density=True)
        else:
            yvals, binedges = np.histogram(data,nbins, range=(xmin,xmax),density=False)

        #lower and upper edges
        binlo = binedges[:-1]
        binhi = binedges[1:]

        return cls(binlo, binhi,yvals)
            
    def xmin(self,binnum=None):
         return self._xmin[binnum] if binnum is not None else self._xmin

    def xmax(self,binnum=None):
         return self._xmax[binnum] if binnum is not None else self._xmax

    def bin(self,binnum=None):
        return self._bins[binnum] if binnum is not None else self._bins

    def bins(self,binnum=None):
        import numpy as np
        return np.asarray(self._bins[binnum]) if binnum is not None else self._bins
    
    def nbins(self):
        return len(self._bins)

    def integral(self,binnum=None):
        return self._bins[binnum:].sum() if binnum is not None else self._bins.sum()

    def binwidth(self):
        return abs(self.xmax(-1) - self.xmin(0))/float(self.nbins() )
    
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

    #store as .dat file
    def to_file(self,path=None):

        path = path if path is not None else "hist.dat"
        import numpy as np
        np.savetxt(path, np.c_[self._xmin, self._xmax, self._bins], fmt=["%.3f","%.3f","%.10f"], delimiter="  ")
        
