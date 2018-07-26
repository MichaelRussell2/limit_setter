# limit-setter

A frequentist hypothesis testing tool, based on the [CLs method](https://indico.cern.ch/event/398949/attachments/799330/1095613/The_CLs_Technique.pdf). Largely a python adaptation of Thomas Junk's [MClimits](https://www-cdf.fnal.gov/~trj/mclimit/production/mclimit.html) code.
Requires numpy.

## Input

   * 1D signal, background and data histograms, either as plain data files in three-column format: x\_min x\_max y\_central or as TH1 objects from a ROOT file
   * Luminosity, scale factors and uncertainties can be tuned in `calc_cls.py`

## Algorithm

   * For each histogram bin, evaluate log-likelihood ratio for signal[bin], background[bin] and data (for setting limits on a test hypothesis, set data = background)
   *  Run number of Monte-Carlo pseudoexperiments, where pseudo-data is drawn from a Poisson distribution with mean = bkg[bin] for null hyp and mean = sig[bin]+bkg[bin]
   * CLb for that bin = fraction of pseudo-experiments with log-likelihood ratio less than data (under null hypothesis)
   * CLsb for that bin = ditto with test hypothesis
   * CLs = CLsb/CLb.

For a histogram, the combined CLs is the product of the CLs for each bin. This assumes un-correlated bins. Correlations can be added in as discussed in the next section.

## Systematic uncertainties

These are assumed to be Gaussian, so are included by convoluting the Poisson PDFs with a Gaussian of std-dev = Delta_bkg  and mean = bkg. This still assumes neighbouring bins to be uncorrelated.


## Output

CLs exclusion limit for signal and signal plus background hypotheses, given data.
