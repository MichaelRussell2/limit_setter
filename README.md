# limit_setter

Adapted from Thomas Junk's [MClimits](https://www-cdf.fnal.gov/~trj/mclimit/production/mclimit.html) code.

Requires [ROOT](https://root.cern.ch/downloading-root) (unfortunately).

## Input

   * 1D signal, background and data histograms, in three-column format: x_min x_max y_central
   * Luminosity, scale factors and uncertainties can be tuned in `hypotheses.cc`

## Output

CLs exclusion limit for signal and signal plus background hypotheses, given data.

## Issues

Can only add flat normalisation uncertainty: shape uncertainties cause errors in minimisation.
