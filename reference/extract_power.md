# Extract power from a wavelet spectra

Extracts the spectral power from a wavelet spectra in the depth domain
using a traced period and boundaries surround the traced period. The
extraction of spectral is useful for cyclostratigraphic studies because
the spectral power of an astronomical cycle is modulated by higher order
astronomical cycles. The spectral power record from an astronomical
cycle can thus be used as a proxy for amplitude modulating cycles The
traced period result from the
[`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
function with boundaries is used to extract spectral power in the depth
domain from a wavelet spectra.

## Usage

``` r
extract_power(
  tracked_period_curve = NULL,
  scalogram = NULL,
  period_up = 1.2,
  period_down = 0.8,
  tracked_cycle_period = NULL,
  extract_cycle_power = NULL,
  ...
)
```

## Arguments

- tracked_period_curve:

  Tracked period of the cycle along the record (matrix or data frame:
  depth/time in column 1, period in column 2, in the depth/time unit of
  the scalogram), e.g. the output of
  [`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
  or
  [`track_period`](https://stratigraphy.github.io/WaverideR/reference/track_period.md)
  completed with
  [`completed_series`](https://stratigraphy.github.io/WaverideR/reference/completed_series.md),
  optionally smoothed with
  [`loess_auto`](https://stratigraphy.github.io/WaverideR/reference/loess_auto.md).
  Rows with `NA` are dropped and the curve is linearly interpolated (and
  extrapolated) to the depth/time axis of the scalogram. Required
  (`Default=NULL`).

- scalogram:

  Wavelet or superlet object (list of class `analyze.wavelet` or
  `analyze.superlet`) created with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  or
  [`analyse_superlet`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md),
  whose power is summed. Required (`Default=NULL`).

- period_up:

  Upper boundary of the band whose power is summed, as a factor of the
  (rescaled) tracked period (numeric, \> `period_down`). `Default=1.2`.

- period_down:

  Lower boundary of the band whose power is summed, as a factor of the
  (rescaled) tracked period (numeric, \> 0). `Default=0.8`.

- tracked_cycle_period:

  Period of the tracked cycle, e.g. 405 for the 405-kyr eccentricity
  cycle (numeric, \> 0), in the same unit as `extract_cycle_power`. The
  tracked curve is multiplied by
  `extract_cycle_power / tracked_cycle_period`. Required
  (`Default=NULL`).

- extract_cycle_power:

  Period of the cycle whose power is extracted (numeric, \> 0), in the
  same unit as `tracked_cycle_period` (both in time or both in depth).
  Use the same value as `tracked_cycle_period` to extract the tracked
  cycle itself. Required (`Default=NULL`).

- ...:

  Old argument names of earlier versions: `wavelet` (used as
  `scalogram`) and `completed_series` (used as `tracked_period_curve`).

## Value

Returns a matrix with 3 columns. The first column is depth/time. The
second column is extracted power. The third column is extracted
power/total power.

## References

Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
Wavelet Analysis. R package version 1.1.
<https://CRAN.R-project.org/package=WaveletComp>

Gouhier TC, Grinsted A, Simko V (2021). R package biwavelet: Conduct
Univariate and Bivariate Wavelet Analyses. (Version 0.20.21),
<https://github.com/tgouhier/biwavelet>

Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet
Analysis. Bulletin of the American Meteorological Society 79:61-78.
<https://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf>

Routines for astrochronologic testing, astronomical time scale
construction, and time series analysis
\<doi:10.1016/j.earscirev.2018.11.015\>

## Author

Code based on the reconstruct function of the 'WaveletComp' R package
which is based on the wavelet 'MATLAB' code written by Christopher
Torrence and Gibert P. Compo. The assignment of the standard deviation
of the uncertainty of the wavelet is based on the work of Gabor (1946)
and Russell et al., (2016) The functionality of this function is is
inspired by the
[integratePower](https://rdrr.io/pkg/astrochron/man/integratePower.html)
function of the 'astrochron' R package.

## Examples

``` r
#Extract the power of the 405 kyr eccentricity cycle from the the magnetic
# susceptibility data set of De pas et al., (2018)
#Perform the CWT on the magnetic susceptibility data set of Pas et al., (2018)

mag_wt <- analyse_wavelet(data = mag,
dj = 1/100,
lowerPeriod = 0.1,
upperPeriod = 254,
verbose = FALSE,
omega_nr = 10)

#Track the 405 kyr eccentricity cycle in a wavelet spectra

#mag_track <- track_period_wavelet(astro_cycle = 405,
#                                   wavelet=mag_wt,
#                                   n.levels = 100,
#                                   periodlab = "Period (metres)",
#                                   x_lab = "depth (metres)")

#Instead of tracking, the tracked solution data set mag_track_solution
#is used
mag_track <- mag_track_solution

mag_track_complete <- completed_series(
  scalogram = mag_wt,
  tracked_curve = mag_track,
  period_up = 1.2,
  period_down = 0.8,
  extrapolate = TRUE,
  genplot = FALSE
)

#Smooth the completed tracking of the 405 kyr eccentricity cycle in the wavelet spectra

mag_track_complete <- loess_auto(time_series = mag_track_complete,
genplot = FALSE, print_span = FALSE)

#extract the spectral power of the 405 kyr eccentricity cycle
mag_power <- extract_power(
tracked_period_curve = mag_track_complete,
scalogram = mag_wt,
period_up = 1.2,
period_down = 0.8,
tracked_cycle_period = 405,
extract_cycle_power = 405
)
```
