# Extract signal from a wavelet spectra using a traced period curve

Extract signal power from the wavelet in the depth domain using the
traced period.

## Usage

``` r
extract_signal(
  tracked_cycle_curve = NULL,
  wavelet = NULL,
  period_up = 1.2,
  period_down = 0.8,
  add_mean = TRUE,
  tracked_cycle_period = NULL,
  extract_cycle = NULL,
  tune = FALSE,
  plot_residual = FALSE
)
```

## Arguments

- tracked_cycle_curve:

  Tracked period of the cycle along the record (matrix or data frame:
  depth/time in column 1, period in column 2, in the depth/time unit of
  the wavelet), e.g. the output of
  [`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
  completed with
  [`completed_series`](https://stratigraphy.github.io/WaverideR/reference/completed_series.md),
  optionally smoothed with
  [`loess_auto`](https://stratigraphy.github.io/WaverideR/reference/loess_auto.md).
  Rows with `NA` are dropped and the curve is linearly interpolated (and
  extrapolated) to the depth/time axis of the wavelet. Required
  (`Default=NULL`).

- wavelet:

  Wavelet object (list) created with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  the signal is reconstructed from its coefficients. Required
  (`Default=NULL`).

- period_up:

  Upper boundary of the band that is reconstructed, as a factor of the
  (rescaled) tracked period (numeric, \> `period_down`). `Default=1.2`.

- period_down:

  Lower boundary of the band that is reconstructed, as a factor of the
  (rescaled) tracked period (numeric, \> 0). `Default=0.8`.

- add_mean:

  Add the mean of the original data series to the extracted cycle
  (logical, TRUE or FALSE); `FALSE` returns it centred on zero.
  `Default=TRUE`.

- tracked_cycle_period:

  Duration in time of the tracked cycle, e.g. 405 (kyr) for the long
  eccentricity cycle (numeric, \> 0), in the same unit as
  `extract_cycle`. The tracked curve is multiplied by
  `extract_cycle / tracked_cycle_period`; also used for the tuning when
  `tune=TRUE`. Required (`Default=NULL`).

- extract_cycle:

  Duration of the cycle to extract (numeric, \> 0), in the same unit as
  `tracked_cycle_period`; use the same value to extract the tracked
  cycle itself. Required (`Default=NULL`).

- tune:

  Convert the extracted cycle from the depth to the time domain with
  [`curve2tune`](https://stratigraphy.github.io/WaverideR/reference/curve2tune.md),
  using `tracked_cycle_curve` and `tracked_cycle_period` (logical, TRUE
  or FALSE). `Default=FALSE`.

- plot_residual:

  Plot the residual (extracted cycle minus data) against depth and its
  histogram (logical, TRUE or FALSE). `Default=FALSE`.

## Value

Returns a matrix with 2 columns The first column is depth/time The
second column is extracted signal

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

## Author

Code based on the reconstruct function of the 'WaveletComp' R package
which is based on the wavelet 'MATLAB' code written by Christopher
Torrence and Gibert P. Compo (1998).

## Examples

``` r
# \donttest{
#Extract the 405 kyr eccentricity cycle from the the magnetic susceptibility \cr
#record of the Sullivan core and use the Gabor uncertainty principle to define \cr
#the mathematical uncertainty of the analysis and use a factor of that standard \cr
#deviation to define boundaries.

#Perform the CWT
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

#Instead of tracking, the tracked solution data set \code{\link{mag_track_solution}} is used \cr
mag_track <- mag_track_solution


mag_track_complete <- completed_series(
  wavelet = mag_wt,
  tracked_curve = mag_track,
  period_up = 1.2,
  period_down = 0.8,
  extrapolate = TRUE,
  genplot = FALSE
)

# smooth the tracking of the 405 kyr eccentricity cycle
mag_track_complete <- loess_auto(time_series = mag_track_complete,
genplot = FALSE, print_span = FALSE)

# extract the 405 kyr eccentricity cycle from the wavelet spectrum and use the \cr
# tracked cycle curve and set factors of the extracted cycle as boundaries

mag_405_ecc  <- extract_signal(
tracked_cycle_curve = mag_track_complete,
wavelet = mag_wt,
period_up = 1.2,
period_down = 0.8,
add_mean = TRUE,
tracked_cycle_period = 405,
extract_cycle = 405,
tune = FALSE,
plot_residual = FALSE
)
# }
```
