# Create an age model using minimal tuning

Create an age model using the minimal tuning technique. This means that
the distance between 2 peaks of an extracted cycle are set to duration
of the interpreted astronomical cycle

## Usage

``` r
minimal_tuning(
  data = NULL,
  pts = 5,
  cycle = 405,
  tune_opt = "max",
  output = 0,
  genplot = FALSE,
  keep_editable = FALSE
)
```

## Arguments

- data:

  Cycle extracted or filtered in the depth domain (e.g. the output of
  [`extract_signal`](https://stratigraphy.github.io/WaverideR/reference/extract_signal.md)
  or columns 1-2 of the element `filtered` of
  [`gaussian_bandpass`](https://stratigraphy.github.io/WaverideR/reference/gaussian_bandpass.md)):
  a two-column matrix or data frame with depth in metres in column 1 (in
  increasing order) and the filtered signal in column 2. Required
  (`Default=NULL`).

- pts:

  Number of points to each side that a point is compared with when the
  maxima and minima of the filtered cycle are detected (integer, \>= 2):
  a point is a maximum (minimum) when it is above (below) the value
  `pts` points further on and the value `pts - 1` points back. Larger
  values give more robust extremes but may miss minor ones (see
  [`max_detect`](https://stratigraphy.github.io/WaverideR/reference/max_detect.md)).
  `Default=5`.

- cycle:

  Duration of the filtered/extracted cycle in kyr (numeric, \> 0, e.g.
  405 for the long eccentricity cycle); the distance between successive
  extremes divided by this duration gives the sedimentation rate.
  `Default=405`.

- tune_opt:

  Which extremes of the cycle are used as tie points of the age model
  (character): `"min"` (minima, one cycle apart), `"max"` (maxima, one
  cycle apart) or `"minmax"` (both, half a cycle apart).
  `Default="max"`.

- output:

  Content of the returned matrix (integer): `0` four columns (depth,
  proxy, sedimentation rate in cm/kyr, time in kyr); `1` two columns
  (depth, sedimentation rate); `2` two columns (depth, time).
  `Default=0`.

- genplot:

  Plot the result belonging to `output` (logical, TRUE or FALSE).
  `Default=FALSE`.

- keep_editable:

  Keep the graphical parameters (`par`) of the plot after the function
  returns, so that features can be added to it (logical, TRUE or FALSE);
  `FALSE` restores the previous settings. `Default=FALSE`.

## Value

The output depends on the output setting If output = 0 output is a
matrix of with 4 columns being (depth,proxy,sedimentation rate and time)
If genplot = TRUE 4 plots are generated; depth vs proxy, depth vs
sedimentation rate, depth vs time and time vs proxy If output = 1 output
is a matrix of with 2 columns being (depth and sedimentation rate ) If
genplot = TRUE a plot of depth vs sedimentation rate is generated If
output = 2 output is a matrix of with 2 columns being (depth and time)
If genplot = TRUE a plot of depth vs time is generated

## References

Routines for astrochronologic testing, astronomical time scale
construction, and time series analysis
\<doi:10.1016/j.earscirev.2018.11.015\>

## Author

Part of the code is based on the
[sedrate2time](https://rdrr.io/pkg/astrochron/man/sedrate2time.html)
function of the 'astrochron' R package

## Examples

``` r
# \donttest{
# Extract the 405kyr eccentricity cycle from the wavelet scalogram
# from the magnetic susceptibility record f the Sullivan core
# of Pas et al., (2018) and then create a age model using minimal tuning
# (e.g.) set the distance between peaks to 405 kyr

mag_wt <- analyse_wavelet(data = mag,
dj = 1/100,
lowerPeriod = 0.1,
upperPeriod = 254,
verbose = FALSE,
omega_nr = 10)


mag_405 <- extract_signal_stable_V2(
 wavelet = mag_wt,
 period_max = 4,
 period_min = 2,
 add_mean = FALSE,
 plot_residual = FALSE,
 keep_editable = FALSE
)

mag_405_min_tuning <- minimal_tuning(data = mag_405,
pts = 5,
cycle = 405,
tune_opt = "max",
output = 0,
genplot = FALSE,
keep_editable = FALSE)


# }
```
