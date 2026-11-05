# Convert the tracked curve to a depth time space

Converts the tracked curve to a depth time space.

## Usage

``` r
curve2time(
  tracked_cycle_curve = NULL,
  tracked_cycle_period = NULL,
  genplot = FALSE,
  keep_editable = FALSE
)
```

## Arguments

- tracked_cycle_curve:

  Tracked period curve of one astronomical cycle (matrix or data frame,
  two columns): depth in metres in column 1 and the period of the
  tracked cycle in metres in column 2, e.g. the output of
  [`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
  (preferably completed with
  [`completed_series`](https://stratigraphy.github.io/WaverideR/reference/completed_series.md)).
  Rows with `NA` are dropped and the curve is interpolated to its median
  depth spacing. Required (`Default=NULL` gives an error).

- tracked_cycle_period:

  Duration of the tracked cycle in kyr (single numeric value, \> 0, e.g.
  405 for the long eccentricity cycle), used to convert the tracked
  period to a sedimentation rate that is integrated to time. Required
  (`Default=NULL` gives an error).

- genplot:

  Plot depth against time (logical, TRUE or FALSE). `Default=FALSE`.

- keep_editable:

  Keep the graphical parameters of the plot so that extra elements can
  be added afterwards (logical, TRUE or FALSE); `FALSE` restores the
  previous `par` settings. `Default=FALSE`.

## Value

The output is a matrix with 2 columns. The first column is depth. The
second column sedimentation rate in cm/kyr. If `genplot=TRUE` then a
depth vs time plot will be plotted.

## References

Routines for astrochronologic testing, astronomical time scale
construction, and time series analysis
\<doi:10.1016/j.earscirev.2018.11.015\>

## Author

Based on the
[sedrate2time](https://rdrr.io/pkg/astrochron/man/sedrate2time.html)
function of the 'astrochron' R package

## Examples

``` r
# \donttest{
#Convert a tracked curve to a depth time space. The examples uses the
#magnetic susceptibility data set of Pas et al., (2018).

#'# perform the CWT
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

#Instead of tracking, the tracked solution data set mag_track_solution is used
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

#convert period in meters to sedrate depth vs time
mag_track_time<- curve2time(tracked_cycle_curve=mag_track_complete,
tracked_cycle_period=405,
genplot=FALSE,
keep_editable=FALSE)

# }
```
