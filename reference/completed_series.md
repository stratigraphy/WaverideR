# Complete the tracking of the period a cycle in a scalogram

Use the traced series and the existing wavelet or superlet scalogram to
complete the tracking of a cycle in scalogram. The selected points using
the
[`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
or
[`track_period`](https://stratigraphy.github.io/WaverideR/reference/track_period.md)
function form a incomplete line unless every point is tracked. However
clicking every individual point along a wavelet ridge is time intensive
and error prone. To avoid errors and save time the `completed_series`
function can be used to complete the tracing of a cycle in a wavelet
spectra. The `completed_series` function interpolates the data points
selected using the
[`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
or
[`track_period`](https://stratigraphy.github.io/WaverideR/reference/track_period.md).
A a search a algorithm then looks up and replaces the interpolated curve
values with the values of the nearest spectral peak in the wavelet
spectra.

## Usage

``` r
completed_series(
  scalogram = NULL,
  tracked_curve = NULL,
  period_up = 1.2,
  period_down = 0.8,
  extrapolate = TRUE,
  genplot = FALSE,
  keep_editable = FALSE,
  ...
)
```

## Arguments

- scalogram:

  Scalogram in which the cycle was tracked: the output of
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  (class `analyze.wavelet`) or of
  [`analyse_superlet`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md)
  (class `analyze.superlet`). Its depth axis sets the positions of the
  completed curve. Required (`Default=NULL` gives an error).

- tracked_curve:

  Incomplete tracked period curve (matrix or data frame, two columns:
  depth in column 1 and the tracked period in column 2, in the depth
  units of the scalogram), e.g. the output of
  [`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
  or
  [`track_period`](https://stratigraphy.github.io/WaverideR/reference/track_period.md);
  rows with `NA` are dropped. It is linearly interpolated onto the depth
  axis of the scalogram before the search for spectral peaks. Required
  (`Default=NULL` gives an error).

- period_up:

  Factor (numeric, \> 1) that sets the upper limit of the search window:
  at every depth the nearest spectral peak (local power maximum along
  the period axis) is used only if its period is below the interpolated
  period times `period_up`; otherwise no value is assigned there (see
  `extrapolate`). `Default=1.2`.

- period_down:

  Factor (numeric, between 0 and 1) that sets the lower limit of the
  search window: the nearest spectral peak is used only if its period is
  above the interpolated period times `period_down`. `Default=0.8`.

- extrapolate:

  Fill the depths where no spectral peak was found within the window by
  linear interpolation between the found values (with the first and last
  found values held constant towards the ends) (logical, TRUE or FALSE);
  `FALSE` leaves `NA` at those depths. `Default=TRUE`.

- genplot:

  Plot the completed curve (red) together with the original tracked
  curve (black) (logical, TRUE or FALSE). `Default=FALSE`.

- keep_editable:

  Keep the graphical parameters of the plot so that extra elements can
  be added afterwards (logical, TRUE or FALSE); `FALSE` restores the
  previous `par` settings. `Default=FALSE`.

- ...:

  Arguments of older versions: `wavelet` is accepted as an alias of
  `scalogram`; other arguments are ignored.

## Value

Returns a matrix with 2 columns The first column is the depth axis The
second column is the completed tracking of the period a cycle of the
wavelet spectra

## Examples

``` r
# \donttest{
#Use the grey_track example points to complete the tracking of the
# precession cycle in the wavelet spectra of the grey scale data set
# of Zeeden et al., (2013).



grey_wt <-
 analyse_wavelet(
   data = grey,
   dj = 1/200,
   lowerPeriod = 0.02,
   upperPeriod = 256,
   verbose = FALSE,
   omega_nr = 8
 )

#The ~22kyr precession cycle is between 0.25 and 1m The grey_track data
#set is a pre-loaded uncompleted tracking of the precession cycle

#grey_track <- track_period_wavelet(
#astro_cycle = 22,
#wavelet = NULL,
#n.levels = 100,
#periodlab = "Period (meters)",
#x_lab = "depth (meters)"
#)



grey_track <- completed_series(
 scalogram = grey_wt,
 tracked_curve  = grey_track,
 period_up  = 1.25,
 period_down  = 0.75,
 extrapolate = TRUE,
genplot = FALSE,
keep_editable=FALSE
)
# }
```
