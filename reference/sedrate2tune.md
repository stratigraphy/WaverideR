# Use a sedimentation curve to convert data to the time domain

Convert a proxy record from the depth to time domain using a
sedimentation rate curve

## Usage

``` r
sedrate2tune(
  data = NULL,
  sed_curve = NULL,
  genplot = FALSE,
  keep_editable = FALSE
)
```

## Arguments

- data:

  Proxy record in the depth domain (required): matrix or data frame with
  depth in metres in column 1 and the proxy value in column 2. The proxy
  values are paired row by row with the times computed on the (regularly
  resampled) depths of `sed_curve`, so `data` must be evenly sampled on
  the same depths as `sed_curve` (e.g. the record from which the
  sedimentation rate curve was made). `Default=NULL`.

- sed_curve:

  Sedimentation rate curve (required): matrix or data frame with depth
  in metres in column 1 and the sedimentation rate in cm/kyr in column
  2, e.g. the output of
  [`minimal_tuning`](https://stratigraphy.github.io/WaverideR/reference/minimal_tuning.md)
  with `output = 1`. It is sorted, cleaned of `NA` rows and linearly
  interpolated to its median depth step before being integrated to time
  (in kyr). `Default=NULL`.

- genplot:

  Plot the proxy record against time (kyr) (logical). `Default=FALSE`.

- keep_editable:

  Keep the graphical parameters of the plot after the function returns,
  so that extra elements can be added (logical); only used with
  `genplot = TRUE`. `Default=FALSE`.

## Value

The output is a matrix with 2 columns. The first column is time The
second column is the proxy value If `genplot=TRUE` then a time vs proxy
value plot will be plotted.

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
# from the magnetic susceptibility record of the Sullivan core
# of Pas et al., (2018) and then create a age model using minimal tuning
# (e.g.) set the distance between peaks to 405 kyr. The age model
# (sedimentation rate curve) is then used to convert the data
# from the depth to the time domain

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
 add_mean = TRUE,
 plot_residual = FALSE,
 keep_editable = FALSE
)

mag_405_min_tuning <- minimal_tuning(data = mag_405,
pts = 5,
cycle = 405,
tune_opt = "max",
output = 1,
genplot = FALSE,
keep_editable = FALSE)

mag_time <- sedrate2tune(
data=mag,
sed_curve=mag_405_min_tuning,
genplot=FALSE,
keep_editable=FALSE)

# }
```
