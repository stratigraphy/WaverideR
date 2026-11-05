# exponentially ramped sedimentation rate model to convert time to stratigraphy

Apply a exponentially increasing increasing sedimentation rate model to
convert time to stratigraphy.

## Usage

``` r
expSedRamp(data, sr_start, sr_end)
```

## Arguments

- data:

  Time series to convert to stratigraphy (matrix or data frame with two
  columns: time in kyr in column 1, increasing and preferably evenly
  sampled, proxy value in column 2). Required, no default.

- sr_start:

  Sedimentation rate at the start (minimum time) of the series, in
  cm/kyr (numeric, \> 0). Required, no default.

- sr_end:

  Sedimentation rate at the end (maximum time) of the series, in cm/kyr
  (numeric, \> 0). The rate changes exponentially (linearly in log2)
  from `sr_start` to `sr_end`. Required, no default.

## Value

Returns a list which contains 10 elements element 1: time series in the
depth (m) domain element 2: sedimentation rate curve

## References

Routines for astrochronologic testing, astronomical time scale
construction, and time series analysis
\<doi:10.1016/j.earscirev.2018.11.015\>

## Author

Based on the [sedRamp](https://rdrr.io/pkg/astrochron/man/sedRamp.html)
function of the 'astrochron' R package.

## Examples

``` r

data_1 <- astrochron::etp(
tmin = 0,
tmax = 4000,
dt = 1,
eWt = 1.5,
oWt = 0.75,
pWt = 1,
esinw = TRUE,
solution = NULL,
standardize = TRUE,
genplot = FALSE,
verbose = FALSE
)
data_ramped_ls <- expSedRamp(data_1, sr_start = 1, sr_end = 5)
```
