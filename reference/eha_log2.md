# Generate a log2 spaced EHA spectra

Compute a log2 spaced Evolutive Harmonic Analysis (EHA) & Evolutive
Power Spectral Analysis This is a wrapper function for the "eha"
function of the 'astrochron' R package

## Usage

``` r
eha_log2(
  data = NULL,
  win = NULL,
  tbw = NULL,
  demean = NULL,
  detrend = NULL,
  upperPeriod = NULL,
  lowerPeriod = NULL,
  pad = NULL,
  padding = "noise"
)
```

## Arguments

- data:

  Input series (matrix or data frame with two columns: depth/time in
  column 1, proxy value in column 2), preferably evenly sampled (the
  padded series is linearly interpolated with
  [`astrochron::linterp`](https://rdrr.io/pkg/astrochron/man/linterp.html)).
  Required (`Default=NULL`).

- win:

  Length of the moving EHA window, in the depth/time unit of `data`
  (numeric, \> 0, shorter than the record). The series is also padded by
  half a window at both ends (see `padding`). Required (`Default=NULL`).

- tbw:

  Time-bandwidth product of the multitaper (MTM) spectra, passed to
  [`eha`](https://rdrr.io/pkg/astrochron/man/eha.html) (numeric, \<= 10,
  e.g. 2, 3 or 4); larger values give smoother spectra with a lower
  frequency resolution. Required (`Default=NULL`).

- demean:

  Remove the mean of each window (logical). Currently not used: the
  windows are always demeaned (`demean=TRUE` is passed to
  [`eha`](https://rdrr.io/pkg/astrochron/man/eha.html)). `Default=NULL`.

- detrend:

  Remove a linear trend from each window (logical). Currently not used:
  the windows are always detrended (`detrend=TRUE` is passed to
  [`eha`](https://rdrr.io/pkg/astrochron/man/eha.html)). `Default=NULL`.

- upperPeriod:

  Longest period to analyse, in the depth/time unit of `data` (numeric,
  \> `lowerPeriod`, at most about `win`); sets `fmin = 1 / upperPeriod`
  of the EHA and the top of the log2 period grid. Required
  (`Default=NULL`).

- lowerPeriod:

  Shortest period to analyse, in the depth/time unit of `data` (numeric,
  \> 0, at least twice the sampling interval); sets
  `fmax = 1 / lowerPeriod` and the bottom of the log2 period grid.
  Required (`Default=NULL`).

- pad:

  Number of points to which each window is zero-padded (integer or
  `NULL`); odd values are increased by one, values below the number of
  points per window are raised to it, and `5 * pad` is passed to
  [`eha`](https://rdrr.io/pkg/astrochron/man/eha.html) to obtain a fine
  low-frequency spacing. `NULL` uses
  `2 * 2^ceiling(log2(points per window))`. `Default=NULL`.

- padding:

  Values used to extend the series by half a window at both ends, so
  that the EHA covers the whole record (character): `"noise"` (AR(1) red
  noise with the mean, standard deviation and lag-1 autocorrelation of
  the end segments; note that this calls
  [`set.seed()`](https://rdrr.io/r/base/Random.html)), `"mean"` (also
  `"average"` or `"avg"`; the mean of the series) or `"zero"` (zeros).
  Any other value is treated as `"mean"`. `Default="noise"`.

## Value

The output is a list (analyze.eha object) which contain 12 objects which
are the result of running the eha. Object 1: depth - depth/time of axis
Object 2: log2_period - log2 scales period Object 3: pwr - Power values
of the EHA run Object 4: amp - Amplitude values of the EHA run Object 5:
prob - Probability values of the EHA run Object 6: f_test - F test
values of the EHA run Scale size Object 7: Power.avg - Average power
values Object 8: amp.avg - Average amplitude values Object 9: prob.avg -
Average probability value Object 10: f_test.avg - Average f test values
Object 11: x - x values of the data set Object 12: y - y values of the
data set Object 13: window size

## References

S.R. Meyers, 2012, Seeing Red in Cyclic Stratigraphy: Spectral Noise
Estimation for \#'Astrochronology: Paleoceanography, 27, PA3228,
\<doi:10.1029/2012PA002307\>

S.R. Meyers, 2019 Cyclostratigraphy and the problem of astrochronologic
testing, Earth-Science Reviews, Volume 190,
\<doi.org/10.1016/j.earscirev.2018.11.015.\>

## Author

Code based on on the "eha" function of the 'astrochron' R package

## Examples

``` r
# \donttest{
#Example 1. Total Solar Irradiance
# data set of Steinhilber et al., (2012)

TSI_eha_log2 <- eha_log2(data = TSI,
win =  8192,
tbw = 4,
demean = TRUE,
detrend = TRUE,
upperPeriod = 8192,
lowerPeriod = 50,
pad = NULL,
padding = "noise")

#Example 2. Magnetic susceptibility data set of Pas et al., (2018)
mag_eha_log2 <- eha_log2(data = mag,
win =  50,
tbw = 4,
demean = TRUE,
detrend = TRUE,
upperPeriod = 50,
lowerPeriod = 2,
pad = NULL,
padding = "noise")

#Example 3. Greyscale data set of Zeeden et al., (2013)
grey_eha_log2 <- eha_log2(data = grey,
win =  20,
tbw = 4,
demean = TRUE,
detrend = TRUE,
upperPeriod = 20,
lowerPeriod = 2,
pad = NULL,
padding = "noise")

# }
```
