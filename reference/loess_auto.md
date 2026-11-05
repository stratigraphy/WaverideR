# Perform an automatically loess based smoothing of a time series

Perform an automatically loess based smoothing of a time series. The
local polynomial regression with automatic smoothing parameter selection
is based on an optimisation using the 'aicc' bias-corrected 'AIC'
criterion and the 'gcv' generalized cross-validation criterion.

## Usage

``` r
loess_auto(
  time_series = NULL,
  genplot = FALSE,
  print_span = FALSE,
  keep_editable = FALSE
)
```

## Arguments

- time_series:

  Series to smooth: a two-column matrix or data frame with depth or time
  in column 1 and the proxy value (e.g. the period of a tracked cycle,
  as from
  [`completed_series`](https://stratigraphy.github.io/WaverideR/reference/completed_series.md))
  in column 2; rows with `NA` are dropped for the span selection.
  Required (`Default=NULL`).

- genplot:

  Plot the original series (blue), the smoothed curve (black) and the
  smoothed curve plus and minus the absolute residual (red, dashed)
  (logical, TRUE or FALSE). `Default=FALSE`.

- print_span:

  Print the selected loess span, a fraction (0-1) of the number of data
  points that is used in each local fit (logical, TRUE or FALSE).
  `Default=FALSE`.

- keep_editable:

  Keep the graphical parameters (`par`) of the plot after the function
  returns, so that features can be added to it (logical, TRUE or FALSE);
  `FALSE` restores the previous settings. `Default=FALSE`.

## Value

A matrix with 3 columns. The first column is depth/time. The second
column is the smoothed curve. The third column is difference between the
original curve and the smoothed curve.

## References

Cleveland, W. S. (1979) Robust locally weighted regression and smoothing
scatter plots. Journal of the American Statistical Association. 74,
829–836. \<doi:10.1080/01621459.1979.10481038\> Hurvich, C.M., Simonoff,
J.S., and Tsai, C.L. (1998), Smoothing Parameter Selection in
Nonparametric Regression Using an Improved Akaike Information Criterion.
Journal of the Royal Statistical Society B. 60, 271–293
\<doi:10.1111/1467-9868.00125\> Golub, G., Heath, M. and Wahba, G.
(1979). Generalized cross validation as a method for choosing a good
ridge parameter. Technometrics. 21, 215–224. \<doi:10.2307/1268518\>

## Author

Based on the the loess.as function of the 'fANCOVA' R package.

## Examples

``` r
# \donttest{
#'smooth the period curve of the 405 kyr eccentricity cycle extracted from
# the magnetic susceptibility data set of Pas et al., (2018)
#perform the CWT on the magnetic susceptibility data set of Pas et al., (2018)

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
  genplot = FALSE,
  keep_editable=FALSE
)

#Smooth the completed tracking of the 405 kyr eccentricity cycle as tracked in the wavelet spectra
mag_track_complete <- loess_auto(time_series = mag_track_complete,
genplot = FALSE, print_span = FALSE,keep_editable=FALSE)
# }
```
