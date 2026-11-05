# Build dynamic-filter control points from a tracked period curve

The `dynamic_filter_controls` function converts a tracked period curve
\[depth, period\] into pass band control points for
[`dynamic_filter`](https://stratigraphy.github.io/WaverideR/reference/dynamic_filter.md).
This replaces the interactive clicking of lower/upper frequency
boundaries on an evolutionary spectrum in Acycle (dynamic_filter_lang.m)
by a reproducible, non-interactive step. The curve can come directly
from the WaverideR functions `track_period_wavelet` (first two columns:
depth, period) or `completed_series` (depth, period), so cycles already
tracked in a wavelet scalogram can be filtered with the Acycle dynamic
filter. The centre frequency is `1 / period`; the pass band is
`f_centre * (1 -/+ rel_bandwidth / 2)`, i.e. `rel_bandwidth` is the
half-power bandwidth relative to the centre frequency (the band is
symmetric in frequency, as the Acycle Gaussian response is). Written for
the R translation of Acycle v3.1 dynamic filtering (Mingsong Li); no
MATLAB equivalent.

## Usage

``` r
dynamic_filter_controls(
  tracked = NULL,
  rel_bandwidth = 0.2,
  n_points = NULL,
  data = NULL
)
```

## Arguments

- tracked:

  Tracked period curve of the cycle (matrix or data frame whose first
  two columns are depth/time and period in the same depth/time unit,
  period \> 0; further columns are ignored), e.g. the output of
  [`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
  or
  [`completed_series`](https://stratigraphy.github.io/WaverideR/reference/completed_series.md).
  Rows with non-finite values are dropped, rows are sorted by depth and
  duplicated depths are removed (first kept). Required (`Default=NULL`).

- rel_bandwidth:

  Half-power bandwidth of the pass band relative to the centre frequency
  `1 / period` (numeric, \> 0 and \< 2): the band is
  `f_centre * (1 -/+ rel_bandwidth / 2)`, so larger values give a wider
  band. `Default=0.2`.

- n_points:

  Maximum number of control points (integer, \>= 1, or `NULL`). When the
  tracked curve has more rows, it is linearly interpolated to `n_points`
  evenly spaced depths spanning the (clipped) curve (`n_points=1` gives
  one point at the middle). `NULL` keeps all rows. `Default=NULL`.

- data:

  Series that will be filtered (two-column matrix or data frame,
  depth/time in column 1, or `NULL`). When given, control points outside
  its depth range are removed, because
  [`dynamic_filter`](https://stratigraphy.github.io/WaverideR/reference/dynamic_filter.md)
  requires control depths inside the data range. `Default=NULL`.

## Value

A matrix with columns `depth`, `f_low`, `f_high` and `f_centre`,
directly usable as the `controls` argument of
[`dynamic_filter`](https://stratigraphy.github.io/WaverideR/reference/dynamic_filter.md).

## References

Arts, M., 2023. WaverideR: Extracting signals from wavelet spectra.
\<https://CRAN.R-project.org/package=WaverideR\>

Li, M., Hinnov, L., Kump, L., 2019. Acycle: Time-series analysis
software for paleoclimate research and education. Computers &
Geosciences 127, 12-22. \<doi:10.1016/j.cageo.2019.02.011\>

## Author

Michiel Arts

## Examples

``` r
# \donttest{
depth <- seq(0, 300, by = 1)
period <- 5 + 10 * depth / 300
ctrl <- dynamic_filter_controls(tracked = cbind(depth, period),
  rel_bandwidth = 0.3, n_points = 20)
head(ctrl)
#>         depth     f_low    f_high  f_centre
#> [1,]  0.00000 0.1700000 0.2300000 0.2000000
#> [2,] 15.78947 0.1538095 0.2080952 0.1809524
#> [3,] 31.57895 0.1404348 0.1900000 0.1652174
#> [4,] 47.36842 0.1292000 0.1748000 0.1520000
#> [5,] 63.15789 0.1196296 0.1618519 0.1407407
#> [6,] 78.94737 0.1113793 0.1506897 0.1310345
# }
```
