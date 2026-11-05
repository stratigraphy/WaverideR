# Legacy power decomposition analysis (Acycle PDA menu / DYNOT core)

The `pda_legacy` function reproduces the power decomposition analysis as
run from the Acycle menu "Power Decomposition Analysis" and as used
inside DYNOT (`pdan.m`; Li et al., 2016). Preprocessing as in the Acycle
menu: rows with NA are removed, data are sorted, duplicated depths are
replaced by their mean, unevenly sampled data (relative tolerance 1e-5)
are linearly interpolated at the median sampling interval, and the
series is optionally padded by half a window at both ends
(`zeropad2.m`). In each window of `trunc(window/dt)` samples the data
are linearly detrended and the adaptive MTM spectrum
(`pmtm(x, nw, pad)`) is computed. Target and total powers are sums of
spectral bins selected with the historical proportional bin mapping
(lower = ceiling(nbin\*f/fNyq), upper = trunc(nbin\*f/fNyq));
overlapping target bands are counted once. Window centers are evenly
spaced between first coordinate + window/2 and last coordinate -
window/2. Translated from Acycle v3.1 pdan.m, zeropad2.m, interpolate.m,
findduplicate.m and AC.m (menu_pda_Callback) (Mingsong Li).

## Usage

``` r
pda_legacy(
  data = NULL,
  target_bands = c(1/45, 1/25),
  window = 500,
  nw = 2,
  f_min = 0,
  f_max = NULL,
  step = 1,
  pad = 5000,
  pad_type = 1,
  genplot = FALSE,
  verbose = FALSE
)
```

## Arguments

- data:

  Two-column numeric matrix or data frame: depth/time in column 1 and
  the proxy value in column 2 (further columns are ignored). Rows with
  `NA` are removed, the data are sorted, duplicated depths are averaged
  and unevenly sampled data are interpolated at the median sampling
  interval. Required (`Default=NULL`).

- target_bands:

  Target frequency bands whose summed power is divided by the power of
  the total band, in cycles per unit of the depth/time column: a numeric
  vector of paired limits `c(lower1, upper1, lower2, upper2, ...)` (any
  order within a pair) or a B x 2 matrix with one band per row;
  overlapping bands are counted once. `Default=c(1/45, 1/25)`.

- window:

  Length of the moving window in the unit of the depth/time column
  (numeric, \> 0); each window holds `trunc(window / dt)` samples.
  `Default=500`.

- nw:

  Time-bandwidth product of the DPSS (Slepian) tapers of the multitaper
  spectrum (numeric, \> 0); larger values smooth the spectrum more.
  `Default=2`.

- f_min:

  Lower cutoff frequency of the total band, in cycles per unit of the
  depth/time column (numeric, \>= 0). `Default=0`.

- f_max:

  Upper cutoff frequency of the total band, in cycles per unit of the
  depth/time column (numeric, \> `f_min`); values above the Nyquist
  frequency are clipped. `NULL` uses the Nyquist frequency.
  `Default=NULL`.

- step:

  Step between successive windows, in number of samples (positive
  integer). `Default=1`.

- pad:

  Length of the FFT used by the multitaper spectrum, i.e. the
  zero-padding number of `pmtm` (positive integer); larger values give a
  finer frequency grid. `Default=5000`.

- pad_type:

  Padding of both ends of the (demeaned) series by half a window before
  the analysis, so that the windows cover the whole record (integer):
  `0` none, `1` zeros, `2` mirror image of the record, `3` mean of the
  first/last half window, `4` random values with the mean and standard
  deviation of the first/last half window. `Default=1`.

- genplot:

  Plot the power ratio against depth/time (logical, TRUE or FALSE).
  `Default=FALSE`.

- verbose:

  Show a progress bar over the windows (logical, TRUE or FALSE).
  `Default=FALSE`.

## Value

A matrix with columns window center, ratio (target / total), target
power and total power.

## References

Li, M., Huang, C., Hinnov, L., Ogg, J., Chen, Z.-Q., Zhang, Y., 2016.
Obliquity-forced climate during the Early Triassic hothouse in China.
Geology 44, 623-626.
[doi:10.1130/G37970.1](https://doi.org/10.1130/G37970.1)

Li, M., Hinnov, L., Kump, L., 2019. Acycle: Time-series analysis
software for paleoclimate research and education. Computers &
Geosciences 127, 12-22.
[doi:10.1016/j.cageo.2019.02.011](https://doi.org/10.1016/j.cageo.2019.02.011)

## Author

Michiel Arts (R translation); original MATLAB code: Mingsong Li and
Linda Hinnov (Acycle); zeropad2 padding options by Nicolas Thibault and
Giovanni Rizzi.

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 2000, by = 2)
x <- sin(2 * pi * t / 40) + 0.5 * sin(2 * pi * t / 100) +
  0.5 * as.numeric(arima.sim(list(ar = 0.7), n = length(t)))
res <- pda_legacy(data = cbind(t, x), target_bands = c(1/45, 1/35),
  window = 300, step = 5, pad = 1000, pad_type = 1)
head(res)
#>      center     ratio target_power total_power
#> [1,]      0 0.3606305     24.11534    66.86996
#> [2,]     10 0.3589282     25.18059    70.15496
#> [3,]     20 0.3444931     25.51678    74.07050
#> [4,]     30 0.3392978     25.92545    76.40913
#> [5,]     40 0.3335794     26.36172    79.02684
#> [6,]     50 0.3324595     27.78380    83.57049
# }
```
