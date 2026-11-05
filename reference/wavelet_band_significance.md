# Significance of scale-averaged wavelet power (Torrence and Compo, section 4d)

Averages the wavelet power over a band of periods (a "scale-averaged
wavelet power" time series, eq. 24 of Torrence and Compo, 1998) and
tests it against a white or AR(1) red noise background. This is the test
to use when one asks "is there significant power in this Milankovitch
band at this depth?" without looking at every scale separately: the band
average is much less noisy than the power at one scale.

## Usage

``` r
wavelet_band_significance(
  object,
  period_range,
  background = c("ar1", "white"),
  siglvl = 0.95,
  method = c("exact", "torrence"),
  lag1 = NULL,
  delta_j0 = NULL,
  genplot = FALSE
)

# S3 method for class 'wr_band_significance'
plot(x, ...)
```

## Arguments

- object:

  Wavelet to test (list of class `analyze.wavelet`), the output of
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  (Morlet wavelet; its `omega_nr`, `dj`, `dt`, `Scale`, `Power` and `y`
  are used). Required, no default.

- period_range:

  Band of periods to average (numeric vector of length 2,
  `c(lowest, highest)`, in the units of the depth/time axis; at least
  one scale of the wavelet must fall inside), for example `c(0.3, 0.5)`
  for an obliquity band. Required, no default.

- background:

  Noise background of the null hypothesis (character, one of `"ar1"`,
  `"white"`): `"ar1"` red noise with the lag-1 autocorrelation `lag1`,
  `"white"` white noise (`lag1` is then set to 0). `Default="ar1"`.

- siglvl:

  Significance level of the test (numeric, between 0 and 1, e.g. 0.95
  for the 95% level). `Default=0.95`.

- method:

  How the null distribution of the band average is computed (character,
  one of `"exact"`, `"torrence"`): `"exact"` the weighted sum of
  exponential variables with a three-moment gamma fit, `"torrence"` the
  scaled chi-square distribution of Torrence and Compo (1998, eq. 28);
  see Details. `Default="exact"`.

- lag1:

  Lag-1 autocorrelation of the AR(1) background (numeric, between 0 and
  1). `NULL` estimates it from the record as in
  [`wavelet_significance`](https://stratigraphy.github.io/WaverideR/reference/wavelet_significance.md)
  (mean of the lag-1 autocorrelation and the square root of the lag-2
  autocorrelation). Ignored for `background = "white"`. `Default=NULL`.

- delta_j0:

  Decorrelation length of the scales, only used by `method = "torrence"`
  (numeric, \> 0); `NULL` uses `0.6 * 6 / omega_nr`, the value 0.6 of
  Torrence and Compo for `omega_nr = 6` rescaled to the `omega_nr` of
  `object`. `Default=NULL`.

- genplot:

  Plot the band-averaged power with the significance level (logical,
  TRUE or FALSE), using the `plot` method. `Default=FALSE`.

- x:

  Result of `wavelet_band_significance` (object of class
  `wr_band_significance`) to plot. Required, no default.

- ...:

  Further graphical arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) (e.g. `xlim`,
  `las`; not `xlab`, `ylab`, `ylim` or `main`, which are set by the
  method).

## Value

An object of class `wr_band_significance`, a list with: `x` (depth/time
axis), `power` (the scale-averaged power, in units of variance when the
record is standardised), `crit` (significance level of the band
average), `pvalue`, `significant` (logical, only points outside the cone
of influence), `valid`, `dof` (equivalent chi-square degrees of
freedom), `background_mean` (expected band average under the null),
`period_range`, `n_scales`, `lag1`, `Cdelta`, `delta_j0`, `method`,
`background`, `siglvl`. There are `print` and `plot` methods.

## Details

Two ways of obtaining the null distribution are offered. Both are
analytic and run instantaneously.

- `method = "exact"`:

  (default) The scale-average of the squared modulus of a linear
  transform of a Gaussian process is a weighted sum of independent
  exponential variables. The weights are the eigenvalues of a small
  (number of scales x number of scales) matrix that is computed from the
  wavelet and the background spectrum with the same Gaussian integrals
  that define the Morlet wavelet. This is exact for any `omega_nr`, any
  `dj`, any band and any background spectrum, and it needs no tabulated
  constants. The distribution of the weighted sum is approximated by a
  three-moment (shifted gamma) fit, which was checked against Monte
  Carlo simulations with the package's own
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  (see the package tests).

- `method = "torrence"`:

  The published closed form: the band average is a scaled chi-square
  variable with \\\nu = 2 n_a (S\_{avg}/S\_{mid}) \sqrt{1 + (n_a \delta
  j / \delta j_0)^2}\\ degrees of freedom (eq. 28), where \\n_a\\ is the
  number of scales, \\S\_{avg}\\ the harmonic-type mean scale and
  \\S\_{mid}\\ the geometric mid-scale of the band, and \\\delta j_0\\
  an empirical decorrelation length. Torrence and Compo only tabulated
  \\\delta j_0 = 0.60\\ for a Morlet wavelet with `omega_nr = 6`.
  Because the Morlet wavelet becomes narrower in log-scale when
  `omega_nr` increases, that constant under-estimates the degrees of
  freedom for larger values, so the test becomes too conservative (at
  `omega_nr = 10` only about 2% instead of 5% of white noise points
  exceed the 95% level). Here the default is rescaled as \\\delta j_0 =
  0.60 \cdot 6 / \omega_0\\, with \\\omega_0\\ the value of `omega_nr`,,
  which was fitted against the exact method for `omega_nr` between 4 and
  12 and is within about 10% of the exact degrees of freedom for bands
  wider than one octave. Use `method = "exact"` when the exact answer
  matters.

Only points whose whole period band lies outside the cone of influence
are tested.

The band average is \\\delta j \delta t / C\_\delta \sum_j \|W_j\|^2 /
s_j\\ (eq. 24). The factor \\C\_\delta\\ is evaluated numerically for
the `omega_nr` and `dj` of the object, normalised so that the average
over all scales of unit-variance white noise is 1. It only sets the
units of `power`; the significance decision does not depend on it.

The test is for a single point in depth. The band average is strongly
autocorrelated along depth, so the number of significant points expected
by chance is 5% of the points, but they occur in runs.

## References

Torrence, C., Compo, G.P., 1998. A practical guide to wavelet analysis.
Bulletin of the American Meteorological Society 79, 61-78.

## See also

[`wavelet_significance`](https://stratigraphy.github.io/WaverideR/reference/wavelet_significance.md)

## Examples

``` r
set.seed(1)
z <- seq(0, 60, by = 0.1)
d <- data.frame(z, sin(2 * pi * z / 4) + as.numeric(arima.sim(list(ar = 0.5), length(z))))
wt <- analyse_wavelet(d, dj = 1 / 10, lowerPeriod = 0.4, upperPeriod = 30,
                       verbose = FALSE, omega_nr = 8)
bs <- wavelet_band_significance(wt, period_range = c(3, 5))
bs
#> Significance of scale-averaged wavelet power (exact, ar1 background, level 0.95)
#>   band: periods 3.2 to 4.85 (7 scales)
#>   equivalent degrees of freedom: 3.84 
#>   expected band average under the null: 0.1003   critical value: 0.243 
#>   significant points: 429 of 513 outside the cone of influence
```
