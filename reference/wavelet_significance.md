# Significance of wavelet and superlet power: analytic and surrogate (Monte Carlo) levels

Calculates the significance level of the spectral power in a wavelet or
superlet scalogram against a noise background, and which cells of the
scalogram are significant.

Two ways are available:

- `method = "analytic"`:

  (wavelets only) the near instantaneous method of Torrence and Compo
  (1998): the theoretical Fourier spectrum of a lag-1 autoregressive
  (red noise) process times a chi-square distribution with 2 degrees of
  freedom for the local wavelet power (their equations 16 and 18), and
  for the global (time averaged) spectrum the chi-square distribution
  with the degrees of freedom of their equation 23. No simulations are
  needed.

- `method = "surrogate"`:

  (wavelets and superlets) noise surrogates of the record are
  transformed in exactly the same way as the data. For every period the
  `siglvl` quantile of the power outside the cone of influence is taken,
  and averaged over the surrogates (scale wise critical level, after
  Weesp et al.). The `siglvl` quantile of the time averaged power of the
  surrogates gives the critical level of the global spectrum. The
  background can be white noise, AR(1) red noise, AR(1) noise estimated
  after removing the periods of the candidate cycles (so that the cycles
  that are tested do not inflate the noise background), or a power law
  (\\1/f^\beta\\) noise.

The result can be added to the plots with the argument `add_signif` of
[`plot_wavelet`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md)
and
[`plot_superlet`](https://stratigraphy.github.io/WaverideR/reference/plot_superlet.md).

## Usage

``` r
wavelet_significance(
  object,
  method = NULL,
  background = c("ar1", "white", "ar1_corrected", "powerlaw"),
  siglvl = 0.95,
  n_sim = 200,
  lag1 = NULL,
  notch_periods = NULL,
  notch_tol = 0.2,
  candidates = NULL,
  candidate_tol = 0.2,
  gamma = NULL,
  run_multicore = FALSE,
  n_cores = NULL,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- object:

  Scalogram to test: a wavelet object (class `analyze.wavelet`, from
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md))
  or a superlet object (class `analyze.superlet`, from
  [`analyse_superlet`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md)).
  Synchrosqueezed transforms are not accepted. Required, no default.

- method:

  How the critical levels are computed (character, `"analytic"` or
  `"surrogate"`): `"analytic"` the Torrence and Compo (1998) levels
  (wavelets only, instantaneous), `"surrogate"` Monte Carlo levels from
  `n_sim` noise surrogates transformed like the data. `NULL` chooses
  `"analytic"` for wavelets and `"surrogate"` for superlets.
  `Default=NULL`.

- background:

  Noise background of the null hypothesis (character, one of `"ar1"`,
  `"white"`, `"ar1_corrected"`, `"powerlaw"`): `"ar1"` red noise with
  lag-1 autocorrelation `lag1`; `"white"` white noise; `"ar1_corrected"`
  AR(1) noise whose coefficient is estimated after removing the bands
  around `notch_periods`; `"powerlaw"` \\1/f^\beta\\ noise with
  \\\beta\\ fitted to the periodogram (outside the `notch_periods`
  bands). The analytic method supports only `"ar1"` and `"white"`.
  `Default="ar1"`.

- siglvl:

  Significance level (numeric, between 0 and 1, e.g. 0.95 for the 95%
  level). `Default=0.95`.

- n_sim:

  Number of noise surrogates (integer, \>= 1), only used with
  `method = "surrogate"`; more surrogates give more stable levels and
  take longer. Use 1000 or more for publication. `Default=200`.

- lag1:

  Lag-1 autocorrelation of the red noise (numeric, between 0 and 1),
  used for `"ar1"` and `"ar1_corrected"`. `NULL` estimates it from the
  record: maximum likelihood (AR(1) fit) for the surrogates, the lag-1
  and lag-2 autocorrelation \\(\alpha_1+\sqrt{\alpha_2})/2\\ of Torrence
  and Compo for the analytic method. `Default=NULL`.

- notch_periods:

  Periods of the candidate cycles (numeric vector, in the depth/time
  unit of the record) that are removed before the AR(1) coefficient
  (`"ar1_corrected"`) or the power law exponent (`"powerlaw"`) is
  estimated, so that the cycles that are tested do not inflate the noise
  background; `NULL` removes nothing. `Default=NULL`.

- notch_tol:

  Relative half width of the removed bands around `notch_periods`
  (numeric, between 0 and 1; 0.2 removes periods within 20 percent of
  every notch period). `Default=0.2`.

- candidates:

  Optional periods of candidate cycles (numeric vector, optionally
  named, in the depth/time unit of the record), only used with
  `method = "surrogate"`. For each, the p-value of the largest time
  averaged power within `candidate_tol` of the period is calculated
  against the surrogates (the largest value in the band is compared with
  the largest value of every surrogate, which corrects for searching a
  band). `Default=NULL`.

- candidate_tol:

  Relative half width of the candidate bands (numeric, between 0 and 1).
  `Default=0.2`.

- gamma:

  Decorrelation factor of the wavelet used in the degrees of freedom of
  the global spectrum of the analytic method (numeric, \> 0; Torrence
  and Compo 1998, eq. 23). `NULL` takes the value of Torrence and Compo
  (1998, table 2) for the wavelet of the object: 2.32 (Morlet), 1.17
  (Paul), 1.43 (DOG of order 2); 2.32 is used for the other wavelets.
  `Default=NULL`.

- run_multicore:

  Run the surrogates in parallel on several cores (logical, TRUE or
  FALSE); with `TRUE` and `n_cores = NULL` a safe number of cores is
  chosen automatically. `Default=FALSE`.

- n_cores:

  Number of cores for the surrogates (integer, \>= 1); `NULL` uses the
  option `waverider.n_cores` or the choice of `run_multicore`. The
  number is capped by the available cores and memory (see
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)).
  `Default=NULL`.

- seed:

  Seed of the random number generator for the surrogates (integer), so
  that the surrogate levels can be reproduced; `NULL` derives it from
  the current random state (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  also gives reproducible results). Results do not depend on the number
  of cores, and the random state of the session is not changed.
  `Default=NULL`.

- verbose:

  Show a progress bar of the surrogate simulations (logical, TRUE or
  FALSE). `Default=FALSE`.

## Value

A list of class `"wr_significance"` with `period`, `crit` (critical
power level per period, in the units of `object$Power`), `gws`,
`gws_crit` (time averaged power outside the cone of influence and its
critical level, per period), `significant` (logical matrix with the same
layout as `object$Power`: power above the critical level and outside the
cone of influence), `valid` (outside the cone of influence),
`candidates` (data frame of p-values), and the settings and estimated
noise parameters (`lag1`, `beta`).

## References

Torrence, C., and G. P. Compo, 1998. A practical guide to wavelet
analysis. Bulletin of the American Meteorological Society 79, 61-78.

Weesp analysis: surrogate significance with white, AR(1), signal
corrected AR(1) and power law backgrounds.

Liu, Y., San Liang, X., Weisberg, R. H., 2007. Rectification of the bias
in the wavelet power spectrum. Journal of Atmospheric and Oceanic
Technology 24, 2093-2102.

## Examples

``` r
set.seed(1)
z <- seq(0, 60, by = 0.1)
d <- data.frame(z, sin(2 * pi * z / 4) + as.numeric(arima.sim(list(ar = 0.5), length(z))))
wt <- analyse_wavelet(d, dj = 1 / 10, lowerPeriod = 0.4, upperPeriod = 30,
                       verbose = FALSE, omega_nr = 6)
sig <- wavelet_significance(wt)                    # analytic, instantaneous
sig$period[which(sig$gws > sig$gws_crit)]          # periods with significant global power
#> [1] 3.429675 3.675835 3.939662 4.222425 4.525483
# \donttest{
sig2 <- wavelet_significance(wt, method = "surrogate", background = "ar1", n_sim = 20, seed = 1)
# }
```
