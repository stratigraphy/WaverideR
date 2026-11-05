# Significance of the lines of a (de-chirped) EHA map against warped noise

The de-chirped map of
[`analyse_chirp_eha`](https://stratigraphy.github.io/WaverideR/reference/analyse_chirp_eha.md)
warps the noise inside every window as well as the cycles, so wavy
structure in it is not by itself evidence. This function runs the same
analysis, with the same windows, gradients (or `steer`) and settings, on
AR1 surrogates of the record and compares the observed line-to-continuum
ratios with those of the surrogates cell by cell. Lines that are only
artefacts of the warping appear in the surrogates as well and are not
significant.

## Usage

``` r
chirp_eha_significance(ce, data, n_sim = 20, rho = NULL, seed = 1)
```

## Arguments

- ce:

  Result of
  [`analyse_chirp_eha`](https://stratigraphy.github.io/WaverideR/reference/analyse_chirp_eha.md)
  (list of class `analyze.wavelet` with its `args` element); its
  settings are reused for the surrogates. Required, no default.

- data:

  The record `ce` was computed from: a matrix or data frame with depth
  in column 1 and the proxy value in column 2. Required, no default.

- n_sim:

  Number of AR(1) surrogates (integer, \> 0). More surrogates give more
  precise p-values and 95th percentiles but take longer (each one
  repeats the full analysis). `Default=20`.

- rho:

  Lag-1 autocorrelation of the AR(1) surrogates (numeric, between -1 and
  1); `NULL` estimates it from the record. `Default=NULL`.

- seed:

  Seed of the random number generator for the surrogates (integer), so
  that the result can be reproduced. `Default=1`.

## Value

`ce` with `Power_sig` (observed `Power_norm` divided by the 95th
percentile of the surrogates: above 1 is significant at that level),
`Power_pval` (fraction of surrogates at or above the observed value) and
`surrogate_q95`.
