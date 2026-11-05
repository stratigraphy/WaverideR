# Windowed timeOpt with precession, short-eccentricity and obliquity interpretations

Runs
[`win_timeOpt`](https://stratigraphy.github.io/WaverideR/reference/win_timeOpt.md)
with fit 1 (precession amplitude modulation), fit 2 (short-eccentricity
amplitude modulation) and fit 3 (obliquity amplitude modulation) on the
same record and compares their best sedimentation rates window by
window.

## Usage

``` r
win_timeOpt_multi(data = NULL, fits = c(1, 2, 3), ratio_tol = 0.15, ...)
```

## Arguments

- data:

  Record to analyse (matrix or data frame with two columns: depth in
  metres in column 1, proxy value in column 2), evenly spaced; passed to
  [`win_timeOpt`](https://stratigraphy.github.io/WaverideR/reference/win_timeOpt.md).
  Must be given (the `NULL` default fails). `Default=NULL`.

- fits:

  The timeOpt fits to run (numeric vector with any of 1, 2, 3): 1
  precession amplitude modulation, 2 short-eccentricity amplitude
  modulation, 3 obliquity amplitude modulation (the `fit` argument of
  [`win_timeOpt`](https://stratigraphy.github.io/WaverideR/reference/win_timeOpt.md)).
  The agreement and ambiguity flags need the corresponding pairs.
  `Default=c(1, 2, 3)`.

- ratio_tol:

  Relative tolerance used to compare the best rates of two fits
  (numeric, \> 0; 0.15 counts rates within a factor 1.15 as equal); used
  for the `agree_*` and `ambiguous_*` flags. `Default=0.15`.

- ...:

  Further arguments passed to
  [`win_timeOpt`](https://stratigraphy.github.io/WaverideR/reference/win_timeOpt.md)
  for every fit (e.g. `window_size`, `sedmin`, `sedmax`, `numsed`,
  `n_null`); `fit` is not allowed (use `fits`).

## Value

A list of class `win_timeOpt_multi.result` with the
[`win_timeOpt`](https://stratigraphy.github.io/WaverideR/reference/win_timeOpt.md)
result of every fit (`fit1`, `fit2`, `fit3`) and `summary`, a data frame
with per window the best rate and significance of each fit and the flags
above.

## Details

A record that holds only one cycle family can be read in two ways:
obliquity cycles read as precession give a peak at about twice the true
rate (fit 1), precession cycles read as obliquity give a peak at about
half the true rate (fit 3). For every window the function flags

- agree_12, agree_13, agree_23:

  both fits are significant and their best rates agree within
  `ratio_tol`: independent support for the same rate

- ambiguous_13:

  both fit 1 and fit 3 are significant and the fit-1 rate is about twice
  the fit-3 rate: precession/obliquity 2:1 ambiguity, which needs other
  evidence (eccentricity cycles, an age tie point) to resolve

- ambiguous_23:

  both fit 2 and fit 3 are significant and the fit-2 rate is about 0.41
  times the fit-3 rate: obliquity (41 kyr with its 173-kyr modulation)
  read as short eccentricity (~100 kyr with the 405-kyr modulation)

## Examples

``` r
# \donttest{
set.seed(1)
mag_multi <- win_timeOpt_multi(data = mag[1:250, ], fits = c(1, 2), window_size = 6,
    sedmin = 0.1, sedmax = 1, numsed = 5, n_null = 6, noise_model = "ar1", step = 25)
head(mag_multi$summary)
#>      x best_fit1 z_fit1 sig_fit1 best_fit2     z_fit2 sig_fit2 agree_12
#> 1 0.00        NA     NA    FALSE 0.3162278 -2.2376431    FALSE    FALSE
#> 2 1.25        NA     NA    FALSE 0.1778279 -1.3857398    FALSE    FALSE
#> 3 2.50        NA     NA    FALSE 0.5623413  0.3560039    FALSE    FALSE
#> 4 3.75        NA     NA    FALSE 0.5623413 -0.2735087    FALSE    FALSE
#> 5 5.00        NA     NA    FALSE 0.3162278 -0.3565558    FALSE    FALSE
#> 6 6.25        NA     NA    FALSE 0.5623413 -0.1424578    FALSE    FALSE
# }
```
