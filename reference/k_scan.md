# Scan k and count how many measured windows a single k can satisfy

Scan k and count how many measured windows a single k can satisfy

## Usage

``` r
k_scan(
  peaks,
  terms = c("k+s3", "k+g5", "k+g2", "k+g4"),
  k_grid = seq(12, 40, by = 0.02),
  tol = 0.04,
  gs = LA2004
)
```

## Arguments

- peaks:

  Measured spectral peak periods in kyr (numeric vector, \> 0). Required
  (no default).

- terms:

  Astronomical terms whose periods are predicted for every trial k
  (character vector), any of `"k+s3"` (obliquity), `"k+g5"`, `"k+g2"`,
  `"k+g4"` and `"k+g3"` (precession).
  `Default=c("k+s3", "k+g5", "k+g2", "k+g4")`.

- k_grid:

  Grid of trial precession constants, expressed as periods in kyr
  (numeric vector, \> 0; the output column `k_kyr`).
  `Default=seq(12, 40, by = 0.02)`.

- tol:

  Relative tolerance for a match between a predicted and a measured
  period (numeric, \> 0; 0.04 is 4%): a term counts as matched when some
  peak lies within a factor `1 + tol` of its predicted period.
  `Default=0.04`.

- gs:

  Secular frequencies in arcsec/yr: a list with elements `g` and `s`,
  each a numeric vector of the nine fundamental frequencies g1-g9 and
  s1-s9. The default `LA2004` is the La2004 table (Laskar et al., 2004)
  defined inside the function. `Default=LA2004`.
