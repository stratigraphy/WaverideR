# Consistency-test inversion (the estimator used in the paper)

After Zeeden et al. (2023): draw a precession constant and a set of
secular frequencies, predict every term, and keep the draw only if all
measured windows are satisfied at once. Two details matter. First, the
secular frequencies are sampled as periods with the Laskar et al. (2011)
uncertainties scaled by sqrt(age / 100 Myr). Second, because a record
tuned to a nominal 405 kyr inherits the error in the true g2-g5 period,
every measured period is rescaled by (sampled 405-kyr period) / 405
before it is compared with the prediction; this is what allows the tuned
time scale itself to stretch slightly.

## Usage

``` r
astro_ctest(
  windows,
  age_ma,
  nsim = 3000,
  k_range = c(10, 34),
  anchor_kyr = 405,
  seed = 1,
  block = 2e+05,
  max_attempts = 5e+07,
  verbose = TRUE
)
```

## Arguments

- windows:

  Measured period windows of the astronomical terms, in kyr on the tuned
  time scale (named list; each element a numeric vector `c(low, high)`).
  Allowed names: `"e405"`, `"k+s3"` (obliquity), `"k+g5"`, `"k+g2"`,
  `"k+g4"` and `"k+g3"` (precession), e.g.
  `list("e405" = c(400, 410), "k+s3" = c(29, 30.4))`. A draw is accepted
  only if every predicted period, rescaled by (sampled 405-kyr
  period)/`anchor_kyr`, falls inside its window. Required, no default.

- age_ma:

  Age of the record in Ma (numeric): either `c(min, max)`, or a single
  value, in which case the range age_ma - 5 to age_ma + 5 Ma is used.
  Ages are sampled uniformly in this range and scale the uncertainty of
  the secular frequencies by sqrt(age/100 Myr). Required, no default.

- nsim:

  Number of accepted (consistent) draws to collect (integer, \> 0).
  `Default=3000`.

- k_range:

  Prior range of the precession constant k in kyr (numeric vector
  `c(min, max)`); k is drawn uniformly in this range.
  `Default=c(10, 34)`.

- anchor_kyr:

  Period (kyr) that was assumed for the 405-kyr eccentricity cycle when
  the record was tuned (numeric, \> 0); the measured windows are
  rescaled by (sampled 405-kyr period)/`anchor_kyr`. `Default=405`.

- seed:

  Seed of the random number generator (integer), passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) so that the Monte
  Carlo results can be reproduced. `Default=1`.

- block:

  Number of random draws made per iteration of the vectorised search
  (integer, \> 0); larger values are faster but use more memory.
  `Default=2e5`.

- max_attempts:

  Maximum total number of draws before the search stops with an error
  because too few consistent draws were found (numeric, \> 0).
  `Default=5e7`.

- verbose:

  Print a message with the number of accepted draws and the acceptance
  rate (logical, TRUE or FALSE). `Default=TRUE`.
