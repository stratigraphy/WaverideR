# Monte Carlo consistency test for tuned (time-domain) or depth-domain cycles

Monte Carlo consistency test for tuned (time-domain) or depth-domain
cycles

## Usage

``` r
astro_mc(
  constraints,
  age,
  domain = c("time", "depth"),
  anchor = 405,
  anchor_term = "e405",
  scale_range = NULL,
  g4g3_range_kyr = NULL,
  s4s3_range_kyr = NULL,
  age_tol = 5,
  nsim = 2000,
  ptest = c(10000, 28000),
  sampleorb = c("U", "N"),
  s3_uncertainty = c("published", "scaled"),
  seed = NULL,
  block = 1e+05,
  max_attempts = 5e+07,
  verbose = TRUE
)
```

## Arguments

- constraints:

  Measured ranges of the cycles you resolved (named list; each element a
  numeric vector `c(min, max)`). Names are astronomical terms: `"e405"`,
  `"e125"` (= e130 and e124), `"e95"` (= e99 and e95), `"k+s3"` (or
  `"obl"`, `"obliquity"`), `"k+g5"`, `"k+g2"`, `"k+g4"`, `"k+g3"`,
  `"prec_short"` (= k+g4 and k+g3), `"prec_all"` (all four precession
  terms), any other term listed by
  [`astro_terms`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md)
  (e.g. `"g4-g3"`, `"s4-s3"`), or a combination such as `"2(g2-g5)"`.
  Several terms in one band: `"k+g4|k+g3"` (all must fit). A peak you
  cannot assign: `"any:k+g5|k+g2"` (at least one must fit). Units: kyr
  when `domain = "time"`, metres when `domain = "depth"`. Required, no
  default.

- age:

  Approximate age of the record in Ma (numeric): `c(min, max)`, or a
  single value, which is widened by `age_tol` Ma on either side.
  Required, no default.

- domain:

  Domain of the measured cycles (character): `"time"` for a tuned record
  (cycles in kyr), `"depth"` for a depth-domain spectrum (cycles in
  metres; `scale_range` must then be given).
  `Default=c("time", "depth")`, i.e. `"time"`.

- anchor:

  Period in kyr that was assumed for the tuning target `anchor_term`
  when the record was tuned (numeric, \> 0); with `domain = "time"`
  every predicted period is divided by (sampled period of
  `anchor_term`)/`anchor`. `NULL` means that the time scale is sampled
  freely within `scale_range`. Ignored when `domain = "depth"`.
  `Default=405`.

- anchor_term:

  Astronomical term the record was tuned to (character, any name
  accepted in `constraints`). `Default="e405"` (= g2-g5).

- scale_range:

  Range sampled uniformly for the scale factor (numeric vector
  `c(min, max)`): the sedimentation rate in m/kyr when
  `domain = "depth"` (then required), or the allowed stretch of the time
  scale when `domain = "time"` and `anchor = NULL` (then `NULL` becomes
  `c(0.9, 1.1)`). `Default=NULL`.

- g4g3_range_kyr, s4s3_range_kyr:

  Not used by the calculation; accepted only for a common interface with
  other functions (numeric vector `c(min, max)` in kyr or `NULL`).
  `Default=NULL`.

- age_tol:

  Half-width in Ma of the age range used when `age` is a single value
  (numeric, \>= 0). `Default=5`.

- nsim:

  Number of accepted (consistent) Monte Carlo samples to collect
  (integer, \> 0). `Default=2000`.

- ptest:

  Range of precession constants k to sample uniformly, in years (numeric
  vector `c(min, max)`). `Default=c(10000, 28000)`.

- sampleorb:

  How the secular frequencies g2-g5 are sampled (character): `"U"`
  uniform within 2 sigma (as in Zeeden et al. 2023) or `"N"` normal with
  1 sigma. Sigma is the Laskar et al. (2011) uncertainty scaled by
  sqrt(age/100 Myr). `Default=c("U", "N")`, i.e. `"U"`.

- s3_uncertainty:

  How the uncertainty of s3 is sampled (only when a constraint involves
  k+s3; character): `"published"` follows the published CTest code, in
  which the relative uncertainty is applied as an absolute one in years
  (s3 is effectively fixed); `"scaled"` applies it relative to the s3
  period, as for g2-g5, giving s3 a realistic spread.
  `Default=c("published", "scaled")`, i.e. `"published"`.

- seed:

  Seed of the random number generator (integer), so that the Monte Carlo
  results can be reproduced; `NULL` leaves the current random state
  unchanged. `Default=NULL`.

- block:

  Number of random draws per iteration of the vectorised search
  (integer; values below 1000 are raised to 1000). Larger blocks are
  faster but use more memory. `Default=1e5`.

- max_attempts:

  Maximum total number of draws before the search stops with an error
  because fewer than `nsim` consistent samples were found (numeric, \>
  0). `Default=5e7`.

- verbose:

  Print a progress message after each block with accepted samples
  (logical, TRUE or FALSE). `Default=TRUE`.

## Value

object of class c("astro_mc", "ctest"), so print(), plot() and
ctest_earthmoon() work on it
