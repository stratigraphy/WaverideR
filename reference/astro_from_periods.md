# Astronomical parameters from a tuned record (cycles already in kyr)

Same output as astro_from_peaks(), but the input is cycle periods in kyr
and the time scale is set by the tuning target rather than by a sampled
sedimentation rate. Give only the terms you actually resolved.

## Usage

``` r
astro_from_periods(
  cycles,
  age,
  anchor = 405,
  anchor_term = "e405",
  nsim = 2000,
  ptest = c(10000, 28000),
  method = c("angular_momentum", "astrogeo22"),
  s3_uncertainty = c("published", "scaled"),
  seed = NULL,
  verbose = TRUE,
  data_dir = NULL,
  ...
)
```

## Arguments

- cycles:

  Measured period ranges of the cycles you resolved in the tuned record,
  in kyr (named list; each element a numeric vector `c(min, max)`).
  Names are astronomical terms as accepted by `constraints` in
  [`astro_mc`](https://stratigraphy.github.io/WaverideR/reference/astro_mc.md)
  (e.g. `"e405"`, `"e125"`, `"e95"`, `"k+s3"` or `"obl"`, `"k+g5"`,
  `"k+g2"`, `"prec_short"`, `"k+g4|k+g3"`, `"any:k+g5|k+g2"`). Passed to
  [`astro_mc`](https://stratigraphy.github.io/WaverideR/reference/astro_mc.md)
  as `constraints`. Required, no default.

- age:

  Approximate age of the record in Ma (numeric): `c(min, max)`, or a
  single value that is widened by `age_tol` (default 5 Ma, see
  [`astro_mc`](https://stratigraphy.github.io/WaverideR/reference/astro_mc.md))
  on either side. Required, no default.

- anchor:

  Period in kyr assumed for the tuning target when the record was tuned
  (numeric, \> 0); every measured period is stretched by (sampled period
  of `anchor_term`)/`anchor`. `NULL` lets the time scale stretch freely
  (see `scale_range` in
  [`astro_mc`](https://stratigraphy.github.io/WaverideR/reference/astro_mc.md)).
  `Default=405`.

- anchor_term:

  Astronomical term the record was tuned to (character), any term name
  accepted by
  [`astro_mc`](https://stratigraphy.github.io/WaverideR/reference/astro_mc.md),
  e.g. `"e405"` (g2-g5). `Default="e405"`.

- nsim:

  Number of accepted Monte Carlo samples to collect (integer, \> 0);
  passed to
  [`astro_mc`](https://stratigraphy.github.io/WaverideR/reference/astro_mc.md).
  `Default=2000`.

- ptest:

  Range of precession constants k to sample, in years (numeric vector
  `c(min, max)`); passed to
  [`astro_mc`](https://stratigraphy.github.io/WaverideR/reference/astro_mc.md).
  `Default=c(10000, 28000)`.

- method:

  How the sampled precession constants are converted into an Earth-Moon
  state (distance, length of day) by
  [`ctest_earthmoon`](https://stratigraphy.github.io/WaverideR/reference/ctest_earthmoon.md)
  (character): `"angular_momentum"` (Earth-Moon angular-momentum
  conservation) or `"astrogeo22"` (inversion of the AstroGeo22 curve).
  `Default=c("angular_momentum", "astrogeo22")`, i.e.
  `"angular_momentum"`.

- s3_uncertainty:

  How the uncertainty of s3 is sampled, passed to
  [`astro_mc`](https://stratigraphy.github.io/WaverideR/reference/astro_mc.md)
  (character): `"published"` follows the published CTest code, in which
  the relative uncertainty is applied as an absolute one in years (s3 is
  effectively fixed); `"scaled"` applies it relative to the s3 period,
  as for g2-g5, giving s3 a realistic spread.
  `Default=c("published", "scaled")`, i.e. `"published"`.

- seed:

  Seed of the random number generator (integer), passed to
  [`astro_mc`](https://stratigraphy.github.io/WaverideR/reference/astro_mc.md)
  so that the Monte Carlo results can be reproduced; `NULL` leaves the
  current random state unchanged. `Default=NULL`.

- verbose:

  Print progress messages of the Monte Carlo sampling (logical, TRUE or
  FALSE). `Default=TRUE`.

- data_dir:

  Folder that contains the AstroGeo22 table `AstroGeo22.txt` (character
  path); `NULL` uses the copy installed with the package. Passed to
  [`ctest_earthmoon`](https://stratigraphy.github.io/WaverideR/reference/ctest_earthmoon.md)
  and
  [`astrogeo22`](https://stratigraphy.github.io/WaverideR/reference/astrogeo22.md).
  `Default=NULL`.

- ...:

  Further arguments passed to
  [`astro_mc`](https://stratigraphy.github.io/WaverideR/reference/astro_mc.md),
  e.g. `age_tol`, `sampleorb`, `scale_range`, `block` or `max_attempts`.

## Value

object of class \`astro_from_peaks\` (see astro_from_peaks())
