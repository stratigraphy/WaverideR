# Astronomical parameters from spectral peaks

Wraps ctest() and ctest_earthmoon(): give the wavelength ranges (m) of
the cycles you identified in a depth-domain spectrum and an approximate
age, and get back the secular frequencies g2-g5 and s3, the precession
constant k, the Earth-Moon distance and the length of day, with
uncertainties.

## Usage

``` r
astro_from_peaks(
  age,
  e405,
  e125,
  e95,
  obl = NULL,
  prec_long = NULL,
  prec_mid = NULL,
  prec_short = NULL,
  age_tol = 5,
  nsim = 2000,
  ptest = c(10000, 28000),
  mode = c("auto", "PT", "T", "P"),
  sampleorb = "U",
  method = c("angular_momentum", "astrogeo22"),
  s3_uncertainty = c("published", "scaled"),
  seed = NULL,
  verbose = TRUE,
  data_dir = NULL
)
```

## Arguments

- age:

  Approximate age of the record in Ma (numeric): `c(min, max)`, or a
  single value, which is widened by `age_tol` Ma on either side.
  Required, no default.

- e405, e125, e95:

  Wavelength ranges in metres of the 405-kyr, ~125-kyr and ~95-kyr
  eccentricity cycles identified in the depth-domain spectrum (each a
  numeric vector `c(min, max)`). All three are required, no default.

- obl:

  Wavelength range in metres of the obliquity cycle (k+s3), if
  identified (numeric vector `c(min, max)`); `NULL` if not identified.
  `Default=NULL`.

- prec_long, prec_mid, prec_short:

  Wavelength ranges in metres of the precession cycles, from longest to
  shortest (k+g5; k+g2; k+g4/k+g3), each a numeric vector `c(min, max)`
  or `NULL`. All three must be given for the precession test to be used.
  `Default=NULL`.

- age_tol:

  Half-width in Ma of the age range used when `age` is a single value
  (numeric, \>= 0). `Default=5`.

- nsim:

  Number of accepted Monte Carlo samples to collect (integer, \> 0);
  passed to
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md).
  `Default=2000`.

- ptest:

  Range of precession constants k to sample, in years (numeric vector
  `c(min, max)`); passed to
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md).
  `Default=c(10000, 28000)`.

- mode:

  Which cycles are used in the consistency test (character): `"auto"`
  picks the mode from the cycles supplied, `"PT"` precession and
  obliquity (tilt), `"T"` obliquity only, `"P"` precession only. Passed
  to
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md).
  `Default=c("auto", "PT", "T", "P")`, i.e. `"auto"`.

- sampleorb:

  How the secular frequencies g2-g5 are sampled, passed to
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md)
  (character): `"U"` uniform within 2 sigma (as in the paper) or `"N"`
  normal. `Default="U"`.

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
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md)
  (character): `"published"` follows the published CTest code, in which
  the relative uncertainty is applied as an absolute one in years (s3 is
  effectively fixed; use it to reproduce the paper); `"scaled"` applies
  it relative to the s3 period, as for g2-g5, giving s3 a realistic
  spread. `Default=c("published", "scaled")`, i.e. `"published"`.

- seed:

  Seed of the random number generator (integer), passed to
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md)
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

## Value

object of class \`astro_from_peaks\`: \`\$frequencies\` (g, s, k in both
arcsec/yr and yr), \`\$earthmoon\` (distance, day length, months, days
per year), \`\$periods\` (the astronomical periods implied),
\`\$reference\` (the AstroGeo22 model at the same age) and the
underlying \`\$ctest\`
