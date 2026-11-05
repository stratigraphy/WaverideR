# Precession constant from one or more single-valued cycles

Each cycle is read as a named astronomical term and inverted directly: k
= 360\*3600 / (f_observed - f_term). Uncertainty is propagated by Monte
Carlo over the cycle's own uncertainty (frequency resolution or \`tol\`)
and over the secular frequencies at that age.

## Usage

``` r
k_from_cycles(
  cycles,
  age_ma,
  record_length_kyr = NULL,
  tol = NULL,
  mult = 1,
  nsim = 20000,
  method = c("angular_momentum", "astrogeo22"),
  combine = c("multiply", "pool"),
  seed = NULL,
  tuned_to = NULL,
  anchor_term = "e405",
  sampleorb = "U",
  g4g3_range_kyr = NULL,
  s4s3_range_kyr = NULL,
  s3_uncertainty = c("published", "scaled"),
  data_dir = NULL
)
```

## Arguments

- cycles:

  Measured cycle periods in kyr (named numeric vector); every name is
  the astronomical term the cycle is read as, e.g.
  `c("k+s3" = 30.8, "k+g5" = 19.9)`. Terms with k (`"k+s3"`, `"k+g5"`,
  `"k+g2"`, `"k+g4"`, `"k+g3"`, ...) give the precession constant; terms
  without k (e.g. `"e405"`, `"g2-g5"`) only constrain the secular
  frequencies (element `secular` of the output).
  [`astro_terms`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md)
  lists the names that can be used. Required (no default).

- age_ma:

  Age of the record in Ma (numeric): one value, or `c(min, max)` to draw
  the age uniformly from that range in the Monte Carlo simulation.
  Required (no default).

- record_length_kyr:

  Length of the analysed record in kyr (numeric, \> 0); with it the
  uncertainty band of each cycle is set by the frequency resolution
  `mult / record_length_kyr` (see
  [`cycle_band`](https://stratigraphy.github.io/WaverideR/reference/cycle_band.md)).
  `NULL` uses the relative half-width `tol` instead. `Default=NULL`.

- tol:

  Relative half-width of the uncertainty band of each cycle, used only
  when `record_length_kyr` is `NULL` (numeric, \> 0; 0.05 is plus or
  minus 5%). `NULL` gives 0.05 (see
  [`cycle_band`](https://stratigraphy.github.io/WaverideR/reference/cycle_band.md)).
  `Default=NULL`.

- mult:

  Multiple of the frequency resolution `1 / record_length_kyr` used as
  half-width of the band (numeric, \> 0; 1 is plus or minus one
  bandwidth). `Default=1`.

- nsim:

  Number of Monte Carlo draws (integer, \>= 0); within its band each
  cycle is drawn uniformly. `0` gives the plain point estimate without
  uncertainty. `Default=20000`.

- method:

  How the precession constant is converted to Earth-Moon distance,
  length of day and obliquity (character, passed to
  [`ctest_earthmoon`](https://stratigraphy.github.io/WaverideR/reference/ctest_earthmoon.md)
  and
  [`earthmoon_state`](https://stratigraphy.github.io/WaverideR/reference/earthmoon_state.md)):
  `"angular_momentum"` (conservation of Earth-Moon angular momentum) or
  `"astrogeo22"` (inversion of the AstroGeo22 curve of Farhat et al.,
  2022). `Default=c("angular_momentum", "astrogeo22")`, i.e.
  `"angular_momentum"`.

- combine:

  How the estimates of several k terms are combined (character):
  `"multiply"` (product of the individual densities, as in Zeeden et
  al., 2023, via
  [`ctest_combine`](https://stratigraphy.github.io/WaverideR/reference/ctest_combine.md))
  or `"pool"` (all samples together). `Default=c("multiply", "pool")`,
  i.e. `"multiply"`.

- seed:

  Seed of the random number generator (integer), so that the Monte Carlo
  results can be reproduced; `NULL` leaves the current random state
  unchanged. `Default=NULL`.

- tuned_to:

  If the record was tuned, the period in kyr that was assumed for the
  tuning target (numeric, e.g. 405). The true duration of that target at
  this age is uncertain, so every measured period is stretched by (true
  target period)/`tuned_to` to propagate that error. `NULL` treats the
  input periods as true durations. `Default=NULL`.

- anchor_term:

  Astronomical term of the tuning target (character; a secular
  difference without k, e.g. `"e405"` or `"g2-g5"`); only used when
  `tuned_to` is given. `Default="e405"`.

- sampleorb:

  Shape of the distribution from which the uncertain secular frequencies
  g2, g3, g4, g5 and s3 are drawn (character): `"U"` uniform or `"N"`
  normal. `Default="U"`.

- g4g3_range_kyr:

  Range `c(min, max)` of the period of the g4-g3 term in kyr (numeric
  vector of length 2), drawn uniformly to allow for a chaotic resonance
  transition (g3 and g4 are moved symmetrically); `NULL` keeps g4-g3 at
  its La2004 value. `Default=NULL`.

- s4s3_range_kyr:

  Range `c(min, max)` of the period of the s4-s3 term in kyr (numeric
  vector of length 2), drawn uniformly as for `g4g3_range_kyr`; `NULL`
  keeps s4-s3 at its La2004 value. `Default=NULL`.

- s3_uncertainty:

  How the uncertainty of s3 is applied (character): `"published"` uses
  the relative uncertainty directly as standard deviation of the s3
  period in years, as in the published code; `"scaled"` multiplies it by
  the s3 period (a much wider spread).
  `Default=c("published", "scaled")`, i.e. `"published"`.

- data_dir:

  Folder that contains the data tables `AstroGeo22.txt` and
  `La04gisi.csv` (character path). `NULL` uses `getOption("cocoR.data")`
  or the copies installed with WaverideR. `Default=NULL`.

## Value

list: \`per_term\` (one row per cycle), \`combined\` (k and Earth-Moon
state), \`samples\`, and \`consistent\` (do the per-term 95% intervals
overlap)
