# CTest: astronomical parameters from interpreted spectral peaks

CTest: astronomical parameters from interpreted spectral peaks

## Usage

``` r
ctest(
  data,
  nsim = 1000,
  mode = c("PT", "P", "T"),
  ptest = c(10000, 28000),
  sampleorb = c("U", "N"),
  seed = NULL,
  max_attempts = 5e+07,
  verbose = TRUE,
  engine = c("vectorized", "loop"),
  block = 1e+05,
  s3_uncertainty = c("published", "scaled")
)
```

## Arguments

- data:

  Interpreted spectral peaks as an 8 x 2 table (matrix, data frame or
  tibble), in the layout of the published spreadsheet, best built with
  [`ctest_input`](https://stratigraphy.github.io/WaverideR/reference/ctest_input.md):
  row 1 the age range of the record in Ma (min, max); rows 2 to 8 the
  frequency ranges in cycles per metre of the 405-kyr, ~125-kyr and
  ~95-kyr eccentricity cycles, the obliquity cycle (k+s3) and the three
  precession cycles (k+g5; k+g2; k+g4 and k+g3), with the lower
  frequency (longer wavelength) in column 1 and the higher frequency in
  column 2. The obliquity row may be `NA` when `mode="P"` and the
  precession rows when `mode="T"`. Required, no default.

- nsim:

  Number of accepted (successful) samples to collect (integer, \> 0).
  `Default=1000`.

- mode:

  Which cycles constrain the precession constant besides eccentricity
  (character): `"PT"` precession and obliquity (tilt), `"P"` precession
  only, `"T"` obliquity only. `Default=c("PT", "P", "T")`, i.e. `"PT"`.

- ptest:

  Range of the precession constant that is sampled (uniformly), given as
  its period in years (numeric vector `c(min, max)`).
  `Default=c(10000, 28000)`.

- sampleorb:

  Distribution from which the secular frequencies g2, g3, g4 and g5 are
  drawn around their Laskar et al. (2011) values, with an uncertainty
  that grows with the square root of age (character): `"U"` uniform
  within 2 sigma (used in the paper) or `"N"` normal.
  `Default=c("U", "N")`, i.e. `"U"`.

- seed:

  Seed of the random number generator (integer), so that the samples can
  be reproduced; `NULL` leaves the current random state unchanged.
  `Default=NULL`.

- max_attempts:

  Maximum number of draws (numeric); the function stops with an error if
  fewer than `nsim` samples have been accepted by then, which usually
  means the interpreted cycle ranges are inconsistent. `Default=5e7`.

- verbose:

  Print progress messages with the number of accepted samples and the
  acceptance rate (logical, TRUE or FALSE). `Default=TRUE`.

- engine:

  Sampling engine (character): `"vectorized"` draws in blocks of `block`
  samples (fast); `"loop"` draws one sample at a time in the published
  draw order and reproduces mmc3/mmc4/mmc5 sample for sample under the
  same seed. `Default=c("vectorized", "loop")`, i.e. `"vectorized"`.

- block:

  Number of draws per block for `engine="vectorized"` (integer; limited
  to 1000 to 5e6); larger blocks are faster but use more memory. Ignored
  by `engine="loop"`. `Default=1e5`.

- s3_uncertainty:

  How the uncertainty of the secular frequency s3 (used for obliquity)
  is applied (character): `"published"` keeps the draw of the published
  code, `rnorm(1, s3, ds3 * sqrt(t / 1e8))`, in which the relative
  uncertainty `ds3` is applied as an absolute one in years, so s3 is
  effectively fixed (use this to reproduce the paper); `"scaled"`
  applies it as for g2 to g5 (`s3 * ds3 * sqrt(t / 1e8)`), which gives
  s3 a realistic spread. `Default=c("published", "scaled")`, i.e.
  `"published"`.

## Value

list with \`summary\` (table as in the paper), \`samples\` (all accepted
samples), \`sampled_range\` of g/s, and \`acceptance\`
