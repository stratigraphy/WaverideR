# Precession constant from the RATIO of two cycles in one section

Scale-free: the ratio of two cycle periods equals the ratio of the two
wavelengths, so this works on an untuned depth series just as well as on
a tuned record, and no sedimentation rate is needed. Solving T1/T2 =
(p + f2)/(p + f1) for the precession frequency p gives p = (f2 - R f1) /
(R - 1), R = T1/T2, T1 \> T2.

## Usage

``` r
k_from_ratio(
  long,
  short,
  term_long = "k+s3",
  term_short = "k+g5",
  age_ma,
  rel_sd = 0.03,
  record_length = NULL,
  mult = 1,
  nsim = 20000,
  method = c("angular_momentum", "astrogeo22"),
  seed = NULL,
  data_dir = NULL
)
```

## Arguments

- long:

  Period of the longer of the two cycles (numeric, \> `short`), in kyr
  for a tuned record or in metres for a depth series; only the ratio
  `long / short` is used. Required (no default).

- short:

  Period of the shorter cycle, in the same unit as `long` (numeric, \>
  0). Required (no default).

- term_long:

  Astronomical term of the longer cycle (character), e.g. `"k+s3"`
  (obliquity), `"k+g5"` (precession) or an eccentricity term such as
  `"e405"`, `"e125"`, `"e95"` (see
  [`astro_terms`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md)).
  At least one of the two terms must contain k. `Default="k+s3"`.

- term_short:

  Astronomical term of the shorter cycle (character), with the same
  options as `term_long`. `Default="k+g5"`.

- age_ma:

  Age of the record in Ma (numeric, one value or `c(min, max)`, of which
  the mean is used); sets the obliquity used by
  [`ctest_earthmoon`](https://stratigraphy.github.io/WaverideR/reference/ctest_earthmoon.md)
  and the AstroGeo22 reference values. Required (no default).

- rel_sd:

  Relative one-sigma uncertainty of each measured cycle (numeric, \>= 0;
  0.03 is 3%); the cycles are drawn from normal distributions. Ignored
  when `record_length` is given. `Default=0.03`.

- record_length:

  Length of the analysed record in the same unit as the cycles (numeric,
  \> 0). If given, each cycle is drawn uniformly from the band set by
  the frequency resolution `mult / record_length` (see
  [`cycle_band`](https://stratigraphy.github.io/WaverideR/reference/cycle_band.md))
  instead of using `rel_sd`. `Default=NULL`.

- mult:

  Multiple of the frequency resolution used as half-width of the band
  when `record_length` is given (numeric, \> 0). `Default=1`.

- nsim:

  Number of Monte Carlo draws (integer, \>= 0); `0` skips the
  uncertainty and returns the plain point estimate. `Default=20000`.

- method:

  How the precession constant is converted to Earth-Moon distance,
  length of day and obliquity (character, passed to
  [`ctest_earthmoon`](https://stratigraphy.github.io/WaverideR/reference/ctest_earthmoon.md)):
  `"angular_momentum"` (conservation of Earth-Moon angular momentum) or
  `"astrogeo22"` (inversion of the AstroGeo22 curve of Farhat et al.,
  2022). `Default=c("angular_momentum", "astrogeo22")`, i.e.
  `"angular_momentum"`.

- seed:

  Seed of the random number generator (integer), so that the Monte Carlo
  results can be reproduced; `NULL` leaves the current random state
  unchanged. `Default=NULL`.

- data_dir:

  Folder that contains the data tables `AstroGeo22.txt` and
  `La04gisi.csv` (character path). `NULL` uses `getOption("cocoR.data")`
  or the copies installed with WaverideR. `Default=NULL`.

## Value

list with the Earth-Moon \`summary\`, \`k_kyr\`, \`k_sd_kyr\`,
\`k_ci95_kyr\`, \`ratio\`, \`samples\`, \`terms\`, \`reference\` and
\`method\`

## Details

Either cycle may be an eccentricity term ("e405", "e125", "e95", ...).
Pairing an eccentricity cycle with a precession or obliquity cycle is
often the most robust option in an untuned section, because the
eccentricity terms carry no k and are the best known frequencies.
