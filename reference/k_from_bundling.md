# Precession constant from bundling: precession cycles per eccentricity bundle

The scale-free way an eccentricity cycle does carry k. Count how many
precession couplets make up one eccentricity bundle (or measure the
ratio of the two periods, or the sideband spacing of the precession peak
relative to its carrier) and invert p = n \* (g_i - g_j) - f_term. No
tuning, no sedimentation rate and no time scale are needed.

## Usage

``` r
k_from_bundling(
  n,
  bundle_term = "e405",
  prec_term = "k+g5",
  age_ma,
  n_sd = 0.5,
  nsim = 20000,
  method = c("angular_momentum", "astrogeo22"),
  sampleorb = "U",
  g4g3_range_kyr = NULL,
  s4s3_range_kyr = NULL,
  seed = NULL,
  data_dir = NULL
)
```

## Arguments

- n:

  Number of precession (or obliquity) cycles counted per eccentricity
  bundle, or the measured ratio of the two periods (numeric, \> 1, may
  be fractional). Required (no default).

- bundle_term:

  Eccentricity term that forms the bundle (character): a named term such
  as `"e405"` (g2-g5), `"e125"` or `"e130"` (g3-g2), `"e124"` (g4-g2),
  `"e99"` (g3-g5), `"e95"` (g4-g5), any other name listed by
  [`astro_terms`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md),
  or an expression of secular frequencies such as `"g2-g5"`.
  `Default="e405"`.

- prec_term:

  Precession or obliquity term that is bundled (character; it must
  contain k), e.g. `"k+g5"` (main precession), `"k+g2"`, `"k+g4"`,
  `"k+s3"` (obliquity), or the aliases `"prec"`/`"precession"` and
  `"obl"`/`"obliquity"`. `Default="k+g5"`.

- age_ma:

  Age of the record in Ma (numeric): one value, or `c(min, max)` to draw
  the age uniformly from that range in the Monte Carlo simulation. It
  sets the uncertainty of the secular frequencies (growing with the
  square root of age) and the obliquity used by
  [`ctest_earthmoon`](https://stratigraphy.github.io/WaverideR/reference/ctest_earthmoon.md).
  Required (no default).

- n_sd:

  One-sigma counting uncertainty of `n`, in cycles (numeric, \>= 0); `n`
  is drawn from a normal distribution with this standard deviation.
  `Default=0.5`.

- nsim:

  Number of Monte Carlo draws (integer, \>= 0); `0` returns the plain
  point estimate without uncertainty. `Default=20000`.

- method:

  How the precession constant is converted to Earth-Moon distance,
  length of day and obliquity (character, passed to
  [`ctest_earthmoon`](https://stratigraphy.github.io/WaverideR/reference/ctest_earthmoon.md)):
  `"angular_momentum"` (conservation of Earth-Moon angular momentum) or
  `"astrogeo22"` (inversion of the AstroGeo22 curve of Farhat et al.,
  2022). `Default=c("angular_momentum", "astrogeo22")`, i.e.
  `"angular_momentum"`.

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
\`k_ci95_kyr\`, \`n\`, \`terms\`, \`samples\`, \`predicted_n_model\`,
\`reference\` and \`method\`
