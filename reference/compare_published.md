# Compare a measured precession constant with published models

Compare a measured precession constant with published models

## Usage

``` r
compare_published(age_ma, k_yr = NULL, obliquity_deg = NULL, data_dir = NULL)
```

## Arguments

- age_ma:

  Age of the record in Ma (single numeric value, 0 to 4429), at which
  the published models (AstroGeo22 and the La2004 extrapolation) are
  evaluated. Required, no default.

- k_yr:

  Measured precession constant of the record, given as its period in
  years (numeric, e.g. about 25700 for today); when given, rows are
  added with the Earth-Moon state inverted from it (by angular momentum
  and on the AstroGeo22 curve, see
  [`earthmoon_from_precession`](https://stratigraphy.github.io/WaverideR/reference/earthmoon_from_precession.md)).
  `NULL` gives only the published models. `Default=NULL`.

- obliquity_deg:

  Obliquity in degrees (numeric) used for the angular-momentum inversion
  of `k_yr`; only used when `k_yr` is given. `NULL` uses the AstroGeo22
  obliquity at `age_ma`. `Default=NULL`.

- data_dir:

  Folder that contains the astronomical data tables `AstroGeo22.txt` and
  `La04gisi.csv` (character path); `NULL` uses the tables installed with
  the package (or the folder set with `options(cocoR.data = ...)`).
  `Default=NULL`.

## Value

data.frame with one row per source: precession period and frequency,
Earth-Moon distance, LOD, obliquity, plus the k+g and k+s periods
