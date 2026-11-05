# Shift the fundamental frequencies and see what every cycle does

Shift the fundamental frequencies and see what every cycle does

## Usage

``` r
chaos_scenario(
  g4_g3_kyr = NULL,
  s4_s3_kyr = NULL,
  hold = c("mean", "g3", "g4"),
  age_ma = 0,
  k_yr = NULL,
  data_dir = NULL
)
```

## Arguments

- g4_g3_kyr:

  Period of the g4-g3 ("2.4 Myr") eccentricity term in the scenario, in
  kyr (single numeric value \> 0, e.g. 1200 for a 2:1 resonance state);
  g3 and g4 are moved (as set by `hold`) so that g4-g3 has this period,
  and every term containing g3 or g4 changes with them. `NULL` keeps the
  La2004 values. `Default=NULL`.

- s4_s3_kyr:

  Period of the s4-s3 ("1.2 Myr") obliquity modulation term in the
  scenario, in kyr (single numeric value \> 0); s3 and s4 are moved
  symmetrically about their mean to reach it. `NULL` keeps the La2004
  values. `Default=NULL`.

- hold:

  How g3 and g4 are moved to reach `g4_g3_kyr` (character): `"mean"`
  moves both symmetrically about their mean, `"g3"` keeps g3 fixed and
  moves g4 only, `"g4"` keeps g4 fixed and moves g3 only.
  `Default=c("mean", "g3", "g4")`, i.e. `"mean"`.

- age_ma:

  Age in Ma (single numeric value, 0 to 4429) at which the precession
  constant for the precession and obliquity (k-bearing) terms is taken
  from AstroGeo22
  ([`astrogeo22`](https://stratigraphy.github.io/WaverideR/reference/astrogeo22.md)).
  Ignored for that purpose when `k_yr` is given. `Default=0`.

- k_yr:

  Your own value of the precession constant, given as its period in
  years (numeric, e.g. about 25700 for today), used instead of the
  AstroGeo22 value; `NULL` uses the model value at `age_ma`.
  `Default=NULL`.

- data_dir:

  Folder that contains the astronomical data tables `AstroGeo22.txt` and
  `La04gisi.csv` (character path); `NULL` uses the tables installed with
  the package (or the folder set with `options(cocoR.data = ...)`).
  `Default=NULL`.

## Value

the astro_terms() table with the scenario period, the La2004 period, and
the change, so it is obvious which cycles move and which do not
