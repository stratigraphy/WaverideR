# The menu of named astronomical terms at a given age

The menu of named astronomical terms at a given age

## Usage

``` r
astro_terms(age_ma = 0, k_yr = NULL, data_dir = NULL)
```

## Arguments

- age_ma:

  Age in Ma (numeric, single value, 0 to 4429) at which the periods are
  computed; the precession constant used for the k-bearing (obliquity
  and precession) terms is taken from
  [`astrogeo22`](https://stratigraphy.github.io/WaverideR/reference/astrogeo22.md)
  at this age unless `k_yr` is given. `Default=0`.

- k_yr:

  Your own precession constant k, as a period in years (numeric, \> 0),
  used instead of the AstroGeo22 value at `age_ma`; `NULL` uses the
  AstroGeo22 value. `Default=NULL`.

- data_dir:

  Folder that contains the data files `AstroGeo22.txt` and
  `La04gisi.csv` (character path); `NULL` uses the copy installed with
  the package. `Default=NULL`.

## Value

data.frame from the longest to the shortest period
