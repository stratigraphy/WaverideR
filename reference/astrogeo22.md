# Earth-Moon state from the AstroGeo22 model at a given age

Earth-Moon state from the AstroGeo22 model at a given age

## Usage

``` r
astrogeo22(age_ma, interpolate = TRUE, data_dir = NULL)
```

## Arguments

- age_ma:

  Age(s) in Ma at which the Earth-Moon state is returned (numeric
  vector, 0 to 4429; one output row per value). Ages outside the table
  give `NA` when `interpolate=TRUE`. Required, no default.

- interpolate:

  How values are taken from the 1-Myr rows of the AstroGeo22 table
  (logical): `TRUE` interpolates linearly between rows, `FALSE` takes
  the nearest row (as calculate_orbit9.m of Acycle does).
  `Default=TRUE`.

- data_dir:

  Folder that contains the file `AstroGeo22.txt` (character path);
  `NULL` uses the copy installed with the package (or the folder set
  with `options(cocoR.data = ...)`). `Default=NULL`.

## Value

data.frame: age_ma, a_RE, a_km, lod_h, obliquity_deg, p_arcsec, p_kyr
and the model's min/max envelope for each
