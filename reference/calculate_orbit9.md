# Nine astronomical periods at a given age

Nine principal orbital periods (calculate_orbit9.m in Acycle): La2004
secular frequencies g_i, s_i and the AstroGeo22 precession constant k(t)
at the age rounded to the nearest Ma.

## Usage

``` r
calculate_orbit9(age_ma, data_dir = NULL)
```

## Arguments

- age_ma:

  Age at which the nine periods are calculated, in Ma (single finite
  numeric value, 0 to 4429); it is rounded to the nearest whole Ma
  because AstroGeo22 is tabulated in 1-Myr steps. Required, no default.

- data_dir:

  Folder that contains the astronomical data tables `AstroGeo22.txt` and
  `La04gisi.csv` (character path); `NULL` uses the tables installed with
  the package (or the folder set with `options(cocoR.data = ...)`).
  `Default=NULL`.

## Value

data.frame with columns term, group, freq_arcsec_yr, period_yr,
period_kyr.
