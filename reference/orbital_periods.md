# All principal astronomical periods for a given precession constant

All principal astronomical periods for a given precession constant

## Usage

``` r
orbital_periods(
  k_yr = NULL,
  p_arcsec = NULL,
  g = NULL,
  s = NULL,
  terms = c("nine", "all"),
  data_dir = NULL
)
```

## Arguments

- k_yr:

  Precession constant k expressed as the period of the axial precession,
  in years (numeric, \> 0), e.g. 20700 for the Devonian estimate of
  Zeeden et al. (2023); ignored when `p_arcsec` is given. One of `k_yr`
  and `p_arcsec` is required. `Default=NULL`.

- p_arcsec:

  Precession constant as a frequency in arcsec/yr (numeric, \> 0; about
  50.5 today); if given it is used instead of `k_yr`. `Default=NULL`.

- g:

  Secular frequencies g1-g9 of the planetary perihelia in arcsec/yr
  (numeric vector of length 9). `NULL` uses the La2004 values of
  `La04gisi.csv`. `Default=NULL`.

- s:

  Secular frequencies s1-s9 of the planetary nodes in arcsec/yr (numeric
  vector of length 9). `NULL` uses the La2004 values of `La04gisi.csv`.
  `Default=NULL`.

- terms:

  Which terms are returned (character): `"nine"` (the Acycle ORBIT9 set:
  g2-g5, g3-g2, g4-g2, g3-g5, g4-g5, k+s3, k+g5, k+g2, k+g4) or `"all"`
  (every gi-gj eccentricity term and every k+si and k+gi term for i = 1
  to 5). `Default=c("nine", "all")`, i.e. `"nine"`.

- data_dir:

  Folder that contains the data tables `AstroGeo22.txt` and
  `La04gisi.csv` (character path). `NULL` uses `getOption("cocoR.data")`
  or the copies installed with WaverideR. `Default=NULL`.

## Value

data.frame: term, group, frequency (arcsec/yr), period (kyr)
