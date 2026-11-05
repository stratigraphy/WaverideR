# Precession constant implied by an observed obliquity or precession cycle

Precession constant implied by an observed obliquity or precession cycle

## Usage

``` r
cycles_to_k(period_kyr, term = "k+s3", g = NULL, s = NULL, data_dir = NULL)
```

## Arguments

- period_kyr:

  Period of the cycle read from the record, in kyr (numeric, \> 0; can
  be a vector of several readings of the same term), i.e. from a tuned
  record or a known sedimentation rate. Required, no default.

- term:

  Astronomical term the cycle represents (character): `"k+s3"`
  (obliquity) or one of the precession terms `"k+g5"`, `"k+g2"`,
  `"k+g4"`, `"k+g3"` (aliases such as `"obl"` and `"prec"` are accepted,
  see
  [`astro_terms`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md));
  a term without k (an eccentricity term) gives an error.
  `Default="k+s3"`.

- g:

  Secular frequencies g1 to g9 in arcsec/yr (numeric vector of length 9)
  used instead of the La2004 values; `NULL` uses the La2004 values from
  `La04gisi.csv`. `Default=NULL`.

- s:

  Secular frequencies s1 to s9 in arcsec/yr (numeric vector of length 9)
  used instead of the La2004 values; `NULL` uses the La2004 values from
  `La04gisi.csv`. `Default=NULL`.

- data_dir:

  Folder that contains the astronomical data tables `AstroGeo22.txt` and
  `La04gisi.csv` (character path); `NULL` uses the tables installed with
  the package (or the folder set with `options(cocoR.data = ...)`).
  `Default=NULL`.

## Value

list: p in arcsec/yr and the precession constant k in years
