# Obliquity at a given age, from AstroGeo22 if available, else today's value

Obliquity at a given age, from AstroGeo22 if available, else today's
value

## Usage

``` r
astro_obliquity(age_ma)
```

## Arguments

- age_ma:

  Age in Ma at which the obliquity is wanted (numeric, single value; if
  a vector is given, only its first value is used). The nearest row of
  the AstroGeo22 table installed with WaverideR is taken (if that copy
  is missing it is downloaded to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html)); if no table is
  available today's obliquity (23.2545 degrees) is returned. Required,
  no default.
