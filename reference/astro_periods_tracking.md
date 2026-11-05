# Astronomical periods and weights for tracking, for any age

Returns the main eccentricity, obliquity and precession periods with
weights proportional to their typical amplitude, for use in
[`track_astro_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md).
The obliquity and precession periods follow from the precession constant
k at the given age (AstroGeo22, via
[`astro_terms`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md))
or from a k given directly (e.g. estimated with
[`test_tuned`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md)),
so that the periods fit deep-time records such as the Silurian and
Devonian. The eccentricity periods are the La2004 values; the 405 kyr
period is stable through the Phanerozoic, the short eccentricity periods
are less certain before about 60 Ma.

## Usage

``` r
astro_periods_tracking(
  age_ma = 0,
  k_yr = NULL,
  include = c("eccentricity", "obliquity", "precession")
)
```

## Arguments

- age_ma:

  Age of the record in Ma (numeric, \>= 0); the precession constant, and
  with it the obliquity and precession periods, are taken from the
  AstroGeo22 solution at this age (passed to
  [`astro_terms`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md)).
  Ignored for the precession constant when `k_yr` is given. `Default=0`
  (present-day periods).

- k_yr:

  Optional precession constant (numeric, single value) used instead of
  the AstroGeo22 value for `age_ma`, given as the period of the axial
  precession in years (passed to
  [`astro_terms`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md),
  which converts it to arcsec/yr as 360 \* 3600 / `k_yr`), e.g.
  estimated with
  [`test_tuned`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md);
  `NULL` uses AstroGeo22. `Default=NULL`.

- include:

  Which families of cycles to return (character vector, any of
  `"eccentricity"`, `"obliquity"`, `"precession"`). When only
  `"eccentricity"` is included, the 405-kyr cycle gets weight 1 (instead
  of 0.5): in records in which only the eccentricity is resolved, the
  stable 405-kyr cycle carries the age model.
  `Default=c("eccentricity", "obliquity", "precession")`.

## Value

A data frame with `term`, `family`, `period` (kyr) and `weight`.

## Examples

``` r
astro_periods_tracking(0)
#>     term       family    period weight
#> 1  g2-g5 eccentricity 405.69120    0.5
#> 2  g4-g2 eccentricity 123.85321    0.3
#> 3  g3-g2 eccentricity 130.69786    0.2
#> 4  g4-g5 eccentricity  94.88563    0.3
#> 5  g3-g5 eccentricity  98.85170    0.2
#> 6   k+s3    obliquity  40.98966    1.0
#> 7   k+s4    obliquity  39.61761    0.4
#> 8   k+s6    obliquity  53.73163    0.3
#> 9   k+g5   precession  23.68197    1.0
#> 10  k+g2   precession  22.37580    0.8
#> 11  k+g4   precession  18.95188    0.6
#> 12  k+g3   precession  19.10498    0.4
astro_periods_tracking(420)        # Silurian
#>     term       family    period weight
#> 1  g2-g5 eccentricity 405.69120    0.5
#> 2  g4-g2 eccentricity 123.85321    0.3
#> 3  g3-g2 eccentricity 130.69786    0.2
#> 4  g4-g5 eccentricity  94.88563    0.3
#> 5  g3-g5 eccentricity  98.85170    0.2
#> 6   k+s3    obliquity  30.49676    1.0
#> 7   k+s4    obliquity  29.73069    0.4
#> 8   k+s6    obliquity  37.03020    0.3
#> 9   k+g5   precession  19.75496    1.0
#> 10  k+g2   precession  18.83767    0.8
#> 11  k+g4   precession  16.35077    0.6
#> 12  k+g3   precession  16.46460    0.4
```
