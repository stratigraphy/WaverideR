# Lunar mean motion and months per year from the Earth-Moon distance (Kepler's third law, n^2 a^3 = G(M_E + M_M), as in Farhat et al. 2022 Table D.2)

Lunar mean motion and months per year from the Earth-Moon distance
(Kepler's third law, n^2 a^3 = G(M_E + M_M), as in Farhat et al. 2022
Table D.2)

## Usage

``` r
months_per_year(a_RE, kind = c("sidereal", "synodic"))
```

## Arguments

- a_RE:

  Earth-Moon distance (semi-major axis of the lunar orbit) in Earth
  radii (numeric vector, \> 0; the present-day value is about 60.14).
  Required (no default).

- kind:

  Type of lunar month that is counted per year (character): `"sidereal"`
  (relative to the stars) or `"synodic"` (new moon to new moon, one
  fewer per year). `Default=c("sidereal", "synodic")`, i.e.
  `"sidereal"`.
