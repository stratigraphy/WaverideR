# Earth-Moon distance from the number of lunar months per year (the observable of tidal rhythmites)

Earth-Moon distance from the number of lunar months per year (the
observable of tidal rhythmites)

## Usage

``` r
moon_distance_from_months(months, kind = c("sidereal", "synodic"))
```

## Arguments

- months:

  Number of lunar months per year (numeric vector, \> 0; about 13.37
  sidereal or 12.37 synodic months today), e.g. counted in tidal
  rhythmites. Required (no default).

- kind:

  Type of lunar month counted in `months` (character): `"sidereal"`
  (relative to the stars) or `"synodic"` (new moon to new moon;
  converted with sidereal = synodic + 1).
  `Default=c("sidereal", "synodic")`, i.e. `"sidereal"`.
