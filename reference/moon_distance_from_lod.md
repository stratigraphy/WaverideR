# Earth-Moon distance implied by a length of day (angular momentum conserved)

Earth-Moon distance implied by a length of day (angular momentum
conserved)

## Usage

``` r
moon_distance_from_lod(lod_h, kind = "solar", L_total = NULL)
```

## Arguments

- lod_h:

  Length of day in hours (numeric vector, \> 0; 24 today), of the type
  set by `kind`. Required (no default).

- kind:

  Type of day given in `lod_h` (character, passed to
  [`earth_omega`](https://stratigraphy.github.io/WaverideR/reference/earth_omega.md)):
  `"solar"` (mean solar day) or `"sidereal"` (one rotation relative to
  the stars). `Default="solar"`.

- L_total:

  Total angular momentum of the Earth-Moon system in kg m^2/s (numeric,
  \> 0), assumed to be conserved. `NULL` uses the present-day value
  (computed from today's Earth-Moon distance and length of day).
  `Default=NULL`.
