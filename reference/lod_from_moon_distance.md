# Length of day implied by an Earth-Moon distance (angular momentum conserved)

Length of day implied by an Earth-Moon distance (angular momentum
conserved)

## Usage

``` r
lod_from_moon_distance(a_RE, kind = "solar", L_total = NULL)
```

## Arguments

- a_RE:

  Earth-Moon distance (semi-major axis of the lunar orbit) in Earth
  radii (numeric vector, \> 0; the present-day value is about 60.14).
  Required (no default).

- kind:

  Type of day that is returned (character, passed to
  [`lod_hours`](https://stratigraphy.github.io/WaverideR/reference/lod_hours.md)):
  `"solar"` (mean solar day) or `"sidereal"` (one rotation relative to
  the stars). `Default="solar"`.

- L_total:

  Total angular momentum of the Earth-Moon system in kg m^2/s (numeric,
  \> 0), assumed to be conserved. `NULL` uses the present-day value
  (computed from today's Earth-Moon distance and length of day).
  `Default=NULL`.
