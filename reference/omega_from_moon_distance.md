# Earth's spin rate implied by an Earth-Moon distance, assuming the Earth-Moon angular momentum is conserved (solar tidal drain neglected).

Earth's spin rate implied by an Earth-Moon distance, assuming the
Earth-Moon angular momentum is conserved (solar tidal drain neglected).

## Usage

``` r
omega_from_moon_distance(a_RE, L_total = NULL)
```

## Arguments

- a_RE:

  Earth-Moon distance (semi-major axis of the lunar orbit) in Earth
  radii (numeric vector, \> 0; the present-day value is about 60.14).
  Required (no default).

- L_total:

  Total angular momentum of the Earth-Moon system in kg m^2/s (numeric,
  \> 0), assumed to be conserved. `NULL` uses the present-day value
  (computed from today's Earth-Moon distance and length of day).
  `Default=NULL`.
