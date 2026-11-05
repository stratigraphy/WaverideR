# Length of day from Earth's angular velocity

Length of day from Earth's angular velocity

## Usage

``` r
lod_hours(omega, kind = c("solar", "sidereal"))
```

## Arguments

- omega:

  Angular velocity of the rotation of the Earth in rad/s (numeric
  vector, \> 0; the present-day value is about 7.29e-5), e.g. the output
  of
  [`omega_from_moon_distance`](https://stratigraphy.github.io/WaverideR/reference/omega_from_moon_distance.md).
  Required (no default).

- kind:

  Type of day that is returned, in hours (character): `"solar"` (mean
  solar day, corrected for the orbital motion of the Earth around the
  Sun) or `"sidereal"` (one rotation relative to the stars).
  `Default=c("solar", "sidereal")`, i.e. `"solar"`.
