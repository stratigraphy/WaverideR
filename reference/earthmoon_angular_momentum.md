# Earth-Moon system angular momentum (spin + lunar orbit), kg m^2 s^-1

Earth-Moon system angular momentum (spin + lunar orbit), kg m^2 s^-1

## Usage

``` r
earthmoon_angular_momentum(a_RE, omega)
```

## Arguments

- a_RE:

  Earth-Moon distance (semi-major axis of the lunar orbit), in Earth
  radii (numeric, scalar or vector, \> 0; present-day value about
  60.14). Required (no default).

- omega:

  Angular velocity (spin rate) of the Earth, in rad/s (numeric, \> 0,
  same length as `a_RE` or length 1), e.g. from
  [`earth_omega`](https://stratigraphy.github.io/WaverideR/reference/earth_omega.md).
  Required (no default).
