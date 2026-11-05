# Invert the AstroGeo22 precession curve: Earth-Moon state for a measured precession frequency / period.

The model curve is not strictly monotonic (tidal resonances), so every
crossing is returned.

## Usage

``` r
earthmoon_from_precession(
  p_kyr = NULL,
  p_arcsec = NULL,
  method = c("astrogeo22", "angular_momentum"),
  obliquity_deg = NULL,
  age_range = NULL,
  data_dir = NULL
)
```

## Arguments

- p_kyr:

  Period of the axial precession, in kyr (numeric, scalar or vector, \>
  0; about 25.7 kyr today), i.e. `360 * 3600 / (p_kyr * 1000)`
  arcsec/yr. Ignored when `p_arcsec` is given; one of the two is
  required. `Default=NULL`.

- p_arcsec:

  Precession constant (rate of axial precession), in arcsec/yr (numeric,
  scalar or vector, \> 0; about 50.47 today). Takes precedence over
  `p_kyr`. `Default=NULL`.

- method:

  How the Earth-Moon state is found (character): `"astrogeo22"` finds
  every age at which the AstroGeo22 precession curve crosses the
  requested value (the curve is not monotonic, so several ages may be
  returned) and returns the AstroGeo22 state there; `"angular_momentum"`
  solves for the Earth-Moon distance (25 to 80 Earth radii) that gives
  the requested precession with
  [`precession_frequency`](https://stratigraphy.github.io/WaverideR/reference/precession_frequency.md),
  assuming conservation of the Earth-Moon angular momentum.
  `Default=c("astrogeo22", "angular_momentum")`, i.e. `"astrogeo22"`.

- obliquity_deg:

  Obliquity of the Earth, in degrees (numeric), used by
  `method="angular_momentum"`; `NULL` uses the present-day value
  23.2545. `Default=NULL`.

- age_range:

  Age range to which the search along the AstroGeo22 curve is
  restricted, in Ma (numeric vector `c(min, max)`, or `NULL` for the
  whole curve). Only used with `method="astrogeo22"`. `Default=NULL`.

- data_dir:

  Directory that contains the file `AstroGeo22.txt` (character, or
  `NULL`); `NULL` uses the data installed with the package (or
  `options(cocoR.data = ...)`). `Default=NULL`.
