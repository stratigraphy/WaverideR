# Luni-solar precession frequency (Farhat et al. 2022, eq. A.7)

Luni-solar precession frequency (Farhat et al. 2022, eq. A.7)

## Usage

``` r
precession_frequency(
  a_RE,
  lod_h = NULL,
  obliquity_deg = .AC$obliquity0_deg,
  Ed = .AC$Ed0,
  lod_kind = "solar"
)
```

## Arguments

- a_RE:

  Earth-Moon distance in Earth radii (numeric, \> 0; the present value
  is about 60.14). Required, no default. It sets the lunar part of the
  precession torque and, when `lod_h = NULL`, the length of day.

- lod_h:

  Length of day in hours (numeric, \> 0), of the type given by
  `lod_kind`. `NULL` derives the Earth's spin rate from `a_RE` by
  conservation of the angular momentum of the Earth-Moon system
  (`omega_from_moon_distance`). `Default=NULL`.

- obliquity_deg:

  Earth's obliquity in degrees (numeric); the precession rate scales
  with its cosine. `Default=.AC$obliquity0_deg`, the present value of
  23.2545 degrees.

- Ed:

  Dynamical ellipticity of the Earth at the present spin rate (numeric,
  dimensionless); it is scaled with the square of the spin rate.
  `Default=.AC$Ed0`, i.e. 0.003243.

- lod_kind:

  Type of day of `lod_h` (character): `"solar"` or `"sidereal"`; also
  used for the length of day returned. `Default="solar"`.

## Value

list: p in arcsec/yr and the precession period in kyr
