# Full Earth-Moon state from a precession constant

Full Earth-Moon state from a precession constant

## Usage

``` r
em_state(k_kyr, obl_deg = AC$obl0_deg)
```

## Arguments

- k_kyr:

  Precession constant expressed as the period of the axial precession,
  in kyr (numeric, scalar or vector, \> 0; about 25.7 today), i.e. a
  precession rate of `360 * 3600 / (k_kyr * 1000)` arcsec/yr. The
  Earth-Moon distance giving this rate is searched between 20 Earth
  radii and the present distance (conserving the Earth-Moon angular
  momentum); values without a solution give `NA`. Required (no default).

- obl_deg:

  Obliquity of the Earth used in the precession equation (Farhat et al.
  2022, eq. A.7), in degrees (numeric). `Default=AC$obl0_deg`, i.e. the
  present-day obliquity, 23.2545 degrees.
