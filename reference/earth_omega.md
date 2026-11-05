# Earth's angular velocity from the length of day

Earth's angular velocity from the length of day

## Usage

``` r
earth_omega(lod_h, kind = c("solar", "sidereal"))
```

## Arguments

- lod_h:

  Length of day, in hours (numeric, scalar or vector, \> 0). Required
  (no default).

- kind:

  Which length of day `lod_h` is (character): `"solar"` (solar day, as
  in the `lod` column of AstroGeo22, converted to the sidereal rotation
  with 1/T_sid = 1/T_sol + 1/yr) or `"sidereal"` (sidereal rotation
  period). `Default=c("solar", "sidereal")`, i.e. `"solar"`.
