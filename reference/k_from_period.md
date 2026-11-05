# Precession constant implied by one measured cycle period

Precession constant implied by one measured cycle period

## Usage

``` r
k_from_period(period_kyr, term = "k+s3", gs = LA2004)
```

## Arguments

- period_kyr:

  Measured period of the cycle in kyr (numeric, \> 0; a vector gives one
  value per period). Required (no default).

- term:

  Astronomical term the cycle is read as (character), one of `"k+s3"`
  (obliquity), `"k+g5"`, `"k+g2"`, `"k+g4"` or `"k+g3"` (precession);
  its secular frequency is subtracted from the observed frequency to
  give the precession constant. `Default="k+s3"`.

- gs:

  Secular frequencies in arcsec/yr: a list with elements `g` and `s`,
  each a numeric vector of the nine fundamental frequencies g1-g9 and
  s1-s9. The default `LA2004` is the La2004 table (Laskar et al., 2004)
  defined inside the function. `Default=LA2004`.
