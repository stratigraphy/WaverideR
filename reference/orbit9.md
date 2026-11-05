# Nine principal astronomical terms for a given precession constant

Nine principal astronomical terms for a given precession constant

## Usage

``` r
orbit9(k_kyr, gs = LA2004)
```

## Arguments

- k_kyr:

  Precession constant k expressed as the period of the axial precession,
  in kyr (numeric, \> 0; about 25.7 kyr today). Required (no default).

- gs:

  Secular frequencies in arcsec/yr: a list with elements `g` and `s`,
  each a numeric vector of the nine fundamental frequencies g1-g9 and
  s1-s9. The default `LA2004` is the La2004 table (Laskar et al., 2004)
  defined inside the function. `Default=LA2004`.

## Value

data frame: term, group, frequency (arcsec/yr), period (kyr)
