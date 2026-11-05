# Which secular resonance state does a record show?

Which secular resonance state does a record show?

## Usage

``` r
resonance_state(ecc_kyr, obl_kyr, max_ratio = 4, data_dir = NULL)
```

## Arguments

- ecc_kyr:

  Measured period of the long eccentricity modulation (the "2.4 Myr"
  cycle, g4-g3), in kyr (numeric, \> 0, e.g. `2400`). Required, no
  default.

- obl_kyr:

  Measured period of the obliquity modulation (the "1.2 Myr" cycle,
  s4-s3), in kyr (numeric, \> 0, e.g. `1200`). Required, no default.

- max_ratio:

  Largest integer p and q considered when searching the simple ratio p:q
  closest to the measured frequency ratio f(obliquity modulation) /
  f(eccentricity modulation) (integer, \>= 1). `Default=4`.

- data_dir:

  Directory that contains the data file `La04gisi.csv` with the La2004
  secular frequencies used for the modern values (character path).
  `NULL` uses the data installed with the package (or
  `options(cocoR.data)`). `Default=NULL`.

## Value

list: the implied frequencies, their ratio, the nearest simple
resonance, the resonance argument (arcsec/yr and as a period), and the
modern values for comparison
