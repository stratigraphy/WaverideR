# Target periodogram from an astronomical time series (e.g. La2004 ETP)

Periodogram of an astronomical series, or of an ETP series built from
eccentricity, obliquity and precession as in gentarget.m.

## Usage

``` r
coco_target_from_series(
  time_kyr,
  value = NULL,
  ecc = NULL,
  obl = NULL,
  prec = NULL,
  ne = 1,
  nt = 0.6,
  np = 0.5,
  f1 = 0,
  f2 = 0.06,
  pad = 10000
)
```

## Arguments

- time_kyr:

  Time axis of the astronomical series, in kyr (numeric vector, evenly
  spaced; the series is sorted by it and its first step is used as the
  sampling interval). Required, no default.

- value:

  Astronomical series to use as target (numeric vector, same length as
  `time_kyr`), e.g. an ETP series; `NULL` builds an ETP series from
  `ecc`, `obl` and `prec`. `Default=NULL`.

- ecc:

  Eccentricity series (numeric vector, same length and order as
  `time_kyr`), used with `obl` and `prec` instead of `value`: the ETP
  series is `ne * z(ecc) + nt * z(obl) - np * z(prec)`, with `z` the
  standardisation to zero mean and unit standard deviation (as in Acycle
  gentarget.m). `Default=NULL`.

- obl:

  Obliquity series (numeric vector, same length and order as
  `time_kyr`), used to build the ETP series when `value` is `NULL`.
  `Default=NULL`.

- prec:

  Precession (index) series (numeric vector, same length and order as
  `time_kyr`), used with a negative sign to build the ETP series when
  `value` is `NULL`. `Default=NULL`.

- ne:

  Weight of the standardised eccentricity in the ETP series (numeric).
  `Default=1`.

- nt:

  Weight of the standardised obliquity in the ETP series (numeric).
  `Default=0.6`.

- np:

  Weight of the standardised precession in the ETP series (numeric).
  `Default=0.5`.

- f1:

  Lowest frequency of the target that is returned, in cycles/kyr
  (numeric, \>= 0). `Default=0`.

- f2:

  Highest frequency of the target that is returned, in cycles/kyr
  (numeric, \> `f1`). `Default=0.06`.

- pad:

  Length (number of points) to which the series is zero-padded for the
  periodogram (integer; at least the series length is used); a larger
  value gives a finer frequency grid. `Default=10000`.

## Value

two-column matrix: frequency (cycles/kyr), power
