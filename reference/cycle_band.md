# Period range implied by the frequency resolution of a record

Period range implied by the frequency resolution of a record

## Usage

``` r
cycle_band(period, record_length = NULL, tol = NULL, mult = 1)
```

## Arguments

- period:

  Measured period of the cycle (single numeric value, \> 0): in kyr for
  a tuned record or in metres for a depth series. Required, no default.

- record_length:

  Length of the analysed record, in the same unit as `period` (numeric,
  \> 0; kyr for a tuned record, metres for a depth series). The
  half-width of the band is `mult / record_length` in frequency (the
  Rayleigh resolution); a warning is given when the record holds fewer
  than two repeats of the cycle. `NULL` uses `tol` instead.
  `Default=NULL`.

- tol:

  Relative half-width of the period range (numeric, between 0 and 1;
  `0.05` means +/-5%), used only when `record_length` is `NULL`; `NULL`
  then uses 0.05. `Default=NULL`.

- mult:

  Number of Rayleigh bandwidths used as half-width when `record_length`
  is given (numeric, \> 0; `1` means +/-1 bandwidth). `Default=1`.

## Value

c(min, max) period range
