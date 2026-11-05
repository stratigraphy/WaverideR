# Frequency uncertainty of the Fourier transform (FFT) for a record length

Calculator for the frequency uncertainty of a Fourier based spectrum.
The practical resolution of an FFT is the bin spacing, the inverse of
the record length (\\\Delta f = 1/T\\). The smallest frequency
uncertainty that is possible according to the Heisenberg-Gabor
uncertainty principle (\\\Delta t \Delta f \ge 1/(4\pi)\\) with \\\Delta
t = T\\ is \\\Delta f = 1/(4\pi T)\\, a factor \\4\pi\\ smaller. The
relative uncertainty at a frequency f is \\\Delta f / f = 1/(f T)\\. For
a record of N samples with sampling frequency fs, \\T = N / f_s\\.

It is the Fourier counterpart of
[`wavelet_uncertainty`](https://stratigraphy.github.io/WaverideR/reference/wavelet_uncertainty.md),
which gives the uncertainty of the wavelet transform. Units are free:
with the record length in metres (or kyr) the frequencies are in cycles
per metre (or per kyr).

## Usage

``` r
fft_frequency_uncertainty(
  T_length = NULL,
  N = NULL,
  fs = 1,
  f = NULL,
  period = NULL
)
```

## Arguments

- T_length:

  Length of the record \\T\\ in time or depth (numeric, \> 0; a vector
  gives one row per value), in any unit (e.g. metres or kyr); the
  frequencies of the output are in cycles per that unit. If `NULL`, it
  is computed as `N / fs`, so either `T_length` or `N` must be given.
  `Default=NULL`.

- N:

  Number of samples of the record (numeric, \> 0); only used when
  `T_length` is `NULL`, together with `fs`, to compute
  `T_length = N / fs`. `Default=NULL`.

- fs:

  Sampling frequency, i.e. samples per unit of time or depth (numeric,
  \> 0; the inverse of the sampling interval); only used with `N`.
  `Default=1`.

- f:

  Frequency or frequencies of interest in cycles per unit of `T_length`
  (numeric vector, \> 0). If given, the output also contains the
  relative uncertainties and the period ranges at these frequencies.
  `NULL` gives only the absolute uncertainties. `Default=NULL`.

- period:

  Period(s) of interest in the unit of `T_length` (numeric vector, \>
  0), an alternative to `f` (converted with `f = 1 / period`); ignored
  when `f` is given. `Default=NULL`.

## Value

A data frame with the record length, the practical uncertainty
`df_practical`, the Heisenberg-Gabor limit `df_gabor`, their ratio and,
when `f` or `period` is given, the frequency, the relative uncertainties
(`rel_practical`, `rel_gabor`) and the period ranges implied by `f - df`
and `f + df` for both definitions.

## References

Gabor, D. (1946). Theory of communication. Part 1: The analysis of
information. Journal of the Institution of Electrical Engineers - Part
III 93, 429-441.

## Examples

``` r
# 2 s record, 100 Hz: practical 0.5 Hz, Heisenberg-Gabor limit 0.040 Hz, 0.5 percent
fft_frequency_uncertainty(T_length = 2, f = 100)
#>   T_length df_practical   df_gabor    ratio   f rel_practical    rel_gabor
#> 1        2          0.5 0.03978874 12.56637 100         0.005 0.0003978874
#>   period period_min_practical period_max_practical period_min_gabor
#> 1   0.01          0.009950249           0.01005025      0.009996023
#>   period_max_gabor
#> 1       0.01000398
# 1000 samples every 0.5 kyr and the 405 kyr eccentricity cycle
fft_frequency_uncertainty(N = 1000, fs = 1 / 0.5, period = 405)
#>   T_length df_practical     df_gabor    ratio           f rel_practical
#> 1      500        0.002 0.0001591549 12.56637 0.002469136          0.81
#>    rel_gabor period period_min_practical period_max_practical period_min_gabor
#> 1 0.06445775    405             223.7569             2131.579         380.4754
#>   period_max_gabor
#> 1          432.904
```
