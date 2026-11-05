# Autoregressive prewhitening of a series for aMTM analysis

The `amtm_prewhiten` function prewhitens an evenly spaced series with an
AR(p) filter whose coefficients are obtained from the autocovariance
sequence of a low-variance aMTM spectrum (NW = 4, K = 7, zero padding to
`2^(ceiling(log2(N)) + nfft_add)`) through the Yule-Walker equations.
The one-step predictions minus the observations form the prewhitened
series; the first `2*order` values are dropped and the mean is removed.
The spectrum of the prewhitened series can be "recoloured" by dividing
it by \\\|1 - \sum_j \phi_j e^{-i 2\pi f j \Delta t}\|^2\\. Translated
from the "Python PreWhiten Filter Example" (prewhiten(),
correct_spec_for_pw()) distributed with amtm_signaldetect (MIT licence),
itself based on an R prewhiten() routine (Wesley Burr / David Riegert).

## Usage

``` r
amtm_prewhiten(data = NULL, order = 2, nw = 4, k = 7, nfft_add = 2)
```

## Arguments

- data:

  Series to prewhiten: a matrix or data frame with depth/time in column
  1 and the proxy value in column 2, evenly spaced (extra columns are
  ignored, rows with non-finite values are dropped and the rows are
  sorted by depth/time). Required (`Default=NULL` gives an error).

- order:

  Order p of the autoregressive (AR) prewhitening filter (integer, \>=
  1; the series must have more than `2*order + 2` values). Higher orders
  remove more of the red-noise background but also more of the signal;
  the first `2*order` samples are lost. `Default=2`.

- nw:

  Time-bandwidth product of the multitaper spectrum from which the
  autocovariance sequence (ACVS) is computed (numeric, \> 0, typically 2
  to 4). `Default=4`.

- k:

  Number of DPSS tapers of that spectrum (integer, \>= 1, at most about
  `2*nw - 1`). `Default=7`.

- nfft_add:

  Zero-padding exponent increment (integer, \>= 0): the spectrum is
  computed on `2^(ceiling(log2(N)) + nfft_add)` frequencies, N being the
  number of samples. Larger values give a finer frequency grid and a
  longer, smoother ACVS. `Default=2`.

## Value

A list with `data` (2-column matrix of the prewhitened series, with the
depth/time values of the retained samples), `ar_coef` (AR coefficients
phi_1..phi_p) and `acvs` (autocovariance sequence used).

## References

Percival, D. B., & Walden, A. T. (2020). Spectral Analysis for
Univariate Time Series. Cambridge University Press.
[doi:10.1017/9781139235723](https://doi.org/10.1017/9781139235723)

## Author

Michiel Arts (R translation); original Python code: Hector A. Salinas et
al.

## Examples

``` r
# \donttest{
set.seed(3)
x <- as.numeric(stats::filter(rnorm(400), 0.9, method = "recursive"))
pw <- amtm_prewhiten(cbind(seq_along(x), x), order = 1)
pw$ar_coef
#> [1] 0.8846462
# }
```
