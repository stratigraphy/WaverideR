# FT and IAAFT surrogates for bispectral significance testing

The `bispectral_surrogate` function generates null-hypothesis surrogate
series. Fourier-transform (FT) phase-randomised surrogates preserve the
periodogram exactly (all Fourier amplitudes, including DC and the
even-length Nyquist bin) while randomising the phases. Iterative
amplitude-adjusted Fourier-transform (IAAFT) surrogates (Schreiber and
Schmitz, 1996) preserve the observed marginal distribution exactly and
iteratively approach the observed Fourier amplitudes; the iteration
stops early when the relative spectral error changes by less than 1e-10.
The reported relative spectral error is evaluated only on the
independent positive frequencies (DC and the even-length Nyquist bin
excluded). Translated from Acycle v3.1 bispectralSurrogate.m (Mingsong
Li).

## Usage

``` r
bispectral_surrogate(
  data = NULL,
  type = "iaaft",
  iaaft_iter = 200,
  n_sim = 1,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  Series to make surrogates of: a numeric vector, or a matrix or data
  frame with depth/time in column 1 and the proxy in column 2 (only
  column 2 is used; a one-column matrix uses that column). All values
  must be finite (no `NA`). Required (`Default=NULL` gives an error).

- type:

  Kind of surrogate (character): `"iaaft"` (iterative amplitude-adjusted
  Fourier transform; keeps the value distribution exactly and the
  Fourier amplitudes approximately) or `"ft"` (Fourier phase
  randomisation; keeps the Fourier amplitudes exactly; `"phase"` and
  `"fourier"` are accepted as synonyms). Case is ignored.
  `Default="iaaft"`.

- iaaft_iter:

  Maximum number of IAAFT iterations per surrogate (positive integer);
  the iteration stops earlier when the relative spectral error no longer
  changes. Ignored when `type="ft"`. `Default=200`.

- n_sim:

  Number of surrogate series to generate (positive integer); each one
  becomes a column of the output matrix `surrogates`. `Default=1`.

- seed:

  Seed of the random number generator (integer), so that the surrogates
  can be reproduced; `NULL` leaves the current random state unchanged.
  `Default=NULL`.

- verbose:

  Show a text progress bar over the surrogates (logical, TRUE or FALSE).
  `Default=FALSE`.

## Value

A list with `surrogates` (matrix, one surrogate per column), `type`,
`iterations` (completed IAAFT iterations per surrogate, 1 for FT) and
`relative_spectral_error` (per surrogate).

## References

Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., Farmer, J.D.,
1992. Testing for nonlinearity in time series: the method of surrogate
data. Physica D 58, 77-94.

Schreiber, T., Schmitz, A., 1996. Improved surrogate data for
nonlinearity tests. Physical Review Letters 77, 635-638.

## Author

Michiel Arts (R translation); original MATLAB code: Mingsong Li (Acycle
bispectral module).

## Examples

``` r
# \donttest{
set.seed(2)
x <- exp(0.35 * as.numeric(arima.sim(list(ar = 0.75), n = 512)))
sur <- bispectral_surrogate(data = x, type = "iaaft", n_sim = 5, seed = 3)
# marginal distribution preserved exactly
max(abs(sort(sur$surrogates[, 1]) - sort(x)))
#> [1] 0
sur$relative_spectral_error
#> [1] 0.012638533 0.013390366 0.011535024 0.011472739 0.009397641
# }
```
