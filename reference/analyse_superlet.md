# Conduct the superlet transform on a time series or signal

Compute the superlet transform using a complex Morlet wavelet following
the approach of Moca et al. (2021). The superlet transform increases
time frequency resolution by combining multiple wavelet orders at each
frequency. Both the frequency sampling and the scaling of superlet order
with frequency can be explicitly controlled.

## Usage

``` r
analyse_superlet(
  data,
  upperPeriod,
  lowerPeriod,
  Nf,
  c1,
  o = NULL,
  mult = TRUE,
  order_scaling = c("none", "log2", "linear", "sqrt", "quadratic", "power"),
  order_alpha = 1,
  verbose = FALSE,
  detrend = c("linear", "mean", "none"),
  standardize = TRUE,
  run_multicore = FALSE,
  n_cores = NULL,
  standardise = NULL,
  wavelet = "morlet",
  wavelet_order = NULL
)

analyze_superlet(
  data,
  upperPeriod,
  lowerPeriod,
  Nf,
  c1,
  o = NULL,
  mult = TRUE,
  order_scaling = c("none", "log2", "linear", "sqrt", "quadratic", "power"),
  order_alpha = 1,
  verbose = FALSE,
  detrend = c("linear", "mean", "none"),
  standardize = TRUE,
  run_multicore = FALSE,
  n_cores = NULL,
  standardise = NULL,
  wavelet = "morlet",
  wavelet_order = NULL
)
```

## Arguments

- data:

  Input record: a matrix or data frame with depth/time in column 1 and
  the proxy value in column 2. Rows with `NA` are removed, the rows are
  sorted by depth/time and the record is linearly interpolated to the
  median sampling interval. Required, no default.

- upperPeriod:

  Longest period to analyse, in the same unit as the depth/time column
  of `data` (numeric, \> `lowerPeriod`); sets the lowest analysed
  frequency. Required, no default.

- lowerPeriod:

  Shortest period to analyse, in the same unit as the depth/time column
  of `data` (numeric, \> 0; at least twice the sampling interval,
  otherwise a warning is given); sets the highest analysed frequency.
  Required, no default.

- Nf:

  Number of frequencies of the superlet spectrum (integer, \>= 1); the
  frequencies are spaced evenly in log2 between `1/upperPeriod` and
  `1/lowerPeriod`. More frequencies give a finer period axis and a
  longer computation time. Required, no default.

- c1:

  Base number of cycles of the Morlet wavelet of order 1 (numeric, \> 0,
  typically 2 to 5): the fundamental wavelet width. The wavelets of
  higher orders have more cycles, see `mult`. Required, no default.

- o:

  Minimum and maximum superlet order (numeric vector of length 2, values
  \>= 1; fractional orders are allowed and give the last wavelet a
  fractional weight). The order changes from `o[1]` at the lowest
  frequency (longest period) to `o[2]` at the highest frequency, as set
  by `order_scaling`. A single number gives the same order at all
  frequencies; `NULL` uses order 1 everywhere (a plain wavelet
  transform). `Default=NULL`.

- mult:

  How the number of cycles grows with the order i within a superlet
  (logical, TRUE or FALSE): `TRUE` multiplicative (\\c_i = i c_1\\);
  `FALSE` additive (\\c_i = c_1 + i - 1\\, Moca et al., 2021).
  `Default=TRUE`.

- order_scaling:

  How the superlet order increases from `o[1]` (lowest frequency) to
  `o[2]` (highest frequency) (character), one of:

  - `"none"` or `"log2"`: linear in log2(frequency);

  - `"linear"`: linear in frequency (as in Moca et al., 2021);

  - `"sqrt"`: with the square root of the normalised frequency (fast
    increase at low frequencies);

  - `"quadratic"`: with the square of the normalised frequency (fast
    increase at high frequencies);

  - `"power"`: with the normalised frequency to the power `order_alpha`.

  `Default="none"` (the first element of
  `c("none", "log2", "linear", "sqrt", "quadratic", "power")`).

- order_alpha:

  Exponent of the `"power"` order scaling (numeric, \> 0); only used
  when `order_scaling="power"`. Values greater than 1 keep the order low
  over most of the frequency range and raise it mainly at the highest
  frequencies (sharpening the short periods); values smaller than 1
  raise it already at low frequencies. `Default=1`.

- verbose:

  Print the sampling interval to which the record was interpolated and
  show a progress bar (logical, TRUE or FALSE). `Default=FALSE`.

- detrend:

  Removal of a trend before the transform (character), one of `"linear"`
  (a least-squares linear trend is subtracted), `"mean"` (the mean is
  subtracted) or `"none"`. A mean or trend that is not removed leaks
  into all periods and gives large edge effects, which is a problem for
  stratigraphic records (for example grey values or counts) that do not
  have a zero mean. `Default="linear"` (the first element of
  `c("linear", "mean", "none")`).

- standardize:

  Divide the (detrended) record by its standard deviation before the
  transform (logical, TRUE or FALSE), as in
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md),
  so that the power does not depend on the units of the proxy.
  `Default=TRUE`.

- run_multicore:

  Compute the frequencies in parallel on several cores (logical, TRUE or
  FALSE). `Default=FALSE`.

- n_cores:

  Number of cores used when `run_multicore=TRUE` (integer, \> 0); `NULL`
  picks a safe automatic number (see
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)).
  Results do not depend on the number of cores. `Default=NULL`.

- standardise:

  British spelling of `standardize` (logical); if given (not `NULL`) it
  overrides `standardize`. `Default=NULL`.

- wavelet:

  Mother wavelet (character string or function), given as a name
  (PyWavelets style names with the order appended, such as `"db4"`,
  `"sym8"`, `"cgau3"`, `"gaus2"` or `"cmor1.5-1"`, also work) or as a
  function `psi_hat(omega)` that returns the Fourier transform of a user
  defined mother wavelet. Analytic (complex) wavelets, which give phase
  and amplitude and can be synchrosqueezed: `"morlet"` (default; aliases
  `"gabor"`, `"amor"`), `"morlet0"` (exactly zero mean Morlet,
  admissible for a small omega_nr; alias `"morlex"`), `"paul"` (alias
  `"cauchy"`), `"morse"` (generalised Morse wavelet of Olhede and
  Walden, 2002; aliases `"gmw"`, `"airy"`), `"bump"`, `"cmor"` (complex
  Morlet with bandwidth and centre frequency as in Matlab/PyWavelets),
  `"cgau"` (analytic Gaussian derivative; `"cmhat"` is the complex
  Mexican hat), `"shannon"` (complex Shannon, box spectrum), `"fbsp"`
  (frequency B-spline), `"meyer"`, `"loggabor"` (alias `"lognormal"`),
  `"gammatone"` (asymmetric in time: sharp onset, slow decay) and
  `"hann"` (raised cosine band). Real wavelets, which give a real
  transform (no phase or envelope): `"dog"` (derivative of Gaussian;
  `"mexh"`, `"ricker"`, `"marr"` and `"sombrero"` are the Mexican hat,
  order 2), `"morl"` (real Morlet), `"hhat"` (Hermitian hat),
  `"poisson"`, and the orthogonal wavelets `"haar"`, `"db"`
  (Daubechies), `"sym"` (Symlets) and `"battle"` (Battle-Lemarie spline
  wavelets), which are computed in the frequency domain from their
  filters. The period axis of the orthogonal wavelets is approximate
  because their spectra are broad and not symmetric. A narrow wavelet
  (Paul, low-order DOG, Shannon, Haar) resolves position well, a wide
  one (Morse with a large beta, Morlet with a high omega_nr, bump)
  resolves period well. `Default="morlet"`.

- wavelet_order:

  Parameter(s) of the wavelet (see
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md));
  when it is NULL the wavelet with the same relative bandwidth as a
  Morlet wavelet with the superlet cycles is used at every frequency
  (cycles are converted as omega_0 = 2 pi cycles / 5). Parameters:
  `paul`: order m (default 4); `dog`: derivative m (default 2); `morse`:
  `c(gamma, beta)` or beta alone with gamma = 3 (default `c(3, 20)`);
  `bump`: `c(mu, sigma)` (default `c(5, 0.6)`); `cmor`:
  `c(bandwidth, centre frequency)` (default `c(1.5, 1)`); `cgau`:
  derivative m (default 4); `shannon`: `c(bandwidth, centre frequency)`
  (default `c(2, 1.5)`); `fbsp`: `c(m, bandwidth, centre frequency)`
  (default `c(2, 2, 1.5)`); `loggabor`:
  `c(centre angular frequency, log-width)` (default `c(6, 0.25)`);
  `meyer` and `hhat` have no parameter; `morlet0`: omega_0 (default
  omega_nr); `morl`: omega_0 (default 5); `poisson`: order n (default
  2); `gammatone`: `c(n, centre angular frequency, bandwidth)` (default
  `c(4, 6, 1.5)`); `hann`: `c(centre angular frequency, half width)`
  (default `c(6, 3)`); `db`: number of vanishing moments 1 to 20
  (default 4, `"haar"` = `"db1"`); `sym`: 2 to 20 (default 8); `battle`:
  spline order 1, 3, 5 or 7 (default 3). For the Morlet wavelet the
  number of cycles is set by `c1`, `o` and `mult`. `Default=NULL`.

## Value

A list with class `"analyze.superlet"` containing Power: time frequency
power spectrum dt: sampling interval after interpolation dj: number of
frequencies Power.avg: average spectral power Period: physical periods
corresponding to frequencies nc: number of columns in the input data nr:
number of frequency levels axis.1: x axis values (time or depth) axis.2:
y axis values (log2 scaled periods) c1: base number of wavelet cycles o:
numeric vector of length two defining the minimum and maximum superlet
order. If NULL, all frequencies are analysed using order one.
freq.scale: scaling of the frequencies order.mode: order.mode ("none",
"log2", "linear", "sqrt", "power") order: superlet order of every period
coi: distance from the start and end of the record within which edge
effects are important (cone of influence), for every period coi.1,
coi.2: the cone of influence as polygons (x and log2 period), drawn by
the plot functions x: interpolated x values y: interpolated signal
values

## References

Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., & Mureșan, R. C. (2021).
Time-frequency super-resolution with superlets. Nature Communications,
12(1), 337.
[doi:10.1038/s41467-020-20539-9](https://doi.org/10.1038/s41467-020-20539-9)

## Author

Adapted from MATLAB implementations by V. V. Moca et al. (2021)

## Examples

``` r
# \donttest{
## Example 1. Using the Total Solar Irradiance data set of Steinhilber et al. (2012)
TSI_sl <-
  analyse_superlet(
    data = TSI,
    upperPeriod = 8192,
    lowerPeriod = 16,
    Nf = 20,
    c1 = 2,
    o = c(1, 5),
    mult = TRUE,
    order_scaling = "log2",
    order_alpha = 1,
    verbose = FALSE
  )
#> Warning: lowerPeriod is shorter than twice the sampling interval (44): periods below the Nyquist period cannot be resolved


## Example 2. Using the magnetic susceptibility data set of Pas et al. (2018)
mag_sl <-
  analyse_superlet(
    data = mag,
    upperPeriod = 254,
    lowerPeriod = 0.1,
    Nf = 20,
    c1 = 2,
    o = c(1, 5),
    mult = TRUE,
    order_scaling = "log2",
    order_alpha = 1,
    verbose = FALSE
  )


## Example 3. Using the greyscale data set of Zeeden et al. (2013)
grey_sl <-
  analyse_superlet(
    data = grey,
    upperPeriod = 256,
    lowerPeriod = 0.02,
    Nf = 20,
    c1 = 2,
    o = c(2,5),
    mult = TRUE,
    order_scaling = "log2",
    order_alpha = 1,
    verbose = FALSE)
# }
```
