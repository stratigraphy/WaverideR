# Conduct the cross wavelet transform on a time series or signal

Compute the cross wavelet transform using a complex Morlet wavelet based
on the "analyze.coherency" function of the WaveletComp package.

## Usage

``` r
analyse_Xwavelet(
  data_1,
  data_2,
  dj = 1/100,
  lowerPeriod = 2,
  upperPeriod = 1024,
  verbose = FALSE,
  omega_nr = 4,
  wavelet = "morlet",
  wavelet_order = NULL
)

analyze_Xwavelet(
  data_1,
  data_2,
  dj = 1/100,
  lowerPeriod = 2,
  upperPeriod = 1024,
  verbose = FALSE,
  omega_nr = 4,
  wavelet = "morlet",
  wavelet_order = NULL
)
```

## Arguments

- data_1:

  First record: a matrix or data frame with depth/time in column 1 and
  the proxy value in column 2 (no `NA` values). Both records are sorted
  and linearly interpolated to a common, regular grid over their
  overlapping interval, with the coarser of the two median sampling
  intervals. Required, no default.

- data_2:

  Second record, in the same format and depth/time unit as `data_1`
  (matrix or data frame: depth/time in column 1, proxy value in column
  2). The cross phase is the phase of `data_1` minus that of `data_2`.
  Required, no default.

- dj:

  Spacing between successive scales in octaves (numeric, \> 0), passed
  to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md):
  the periods increase by powers of 2 and `dj=1/x` gives x periods per
  doubling of the period. Smaller values give a finer period axis but
  increase the computation time. `Default=1/100`.

- lowerPeriod:

  Shortest period to analyse, in the same unit as the depth/time column
  of the records (numeric, \> 0; at least twice the common sampling
  interval). `Default=2`.

- upperPeriod:

  Longest period to analyse, in the same unit as the depth/time column
  of the records (numeric, \> `lowerPeriod`). `Default=1024`.

- verbose:

  Print the sampling interval to which the records were interpolated,
  passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  (logical, TRUE or FALSE). `Default=FALSE`.

- omega_nr:

  Number of cycles (central angular frequency omega_0) of the Morlet
  wavelet, passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  (numeric, \> 0). Larger values give a better period resolution and a
  worse depth/time resolution. `Default=4`.

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

  Parameter(s) of the wavelet (numeric scalar or vector; `NULL` uses the
  defaults listed here): `paul`: order m (default 4); `dog`: derivative
  m (default 2); `morse`: `c(gamma, beta)` or beta alone with gamma = 3
  (default `c(3, 20)`); `bump`: `c(mu, sigma)` (default `c(5, 0.6)`);
  `cmor`: `c(bandwidth, centre frequency)` (default `c(1.5, 1)`);
  `cgau`: derivative m (default 4); `shannon`:
  `c(bandwidth, centre frequency)` (default `c(2, 1.5)`); `fbsp`:
  `c(m, bandwidth, centre frequency)` (default `c(2, 2, 1.5)`);
  `loggabor`: `c(centre angular frequency, log-width)` (default
  `c(6, 0.25)`); `meyer` and `hhat` have no parameter; `morlet0`:
  omega_0 (default omega_nr); `morl`: omega_0 (default 5); `poisson`:
  order n (default 2); `gammatone`:
  `c(n, centre angular frequency, bandwidth)` (default `c(4, 6, 1.5)`);
  `hann`: `c(centre angular frequency, half width)` (default `c(6, 3)`);
  `db`: number of vanishing moments 1 to 20 (default 4, `"haar"` =
  `"db1"`); `sym`: 2 to 20 (default 8); `battle`: spline order 1, 3, 5
  or 7 (default 3). Passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  for the Morlet wavelet use `omega_nr`. `Default=NULL`.

## Value

A list with class `"analyze.Xwavelet"` containing Wave: complex wavelet
coefficients Power: time frequency power spectrum dt: sampling interval
after interpolation Phase: instantaneous phase of the signal dj: number
of frequencies Power.avg: average spectral power Period: physical
periods corresponding to frequencies nc: number of columns in the input
data nr: number of frequency levels axis.1: x axis values (time or
depth) axis.2: y axis values (log2 scaled periods) c1: base number of
wavelet cycles o: numeric vector of length two defining the minimum and
maximum superlet order. If NULL, all frequencies are analysed using
order one. x1: interpolated x values first data set y1: interpolated
signal values first data set x2: interpolated x values second data set
y2: interpolated signal values second data set

## References

Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., & Mureșan, R. C. (2021).
Time-frequency super-resolution with superlets. Nature Communications,
12(1), 337.
[doi:10.1038/s41467-020-20539-9](https://doi.org/10.1038/s41467-020-20539-9)

## Author

The "analyse_Xsuperlet" that generates the input for the plotting
function is based on the matlab code in Moca et al. (2021) and the the
"analyze.coherency" function of the 'WaveletComp' R package

## Examples

``` r
# \donttest{
#Example 1. A cross superlet of two etp solutions with noise overprint
etp_1 <- astrochron::etp(
 tmin = 0,
 tmax = 500,
 dt = 2,
 eWt = 1.5,
 oWt = 0.75,
 pWt = 1,
 esinw = TRUE,
 standardize = TRUE,
 genplot = FALSE,
 verbose = FALSE
)

etp_2 <- astrochron::etp(
 tmin = 0,
 tmax = 500,
 dt = 2,
 eWt = 1,
 oWt = 0.5,
 pWt = 1.5,
 esinw = TRUE,
 standardize = TRUE,
 genplot = FALSE,
 verbose = FALSE
)

etp_1[, 2] <- etp_1[, 2] + wr_coloured_noise(
 nrow(etp_1),
 sd = sd(etp_1[, 2]) / 1.5,
 mean = mean(etp_1[, 2]),
 phi = 0.9
)
etp_2[, 2] <- etp_2[, 2] + wr_coloured_noise(
 nrow(etp_2),
 sd = sd(etp_2[, 2]) / 1.5,
 mean = mean(etp_2[, 2]),
 phi = 0.9
)

Xetp <- analyse_Xwavelet(
 data_1 = etp_1,
 data_2  = etp_2,
 upperPeriod = 512,
     dj = 1/20,
 lowerPeriod = 4,
 verbose = FALSE,
 omega_nr = 8
)
# }
```
