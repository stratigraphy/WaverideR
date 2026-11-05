# Cross-wavelet coherence analysis

Computes the cross-wavelet transform and wavelet coherence between two
time/depth series using a superlet-based wavelet approach. The function
internally interpolates both input series onto a common grid, computes
individual wavelet transforms, and derives cross-wavelet power, phase,
and coherence with configurable smoothing.

## Usage

``` r
analyse_wavelet_coherence(
  data_1,
  data_2,
  dj = 1/100,
  lowerPeriod = 2,
  upperPeriod = 1024,
  verbose = FALSE,
  omega_nr = 8,
  wavelet = "morlet",
  wavelet_order = NULL,
  n_simulations = 10,
  window.type.t = 4,
  window.type.s = 4,
  abs_window_t = 500,
  abs_window_s = 0.5
)

analyze_wavelet_coherence(
  data_1,
  data_2,
  dj = 1/100,
  lowerPeriod = 2,
  upperPeriod = 1024,
  verbose = FALSE,
  omega_nr = 8,
  wavelet = "morlet",
  wavelet_order = NULL,
  n_simulations = 10,
  window.type.t = 4,
  window.type.s = 4,
  abs_window_t = 500,
  abs_window_s = 0.5
)
```

## Arguments

- data_1:

  First record: a two-column matrix or data frame with depth/time in
  column 1 and the proxy value in column 2 (no `NA` values). Both
  records are sorted and linearly interpolated to a common, regular grid
  over their overlapping interval, with the coarser of the two median
  sampling intervals. Required, no default.

- data_2:

  Second record, in the same format and depth/time unit as `data_1`
  (two-column matrix or data frame: depth/time, proxy). The phase is the
  phase of `data_1` minus that of `data_2`. Required, no default.

- dj:

  Spacing between successive scales in octaves (numeric, \> 0), passed
  to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  `dj=1/x` gives x periods per doubling of the period, smaller values
  give a finer period axis. `Default=1/100`.

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
  (numeric, \> 0); controls the trade-off between period and depth/time
  resolution. `Default=8`.

- wavelet:

  Mother wavelet (character string or function), passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md),
  given as a name (PyWavelets style names with the order appended, such
  as `"db4"`, `"sym8"`, `"cgau3"`, `"gaus2"` or `"cmor1.5-1"`, also
  work) or as a function `psi_hat(omega)` that returns the Fourier
  transform of a user defined mother wavelet. Analytic (complex)
  wavelets, which give phase and amplitude and can be synchrosqueezed:
  `"morlet"` (default; aliases `"gabor"`, `"amor"`), `"morlet0"`
  (exactly zero mean Morlet, admissible for a small omega_nr; alias
  `"morlex"`), `"paul"` (alias `"cauchy"`), `"morse"` (generalised Morse
  wavelet of Olhede and Walden, 2002; aliases `"gmw"`, `"airy"`),
  `"bump"`, `"cmor"` (complex Morlet with bandwidth and centre frequency
  as in Matlab/PyWavelets), `"cgau"` (analytic Gaussian derivative;
  `"cmhat"` is the complex Mexican hat), `"shannon"` (complex Shannon,
  box spectrum), `"fbsp"` (frequency B-spline), `"meyer"`, `"loggabor"`
  (alias `"lognormal"`), `"gammatone"` (asymmetric in time: sharp onset,
  slow decay) and `"hann"` (raised cosine band). Real wavelets, which
  give a real transform (no phase or envelope): `"dog"` (derivative of
  Gaussian; `"mexh"`, `"ricker"`, `"marr"` and `"sombrero"` are the
  Mexican hat, order 2), `"morl"` (real Morlet), `"hhat"` (Hermitian
  hat), `"poisson"`, and the orthogonal wavelets `"haar"`, `"db"`
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
  or 7 (default 3). For the Morlet wavelet use `omega_nr`. Passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md).
  `Default=NULL`.

- n_simulations:

  Number of simulations (integer). Not used: reserved for a future
  significance test. `Default=10`.

- window.type.t:

  Type of the smoothing window applied in the depth/time direction
  before the coherence is computed, given as a number or a character
  code: `0` or `"none"` (no smoothing), `1` or `"bar"` (Bartlett), `2`
  or `"tri"` (triangular), `3` or `"box"` (boxcar), `4` or `"han"`
  (Hanning), `5` or `"ham"` (Hamming), `6` or `"bla"` (Blackman).
  `Default=4` (Hanning).

- window.type.s:

  Type of the smoothing window applied in the period (scale) direction;
  same options as `window.type.t`. `Default=4` (Hanning).

- abs_window_t:

  Length of the smoothing window in the depth/time direction, in the
  same unit as the depth/time column of the records (numeric, \> 0);
  converted to an odd number of samples and capped at about a quarter of
  the record length (with a warning). Longer windows give smoother but
  less local coherence. `Default=500`.

- abs_window_s:

  Length of the smoothing window in the period (scale) direction, in
  octaves (numeric, \> 0; converted to an odd number of scale steps of
  size `dj`). `Default=0.5`.

## Value

A list of class containing: Wave Smoothed cross-wavelet transform
(complex) Coherence Wavelet coherence matrix Phase Phase difference
between the two signals Coh.avg Average coherence over time Period
(m)vector corresponding to scales Scale Wavelet scales dt Time step of
interpolated grid dj Scale resolution axis.1 Time/depth axis axis.2
Period axis nc, nr Dimensions of the wavelet matrices x1, y1
Interpolated first dataset x2, y2 Interpolated second dataset

## References

Torrence, C., and G. P. Compo (1998). A practical guide to wavelet
analysis. *Bulletin of the American Meteorological Society*, 79(1),
61–78.

## Author

plotting code based on the "analyze.coherency" function of the
'WaveletComp' R package

## Examples

``` r
# \donttest{
# Generate two synthetic signals
t <- seq(0, 1000, by = 1)
x1 <- sin(2 * pi * t / 100) + rnorm(length(t), 0, 0.5)
x2 <- sin(2 * pi * t / 100 + pi/4) + rnorm(length(t), 0, 0.5)

data_1 <- cbind(t, x1)
data_2 <- cbind(t, x2)

coh <- analyse_wavelet_coherence(
  data_1 = data_1,
  data_2 = data_2,
  lowerPeriod = 2,
  upperPeriod = 256
)

# }
```
