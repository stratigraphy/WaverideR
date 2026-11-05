# Conduct the adaptive continuous wavelet transform on a time series or signal

Compute the continuous wavelet transform (CWT) using a complex Morlet
wavelet with scale dependent wavelet width. The number of oscillatory
cycles in the Morlet wavelet is allowed to vary smoothly with scale,
enabling adaptive time frequency resolution similar in spirit to
superlet based approaches.

## Usage

``` r
analyse_Awavelet(
  data = NULL,
  dj = 1/100,
  lowerPeriod = 2,
  upperPeriod = 1024,
  verbose = FALSE,
  omega_min = 6,
  omega_max = 12,
  scaling = c("none", "log2", "linear", "sqrt", "quadratic", "power"),
  alpha = 1,
  n_simulations = 10,
  run_multicore = FALSE,
  wavelet = "morlet",
  wavelet_order = NULL
)

analyze_Awavelet(
  data = NULL,
  dj = 1/100,
  lowerPeriod = 2,
  upperPeriod = 1024,
  verbose = FALSE,
  omega_min = 6,
  omega_max = 12,
  scaling = c("none", "log2", "linear", "sqrt", "quadratic", "power"),
  alpha = 1,
  n_simulations = 10,
  run_multicore = FALSE,
  wavelet = "morlet",
  wavelet_order = NULL
)
```

## Arguments

- data:

  Input record: a matrix or data frame with depth/time in column 1 and
  the proxy value in column 2. Rows with `NA` are removed, the rows are
  sorted by depth/time and the record is linearly interpolated to an
  even spacing equal to the median sampling interval. Required
  (`Default=NULL` gives an error).

- dj:

  Spacing between successive scales in octaves (numeric, \> 0): scales
  increase as `2^(j * dj)`, so `1/100` gives 100 scales per doubling of
  the period. Smaller values give a finer period axis at the cost of
  computation time. `Default=1/100`.

- lowerPeriod:

  Shortest period to analyse, in the same unit as the depth/time column
  of `data` (numeric, \> 0; at least twice the sampling interval).
  `Default=2`.

- upperPeriod:

  Longest period to analyse, in the same unit as the depth/time column
  of `data` (numeric, \> `lowerPeriod`; normally not more than the
  length of the record). `Default=1024`.

- verbose:

  Print the sampling interval to which the record was interpolated
  (logical, TRUE or FALSE). `Default=FALSE`.

- omega_min:

  Number of oscillatory cycles of the Morlet wavelet (omega_0) at the
  shortest period (numeric, \> 0, typically 5 to 6). Fewer cycles give a
  better depth/time resolution, more cycles a better period resolution.
  `Default=6`.

- omega_max:

  Number of oscillatory cycles of the Morlet wavelet at the longest
  period (numeric, \> 0). The number of cycles changes from `omega_min`
  to `omega_max` along the period axis as set by `scaling`.
  `Default=12`.

- scaling:

  How the number of wavelet cycles changes from `omega_min` (shortest
  period) to `omega_max` (longest period) (character), one of:

  - `"none"`: linear in the scale index, which is equivalent to linear
    in log2(period);

  - `"log2"`: linear in log2(period);

  - `"linear"`: linear in period;

  - `"sqrt"`: with the square root of the normalised log2(period) (fast
    change at short periods);

  - `"quadratic"`: with the square of the normalised log2(period) (fast
    change at long periods);

  - `"power"`: with the normalised log2(period) to the power `alpha`.

  `Default="none"` (the first element of
  `c("none", "log2", "linear", "sqrt", "quadratic", "power")`).

- alpha:

  Exponent of the `"power"` scaling (numeric, \> 0); only used when
  `scaling="power"`. Values greater than 1 keep the number of cycles
  close to `omega_min` over most of the period range and increase it
  mainly at the longest periods; values smaller than 1 do the opposite.
  `Default=1`.

- n_simulations:

  Number of Monte Carlo simulations (integer). Not used: kept only for
  compatibility with the interface of the other wavelet functions (no
  significance test is made). `Default=10`.

- run_multicore:

  Run in parallel on several cores (logical, TRUE or FALSE). Not used:
  kept only for compatibility with the interface of the other wavelet
  functions. `Default=FALSE`.

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
  one (Morse with a large beta, Morlet with a high number of cycles,
  bump) resolves period well. For wavelets other than the Morlet, the
  wavelet at each scale is given the relative bandwidth of a Morlet
  wavelet with the local number of cycles (between `omega_min` and
  `omega_max`). `Default="morlet"`.

- wavelet_order:

  Parameter(s) of the wavelet (numeric scalar or vector, see
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md));
  when it is `NULL` the wavelet whose relative bandwidth equals that of
  a Morlet wavelet with `omega_min` to `omega_max` cycles is used at
  every scale. Parameters of the wavelet: `paul`: order m (default 4);
  `dog`: derivative m (default 2); `morse`: `c(gamma, beta)` or beta
  alone with gamma = 3 (default `c(3, 20)`); `bump`: `c(mu, sigma)`
  (default `c(5, 0.6)`); `cmor`: `c(bandwidth, centre frequency)`
  (default `c(1.5, 1)`); `cgau`: derivative m (default 4); `shannon`:
  `c(bandwidth, centre frequency)` (default `c(2, 1.5)`); `fbsp`:
  `c(m, bandwidth, centre frequency)` (default `c(2, 2, 1.5)`);
  `loggabor`: `c(centre angular frequency, log-width)` (default
  `c(6, 0.25)`); `meyer` and `hhat` have no parameter; `morlet0`:
  omega_0 (default: the local number of cycles); `morl`: omega_0
  (default 5); `poisson`: order n (default 2); `gammatone`:
  `c(n, centre angular frequency, bandwidth)` (default `c(4, 6, 1.5)`);
  `hann`: `c(centre angular frequency, half width)` (default `c(6, 3)`);
  `db`: number of vanishing moments 1 to 20 (default 4, `"haar"` =
  `"db1"`); `sym`: 2 to 20 (default 8); `battle`: spline order 1, 3, 5
  or 7 (default 3). For the Morlet wavelet use `omega_min` and
  `omega_max`; a non-NULL `wavelet_order` uses one fixed wavelet at all
  scales. `Default=NULL`.

## Value

A list with class `"analyze.Awavelet"` containing:

- Wave: complex wavelet coefficients

- Phase: instantaneous phase

- Ampl: wavelet amplitude

- Power: wavelet power spectrum

- dt: sampling interval

- dj: scale spacing

- Power.avg: average power per period

- Period: physical periods

- Scale: wavelet scales

- coi.1: cone of influence x coordinates

- coi.2: cone of influence y coordinates

- nc: number of samples

- nr: number of scales

- axis.1: x axis values

- axis.2: log2 scaled periods

- omega_min: minimum wavelet cycles

- omega_max: maximum wavelet cycles

- omega: scale dependent wavelet cycles

- scaling: scaling method used

- alpha: power law exponent

- x: interpolated x values

- y: interpolated signal values

## References

Torrence, C., and G. P. Compo (1998). A practical guide to wavelet
analysis. Bulletin of the American Meteorological Society, 79, 61–78.

Morlet, J., G. Arens, E. Fourgeau, and D. Giard (1982). Wave propagation
and sampling theory. Geophysics, 47, 203–236.

## Author

Adapted from the WaveletComp and biwavelet R packages, which are based
on the MATLAB wavelet code by Torrence and Compo, with extensions for
scale dependent wavelet width.

## Examples

``` r
# \donttest{
#Example 1. Using the Total Solar Irradiance data set of Steinhilber et al., (2012)
TSI_wt <-
 analyse_Awavelet(
   data = TSI,
   dj = 1/200,
   lowerPeriod = 16,
   upperPeriod = 8192,
   verbose = FALSE,
omega_min = 6,
omega_max = 12,
scaling = "log2",
alpha = 1
 )


#Example 2. Using the magnetic susceptibility data set of Pas et al., (2018)
mag_wt <-
analyse_Awavelet(
data = mag,
dj = 1/100,
lowerPeriod = 0.1,
upperPeriod = 254,
verbose = FALSE,
omega_min = 6,
omega_max = 12,
scaling = "log2",
alpha = 1
)

#Example 3. Using the greyscale data set of Zeeden et al., (2013)
grey_wt <-
 analyse_Awavelet(
   data = grey,
   dj = 1/200,
   lowerPeriod = 0.02,
   upperPeriod = 256,
   verbose = FALSE,
omega_min = 6,
omega_max = 12,
scaling = "log2",
alpha = 1,
 )

# }
```
