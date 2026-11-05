# Conduct the continuous wavelet transform on a time series/signal

Compute the continuous wavelet transform (CWT) using a Morlet wavelet

## Usage

``` r
analyse_wavelet(
  data = NULL,
  dj = 1/100,
  lowerPeriod = 2,
  upperPeriod = 1024,
  verbose = FALSE,
  omega_nr = 8,
  wavelet = "morlet",
  wavelet_order = NULL,
  pval = FALSE,
  n_simulations = 10,
  run_multicore = FALSE,
  n_cores = NULL
)

analyze_wavelet(
  data = NULL,
  dj = 1/100,
  lowerPeriod = 2,
  upperPeriod = 1024,
  verbose = FALSE,
  omega_nr = 8,
  wavelet = "morlet",
  wavelet_order = NULL,
  pval = FALSE,
  n_simulations = 10,
  run_multicore = FALSE,
  n_cores = NULL
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

  Spacing between successive scales in octaves (numeric, \> 0). The
  periods analysed increase by powers of 2 (1, 2, 4, 8, ...), and `dj`
  sets how many steps are made per doubling of the period: `dj=1/x`
  gives x periods per octave. Smaller values give a finer period axis
  but increase the computation time. `Default=1/100`.

- lowerPeriod:

  Shortest period to analyse, in the same unit as the depth/time column
  of `data` (numeric, \> 0; at least twice the sampling interval).
  Together with `upperPeriod` it sets the range of cycles analysed;
  because the period axis is log2-scaled, a power of 2 gives the neatest
  plots. `Default=2`.

- upperPeriod:

  Longest period to analyse, in the same unit as the depth/time column
  of `data` (numeric, \> `lowerPeriod`; normally not more than the
  length of the record); a power of 2 gives the neatest plots.
  `Default=1024`.

- verbose:

  Print the sampling interval to which the record was interpolated and,
  when `pval=TRUE`, show a progress bar of the simulations (logical,
  TRUE or FALSE). `Default=FALSE`.

- omega_nr:

  Number of cycles of the Morlet wavelet, i.e. its central angular
  frequency omega_0 (numeric, \> 0, typically 5 to 15). A larger value
  gives a longer wavelet with more oscillations: a better period
  resolution and a worse depth/time resolution. Also used as the default
  parameter of `"morlet0"` and to match the bandwidth of some other
  wavelets. `Default=8`.

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
  or 7 (default 3). For the Morlet wavelet use `omega_nr`. The period
  axis is always the Fourier period, i.e. the period of the sinusoid at
  which the wavelet has its maximum response, so the results of
  different wavelets are comparable. `Default=NULL`.

- pval:

  Calculate p-values of the wavelet power by Monte Carlo simulation
  (logical, TRUE or FALSE). `Default=FALSE`. Each simulation analyses a
  surrogate record made of the average of a shuffled copy of the data
  (white noise) and an AR(1) red-noise series with the mean and standard
  deviation of the data and a lag-1 autocorrelation drawn around the
  value estimated in a moving window of half the record length. The
  shuffled part keeps high-amplitude short-period cycles from becoming
  significant; the red-noise part gives a spectrum similar to that of
  the data without being biased by possible astronomical peaks in it.
  The p-value is the fraction of simulations whose power is at least the
  observed power (returned as `Power.pval` and `Power.avg.pval`).

- n_simulations:

  Number of Monte Carlo simulations used for the p-values (integer, \>
  0); only used when `pval=TRUE`. More simulations give more precise
  p-values but take longer. `Default=10`.

- run_multicore:

  Run the Monte Carlo simulations in parallel on several cores (logical,
  TRUE or FALSE); only used when `pval=TRUE`. `Default=FALSE`.

- n_cores:

  Number of cores used when `run_multicore=TRUE` (integer, \> 0); `NULL`
  picks an automatic, safe number: physical cores minus one, at most 4,
  and limited by the free memory. Set `options(waverider.n_cores = ...)`
  to change that default. A value \> 1 also switches on parallel
  computing. Results do not depend on the number of cores.
  `Default=NULL`.

## Value

The output is a list (wavelet object) which contain 20 objects which are
the result of the continuous wavelet transform (CWT). Object 1: Wave -
Wave values of the wavelet Object 2: Phase - Phase of the wavelet Object
3: Ampl - Amplitude values of the wavelet Object 4: Power - Power values
of the wavelet Object 5: dt - Step size Object 6: dj - Scale size Object
7: Power.avg - Average power values Object 8: Period - Period values
Object 9: Scale - Scale value Object 10: coi.1 - Cone of influence
values 1 Object 11: coi.2 - Cone of influence values 2 Object 12: nc -
Number of columns Object 13: nr - Number of rows Object 14: axis.1 -
axis values 1 Object 15: axis.2 - axis values 2 Object 16: omega_nr -
Number of cycles in the wavelet Object 17: x - x values of the data set
Object 18: y - y values of the data set Object 19: average p value of
the spectral power Object 20: p value of spectral power

## References

Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
Wavelet Analysis. R package version 1.1.
<https://CRAN.R-project.org/package=WaveletComp>

Gouhier TC, Grinsted A, Simko V (2021). R package biwavelet: Conduct
Univariate and Bivariate Wavelet Analyses. (Version 0.20.21),
<https://github.com/tgouhier/biwavelet>

Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet
Analysis. Bulletin of the American Meteorological Society 79:61-78.
<https://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf>

Morlet, Jean, Georges Arens, Eliane Fourgeau, and Dominique Glard. "Wave
propagation and sampling theory—Part I: Complex signal and scattering in
multilayered media. " Geophysics 47, no. 2 (1982): 203-221.

J. Morlet, G. Arens, E. Fourgeau, D. Giard; Wave propagation and
sampling theory; Part II, Sampling theory and complex waves. Geophysics
1982 47 (2): 222–236.

## Author

Code based on on the "WaveletComp" function of the 'WaveletComp' R
package and "wt" function of the 'biwavelet' R package which are based
on the wavelet MATLAB code written by Christopher Torrence and Gibert P.
Compo.

## Examples

``` r
# \donttest{
#Example 1. Using the Total Solar Irradiance data set of Steinhilver et al., (2012)
TSI_wt <-
 analyse_wavelet(
   data = TSI,
   dj = 1/200,
   lowerPeriod = 16,
   upperPeriod = 8192,
   verbose = FALSE,
   omega_nr = 6,
   pval=FALSE,
   n_simulations=10,
   run_multicore = FALSE
 )


#Example 2. Using the magnetic susceptibility data set of Pas et al., (2018)
mag_wt <-
analyse_wavelet(
data = mag,
dj = 1/100,
lowerPeriod = 0.1,
upperPeriod = 254,
verbose = FALSE,
omega_nr = 10,
pval=FALSE,
n_simulations=10,
run_multicore = FALSE
)

#Example 3. Using the greyscale data set of Zeeden et al., (2013)
grey_wt <-
 analyse_wavelet(
   data = grey,
   dj = 1/200,
   lowerPeriod = 0.02,
   upperPeriod = 256,
   verbose = FALSE,
   omega_nr = 8,
   pval=FALSE,
   n_simulations=10,
   run_multicore = FALSE
 )

# }
```
