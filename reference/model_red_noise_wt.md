# Models average spectral power based curves based on a red-noise signal generated using the characteristics of an input signal.

The `model_red_noise_wt` function is used to generate average spectral
power curves based on and input signal and set wavelet settings.

## Usage

``` r
model_red_noise_wt(
  wavelet = NULL,
  n_simulations = NULL,
  run_multicore = FALSE,
  n_cores = NULL,
  verbose = FALSE
)
```

## Arguments

- wavelet:

  Wavelet object, the output of
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md).
  Its record (`x`, `y`), period range, `dj` and `omega_nr` are used:
  every simulation is red (AR1) noise with the mean, standard deviation,
  length and (randomly perturbed) lag-1 autocorrelation of the record,
  analysed with the same wavelet settings. Required (`Default=NULL`).

- n_simulations:

  Number of red-noise simulations (integer, \> 0; e.g. 1000 for stable
  percentiles); one output column per simulation. Required
  (`Default=NULL`).

- run_multicore:

  Run the simulations in parallel on several cores (logical, TRUE or
  FALSE). Results do not depend on it. `Default=FALSE`.

- n_cores:

  Number of cores used when `run_multicore = TRUE` (integer or `NULL`).
  `NULL` gives an automatic, safe number (physical cores minus one, at
  most 4, limited by the free memory); change this default with
  `options(waverider.n_cores = ...)`. A value \> 1 also switches on
  parallel computing. Results do not depend on the number of cores.
  `Default=NULL`.

- verbose:

  Print the interpolation interval of the record and show a progress bar
  of the simulations (logical, TRUE or FALSE). `Default=FALSE`.

## Value

Returns a matrix in which each column represents the average spectral
power resulting from a red-noise run.

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

Code based on the "analyze.wavelet" function of the 'WaveletComp' R
package and "wt" function of the 'biwavelet' R package which are based
on the wavelet 'MATLAB' code written by Christopher Torrence and Gibert
P. Compo (1998).

## Examples

``` r
# \donttest{
#'#generate average spectral power curves based on red noise curves which are
# based on the magnetic susceptibility record of the Sullivan core of Pas et al., (2018)

mag_wt <- analyse_wavelet(data = mag,
dj = 1/100,
lowerPeriod = 0.1,
upperPeriod = 254,
verbose = FALSE,
omega_nr = 10)

#increase n_simulations to better define the red noise spectral power curve
mag_wt_red_noise <- model_red_noise_wt(wavelet=mag_wt,
n_simulations=10, # increase number for better constrained results
run_multicore=FALSE,
verbose=FALSE)
# }

```
