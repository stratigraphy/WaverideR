# Extract a signal/cycle from a wavelet spectra using a set period and boundaries

Extracts a cycle from the wavelet object created using the
[`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
function using a fixed period and fixed period boundaries defined as
factors of the original cycle

## Usage

``` r
extract_signal_stable(
  wavelet = NULL,
  cycle = NULL,
  period_up = 1.2,
  period_down = 0.8,
  add_mean = TRUE,
  plot_residual = FALSE,
  keep_editable = FALSE
)
```

## Arguments

- wavelet:

  Wavelet object (list) created with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  the signal is reconstructed from its coefficients. Required
  (`Default=NULL`).

- cycle:

  Constant period of the cycle to extract, in the depth/time unit of the
  wavelet (numeric, \> 0). Required (`Default=NULL`).

- period_up:

  Upper boundary of the band that is reconstructed, as a factor of
  `cycle` (numeric, \> `period_down`). `Default=1.2`.

- period_down:

  Lower boundary of the band that is reconstructed, as a factor of
  `cycle` (numeric, \> 0). `Default=0.8`.

- add_mean:

  Add the mean of the original data series to the extracted cycle
  (logical, TRUE or FALSE); `FALSE` returns it centred on zero.
  `Default=TRUE`.

- plot_residual:

  Plot the residual (extracted cycle minus data) against depth/time and
  its histogram (logical, TRUE or FALSE). `Default=FALSE`.

- keep_editable:

  Keep the plot layout active after the call so that extra features can
  be added to the residual plot (logical, TRUE or FALSE); `FALSE`
  restores the previous graphical parameters. Only used when
  `plot_residual=TRUE`. `Default=FALSE`.

## Value

\#'Returns a matrix with 2 columns. The first column is time/depth. The
second column is the extracted signal/cycle.

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

## Author

Code based on the reconstruct function of the 'WaveletComp' R package
which is based on the wavelet 'MATLAB' code written by Christopher
Torrence and Gibert P. Compo (1998).

## Examples

``` r
#Example in which the ~210yr de Vries cycle is extracted from the Total Solar
#Irradiance data set of Steinhilber et al., (2012)

#Perform the CWT
TSI_wt <-
analyse_wavelet(
data = TSI,
dj = 1/200,
lowerPeriod = 16,
upperPeriod = 8192,
   verbose = FALSE,
   omega_nr = 6
 )

#Extract the 210 yr de Vries cycle from the wavelet spectra
de_Vries_cycle <- extract_signal_stable(wavelet=TSI_wt,
cycle=210,
period_up =1.25,
period_down = 0.75,
add_mean=TRUE,
plot_residual=FALSE,
keep_editable=FALSE)
```
