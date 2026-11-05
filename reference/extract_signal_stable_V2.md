# Extract signal from a wavelet spectrum using a upper and lower period boundary

Extract a signal from the wavelet using a upper and lower period
boundary

## Usage

``` r
extract_signal_stable_V2(
  wavelet = NULL,
  period_max = NULL,
  period_min = NULL,
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

- period_max:

  Upper period boundary of the band that is reconstructed, in the
  depth/time unit of the wavelet (numeric, \> `period_min`; the nearest
  wavelet period is used). Required (`Default=NULL`).

- period_min:

  Lower period boundary of the band that is reconstructed, in the
  depth/time unit of the wavelet (numeric, \> 0; the nearest wavelet
  period is used). Required (`Default=NULL`).

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

Signal extracted from the wavelet spectra. Output is a matrix with the
first column being depth/time and the second column is the cycle
extracted from the proxy record.

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
# Irradiance data set of Steinhilber et al., (2012)

TSI_wt <-
analyse_wavelet(
data = TSI,
dj = 1/200,
lowerPeriod = 16,
upperPeriod = 8192,
   verbose = FALSE,
   omega_nr = 6
 )

de_Vries_cycle <- extract_signal_stable_V2(wavelet=TSI_wt,
period_max = 240,
period_min = 180,
add_mean=TRUE,
plot_residual=FALSE,
keep_editable=FALSE)

```
