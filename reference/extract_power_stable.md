# Extract power from a wavelet spectra by using a constant period/duration

Extract spectral power from the wavelet using a constant period/duration
and boundaries as selection criteria. The extraction of spectral is
useful for cyclostratigraphic studies because the spectral power of an
astronomical cycle is modulated by higher order astronomical cycles. The
spectral power record from an astronomical cycle can thus be used as a
proxy for amplitude modulating cycles. The spectral power is extracted
from a wavelet spectra which was created using the
[`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
function for a given, `cycle`, `period_up` and `period_down`

## Usage

``` r
extract_power_stable(
  scalogram = NULL,
  cycle = NULL,
  period_up = 1.2,
  period_down = 0.8,
  ...
)
```

## Arguments

- scalogram:

  Wavelet or superlet object (list of class `analyze.wavelet` or
  `analyze.superlet`) created with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  or
  [`analyse_superlet`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md),
  whose power is summed. Required (`Default=NULL`).

- cycle:

  Constant period of the cycle whose power is extracted, in the
  depth/time unit of the scalogram (numeric, \> 0). Required
  (`Default=NULL`).

- period_up:

  Upper boundary of the band whose power is summed, as a factor of
  `cycle` (numeric, \> `period_down`). `Default=1.2`.

- period_down:

  Lower boundary of the band whose power is summed, as a factor of
  `cycle` (numeric, \> 0). `Default=0.8`.

- ...:

  Old argument names of earlier versions: `wavelet` (used as
  `scalogram`).

## Value

Returns a matrix with 3 columns. The first column is depth/time. The
second column is extracted power. The third column is extracted
power/total power.

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

Routines for astrochronologic testing, astronomical time scale
construction, and time series analysis
\<doi:10.1016/j.earscirev.2018.11.015\>

## Author

Code based on the reconstruct function of the 'WaveletComp' R package
which is based on the wavelet 'MATLAB' code written by Christopher
Torrence and Gibert P. Compo (1998). The functionality of this function
is is inspired by the
[integratePower](https://rdrr.io/pkg/astrochron/man/integratePower.html)
function of the 'astrochron' R package

## Examples

``` r
#Extract the spectral power of the 210 yr de Vries cycle from the Total Solar
#Irradiance data set of Steinhilber et al., (2012).

TSI_wt <-
 analyse_wavelet(
   data = TSI,
   dj = 1/200,
   lowerPeriod = 16,
   upperPeriod = 8192,
   verbose = FALSE,
   omega_nr = 6
 )
TSI_wt_pwr_de_Vries_cycle <-  extract_power_stable(
 scalogram = TSI_wt,
 cycle = 210,
 period_up = 1.2,
 period_down = 0.8
)

```
