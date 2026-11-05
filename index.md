# WaverideR

Extracts and analyses cyclic signals in stratigraphic and other time
series using the continuous wavelet transform, superlets and
synchrosqueezed transforms: spectra, significance tests,
sedimentation-rate estimation, tracking of astronomical cycles and
periods, age models with uncertainty and de-chirped spectra, together
with the astrochronology functions of acycle2R and the Bayesian
age-model code of astroMCBchron.

Website with the manual, the function reference and the tutorial:
<https://stratigraphy.github.io/WaverideR/>

Version 0.6.0 (see `NEWS.md` for all changes).

## Installation

``` r

# development version
# install.packages("remotes")
remotes::install_github("stratigraphy/WaverideR")
```

## A short example

``` r

library(WaverideR)

wt <- analyse_wavelet(data = mag[1:300, ], dj = 1/12, lowerPeriod = 0.2, upperPeriod = 32,
                      omega_nr = 8, verbose = FALSE)
sig <- wavelet_significance(wt)                       # red-noise significance (Torrence & Compo)
plot_wavelet(wt, add_avg = TRUE, add_signif = sig)

# significance of the power in a period band, valid for any omega_nr
wavelet_band_significance(wt, period_range = c(0.8, 1.6))
```

[`vignette("WaverideR-intro")`](https://stratigraphy.github.io/WaverideR/articles/WaverideR-intro.md)
walks through the workflow.

## What is in the package

- **Transforms**:
  [`analyse_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md),
  [`analyse_superlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md),
  cross spectra and coherence
  ([`analyse_Xwavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_Xwavelet.md),
  [`analyse_Xsuperlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_Xsuperlet.md),
  [`analyse_wavelet_coherence()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet_coherence.md)),
  evolutive harmonic analysis
  ([`eha_log2()`](https://stratigraphy.github.io/WaverideR/reference/eha_log2.md)),
  and multitaper and bispectral tools.
- **Significance**:
  [`wavelet_significance()`](https://stratigraphy.github.io/WaverideR/reference/wavelet_significance.md)
  (analytic and surrogate),
  [`wavelet_band_significance()`](https://stratigraphy.github.io/WaverideR/reference/wavelet_band_significance.md)
  (scale-averaged, any `omega_nr`),
  [`wavelet_uncertainty()`](https://stratigraphy.github.io/WaverideR/reference/wavelet_uncertainty.md).
- **Sedimentation rate**:
  [`flmw()`](https://stratigraphy.github.io/WaverideR/reference/flmw.md),
  [`flmw_path()`](https://stratigraphy.github.io/WaverideR/reference/flmw_path.md),
  [`sum_power_sedrate()`](https://stratigraphy.github.io/WaverideR/reference/sum_power_sedrate.md),
  [`model_red_noise_wt()`](https://stratigraphy.github.io/WaverideR/reference/model_red_noise_wt.md),
  COCO and eCOCO.
- **Tying to astronomical solutions**:
  [`astro_anchor()`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md)
  and the interactive, random-order
  [`astro_anchor_interactive()`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor_interactive.md)
  (tie, untie, tune and inspect),
  [`plot_astro_anchor()`](https://stratigraphy.github.io/WaverideR/reference/plot_astro_anchor.md),
  [`anchor2time()`](https://stratigraphy.github.io/WaverideR/reference/anchor2time.md).
- **Tracking and age models**:
  [`track_period()`](https://stratigraphy.github.io/WaverideR/reference/track_period.md),
  [`track_ecoco()`](https://stratigraphy.github.io/WaverideR/reference/track_ecoco.md),
  [`curve2time()`](https://stratigraphy.github.io/WaverideR/reference/curve2time.md),
  [`curve2time_unc()`](https://stratigraphy.github.io/WaverideR/reference/curve2time_unc.md),
  [`astroMCBchron()`](https://stratigraphy.github.io/WaverideR/reference/astroMCBchron.md),
  [`modBchron()`](https://stratigraphy.github.io/WaverideR/reference/modBchron.md).
- **Plots**: all scalograms share one set of panel functions, so
  layouts, palettes and options behave the same.

## Parallel computing

Everything is sequential by default. `run_multicore = TRUE` and
`n_cores` switch on base-R `parallel` (forked workers on Linux and
macOS, PSOCK on Windows and in RStudio). The number of cores is capped
by the physical cores, the free memory and the CRAN limit, and results
do not depend on it.

## Citation

When using WaverideR please cite both the paper that introduced the
package and the package itself, for example “… using the WaverideR R
package (Arts et al., 2024; Arts, 2023)”:

- Arts, M., Corradini, C., Pondrelli, M., Pas, D., Da Silva, A.-C.
  (2024). Age and orbital forcing in the upper Silurian Cellon section
  (Carnic Alps, Austria) uncovered using the WaverideR R package.
  *Frontiers in Earth Science* 12, 1357751.
  <https://doi.org/10.3389/feart.2024.1357751>
- Arts, M.C.M. (2023). WaverideR: Extracting Signals from Wavelet
  Spectra. R package. <https://CRAN.R-project.org/package=WaverideR>

`citation("WaverideR")` prints both entries.
