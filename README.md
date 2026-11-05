# WaverideR

Extracts and analyses cyclic signals in stratigraphic and other time series using the continuous
wavelet transform, superlets and synchrosqueezed transforms: spectra, significance tests,
sedimentation-rate estimation, tracking of astronomical cycles and periods, age models with
uncertainty and de-chirped spectra, together with the astrochronology functions of acycle2R and the
Bayesian age-model code of astroMCBchron.

Website with the manual, the function reference and the tutorial: <https://stratigraphy.github.io/WaverideR/>

Version 0.6.0 (see `NEWS.md` for all changes).

## Installation

```r
# development version
# install.packages("remotes")
remotes::install_github("stratigraphy/WaverideR")
```

## A short example

```r
library(WaverideR)

wt <- analyse_wavelet(data = mag[1:300, ], dj = 1/12, lowerPeriod = 0.2, upperPeriod = 32,
                      omega_nr = 8, verbose = FALSE)
sig <- wavelet_significance(wt)                       # red-noise significance (Torrence & Compo)
plot_wavelet(wt, add_avg = TRUE, add_signif = sig)

# significance of the power in a period band, valid for any omega_nr
wavelet_band_significance(wt, period_range = c(0.8, 1.6))
```

`vignette("WaverideR-intro")` walks through the workflow.

## What is in the package

* **Transforms**: `analyse_wavelet()`, `analyse_superlet()`, cross spectra and coherence
  (`analyse_Xwavelet()`, `analyse_Xsuperlet()`, `analyse_wavelet_coherence()`), evolutive harmonic
  analysis (`eha_log2()`), and multitaper and bispectral tools.
* **Significance**: `wavelet_significance()` (analytic and surrogate), `wavelet_band_significance()`
  (scale-averaged, any `omega_nr`), `wavelet_uncertainty()`.
* **Sedimentation rate**: `flmw()`, `flmw_path()`, `sum_power_sedrate()`, `model_red_noise_wt()`, COCO and eCOCO.
* **Tying to astronomical solutions**: `astro_anchor()` and the interactive, random-order `astro_anchor_interactive()` (tie, untie, tune and inspect), `plot_astro_anchor()`, `anchor2time()`.
* **Tracking and age models**: `track_period()`, `track_ecoco()`, `curve2time()`, `curve2time_unc()`,
  `astroMCBchron()`, `modBchron()`.
* **Plots**: all scalograms share one set of panel functions, so layouts, palettes and options behave the same.

## Parallel computing

Everything is sequential by default. `run_multicore = TRUE` and `n_cores` switch on base-R `parallel`
(forked workers on Linux and macOS, PSOCK on Windows and in RStudio). The number of cores is capped by
the physical cores, the free memory and the CRAN limit, and results do not depend on it.

## Citation

When using WaverideR please cite both the paper that introduced the package and the package itself, for example
"... using the WaverideR R package (Arts et al., 2024; Arts, 2023)":

* Arts, M., Corradini, C., Pondrelli, M., Pas, D., Da Silva, A.-C. (2024). Age and orbital forcing in the upper
  Silurian Cellon section (Carnic Alps, Austria) uncovered using the WaverideR R package. *Frontiers in Earth
  Science* 12, 1357751. https://doi.org/10.3389/feart.2024.1357751
* Arts, M.C.M. (2023). WaverideR: Extracting Signals from Wavelet Spectra. R package.
  https://CRAN.R-project.org/package=WaverideR

`citation("WaverideR")` prints both entries.
