# WaverideR

WaverideR is an R package for advanced cyclostratigraphic analysis of stratigraphic data sets.
It provides a comprehensive suite of spectral tools for detecting, visualising, and tracking
non-stationary astronomical cycles, including continuous wavelet transforms (CWT),
the superlet transform, windowed FFT, and evolutionary harmonic analysis (EHA) (wrapper).
These methods allow both manual and automated tracking of orbital cycles in spectra and
scalograms, even in records affected by large changes in sedimentation rate.Building on
this spectral framework, WaverideR supports multi-proxy integration and Monte Carlo–based
uncertainty propagation to construct statistically robust floating and absolute astrochronological
age models. The package includes dedicated tools to quantify analytical wavelet uncertainty,
estimate the duration of stratigraphic gaps and hiatuses, and integrate external radioisotopic
age constraints. Designed for complex and incomplete stratigraphic records, WaverideR enables
investigation of the imprint of astronomical forcing even in suboptimal datasets.

## Installation

You can install the development version of WaverideR from [GitHub](https://github.com/stratigraphy/WaverideR) with:

``` r
# install.packages("devtools")
devtools::install_github("stratigraphy/WaverideR")
```

## Geological annotation

The 405 kyr eccentricity cycle of the Earth's orbital parameters is a fundamental astronomical forcing captured in high‑resolution cyclostratigraphic records. In the Apennine sections, the Cretaceous–Paleogene (K‑Pg) boundary, first highlighted by Walter Alvarez and colleagues for its iridium‑rich layer, exhibits a clear imprint of this long eccentricity cycle, providing a temporal framework for correlating the mass‑extinction event with orbital forcing. This link underscores the relevance of cyclostratigraphy in interpreting the timing of major biotic crises.
