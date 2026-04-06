# WaverideR

WaverideR is an R package for advanced cyclostratigraphic analysis of stratigraphic data sets.
It provides a comprehensive suite of spectral tools for detecting, visualising, and tracking
non-stationary astronomical cycles, including continuous wavelet transforms (CWT),
the superlet transform, windowed FFT, and evolutionary harmonic analysis (EHA) (wrapper).
These methods allow both manual and automated tracking of orbital cycles in spectra and
scalograms, even in records affected by large changes in sedimentation rate.Building on
this spectral framework, WaverideR supports multi-proxy integration and Monte Carlo–based
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

## Geological Annotation

This repository now includes a note on cyclostratigraphic signals relevant to the WaverideR data. The 405 kyr eccentricity cycle, a fundamental component of the Milankovitch orbital forcing, is recognised as the longest and most stable eccentricity term in the Earth's pre‑cessional climate record. Recent cyclostratigraphic studies of the Apennine sections have linked this 405 kyr signal to the timing of the Cretaceous–Paleogene (K‑Pg) boundary, building on the pioneering work of Walter Alvarez and collaborators, who identified the iridium‑rich layer marking the mass‑extinction event. Integrating the 405 kyr eccentricity modulation with the Alvarez K‑Pg boundary chronostratigraphy helps refine the astronomical calibration of these key intervals.

