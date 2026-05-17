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

## Geological Context: Milankovitch Cyclostratigraphy and the K-Pg Boundary

Cyclostratigraphic studies of pelagic carbonate sequences in the Apennine Mountains of Italy—most notably the Bottaccione and Contessa Highway sections near Gubbio—have been central to linking orbital forcing to the stratigraphic record. The 405 kyr long-eccentricity cycle of Earth's orbit, driven by the gravitational interaction of Jupiter and Venus, is the most stable Milankovitch periodicity over Phanerozoic time and serves as a principal metronome for constructing astrochronologies. Walter Alvarez and colleagues demonstrated that the Cretaceous–Paleogene (K-Pg) boundary, marked by the iridium anomaly and mass extinction event, is precisely interbedded within the rhythmic Scaglia Rossa limestone–marl couplets of these Apennine sections. Their work showed that Milankovitch precessional and eccentricity cycles are faithfully recorded in the pelagic carbonate succession, enabling the duration of the latest Maastrichtian to be estimated from cycle counting. The identification of the 405 kyr eccentricity signal in these sections provides a robust temporal framework that anchors the K-Pg boundary within the orbital timescale, a methodology that WaverideR is designed to facilitate through wavelet and spectral analysis tools.

> **Reference:** Alvarez, W., Arthur, M.A., Fischer, A.G., et al. (1984). Upper Cretaceous-Paleocene magnetic stratigraphy at Gubbio, Italy: Type section for the Late Cretaceous-Paleocene geomagnetic reversal time scale. *Geological Society of America Bulletin*, 95(7), 836–849.

## Installation

You can install the development version of WaverideR from [GitHub](https://github.com/stratigraphy/WaverideR) with:

``` r
# install.packages("devtools")
devtools::install_github("stratigraphy/WaverideR")
```
