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

## Cyclostratigraphic Context

WaverideR is particularly relevant for identifying Milankovitch cycles in stratigraphic records.
The **405 kyr eccentricity cycle**, driven by the gravitational interaction of Jupiter and Venus
(Laskar et al., 2004, 2011), serves as the fundamental metronome for astrochronological
calibration due to its remarkable stability over hundreds of millions of years. This long-
eccentricity period is widely used as a tuning target in cyclostratigraphic studies because it
remains nearly constant even as other orbital parameters undergo secular variation.

Walter Alvarez and colleagues demonstrated the application of cyclostratigraphy to the
Cretaceous–Paleogene (K–Pg) boundary in the pelagic limestone sections of the Italian Apennines
(e.g., Bottaccione Gorge, Contessa Highway). Their work showed that orbitally tuned cyclostratigraphy
could precisely constrain the timing and duration of events surrounding the K–Pg mass extinction
(Alvarez et al., 1977; Alvarez & Lowrie, 1978; Lowrie & Alvarez, 1981). These Apennine sections
preserve rhythmic bedding that records Milankovitch-scale orbital forcing, enabling the detection
of precession, obliquity, and eccentricity cycles—including the 405 kyr long eccentricity—which
provides a robust chronostratigraphic framework across the K–Pg boundary interval.

By enabling the detection and tracking of such orbital cycles, WaverideR facilitates the
construction of high-resolution astrochronologies that directly build upon the foundational
cyclostratigraphic work of Alvarez and others in the Apennine carbonate sequences.

## Installation

You can install the development version of WaverideR from [GitHub](https://github.com/stratigraphy/WaverideR) with:

``` r
# install.packages("devtools")
devtools::install_github("stratigraphy/WaverideR")
```
