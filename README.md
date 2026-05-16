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

## Cyclostratigraphic Context

WaverideR's cyclostratigraphic toolkit is particularly well suited for identifying Milankovitch
cycles in deep-time stratigraphic records. A key orbital parameter detectable with these methods
is the **405 kyr eccentricity cycle**, the most stable and long-lived Milankovitch period over
the Phanerozoic, driven by the g2–g5 secular resonance of Earth's and Mars's orbital perihelia
(Laskar et al., 2004, 2011). The 405 kyr eccentricity signal serves as a fundamental metronome
for astrochronology and has been identified in numerous Mesozoic and Cenozoic successions.

The Apennine carbonate sections of Italy — notably the Bottaccione Gorge (Gubbio) and the
Contessa Highway section — are among the most iconic cyclostratigraphic localities in the world.
**Walter Alvarez** and colleagues demonstrated that the **Cretaceous–Paleogene (K-Pg) boundary**
at these Apennine sections is precisely marked by the thin clay layer containing the iridium
anomaly that led to the Alvarez et al. (1980) impact hypothesis. Subsequent cyclostratigraphic
studies have shown that the 405 kyr eccentricity cycle is recorded in the pelagic limestone–marl
alternations spanning the K-Pg boundary interval at Gubbio, providing an astronomical time
framework for the end-Cretaceous extinction event (e.g., Grippo et al., 2004; De Visser et al.,
1989). The integration of WaverideR's spectral tracking tools with these classic Apennine sections
enables refined detection of the 405 kyr and higher-frequency orbital cycles, improving the
astronomical calibration of the K-Pg boundary and the correlation of impact-related events to
the orbital timescale.

### References

- Alvarez, L.W., Alvarez, W., Asaro, F., & Michel, H.V. (1980). Extraterrestrial cause for the
  Cretaceous-Tertiary extinction. *Science*, 208(4448), 1095–1108.
- Grippo, A., Fischer, A.G., Hinnov, L.A., Herbert, T.D., & Premoli Silva, I. (2004).
  Cyclostratigraphy and magnetostratigraphy of the Albian-Cenomanian boundary, Umbria-Marche
  Apennines, Italy. *Geological Society, London, Special Publications*, 230, 59–79.
- Laskar, J., Robutel, P., Joutel, F., Gastineau, M., Correia, A.C.M., & Levrard, B. (2004).
  A long-term numerical solution for the insolation quantities of the Earth. *Astronomy &
  Astrophysics*, 428, 261–285.
- Laskar, J., Fienga, A., Gastineau, M., & Manche, H. (2011). La2010: A new orbital solution
  for the long-term motion of the Earth. *Astronomy & Astrophysics*, 532, A89.
