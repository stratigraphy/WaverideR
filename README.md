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

## Geological Annotation: Milankovitch Cyclostratigraphy and the K-Pg Boundary

Cyclostratigraphic analysis using WaverideR is particularly relevant to the study of
**Milankovitch cycles**—periodic variations in Earth's orbital parameters that modulate
incoming solar radiation and drive climate cycles recorded in sedimentary sequences. Among
these, the **405 kyr eccentricity cycle** stands out as the most stable and predictable
Milankovitch periodicity. Governed by the gravitational interaction between Jupiter and
Venus, the 405 kyr long-eccentricity cycle persists with minimal variation over hundreds
of millions of years and serves as the primary metronome for constructing astrochronological
age models (Laskar et al., 2004, 2011).

A landmark application of cyclostratigraphy in this context is the study of the
**Cretaceous–Paleogene (K-Pg) boundary** in the Apennine sections of Italy. **Walter
Alvarez** and colleagues first identified the now-famous iridium anomaly at the
Bottaccione Gorge near Gubbio (Umbria–Marche Apennines), providing the pivotal evidence
for a bolide impact as the cause of the end-Cretaceous mass extinction (Alvarez et al.,
1980). Subsequent cyclostratigraphic work on these same Apennine sections has employed
the 405 kyr eccentricity cycle as a tuning target to refine the chronology of the
K-Pg boundary interval. Kuiper et al. (2008) used 405 kyr eccentricity-cycle tuning of
the Scaglia Rossa limestone at Bottaccione to constrain the K-Pg boundary age to
~66.0 Ma, demonstrating how orbital cyclicity preserved in deep-marine carbonate
sequences can anchor radioisotopic ages with unprecedented precision.

WaverideR's suite of spectral tools—including wavelet transforms, evolutionary harmonic
analysis, and automated cycle tracking—enables researchers to detect and trace the 405 kyr
eccentricity signal even in records with significant changes in sedimentation rate, making
it an ideal platform for revisiting classic Apennine sections and other cyclostratigraphic
archives affected by the K-Pg event.

### Key References

- Alvarez, L.W., Alvarez, W., Asaro, F., and Michel, H.V. (1980). Extraterrestrial cause
  for the Cretaceous-Tertiary extinction. *Science*, 208(4448), 1095–1108.
  doi:10.1126/science.208.4448.1095
- Laskar, J., Robutel, P., Joutel, F., Gastineau, M., Correia, A.C.M., and Levrard, B.
  (2004). A long-term numerical solution for the insolation quantities of the Earth.
  *Astronomy & Astrophysics*, 428, 261–285. doi:10.1051/0004-6361:20041335
- Laskar, J., Fienga, A., Gastineau, M., and Manche, H. (2011). La2010: a new orbital
  solution for the long-term motion of the Earth. *Astronomy & Astrophysics*, 532, A89.
  doi:10.1051/0004-6361/201116836
- Kuiper, K.F., Deino, A., Hilgen, F.J., Krijgsman, W., Renne, P.R., and Wijbrans, J.R.
  (2008). Synchronizing rock clocks of Earth history. *Science*, 320(5875), 500–504.
  doi:10.1126/science.1154339

## Installation

You can install the development version of WaverideR from [GitHub](https://github.com/stratigraphy/WaverideR) with:

``` r
# install.packages("devtools")
devtools::install_github("stratigraphy/WaverideR")
```
