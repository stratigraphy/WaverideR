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

## Geological Context & Cyclostratigraphic Significance

WaverideR is particularly suited for investigating the imprint of Milankovitch orbital forcing
in pelagic carbonate sequences. A prominent application is the detection and tracking of the
**405 kyr eccentricity cycle** — the most stable and predictable Milankovitch frequency over
geological time, driven by the gravitational interaction between Jupiter and Venus. The 405 kyr
long-eccentricity term serves as a fundamental metronome for astrochronological calibration and
has been identified in numerous deep-marine and pelagic records worldwide.

The pelagic limestone–marl successions of the **Apennine sections** (central Italy) — including
the classic outcrops at Gubbio (Bottaccione Gorge), Furlo, and the Contessa Highway — provide
some of the most celebrated cyclostratigraphic records. These sections were central to
**Walter Alvarez's** pioneering work on the **Cretaceous–Paleogene (K–Pg) boundary**, where the
now-famous iridium anomaly was first discovered in the Scaglia Rossa Formation (Alvarez et al.,
1980, *Science* 208: 1095–1108). The rhythmic bedding in these Apennine sections reflects
Milankovitch-driven oscillations in carbonate productivity and terrigenous dilution, and spectral
analysis of these records has resolved precession, obliquity, and eccentricity components —
notably the 405 kyr eccentricity cycle — enabling the construction of floating astrochronologies
anchored to the K–Pg boundary event.

WaverideR's wavelet-based and evolutionary spectral tools are ideally suited for tracking such
non-stationary cyclostratigraphic signals through intervals of variable sedimentation rate,
including across major stratigraphic boundaries like the K–Pg, where abrupt changes in
depositional regime and hiatuses are common.

### Key References

- Alvarez, L.W., Alvarez, W., Asaro, F., & Michel, H.V. (1980). Extraterrestrial cause for the
  Cretaceous-Tertiary extinction. *Science*, 208(4448), 1095–1108.
- Hinnov, L.A. (2018). Cyclostratigraphy and astrochronology in 2018. *Stratigraphy & Timescales*, 3, 1–80.
- Laskar, J., Fienga, A., Gastineau, M., & Manche, H. (2011). La2010: A new orbital solution for
  the long-term motion of the Earth. *Astronomy & Astrophysics*, 532, A89.
- Hilgen, F.J., Kuiper, K.F., & Lourens, L.J. (2010). Evaluation of the astronomical time scale
  for the Paleocene and earliest Eocene. *Earth and Planetary Science Letters*, 300(1–2), 55–66.

## Installation

You can install the development version of WaverideR from [GitHub](https://github.com/stratigraphy/WaverideR) with:

``` r
# install.packages("devtools")
devtools::install_github("stratigraphy/WaverideR")
```