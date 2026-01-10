#' @name WaverideR
#' @title An R package for cyclostratigraphic analysis
#' @description WaverideR is an R package for advanced cyclostratigraphic
#'analysis of stratigraphic data sets. It provides a comprehensive suite
#'of spectral tools for detecting, visualising, and tracking non-stationary
#'astronomical cycles, including continuous wavelet and superlet transform (CWT),
#'the superlet transform, windowed FFT, and evolutionary harmonic
#'analysis (EHA) (wrapper). These methods allow both manual and automated
#'tracking of orbital cycles in spectra and scalograms, even in records
#'affected by large changes in sedimentation rate.Building on this
#'spectral framework, WaverideR supports multi-proxy integration and
#'Monte Carlo–based uncertainty propagation to construct statistically
#'robust floating and absolute astrochronological age models. The package
#'includes dedicated tools to quantify analytical wavelet uncertainty,
#'estimate the duration of stratigraphic gaps and hiatuses, and integrate
#'external radioisotopic age constraints. Designed for complex and
#'incomplete stratigraphic records, WaverideR enables investigation
#'of the imprint of astronomical forcing even in suboptimal datasets.
#'
#'@references
#'The 'WaverideR' package builds upon existing literature and existing codebase.
#'The following list of articles is relevant for the 'WaverideR' R package and
#'its functions. Individual articles are also cited in the descriptions of
#'function when relative for set function. The articles in the list below can
#'be grouped in four subjects: (1) Cyclostratigraphic data analysis, (2)
#'example data sets, (3) the (continuous) wavelet transform and (4)
#'astronomical solutions). For each of these categories the
#'relevance of set articles will be explained in the framework of
#'the 'WaverideR' R package. \cr
#'
#'# 1. Cyclostratigraphic data analysis \cr
#'
#' Meyers, S. R. (2019). Cyclostratigraphy and the problem of astrochronologic
#' testing. Earth-Science Reviews, 190, 190--223.
#' \doi{10.1016/j.earscirev.2018.11.015} \cr
#'
#' Meyers, S. R. (2012). Seeing red in cyclic stratigraphy: Spectral noise
#' estimation for astrochronology. Paleoceanography, 27, PA3228.
#' \doi{10.1029/2012PA002307} \cr
#'
#' Li, M., Hinnov, L. A., and Kump, L. R. (2019). Acycle: Time-series analysis
#' software for paleoclimate research and education. Computers and Geosciences,
#' 127, 12--22. \doi{10.1016/j.cageo.2019.02.011} \cr
#'
#' Li, M., Kump, L. R., Hinnov, L. A., and Mann, M. E. (2018). Tracking variable
#' sedimentation rates and astronomical forcing in Phanerozoic paleoclimate proxy
#' series with evolutionary correlation coefficients and hypothesis testing.
#' Earth and Planetary Science Letters, 501, 165--179.
#' \doi{10.1016/j.epsl.2018.08.041} \cr
#'
#' Wouters, S., Crucifix, M., Sinnesael, M., Da Silva, A.-C., Zeeden, C.,
#' Zivanovic, M., Boulvain, F., and Devleeschouwer, X. (2022). A decomposition
#' approach to cyclostratigraphic signal processing. Earth-Science Reviews, 225,
#' 103894. \doi{10.1016/j.earscirev.2021.103894} \cr
#'
#' Wouters, S., Da Silva, A.-C., Boulvain, F., and Devleeschouwer, X. (2021).
#' StratigrapheR: Concepts for litholog generation in R. The R Journal, 13.
#' \doi{10.32614/RJ-2021-039} \cr
#'
#' Huang, N. E., Wu, Z., Long, S. R., Arnold, K. C., Chen, X., and Blank, K.
#' (2009).On instantaneous frequency. Advances in Adaptive Data Analysis,
#' 1, 177--229. \doi{10.1142/S1793536909000096} \cr
#'
#' Cleveland, W. S. (1979). Robust locally weighted regression and smoothing
#' scatterplots. Journal of the American Statistical Association, 74, 829--836.
#' \doi{10.1080/01621459.1979.10481038} \cr
#'
#' Hurvich, C. M., Simonoff, J. S., and Tsai, C.-L. (1998). Smoothing parameter
#' selection in nonparametric regression using an improved Akaike information
#' criterion. Journal of the Royal Statistical Society: Series B, 60, 271--293.
#' \doi{10.1111/1467-9868.00125} \cr
#'
#' Golub, G. H., Heath, M., and Wahba, G. (1979). Generalized cross-validation
#' as a method for choosing a good ridge parameter. Technometrics, 21, 215--224.
#' \doi{10.2307/1268518} \cr
#'
#' #2. Example data sets \cr
#'
#' Pas, D., Hinnov, L. A., Day, J. E., Kodama, K., Sinnesael, M., and Liu, W.
#' (2018). Cyclostratigraphic calibration of the Famennian Stage
#' (Late Devonian, Illinois Basin, USA). Earth and Planetary Science Letters,
#' 488, 102--114. \doi{10.1016/j.epsl.2018.02.010} \cr
#'
#' Steinhilber, F., Abreu, J., Beer, J., Brunner, I., Christl, M., Fischer, H.,
#' Heikkilä, U., Kubik, P. W., Mann, M. E., McCracken, K. G., Miller, H., Miyahara,
#' H., Oerter, H., and Wilhelms, F. (2012). 9400 years of cosmic radiation and solar
#' activity from ice cores and tree rings. Proceedings of the National Academy of
#' Sciences of the United States of America, 109, 5967--5971.
#' \doi{10.1073/pnas.1118965109} \cr
#'
#' Zeeden, C., Hilgen, F. J., Westerhold, T., Lourens, L. J., Röhl, U., and
#' Bickert, T. (2013). Revised Miocene splice, astronomical tuning and calcareous
#' plankton biochronology of ODP Site 926 between 5 and 14.4 Ma.
#' Palaeogeography, Palaeoclimatology, Palaeoecology, 369, 430--451.
#' \doi{10.1016/j.palaeo.2012.11.009} \cr
#'
#' #3. Continuous wavelet and superlet transform\cr
#'
#' Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., and Mureșan, R. C. (2021).
#' Time-frequency super-resolution with superlets. Nature Communications, 12, 337.
#' \doi{10.1038/s41467-020-20539-9} \cr
#'
#' Morlet, J., Arens, G., Fourgeau, E., and Giard, D. (1982a). Wave propagation and
#' sampling theory. Part I: Complex signal and scattering in multilayered media.
#' Geophysics, 47, 203--221. \cr
#'
#' Morlet, J., Arens, G., Fourgeau, E., and Giard, D. (1982b). Wave propagation and
#' sampling theory. Part II: Sampling theory and complex waves. Geophysics, 47,
#' 222--236. \cr
#'
#' Torrence, C., and Compo, G. P. (1998). A practical guide to wavelet analysis.
#' Bulletin of the American Meteorological Society, 79, 61--78. \cr
#'
#' Gouhier, T. C., Grinsted, A., and Simko, V. (2021). biwavelet: Conduct univariate
#' and bivariate wavelet analyses. R package version 0.20.21. \cr
#'
#' Roesch, A., and Schmidbauer, H. (2018). WaveletComp: Computational wavelet
#' analysis. R package version 1.1. \cr
#'
#' Gabor, D. (1946). Theory of communication. Part I: The analysis of information.
#' Journal of the Institution of Electrical Engineers, Part III, 93, 429--441. \cr
#'
#' #4. Astronomical solutions \cr
#'
#' Laskar, J., Robutel, P., Joutel, F., Gastineau, M., Correia, A. C. M., and
#' Levrard, B. (2004). A long-term numerical solution for the insolation quantities
#' of the Earth. Astronomy and Astrophysics, 428, 261--285.
#' \doi{10.1051/0004-6361:20041335} \cr
#'
#' Laskar, J., Fienga, A., Gastineau, M., and Manche, H. (2011a). La2010: A new
#' orbital solution for the long-term motion of the Earth. Astronomy and
#' Astrophysics, 532, A89. \doi{10.1051/0004-6361/201116836} \cr
#'
#' Zeebe, R. E., and Lourens, L. J. (2019). Solar System chaos and the
#' Paleocene--Eocene boundary age constrained by geology and astronomy.
#' Science, 365, 926--929. \doi{10.1126/science.aax0612}
#'
#' @details Package: 'WaverideR'
#'
#' Type: R package
#'
#' Version: 0.5.0 (1st quarter 2026)
#'
#' License: GPL (= 2)
#'
#' @note
#' If you want to use this package for publication or research
#' purposes, please cite:
#'
#' Arts, M.C.M (2023).
#' WaverideR: Extracting Signals from Wavelet Spectra.
#' https://CRAN.R-project.org/package=WaverideR
#'
#' @author Michiel Arts
#'
#' Maintainer: Michiel Arts \email{michiel.arts@stratigraphy.eu}
NULL
