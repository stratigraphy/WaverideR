#' @name WaverideR
#' @title Extracting Signals from Wavelet Spectra
#' @description The continuous wavelet transform enables the observation of transient/non-stationary
#' cyclicity in time-series. The goal of cyclostratigraphic studies is to define frequency/period
#' in the depth/time domain. By conducting the continuous wavelet transform on cyclostratigraphic
#' data series one can observe and extract cyclic signals/signatures from signals.
#' These results can then be visualized and interpreted enabling one to identify/interpret
#' cyclicity in the geological record, which can be used to construct astrochronological
#' age-models and identify and interpret cyclicity in past en present climate systems.
#'
#'@references
#'The 'WaverideR' package builds upon existing literature and an existing codebase.
#'The following list of articles are articles which are relevant for the functions of the 'WaverideR'
#'package but also for the 'WaverideR' R package as a whole.
#'Individual articles are (re)cited in the descriptions of function when relative for set function.
#'The list of Articles is can be grouped in four subjects (Cyclostratigraphic data analysis,
#'example data sets,the (continuous) wavelet transform and astronomical solutions). For each of
#'these categories the relevance of set articles will be given and then explained why set article is
#'important for the 'WaverideR' R package.
#'
#'# 1. Cyclostratigraphic data analysis
#'
#'Stephen R. Meyers,Cyclostratigraphy and the problem of astrochronologic testing,
#'Earth-Science Reviews,Volume 190,2019,Pages 190-223,ISSN 0012-8252
#'\doi{<doi:10.1016/j.earscirev.2018.11.015>} \cr
#'The article of Meyers (2019) is the article which needs to cited when one wants to cite the
#'astrochron R package. The 'astrochron' R package is the most extensive
#' R package with regards to cyclostratigraphic analysis. As such many of the functionalities of the 'WaverideR' R package
#' are inspired/based on the 'astrochron' R package.The major difference between the 'astrochron' R package
#' and the 'WaverideR' package is that the 'astrochron' R package relies on the Fast Fourier Transform
#' whereas the 'WaverideR' R package relies on the continuous wavelet transform (CWT).
#'
#'S.R. Meyers, 2012, Seeing Red in Cyclic Stratigraphy: Spectral Noise Estimation for
#'Astrochronology: Paleoceanography, 27, PA3228, \doi{<doi:10.1029/2012PA002307>}
#'The article of Meyers (2012)  explain the (Multitaper method) MTM technique implemented into The 'astrochron' R package
#'The MTM method can be used to assign confidence levels to spectral peaks and distinguish spectral peaks
#'from harmonic spectral peaks
#'
#'
#'Acycle: Time-series analysis software for paleoclimate research and education,
#'Mingsong Li, Linda Hinnov, Lee Kump,
#'Computers & Geosciences,Volume 127,2019,Pages 12-22,ISSN 0098-3004,
#'\doi{<doi:10.1016/j.cageo.2019.02.011>}\cr
#'The 'Acycle' software package is a 'Matlab' based program, which is used for cyclostratigraphic studies.
#'Acycle has many relies mostly on the Fast Fourier Transform. Many of the Fast Fourier Transform based functions
#'formed as inspiration for the continuous wavelet transform (CWT) based functions of the 'WaverideR' R package.
#'
#'Tracking variable sedimentation rates and astronomical forcing in Phanerozoic paleoclimate proxy series with evolutionary correlation coefficients and hypothesis testing,
#'Mingsong Li, Lee R. Kump, Linda A. Hinnov, Michael E. Mann,
#'Earth and Planetary Science Letters,Volume 501,
#'2018,Pages 165-179,ISSN 0012-821X,\doi{<doi:10.1016/j.epsl.2018.08.041>}\cr
#'Li et al., (2019) introduces the Coco and eCoco functions of the Acycle software package.
#'the 'Coco' and 'eCoco' function of the 'Acycle' software are able to estimate the sedimentation rate based on
#'spectral characteristics of astronomical cycles. The 'Coco' and 'eCoco' function and form the inspiration for the inspiration
#'for the \link[WaverideR]{flmw} \link[WaverideR]{sum_power_sedrate} functions.
#'
#'Wouters, S., Crucifix, M., Sinnesael, M., Da Silva, A.C., Zeeden, C., Zivanovic, M., Boulvain, F.,
#'Devleeschouwer, X., 2022, "A decomposition approach to cyclostratigraphic signal processing".
#'Earth-Science Reviews 225 (103894).\doi{<doi:10.1016/j.earscirev.2021.103894>}\cr
#'Wouters  et al., (2022) introduces the Empirical Mode Decomposition (EMD) as
#'part of the 'DecomposeR' R package. EMD is a non-Fast Fourier Transform based spectral
#'analysis technique. The Hilbert transform function (\link[DecomposeR]{inst.pulse}) of set package is used as
#'a control for the estimation of the amplitude extraction using wavelets
#'
#'Wouters, S., Da Silva, A.-C., Boulvain, F., and Devleeschouwer, X.. 2021.
#' StratigrapheR: Concepts for Litholog Generation in R.
#' The R Journal. \doi{doi:10.32614/RJ-2021-039>}\cr
#'Wouters et al., (2021) introduces the  \link[StratigrapheR]{StratigrapheR} R package. This package
#'contains functions which format, process, and plot lithologs. The litholog format of Wouters et al., (2021)
#'is used as the standardized input format which is used to convert lithologs to format on which the
#' continuous wavelet transform (CWT) can be conducted.
#'
#'Huang, Norden E., Zhaohua Wu, Steven R. Long, Kenneth C. Arnold, Xianyao Chen, and Karin Blank. 2009.
#' "On Instantaneous Frequency". Advances in Adaptive Data Analysis 01 (02): 177–229. \doi{<doi:10.1142/S1793536909000096>}\cr
#'The Hilbert transform function (\link[DecomposeR]{inst.pulse}) of the 'DecomposeR' R package is based on
#'the work of Huang et al., (2009).
#'
#'Cleveland, W. S. (1979) Robust locally weighted regression and smoothing scatter plots. Journal of the American Statistical Association. 74, 829–836. \doi{doi:10.1080/01621459.1979.10481038}
#'Cleveland (1979) explains how the robust locally weighted regression works and how it can be used to smooth data sets
#'
#'#'Hurvich, C.M., Simonoff, J.S., and Tsai, C.L. (1998), Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion. Journal of the Royal Statistical Society B. 60, 271–293 \doi{<doi:10.1111/1467-9868.00125>}
#'Hurvich et al., (1998) explains how the Improved Akaike Information Criterion can be used to optimally smooth data sets
#'
#'#'Golub, G., Heath, M. and Wahba, G. (1979). Generalized cross validation as a method for choosing a good ridge parameter. Technometrics. 21, 215–224. \doi{<doi:10.2307/1268518>}
#'Golub et al., (1979) explains how the Generalized cross validation can be used to optimally smooth data sets
#'
#'# 2. Example data sets
#'
#'Damien Pas, Linda Hinnov, James E. (Jed) Day, Kenneth Kodama, Matthias Sinnesael, Wei Liu,
#'Cyclostratigraphic calibration of the Famennian stage (Late Devonian, Illinois Basin, USA),
#'Earth and Planetary Science Letters,
#'Volume 488,2018,Pages 102-114,ISSN 0012-821X,
#' \doi{<doi:10.1016/j.epsl.2018.02.010>}\cr
#'The data set of Pas et al, (2018) is a magnetic susceptibility data measured on the Fammennian aged shales of
#'the from the Illinois basin in the USA. The data set contains the imprint of astronomical cycles in the a
#'Paleozoic succession making it a good example for times (>250Ma) when no astronomical solutions are available.
#'
#'
#'Steinhilber, Friedhelm & Abreu, Jacksiel & Beer, Juerg & Brunner,
#'Irene & Christl, Marcus & Fischer, Hubertus & Heikkilä, U. & Kubik,
#' Peter & Mann, Mathias & Mccracken, K. & Miller, Heinrich & Miyahara,
#' Hiroko & Oerter, Hans & Wilhelms, Frank. (2012).
#' 9,400 Years of cosmic radiation and solar activity from ice cores and tree rings.
#' Proceedings of the National Academy of Sciences of the United States of America.
#' 109. 5967-71. 10.1073/pnas.1118965109.
#' \doi{<doi:10.1073/pnas.1118965109>}\cr
#'The Total Solar Irradiance record of Steinhilber et al., (2012) is a Holocene record of normalized
#'Total Solar Irradiance in the time domain. The data set is a good example for
#'studying/extracting sub-Milankovitch <5000yr from a relatively (geologically) speaking young record.
#'
#' Christian Zeeden, Frederik Hilgen, Thomas Westerhold, Lucas Lourens, Ursula Röhl, Torsten Bickert,
#' Revised Miocene splice, astronomical tuning and calcareous plankton biochronology of ODP Site 926 between 5 and 14.4Ma,
#' Palaeogeography, Palaeoclimatology, Palaeoecology,Volume 369,2013,Pages 430-451,ISSN 0031-0182,
#' \doi{<doi:10.1016/j.palaeo.2012.11.009>}\cr
#' The record of Zeeden et al., (2013) consists of a grey scale record from Miocene sediment cores from offshore
#' Brazil. The record contains a clear imprint of astronomical cycles as such it is a good Neogene example data set
#' to demonstrate the functionalities of the 'WaverideR' R package
#'
#'# 3. The (continuous) wavelet transform
#'
#'Morlet, Jean, Georges Arens, Eliane Fourgeau, and Dominique Glard.
#'"Wave propagation and sampling theory—Part I: Complex signal and scattering in multilayered media.
#'" Geophysics 47, no. 2 (1982): 203-221.
#' \doi{<doi.org/10.1190/1.1441328>}\cr
#' Morlet et al., (1982a) together with Morlet et al., (1982b) are the original publications which explain the
#' use of the wavelet to analyse signal.
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#' Geophysics 1982 47 (2): 222–236. \doi{<doi.org/10.1190/1.1441329>}\cr
#' Morlet et al., (1982a) together with Morlet et al., (1982b) are the original publications which explain the
#' use of the wavelet to analyse signal.
#'
#'Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#'Bulletin of the American Meteorological Society 79:61-78.
#'\url{https://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf}
#''Torrence and Compo (1998) shows how the continuous wavelet transform can be used to analyse
#'cyclicity in paleo-climatic data-sets. The equations in this publication forms the basis for many
#'wavelet based packages/software applications.
#'
#'
#'Gouhier TC, Grinsted A, Simko V (2021).
#'R package biwavelet: Conduct Univariate and Bivariate Wavelet Analyses. (Version 0.20.21),
#'\url{https://github.com/tgouhier/biwavelet}\cr
#'Gouhier et al., (2021) is the implementation of equations of Torrence and Compo (1998) in the form of the
#''biwavelet' R package
#'
#'Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
#'Wavelet Analysis. R package version 1.1.
#'\url{https://CRAN.R-project.org/package=WaveletComp}\cr
#'Roesch and Schmidbauer et al., (2018) is the article of the 'WaveletComp' R package which
#'is a built upon the functionalities of the 'biwavelet' R package
#'
#'Russell, Brian, and Jiajun Han. "Jean Morlet and the continuous wavelet transform.
#'" CREWES Res. Rep 28 (2016): 115. \url{https://www.crewes.org/Documents/ResearchReports/2016/CRR201668.pdf}\cr
#'Russell and Han (2016) gives a concise summary of the work of Morlet et al., (1982a) and Morlet et al., (1982b) and the
#'developments since then. The publication also describes how the Gabor uncertainty principle (Gabor 1946) affects the frequency
#'uncertainty of the wavelet which can be used to calculate the analytical uncertainty of a given wavelet spectra.
#'
#'Gabor, Dennis. "Theory of communication. Part 1: The analysis of information."
#' Journal of the Institution of Electrical Engineers-part III: radio and
#' communication engineering 93, no. 26 (1946): 429-441. \url{http://genesis.eecg.toronto.edu/gabor1946.pdf}\cr
#'Gabor (1946) describes the Gabor uncertainty principle which states how the uncertainty in time and frequency are related
#'in time series analysis.
#'
#'# 4. Astronomical solutions
#'
#'J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004,
#'A long term numerical solution for the insolation quantities of the Earth: Astron. Astrophys.,
#' Volume 428, 261-285. \doi{<doi:10.1051/0004-6361:20041335>}\cr
#' Laskar et al., (2004) is an astronomical solution which can be used to anchor geological data to absolute ages.
#'
#'Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011a,
#' La2010: A new orbital solution for the long-term motion of the Earth: Astron. Astrophys.,
#' Volume 532, A89 \doi{<doi:10.1051/0004-6361/201116836>}\cr
#'Laskar et al., (2011a) is an astronomical solution which can be used to anchor geological data to absolute ages.
#'
#'Laskar, J., Gastineau, M., Delisle, J.-B., Farres, A., Fienga, A.:
#'2011b, Strong chaos induced by close encounters with Ceres and Vesta, Astron: Astrophys.,
#'Volume 532, L4.  \doi{<doi:10.1051/0004-6361/201117504>}\cr
#'Laskar et al., (2011b) is an astronomical solution which can be used to anchor geological data to absolute ages.
#'
#'J. Laskar,Chapter 4 - Astrochronology,Editor(s): Felix M. Gradstein, James G. Ogg, Mark D. Schmitz, Gabi M. Ogg,Geologic Time Scale 2020,Elsevier,2020,Pages 139-158,ISBN 9780128243602,
#' '\doi{<doi:10.1016/B978-0-12-824360-2.00004-8>}\cr
#'Laskar et al., (2019) explains how astronomical solutions are created and how they should/can be used
#'
#'Zeebe, Richard E. "Numerical solutions for the orbital motion of the Solar System over the past 100 Myr: limits and new results."
#'The Astronomical Journal 154, no. 5 (2017): 193. \doi{<doi:10.3847/1538-3881/aa8cce>} \cr
#'Zeebe (2017) is an astronomical solution which can be used to anchor geological data to absolute ages.
#'
#'Richard E. Zeebe Lucas J. Lourens ,Solar System chaos and the Paleocene–Eocene boundary age constrained by geology and astronomy.Science365,926-929(2019)
#'\doi{<doi:10.1126/science.aax0612>}\cr
#'Zeebe and Lourens (2019) is an astronomical solution which can be used to anchor geological data to absolute ages.
#'
#'Zeebe, R. E. and Lourens, L. J.
#'Geologically constrained astronomical solutions for the Cenozoic era,
#'Earth and Planetary Science Letters, 2022 \doi{<doi:10.1016/j.epsl.2022.117595>}\cr
#'Zeebe and Lourens (2022) is an astronomical solution which can be used to anchor geological data to absolute ages.
#'
#'
#' @details Package: 'WaverideR'
#'
#' Type: R package
#'
#' Version: 0.3.1 (begin of 2023)
#'
#' License: GPL (>= 2)
#'
#' @note
#' If you want to use this package for publication or research
#' purposes, please cite: to be submitted
#'
#' @author Michiel Arts
#'
#' Maintainer: Michiel Arts \email{michiel.arts@stratigraphy.eu}
NULL
