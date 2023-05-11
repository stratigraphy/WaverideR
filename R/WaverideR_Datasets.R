#' @title Example data sets for the 'WaverideR' package
#' @name WaverideR_Datasets
#' @description Data sets for testing the 'WaverideR' R package:\cr
#' The \code{age_model_zeeden} data set is and age model (anchor points) for\cr
#' the IODP 926 grey scale (154-174m) record of Zeeden et al., (2013)\cr
#' \cr
#' The \code{astrosignal_example} data set consists of pre-generated ETP (eccentricity-tilt-precession)\cr
#' data set based on the p-0.5t  la2004 solution and was generated using \cr
#' the \link[astrochron]{etp} function of the 'astrochron' R package \cr
#' \cr
#' The \code{depth_rank_example} data set is synthetic succession of sedimentary\cr
#' The \code{grey} data set is the grey scale record of IODP 926 for the interval (154-174m) which originates\cr from
#'  Zeeden et al., (2013) \cr
#'  \cr
#' The \code{grey_track} data set consists of tracking points of the\cr
#'  precession (22 kyr cycle) in the IODP 926 grey scale (154-174m) record of Zeeden et al., (2013)\cr
#'  \cr
#' The \code{mag} data set is the magnetic susceptibility record of Pas et al., (2018)\cr
#' \cr
#' The \code{mag_track_solution} is the period of the 405 kyr eccentricity cycle in\cr
#' the magnetic susceptibility record of from Pas et al., (2018)\cr
#' \cr
#' The \code{TSI} data set is the Total Solar Irradiance record of Steinhilber et al., (2012)\cr
#'
#' @references
#'Damien Pas, Linda Hinnov, James E. (Jed) Day, Kenneth Kodama, Matthias Sinnesael, Wei Liu,
#'Cyclostratigraphic calibration of the Famennian stage (Late Devonian, Illinois Basin, USA),
#'Earth and Planetary Science Letters,
#'Volume 488,2018,Pages 102-114,ISSN 0012-821X,
#'\doi{<doi:10.1016/j.epsl.2018.02.010>}
#'
#'Steinhilber, Friedhelm & Abreu, Jacksiel & Beer, Juerg & Brunner,
#'Irene & Christl, Marcus & Fischer, Hubertus & Heikkilä, U. & Kubik,
#' Peter & Mann, Mathias & Mccracken, K. & Miller, Heinrich & Miyahara,
#' Hiroko & Oerter, Hans & Wilhelms, Frank. (2012).
#' 9,400 Years of cosmic radiation and solar activity from ice cores and tree rings.
#' Proceedings of the National Academy of Sciences of the United States of America.
#' 109. 5967-71. 10.1073/pnas.1118965109.
#' \doi{<doi:10.1073/pnas.1118965109>}
#'
#' Christian Zeeden, Frederik Hilgen, Thomas Westerhold, Lucas Lourens, Ursula Röhl, Torsten Bickert,
#' Revised Miocene splice, astronomical tuning and calcareous plankton biochronology of ODP Site 926 between 5 and 14.4Ma,
#' Palaeogeography, Palaeoclimatology, Palaeoecology,Volume 369,2013,Pages 430-451,ISSN 0031-0182,
#' \doi{<doi:10.1016/j.palaeo.2012.11.009>}
#'
#'Stephen R. Meyers,Cyclostratigraphy and the problem of astrochronologic testing,
#'Earth-Science Reviews,Volume 190,2019,Pages 190-223,ISSN 0012-8252
#'\doi{<doi:10.1016/j.earscirev.2018.11.015>}
#'
#'J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004,
#'A long term numerical solution for the insolation quantities of the Earth: Astron. Astrophys.,
#' Volume 428, 261-285. \doi{<doi:10.1051/0004-6361:20041335>} \cr
#'
NULL

#' @title Magnetic susceptibility data of the Sullivan core of Pas et al., (2018)
#' @rdname mag
#' @name mag
#' @description The magnetic susceptibility data set consists
#' of the magnetic susceptibility measurements of Pas et al., (2018), which measured the magnetic
#'  susceptibility on the Sullivan core which is of Famennian age.
#'
#' @details
#'Column 1: depth value (meters depoth)\cr
#'Column 2: magnetic susceptibility  value\cr
#'
#' @references
#'Damien Pas, Linda Hinnov, James E. (Jed) Day, Kenneth Kodama, Matthias Sinnesael, Wei Liu,
#'Cyclostratigraphic calibration of the Famennian stage (Late Devonian, Illinois Basin, USA),
#'Earth and Planetary Science Letters,
#'Volume 488,
#'2018,
#'Pages 102-114,
#'ISSN 0012-821X,
#' \doi{<doi:1016/j.epsl.2018.02.010>}
NULL

#' @title Total solar irradiation data (0-9400ka) of steinhilber et al., (2012)
#' @rdname TSI
#' @name TSI
#'
#' @description The Total solar irradiation data set consists of the TSI values of Steinhilber et al., (2012)
#'
#' @details
#'Column 1: Age (kyr)\cr
#'Column 2: Total solar Irradiation (TSI)\cr
#'
#' @references
#'Steinhilber, Friedhelm & Abreu, Jacksiel & Beer, Juerg & Brunner,
#'Irene & Christl, Marcus & Fischer, Hubertus & Heikkilä, U. & Kubik,
#' Peter & Mann, Mathias & Mccracken, K. & Miller, Heinrich & Miyahara,
#' Hiroko & Oerter, Hans & Wilhelms, Frank. (2012).
#' 9,400 Years of cosmic radiation and solar activity from ice cores and tree rings.
#' Proceedings of the National Academy of Sciences of the United States of America.
#' 109. 5967-71. 10.1073/pnas.1118965109.
#' \doi{<doi:10.1073/pnas.1118965109>}
NULL


#' @title Grey scale record IODP 926 of Zeeden et al., (2013)
#' @rdname grey
#' @name grey
#'
#' @description IODP 926 grey scale record of Zeeden et al., (2013) for the (154-174m) interval.
#' The (154-174m) interval spans the Miocene.
#'
#' @details
#'Column 1: depth (meters)\cr
#'Column 2: greyscale value\cr
#'
#' @references
#' Christian Zeeden, Frederik Hilgen, Thomas Westerhold, Lucas Lourens, Ursula Röhl, Torsten Bickert,
#' Revised Miocene splice, astronomical tuning and calcareous plankton biochronology of ODP Site 926 between 5 and 14.4Ma,
#' Palaeogeography, Palaeoclimatology, Palaeoecology,Volume 369,2013,Pages 430-451,ISSN 0031-0182,
#' \doi{<doi:10.1016/j.palaeo.2012.11.009>}
NULL


#' @title Age model of Zeeden et al., (2013) for the (154-174m) interval of the IODP 926 grey scale record
#' @rdname age_model_zeeden
#' @name age_model_zeeden
#'
#' @description
#' Age model (anchor points) of the IODP 926 grey scale (154-174m) record of Zeeden et al., (2013) \cr
#' Anchored to the eccentricity-tilt-precession model p-0.5t of la 2004.
#'
#' @details
#'Column 1: Depth (meters)\cr
#'Column 2: Age (kyr)\cr
#'
#' @references
#' Christian Zeeden, Frederik Hilgen, Thomas Westerhold, Lucas Lourens, Ursula Röhl, Torsten Bickert,
#' Revised Miocene splice, astronomical tuning and calcareous plankton biochronology of ODP Site 926 between 5 and 14.4Ma,
#' Palaeogeography, Palaeoclimatology, Palaeoecology,Volume 369,2013,Pages 430-451,ISSN 0031-0182,
#' \doi{<doi:10.1016/j.palaeo.2012.11.009>}
#'
#'
#'J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004,
#'A long term numerical solution for the insolation quantities of the Earth: Astron. Astrophys.,
#' Volume 428, 261-285. \doi{<doi:10.1051/0004-6361:20041335>}
NULL

#' @title Tracking points of the precession (22 kyr cycle) IODP 926 grey scale (154-174m) record of Zeeden et al., (2013)
#' @rdname grey_track
#' @name grey_track
#'
#' @description
#' Example data which consists of tracking points of the precession (22 kyr cycle) in the wavelet
#' spectra of the IODP 926 grey scale (154-174m) record of Zeeden et al., (2013)
#'
#' @details
#'Column 1: Depth (meters)\cr
#'Column 2: period (meters)\cr
#'
#' @references
#' Christian Zeeden, Frederik Hilgen, Thomas Westerhold, Lucas Lourens, Ursula Röhl, Torsten Bickert,
#' Revised Miocene splice, astronomical tuning and calcareous plankton biochronology of ODP Site 926 between 5 and 14.4Ma,
#' Palaeogeography, Palaeoclimatology, Palaeoecology,Volume 369,2013,Pages 430-451,ISSN 0031-0182,
#' \doi{<doi:10.1016/j.palaeo.2012.11.009>}
NULL

#' @title Period of the 405 kyr ecc cycle in the magnetic susceptibility record of the Sullivan core
#' @rdname mag_track_solution
#' @name mag_track_solution
#'
#' @description
#' Data points which give the period (in meters) of the 405 kyr eccentricity cycle tracked
#' in the wavelet spectra of the magnetic susceptibility record of the Sullivan core \cr
#' The period was tracked using the \code{\link{track_period_wavelet}} function\cr
#' The tracking is based on the original age model of Pas et al., (2018)\cr
#'
#'@details
#'Column 1: Depth (meters)\cr
#'Column 2: tracked period of 405 kyr eccentricity cycle (meters)\cr
#'
#' @references
#' Damien Pas, Linda Hinnov, James E. (Jed) Day, Kenneth Kodama, Matthias Sinnesael, Wei Liu,
#'Cyclostratigraphic calibration of the Famennian stage (Late Devonian, Illinois Basin, USA),
#'Earth and Planetary Science Letters,
#'Volume 488,
#'2018,
#'Pages 102-114,
#'ISSN 0012-821X,
#' \doi{<doi:10.1016/j.epsl.2018.02.010>}
NULL


#' @title An example depth rank series
#' @rdname depth_rank_example
#' @name depth_rank_example
#'
#' @description
#' The \code{\link{depth_rank_example}} example data set is a depth rank series which \cr
#' can be used as input for the \code{\link{lithlog_disc}} function which creates a \cr
#' discritzed record which can then be used as input in the \code{\link{analyze_wavelet}}\cr
#' function
#'
#'@details
#'Column 1: depth (meters)\cr
#'Column 2: depth rank \cr
NULL

#' @title An ETP astronomical solution
#' @rdname astrosignal_example
#' @name astrosignal_example
#'
#' @description
#'The \code{\link{astrosignal_example}} is a pre-generated ETP \cr
#'(eccentricity-tilt-precession) (p-0.5t based on the la2004 solution) \cr
#'the \code{\link{astrosignal_example}} can be used to anchor the \code{\link{grey}}\cr
#' data set to an astronomical solution eg. \code{\link{astrosignal_example}} \cr
#' using the  \code{\link{astro_anchor}} function. the data set was generated using the \cr
#' \link[astrochron]{etp} function of the 'astrochron' R package.
#'  The pre-generated ETP spans 5000 to 6000kyr.
#'
#' #'@details
#'Column 1: time (kyr)\cr
#'Column 2: ETP  \cr
#'
#' @author
#'Generated using the \link[astrochron]{etp}
#'function of the \link[astrochron]{astrochron-package}.
#'
#'@references
#'Stephen R. Meyers,Cyclostratigraphy and the problem of astrochronologic testing,
#'Earth-Science Reviews,Volume 190,2019,Pages 190-223,ISSN 0012-8252
#'\doi{<doi:10.1016/j.earscirev.2018.11.015>}
#'
#'J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004,
#'A long term numerical solution for the insolation quantities of the Earth: Astron. Astrophys.,
#' Volume 428, 261-285. \doi{<doi:10.1051/0004-6361:20041335>}
NULL
