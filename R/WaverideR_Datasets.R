#' @title Example data sets for the WaverideR package
#' @name WaverideR_Datasets
#' @description Data sets for testing the WaverideR package:\cr
#' The \code{mag} data set is the magnetic susceptibility record of De Pas et al., (2018)\cr
#' The \code{TSI} data set is the Total Solar Irradiance record of Steinhilber et al., (2012)\cr
#' The \code{grey} data set is the grey scale record of IODP 926 for the interval (154-174m) which originates from
#'  Zeeden et al., (2013) \cr
#' The \code{mag_track_solution} is the period of the 405 kyr ecc cycle in
#' the magnetic susceptibility record of rom De Pas et al., (2018)\cr
#' the \code{age_model_zeeden} data set is and age model (anchor points) for
#' the IODP 926 grey scale (154-174m) record of Zeeden et al., (2013)\cr
#' the \code{grey_track} data set consists of tracking points of the
#'  precession (22 kyr cycle) in the IODP 926 grey scale (154-174m) record of Zeeden et al., (2013)
#'
#' @references
#'Damien Pas, Linda Hinnov, James E. (Jed) Day, Kenneth Kodama, Matthias Sinnesael, Wei Liu,
#'Cyclostratigraphic calibration of the Famennian stage (Late Devonian, Illinois Basin, USA),
#'Earth and Planetary Science Letters,
#'Volume 488,2018,Pages 102-114,ISSN 0012-821X,
#' <doi:10.1016/j.epsl.2018.02.010>
#'
#'Steinhilber, Friedhelm & Abreu, Jacksiel & Beer, Juerg & Brunner,
#'Irene & Christl, Marcus & Fischer, Hubertus & Heikkilä, U. & Kubik,
#' Peter & Mann, Mathias & Mccracken, K. & Miller, Heinrich & Miyahara,
#' Hiroko & Oerter, Hans & Wilhelms, Frank. (2012).
#' 9,400 Years of cosmic radiation and solar activity from ice cores and tree rings.
#' Proceedings of the National Academy of Sciences of the United States of America.
#' 109. 5967-71. 10.1073/pnas.1118965109.
#' <doi:10.1073/pnas.1118965109>
#'
#' Christian Zeeden, Frederik Hilgen, Thomas Westerhold, Lucas Lourens, Ursula Röhl, Torsten Bickert,
#' Revised Miocene splice, astronomical tuning and calcareous plankton biochronology of ODP Site 926 between 5 and 14.4Ma,
#' Palaeogeography, Palaeoclimatology, Palaeoecology,Volume 369,2013,Pages 430-451,ISSN 0031-0182,
#' <doi:10.1016/j.palaeo.2012.11.009>
#'
NULL

#' @title Magnetic susceptibility data of the Sullivan core of De Pas et al., (2018)
#' @rdname mag
#' @name mag
#' @description The magnetic susceptibility data set consists
#' of the magnetic susceptibility measurements of De Pas et al., (2018), which measured the magnetic
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
#' <doi:1016/j.epsl.2018.02.010>
NULL

#' @title Total solar irradiation data (0-9400ka) of steinhilber et al., (2012)
#' @rdname TSI
#' @name TSI
#'
#' @description The Total solar irradiation data set consists of the TSI values of steinhilber et al., (2012)
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
#' <doi:10.1073/pnas.1118965109>
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
#' <doi:10.1016/j.palaeo.2012.11.009>
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
#' <doi:10.1016/j.palaeo.2012.11.009>
#'
#'
#' Laskar, Jacques, Philippe Robutel, Frédéric Joutel, Mickael Gastineau, Alexandre CM Correia,
#'  and Benjamin Levrard. "A long-term numerical solution for the insolation quantities of the Earth.
#'  " Astronomy & Astrophysics 428, no. 1 (2004): 261-285.
#' \url{https://www.aanda.org/articles/aa/pdf/2004/46/aa1335.pdf}
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
#' <doi:10.1016/j.palaeo.2012.11.009>
NULL

#' @title Period of the 405 kyr ecc cycle in the magnetic susceptibility record of the Sullivan core
#' @rdname mag_track_solution
#' @name mag_track_solution
#'
#' @description
#' Data points which give the period (in meters) of the 405 kyr eccentricity cycle tracked
#' in the wavelet spectra of the magnetic susceptibility record of the Sullivan core \cr
#' The period was tracked using the \code{\link{track_period_wavelet}} function\cr
#' The tracking is based on the original age model of De Pas et al., (2018)\cr
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
#' <doi:10.1016/j.epsl.2018.02.010>
NULL







