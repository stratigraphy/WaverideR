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
#'\cr
#' The \code{Bisciaro_Mg_wt_track} data set is the 110-kyr (short eccentricity) \cr
#' cycle tracked in the wavelet scalogram of the Magnesium (XRF) record of Arts (2014)\cr
#'\cr
#' The \code{Bisciaro_Mn_wt_track} data set is the 110-kyr (short eccentricity) \cr
#' cycle tracked in the wavelet scalogram of the Manganese (XRF)record of Arts (2014)\cr
#'\cr
#' The \code{Bisciaro_al_wt_track} data set is the 110-kyr (short eccentricity) \cr
#' cycle tracked in the wavelet scalogram of the Aluminum (XRF) record of Arts (2014)\cr
#'\cr
#' The \code{Bisciaro_ca_wt_track} data set is the 110-kyr (short eccentricity) \cr
#' cycle tracked in the wavelet scalogram of the Calcium (XRF) record of Arts (2014)\cr
#'\cr
#' The \code{Bisciaro_sial_wt_track} data set is the 110-kyr (short eccentricity) \cr
#' cycle tracked in the wavelet scalogram of the Silicon/Aluminum (XRF) record of Arts (2014)\cr
#'\cr
#' The \code{Bisciaro_XRF} is the XRF data set of Arts (2014)\cr
#'\cr
#'
#' The \code{anchor_points_Bisciaro_al} data set consist of the tie points between the
#' Bisciaro_al record of Arts (2014) and the la2011 solution of laskar et al., (20111)\cr
#' \cr
#'
#' The \code{GTS_info} data set contains the color coding and ages and uncertainties
#' of Geologic Time Scale 2020 of Ogg (et al., 2021)\cr
#' \cr
#'
#' @references
#'Damien Pas, Linda Hinnov, James E. (Jed) Day, Kenneth Kodama, Matthias Sinnesael, Wei Liu,
#'Cyclostratigraphic calibration of the Famennian stage (Late Devonian, Illinois Basin, USA),
#'Earth and Planetary Science Letters,
#'Volume 488,2018,Pages 102-114,ISSN 0012-821X,
#'<doi:10.1016/j.epsl.2018.02.010>
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
#'Stephen R. Meyers,Cyclostratigraphy and the problem of astrochronologic testing,
#'Earth-Science Reviews,Volume 190,2019,Pages 190-223,ISSN 0012-8252
#'<doi:10.1016/j.earscirev.2018.11.015>
#'
#'J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004,
#'A long term numerical solution for the insolation quantities of the Earth: Astron. Astrophys.,
#' Volume 428, 261-285. <doi:10.1051/0004-6361:20041335> \cr
#'
#'Laskar, J., M. Gastineau, J. B. Delisle, A. Farrés, and A. Fienga (2011b),
#'Strong chaos induced by close encounters with Ceres and Vesta, Astron. Astrophys.,
#'532, L4,<doi:10.1051/0004-6361/201117504> \cr
#'
#'M.C.M. Arts, 2014,
#'Magnetostratigrpahy and geochemical analysis of the early Miocene Bisciaro Formation
#'in the Contessa Valley (Northern Italy). Unpublished Bsc. thesis \cr
#'
#'Ogg, Gabi & Ogg, James & Gradstein, Felix. (2021).
#'Recommended color coding of stages - Appendix 1
#'from Geologic Time Scale 2020.
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
#' <doi:1016/j.epsl.2018.02.010>
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
#'J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004,
#'A long term numerical solution for the insolation quantities of the Earth: Astron. Astrophys.,
#' Volume 428, 261-285. <doi:10.1051/0004-6361:20041335>
NULL

#' @title Tracking points of the precession (22 kyr cycle) IODP 926 grey scale (154-174m) record of Zeeden et al., (2013)
#' @rdname grey_track
#' @name grey_track
#'
#' @description
#' Example data which consists of tracking points of the precession (22 kyr cycle) in the wavelet
#' scalogram of the IODP 926 grey scale (154-174m) record of Zeeden et al., (2013)
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
#' in the wavelet scalogram of the magnetic susceptibility record of the Sullivan core \cr
#' The period was tracked using the \code{\link{track_period_wavelet}} function\cr
#' The tracking is based on the original age model of Pas et al., (2018)\cr
#'
#'@details
#'Column 1: Depth (meters)\cr
#'Column 2: tracked period of 405 kyr eccentricity cycle (meters)\cr
#'
#' @references
#'Damien Pas, Linda Hinnov, James E. (Jed) Day, Kenneth Kodama, Matthias Sinnesael, Wei Liu,
#'Cyclostratigraphic calibration of the Famennian stage (Late Devonian, Illinois Basin, USA),
#'Earth and Planetary Science Letters,
#'Volume 488,
#'2018,
#'Pages 102-114,
#'ISSN 0012-821X,
#' <doi:10.1016/j.epsl.2018.02.010>
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
#' @details
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
#'<doi:10.1016/j.earscirev.2018.11.015>
#'
#'J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004,
#'A long term numerical solution for the insolation quantities of the Earth: Astron. Astrophys.,
#' Volume 428, 261-285. <doi:10.1051/0004-6361:20041335>
NULL


#' @title Example anchor points for the grey scale data set of Zeeden et al., (2013)
#' @rdname anchor_points_grey
#' @name anchor_points_grey
#'
#' @description
#'An example of anchor points generated using \code{\link{astro_anchor}} function \cr
#'The anchor points were generated for the \code{\link{grey}} grey data set of \cr
#'Zeeden et al., (2013) and anchored to the code{\link{astrosignal_example}} \cr
#'astronomical solution which is a pre-generated ETP (eccentricity-tilt-precession) \cr
#'solution(p-0.5t based on the la2004 solution) based on Laskar et al., (20004) \cr
#'astronomical solution.
#'
#' @details
#'Column 1: depth proxy record\cr
#'Column 2: time astronomical solution  \cr
#'Column 3: y-scale value proxy record\cr
#'Column 4: y-scale value astronomical solution  \cr
#'
#'
#'@references
#' Christian Zeeden, Frederik Hilgen, Thomas Westerhold, Lucas Lourens, Ursula Röhl, Torsten Bickert,
#' Revised Miocene splice, astronomical tuning and calcareous plankton biochronology of ODP Site 926 between 5 and 14.4Ma,
#' Palaeogeography, Palaeoclimatology, Palaeoecology,Volume 369,2013,Pages 430-451,ISSN 0031-0182,
#' <doi:10.1016/j.palaeo.2012.11.009>
#'
#'J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004,
#'A long term numerical solution for the insolation quantities of the Earth: Astron. Astrophys.,
#' Volume 428, 261-285. <doi:10.1051/0004-6361:20041335>
NULL


#' @title Period of the short kyr ecc cycle in the Mg record of the Bisciaro Fm
#' @rdname Bisciaro_Mg_wt_track
#' @name Bisciaro_Mg_wt_track
#'
#' @description
#' Data points which give the period (in meters) of the short kyr eccentricity cycle tracked \cr
#' in the wavelet scalogram of the magnesium (XRF) record of the Bisciaro Formation \cr
#' The period was tracked using the \code{\link{track_period_wavelet}} function\cr
#' The tracking is based on a reinterpretation of Arts (2014)\cr
#'
#' @details
#'Column 1: depth proxy record\cr
#'Column 2: period tracked in the wavelet scalogram of the Magnesium (XRF) record
#'
#'
#'@references
#'M.C.M. Arts, 2014,
#'Magnetostratigrpahy and geochemical analysis of the early Miocene Bisciaro Formation
#'in the Contessa Valley (Northern Italy). Unpublished Bsc. thesis \cr
NULL

#' @title Period of the short kyr ecc cycle in the Mn record of the Bisciaro Fm
#' @rdname Bisciaro_Mn_wt_track
#' @name Bisciaro_Mn_wt_track
#'
#' @description
#' Data points which give the period (in meters) of the short kyr eccentricity cycle tracked \cr
#' in the wavelet scalogram of the manganese (XRF) record of the Bisciaro Formation \cr
#' The period was tracked using the \code{\link{track_period_wavelet}} function\cr
#' The tracking is based on a reinterpretation of Arts (2014)\cr
#'
#' @details
#'Column 1: depth proxy record\cr
#'Column 2: period tracked in the wavelet scalogram of the manganese (XRF) record
#'
#'
#'@references
#'M.C.M. Arts, 2014,
#'Magnetostratigrpahy and geochemical analysis of the early Miocene Bisciaro Formation
#'in the Contessa Valley (Northern Italy). Unpublished Bsc. thesis \cr
NULL

#' @title Period of the short kyr ecc cycle in the Al record of the Bisciaro Fm
#' @rdname Bisciaro_al_wt_track
#' @name Bisciaro_al_wt_track
#'
#' @description
#' Data points which give the period (in meters) of the short kyr eccentricity cycle tracked \cr
#' in the wavelet scalogram of the aluminium (XRF) record of the Bisciaro Formation \cr
#' The period was tracked using the \code{\link{track_period_wavelet}} function\cr
#' The tracking is based on a reinterpretation of Arts (2014)\cr
#'
#' @details
#'Column 1: depth proxy record\cr
#'Column 2: period tracked in the wavelet scalogram of the Aluminium (XRF) record
#'
#'
#'@references
#'M.C.M. Arts, 2014,
#'Magnetostratigrpahy and geochemical analysis of the early Miocene Bisciaro Formation
#'in the Contessa Valley (Northern Italy). Unpublished Bsc. thesis \cr
NULL



#' @title Period of the short kyr ecc cycle in the Ca record of the Bisciaro Fm
#' @rdname Bisciaro_ca_wt_track
#' @name Bisciaro_ca_wt_track
#'
#' @description
#' Data points which give the period (in meters) of the short kyr eccentricity cycle tracked \cr
#' in the wavelet scalogram of the calcium (XRF) record of the Bisciaro Formation \cr
#' The period was tracked using the \code{\link{track_period_wavelet}} function\cr
#' The tracking is based on a reinterpretation of Arts (2014)\cr
#'
#' @details
#'Column 1: depth proxy record\cr
#'Column 2: period tracked in the wavelet scalogram of the calcium (XRF) record
#'
#'
#'@references
#'M.C.M. Arts, 2014,
#'Magnetostratigrpahy and geochemical analysis of the early Miocene Bisciaro Formation
#'in the Contessa Valley (Northern Italy). Unpublished Bsc. thesis \cr
NULL


#' @title Period of the short kyr ecc cycle in the si/Al record of the Bisciaro Fm
#' @rdname Bisciaro_sial_wt_track
#' @name Bisciaro_sial_wt_track
#'
#' @description
#' Data points which give the period (in meters) of the short kyr eccentricity cycle tracked \cr
#' in the wavelet scalogram of the silicon/aluminium (XRF) record of the Bisciaro Formation \cr
#' The period was tracked using the \code{\link{track_period_wavelet}} function\cr
#' The tracking is based on a reinterpretation of Arts (2014)\cr
#'
#' @details
#'Column 1: depth proxy record\cr
#'Column 2: period tracked in the wavelet scalogram of the silicon/aluminium (XRF) record
#'
#'
#'@references
#'M.C.M. Arts, 2014,
#'Magnetostratigrpahy and geochemical analysis of the early Miocene Bisciaro Formation
#'in the Contessa Valley (Northern Italy). Unpublished Bsc. thesis \cr
NULL


#' @title XRF records of the Bisciaro Fm
#' @rdname Bisciaro_XRF
#' @name Bisciaro_XRF
#'
#' @description
#'XRF proxy records from the early Miocene Bisciaro Formation
#'in the Contessa Valley (Northern Italy)
#'
#' @details
#'Column 1: depth proxy record\cr
#'Column 2-71: XRF proxy records
#'
#'
#'@references
#'M.C.M. Arts, 2014,
#'Magnetostratigrpahy and geochemical analysis of the early Miocene Bisciaro Formation
#'in the Contessa Valley (Northern Italy). Unpublished Bsc. thesis \cr
NULL

#' @title XRF records of the Bisciaro Fm
#' @rdname anchor_points_Bisciaro_al
#' @name anchor_points_Bisciaro_al
#'
#' @description
#' data set consist of the tie points between the Bisciaro_al record
#' of Arts (2014) and the la2011 solution of laskar et al., (20111)\cr
#'
#' @details
#'The data set is a matrix with the 4 columns.
#'The first column is the depth/time of the al proxy record tie-points.
#'The second column is the time value of the la2011 astronomical solution tie-points.
#'The third column is the Al value of the a; tie-point.
#'The fourth column is the eccentricity value of the la2011 astronomical solution tie-point.
#'
#'@references
#'M.C.M. Arts, 2014,
#'Magnetostratigrpahy and geochemical analysis of the early Miocene Bisciaro Formation
#'in the Contessa Valley (Northern Italy). Unpublished Bsc. thesis \cr
#'
#'Laskar, J., M. Gastineau, J. B. Delisle, A. Farrés, and A. Fienga (2011b),
#'Strong chaos induced by close encounters with Ceres and Vesta, Astron. Astrophys.,
#'532, L4,<doi:10.1051/0004-6361/201117504> \cr
#'
NULL

#' @title Information of the Geological timescale 2020
#' @rdname GTS_info
#' @name GTS_info
#' @description GTS_info data set consists the information of the Geological
#' timescale 2020 including the color data of Ogg et al., (2021)
#' The ages, durations, uncertainties and colors of the Geological
#' timescale 2020 are included in the data set
#'
#' @details
#' Column 1: name	\cr
#' Column 2: type	\cr
#' Column 1: top age\cr
#' Column 1: top error\cr
#' Column 1: bottom age\cr
#' Column 1: bottom error\cr
#' Column 1: Cyan value\cr
#' Column 1: Magenta value\cr
#' Column 1: Yellow value\cr
#' Column 1: Key  value\cr
#' Column 1: Red Value\cr
#' Column 1: Green value\cr
#' Column 1: Blue value	\cr
#' Column 1: font style	\cr
#' Column 1: font color\cr
#'
#' @references
#'Ogg, Gabi & Ogg, James & Gradstein, Felix. (2021).
#'Recommended color coding of stages - Appendix 1
#'from Geologic Time Scale 2020.
NULL




