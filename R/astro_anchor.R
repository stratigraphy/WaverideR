#'@title Anchor proxy record to an astronomical solution
#'
#' @description Anchor the extracted signal to an astronomical solution using a GUI.
#' The \code{\link{astro_anchor}} function allows one to tie minima or maxima in the proxy record
#' to minima or maxima in an astronomical solution.
#' By tying the proxy record to an astronomical solution one will generate tie-points which
#' can be used to generate a astrochronological age-model
#' As minima or maxima in the proxy record are tied to minima or maxima in an astronomical solution it is
#' important to provide input which has clearly definable minima and maxima.
#' As such input should be of a "sinusoidal" nature otherwise the \code{extract_astrosolution=TRUE}
#' and/or \code{extract_proxy_signal=TRUE} options need to be set to TRUE to create sinusoidal signals.
#'
#'
#'
#'
#' @description Astronomical solutions option are:
#'\itemize{
#'     \item La2004 Eccentricity solution available via the \link[astrochron]{getLaskar} function or downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2004 Obliquity solution available via the \link[astrochron]{getLaskar} function or downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2004 Precession solution available via the \link[astrochron]{getLaskar} function or downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010a Eccentricity solution available via the \link[astrochron]{getLaskar}function or downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010a Obliquity solution downloadable via the \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010a Precession solution downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010b Eccentricity solution available via the \link[astrochron]{getLaskar} function or downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010b Obliquity solution downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010b Precession solution downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010c Eccentricity solution available via the \link[astrochron]{getLaskar} function or downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010c Obliquity solution downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010c Precession solution downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010d Eccentricity solution available via the \link[astrochron]{getLaskar} function or downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010d Obliquity solution downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2010d Precession solution downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item La2011 Eccentricity solution available via the \link[astrochron]{getLaskar} function or downloadable via \url{http://vo.imcce.fr/insola/earth/online/earth/earth.html}
#'	   \item ZB17a Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17a Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17b Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17b Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17c Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17c Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17d Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17d Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17e Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17e Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17f Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17f Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17h Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17h Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17i Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17i Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17j Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17j Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17k Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17k Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17p Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB17p Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB18a Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB18a Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB20a Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB20a Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB20b Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB20b Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB20c Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB20c Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB20d Eccentricity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item ZB20d Obliquity solution downloadable via \url{https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html}
#'	   \item 405kyr eccentricity 405 metronome can be generated using the formula:\cr
#'	     e405=0.027558-0.010739*cos(0.0118+2(pi)*(t/405000)) (laskar et al., 2004 & laskar 2020)
#'	   \item 173kyr obliquity  metronome can be generated using using the formula:\cr
#'	    es3-s6(t) = 0.144*cos(1.961+2(pi)*(t/172800) (laskar et al., 2004 & laskar 2020)
#'	   \item An etp model using the \link[astrochron]{etp} function of the 'astrochron' R package
#'	   }
#'
#'@param astro_solution Input is an astronomical solution which the proxy record will be anchored to,
#' the input should be a matrix or data frame with the first column being
#' age and the second column should be a insolation/angle/value
#'@param proxy_signal Input is the proxy data set which will
#'be anchored to an astronomical solution, the input should be a matrix or
#' data frame with the first column being  depth/time and the second column should be a proxy value.
#' For the best results either the astronomical components need to be pre-extracted
#' before anchoring. This means that a filtering/cycle extracting need to be applied to
#' the input data or the  extract_proxy_signal option needs to be set to TRUE.
#'@param proxy_min_or_max Tune proxy maxima or minima to the astronomical solution \code{Default="max"}.
#'@param clip_astrosolution Clip the astronomical solution \code{Default=FALSE}.
#'@param astrosolution_min_or_max Tune to maximum or minimum values of
#'the astronomical solution \code{Default="max"}
#'@param clip_high Upper value to clip to.
#'@param clip_low Lower value to clip to.
#'@param extract_astrosolution Extract a certain astronomical cycle/component from a
#' astronomical solution prior to anchoring \code{Default=FALSE}.
#'@param astro_period_up Specifies the upper period of the astronomical cycle which is
#'extracted from an astronomical solution. The astro_period_up is a
#'factor with which the astronomical component is multiplied by.  \code{Default=1.2}
#'@param astro_period_down Specified the lower period of the astronomical cycle which is
#'extracted from an astronomical solution. The astro_period_down value is a
#'factor with which the astronomical component is multiplied by. \code{Default=0.8}
#'@param astro_period_cycle Period (in kyr) of the to be extracted astronomical component
#'from the astronomical solution.
#'@param extract_proxy_signal Extract a certain astronomical cycle/component from a
#'proxy signal  \code{Default=FALSE}.
#'@param proxy_period_up Specifies the upper period of the astronomical cycle to be extracted
#' from the proxy record. The upper bound value is a factor with which the
#' (proxy_period_cycle) value is multiplied by. \code{Default=1.2}.
#'@param proxy_period_down Specifies the upper period of the astronomical cycle to be
#'extracted from the proxy record. The lower bound value is a factor with
#' which the astronomical component (proxy_period_cycle) value is multiplied by. \code{Default=0.8}.
#'@param proxy_period_cycle Period in kyr of the astronomical cycle/component which is to be extracted
#' from the proxy record.
#' @param pts The pts parameter specifies how many points to the left/right up/down the peak detect algorithm goes in detecting
#'a peak. The peak detecting algorithm works by comparing the values left/right up/down of it, if the values are both higher or lower
#'then the value a peak. To deal with error produced by this algorithm the pts parameter can be changed which can
#'aid in peak detection. Usually increasing the pts parameter means more peak certainty, however it also means that minor peaks might not be
#'picked up by the algorithm \code{Default=3}
#' @param verbose print text \code{Default=FALSE} set verbose to TRUE to allow for anchoring using text feedback commands
#'@param time_dir The direction of the proxy record which is assumed for tuning if time increases with increasing depth/time values
#'(e.g. bore hole data which gets older with increasing depth ) then time_dir should be set to TRUE
#'if time decreases with depth/time values (eg stratigraphic logs where 0m is the bottom of the section)
#'then time_dir should be set to FALSE \code{time_dir=TRUE}
#'@param genplot Generate plot \code{Default="FALSE"}
#'
#'@references
#'J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004,
#'A long term numerical solution for the insolation quantities of the Earth: Astron. Astrophys.,
#' Volume 428, 261-285. <doi:10.1051/0004-6361:20041335> \cr
#'
#'Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011a,
#' La2010: A new orbital solution for the long-term motion of the Earth: Astron. Astrophys.,
#' Volume 532, A89 <doi:10.1051/0004-6361/201116836> \cr
#'
#'Laskar, J., Gastineau, M., Delisle, J.-B., Farres, A., Fienga, A.:
#'2011b, Strong chaos induced by close encounters with Ceres and Vesta, Astron: Astrophys.,
#'Volume 532, L4.  <doi:10.1051/0004-6361/201117504>\cr
#'
#'J. Laskar,Chapter 4 - Astrochronology,Editor(s): Felix M. Gradstein, James G. Ogg, Mark D. Schmitz, Gabi M. Ogg,Geologic Time Scale 2020,Elsevier,2020,Pages 139-158,ISBN 9780128243602,
#' <doi:10.1016/B978-0-12-824360-2.00004-8>\cr
#'
#'Zeebe, R. E. and Lourens, L. J.
#'Geologically constrained astronomical solutions for the Cenozoic era,
#'Earth and Planetary Science Letters, 2022 <\doi{doi:10.1016/j.epsl.2022.117595}>\cr
#'
#'Richard E. Zeebe Lucas J. Lourens ,Solar System chaos and the Paleocene–Eocene boundary age constrained by geology and astronomy.Science365,926-929(2019)
#'<doi:10.1126/science.aax0612>\cr
#'
#'Zeebe, Richard E. "Numerical solutions for the orbital motion of the Solar System over the past 100 Myr: limits and new results."
#'The Astronomical Journal 154, no. 5 (2017): 193. <doi:10.3847/1538-3881/aa8cce> \cr
#'
#'Stephen R. Meyers,Cyclostratigraphy and the problem of astrochronologic testing,
#'Earth-Science Reviews,Volume 190,2019,Pages 190-223,ISSN 0012-8252
#'<doi:10.1016/j.earscirev.2018.11.015>
#'
#' @return
#'The output is a matrix with the 4 columns.
#'The first column is the depth/time of the proxy tie-point.
#'The second column is the time value of the astronomical solution tie-point.
#'The third column is the proxy value of the proxy tie-point.
#'The fourth column is the proxy/insolation value of the astronomical solution  tie-point.
#'If genplot is set to true then at plot of the of the achored points will be plotted
#'
#'
#'@examples
#'\donttest{
#'# Use the grey_track example tracking points to anchor the grey scale data set
#'# of Zeeden et al., (2013) to the p-0.5t la2004 solution
#'
#'grey_wt <-
#'  analyze_wavelet(
#'    data = grey,
#'    dj = 1/200,
#'    lowerPeriod = 0.02,
#'    upperPeriod = 256,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'#Use the pre-tracked grey_track curve which traced the precession cycle
#'grey_track <- completed_series(
#'  wavelet = grey_wt,
#'  tracked_curve  = grey_track,
#'  period_up  = 1.25,
#'  period_down  = 0.75,
#'  extrapolate = TRUE,
#' genplot = FALSE
#')
#
#'# Extract precession, obliquity and eccentricity to create a synthetic insolation curve
#'
#'grey_prec <- extract_signal(
#'tracked_cycle_curve = grey_track[,c(1,2)],
#'wavelet = grey_wt,
#'period_up = 1.2,
#'period_down = 0.8,
#'add_mean = FALSE,
#'tracked_cycle_period = 22,
#'extract_cycle = 22,
#'tune = FALSE,
#'plot_residual = FALSE
#')
#'
#'grey_obl <- extract_signal(
#'  tracked_cycle_curve = grey_track[,c(1,2)],
#'  wavelet = grey_wt,
#'  period_up = 1.2,
#'  period_down = 0.8,
#'  add_mean = FALSE,
#'  tracked_cycle_period = 22,
#'  extract_cycle = 110,
#'  tune = FALSE,
#'  plot_residual = FALSE
#')
#'
#'grey_ecc <- extract_signal(
#'  tracked_cycle_curve = grey_track[,c(1,2)],
#'  wavelet = grey_wt,
#'  period_up = 1.25,
#'  period_down = 0.75,
#'  add_mean = FALSE,
#'  tracked_cycle_period = 22,
#'  extract_cycle = 40.8,
#'  tune = FALSE,
#'  plot_residual = FALSE
#')
#'
#'insolation_extract <- cbind(grey_ecc[,1],grey_prec[,2]+grey_obl[,2]+grey_ecc[,2]+mean(grey[,2]))
#'insolation_extract <- as.data.frame(insolation_extract)
#'insolation_extract_mins <- min_detect(insolation_extract,pts=3)
#'
#'#use the astrosignal_example to tune to which is an \cr
#'# ETP solution (p-0.5t la2004 solution)
#'astrosignal_example <- na.omit(astrosignal_example)
#'astrosignal_example[,2] <- -1*astrosignal_example[,2]
#'astrosignal <- as.data.frame(astrosignal_example)
#'
#'#anchor the synthetic insolation curve extracted from the grey scale record to the insolation curve.
#'
#'anchor_pts <- astro_anchor(
#'astro_solution = astrosignal,
#'proxy_signal = insolation_extract,
#'proxy_min_or_max = "min",
#'clip_astrosolution = FALSE,
#'astrosolution_min_or_max = "min",
#'clip_high = NULL,
#'clip_low = NULL,
#'extract_astrosolution  = FALSE,
#'astro_period_up  = NULL,
#'astro_period_down  = NULL,
#'astro_period_cycle  = NULL,
#'extract_proxy_signal  = FALSE,
#'proxy_period_up  = NULL,
#'proxy_period_down  = NULL,
#'proxy_period_cycle  = NULL,
#'pts=3,
#'verbose=FALSE, #set verbose to TRUE to allow for anchoring using text feedback commands
#'genplot=FALSE
#')}
#'
#' @export
#' @importFrom graphics points
#' @importFrom grDevices dev.new
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom grDevices graphics.off
#' @importFrom astrochron etp


astro_anchor <- function(astro_solution = NULL,
                         proxy_signal = NULL,
                         proxy_min_or_max = "max",
                         clip_astrosolution = FALSE,
                         astrosolution_min_or_max = "max",
                         clip_high = NULL,
                         clip_low = NULL,
                         extract_astrosolution = FALSE,
                         astro_period_up = 1.2,
                         astro_period_down = 0.8,
                         astro_period_cycle = NULL,
                         extract_proxy_signal = FALSE,
                         proxy_period_up = 1.2,
                         proxy_period_down = 0.8,
                         proxy_period_cycle = NULL,
                         pts=3,
                         verbose = FALSE,
                         time_dir = TRUE,
                         genplot=FALSE) {

  if (clip_astrosolution == TRUE) {
    astro_solution <- astro_solution[astro_solution[, 1]  >= clip_low,]
    astro_solution <-
      astro_solution[astro_solution[, 1]  <= clip_high,]
  }

  if (extract_astrosolution == TRUE) {
    extract_astrosolution_wt <- analyze_wavelet(
      data = astro_solution,
      dj = 1 / 200,
      lowerPeriod = (astro_solution[2, 1] - astro_solution[1, 1]),
      upperPeriod = (astro_solution[nrow(astro_solution), 1] - astro_solution[1, 1]),
      verbose = FALSE,
      omega_nr = 6
    )

    astro_solution <- extract_signal_stable(
      wavelet = extract_astrosolution_wt,
      cycle = astro_period_cycle,
      period_up = astro_period_up,
      period_down = astro_period_down,
      add_mean = TRUE,
      plot_residual = FALSE
    )

  }

  if (extract_proxy_signal == TRUE) {
    proxy_signal_wt <- analyze_wavelet(
      data = astro_solution,
      dj = 1 / 200,
      lowerPeriod = (astro_solution[2, 1] - astro_solution[1, 1]),
      upperPeriod = (astro_solution[nrow(astro_solution), 1] - astro_solution[1, 1]),
      verbose = FALSE,
      omega_nr = 6
    )

    proxy_signal <- extract_signal_stable(
      wavelet = proxy_signal_wt,
      cycle = proxy_period_cycle,
      period_up = proxy_period_up,
      period_down = proxy_period_down,
      add_mean = TRUE,
      plot_residual = FALSE
    )
  }


  astro_solution <- as.data.frame(astro_solution)
  if (astrosolution_min_or_max == "max") {
    astro_maxdetect_error_corr <- max_detect(data=astro_solution,pts=pts)
  }


  if (astrosolution_min_or_max == "min") {
    astro_maxdetect_error_corr <- min_detect(data=astro_solution,pts=pts)
  }

  all_data_pre_max <- as.data.frame(proxy_signal)

  if (proxy_min_or_max == "max") {
    proxy_signal_max <- max_detect(all_data_pre_max,pts=pts)
  }

  if (proxy_min_or_max == "min") {
    proxy_signal_max <- min_detect(all_data_pre_max,pts=pts)
  }


  if (time_dir == TRUE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1, 1))
    proxy_signal_max <- proxy_signal_max[, c(1, 2)]
    plot(proxy_signal,
         type = "l",
         xlab = "depth/time",
         ylab = "proxy")
    points(proxy_signal_max,
           xlab = "depth/time",
           ylab = "proxy")


    tie_points <- matrix(data = NA, ncol = 4)
    colnames(tie_points) <-
      c("data_x", "data_y", "insolation_x", "insolation_y")
    tie_points <- tie_points[-c(1), ]

    for (i in 1:nrow(proxy_signal_max)) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      next_step <- "NO"
      while (next_step == "NO") {
        next_step <- "NO"
        dev.new(width = 20,
                height = 10,
                noRStudioGD = TRUE)
        par(mfrow = c(2, 1))
        plot(
          all_data_pre_max,
          type = "l",
          xlab = "depth/time",
          ylab = "proxy"
        )
        points(
          proxy_signal_max,
          type = "p",
          pch = 1,
          col = "black",
          lwd = "1",
          xlab = "depth/time",
          ylab = "proxy"
        )
        par(new = TRUE)
        plot(
          x = proxy_signal_max[i, 1],
          y = proxy_signal_max[i, 2],
          type = "p",
          pch = 19,
          cex = 1.5,
          col = "red",
          lwd = "1",
          xlim = c(min(all_data_pre_max[, 1]), max(all_data_pre_max[, 1])),
          ylim = c(min(all_data_pre_max[, 2]), max(all_data_pre_max[, 2])),
          xlab = "depth/time",
          ylab = "proxy"
        )
        par(new = TRUE)
        plot(
          x = tie_points[, 1],
          y = tie_points[, 2],
          type = "p",
          pch = 19,
          cex = 1.5,
          col = "green",
          lwd = "1",
          xlim = c(min(all_data_pre_max[, 1]), max(all_data_pre_max[, 1])),
          ylim = c(min(all_data_pre_max[, 2]), max(all_data_pre_max[, 2])),
          xlab = "depth/time",
          ylab = "proxy"
        )
        plot(
          astro_solution,
          type = "l",
          xlim = c(min(astro_solution[, 1]), max(astro_solution[, 1])),
          ylim = c(min(astro_solution[, 2]), max(astro_solution[, 2])),
          xlab = "depth/time",
          ylab = "tie_point value"
        )
        par(new = TRUE)
        plot(
          astro_maxdetect_error_corr[, c(1, 2)],
          type = "p",
          xlim = c(min(astro_solution[, 1]), max(astro_solution[, 1])),
          ylim = c(min(astro_solution[, 2]), max(astro_solution[, 2])),
          xlab = "depth/time",
          ylab = "tie_point value"
        )
        par(new = TRUE)
        plot(
          tie_points[, c(3, 4)],
          type = "p",
          pch = 19,
          cex = 1.5,
          col = "green",
          lwd = "1",
          xlim = c(min(astro_solution[, 1]), max(astro_solution[, 1])),
          ylim = c(min(astro_solution[, 2]), max(astro_solution[, 2])),
          xlab = "depth/time",
          ylab = "tie_point value"
        )
        pts <- tuning_pts(x = astro_maxdetect_error_corr[, 1],
                          y = astro_maxdetect_error_corr[, 2])


        if (verbose == TRUE) {
          var = readline(prompt <-
                           "did you select the right point Y/N : ")
        }
        else
          (var = "Y")

        if (var == "Y" |
            var == "Yes" | var == "YES" | var == "yes" | var == "y")
        {
          if (verbose == TRUE) {
            cat("okay next point")
          }
          if (identical(pts, integer(0))) {
            if (verbose == TRUE) {
              cat("no points selected on to the next point")
            }
            sel_astro_maxdetect_error_corr <-
              matrix(data = NA, ncol = 2)
            sel_proxy_signal_max <-
              (as.data.frame(proxy_signal_max[i,]))
            tie_row <-
              cbind(sel_proxy_signal_max,
                    sel_astro_maxdetect_error_corr)
            colnames(tie_row) <-
              c("data_x",
                "data_y",
                "insolation_x",
                "insolation_y")
            tie_points <- rbind(tie_points, tie_row)
            graphics.off()
            next_step <- "YES"
          }
          else
          {
            sel_astro_maxdetect_error_corr <-
              as.data.frame(astro_maxdetect_error_corr[pts, c(1, 2)])
            sel_proxy_signal_max <-
              (as.data.frame(proxy_signal_max[i,]))
            tie_row <-
              cbind(sel_proxy_signal_max,
                    sel_astro_maxdetect_error_corr)
            colnames(tie_row) <-
              c("data_x",
                "data_y",
                "insolation_x",
                "insolation_y")
            tie_points <- rbind(tie_points, tie_row)
            graphics.off()
            next_step <- "YES"
          }
        }
        else if (var == "N" |
                 var == "NO" |
                 var == "No" | var == "no" | var == "n")
        {
          if (verbose == TRUE) {
            cat("okay selection of point will be redone")
          }
          graphics.off()
          next_step <- "NO"

        }
        else
        {
          if (verbose == TRUE) {
            cat("You did not type yes or no quiting of selection process")
          }
          graphics.off()
          next_step <- "Finished"

        }

      }
    }
  }
  else{
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1, 1))
    proxy_signal_max <- proxy_signal_max[, c(1, 2)]
    plot(
      proxy_signal,
      type = "l",
      xlab = "depth/time",
      ylab = "proxy",
      xlim = rev(c(
        min(proxy_signal[, 1]), max(proxy_signal[, 1])
      ))
    )
    points(proxy_signal_max,
           xlab = "depth/time",
           ylab = "proxy")

    tie_points <- matrix(data = NA, ncol = 4)
    colnames(tie_points) <-
      c("data_x", "data_y", "insolation_x", "insolation_y")
    tie_points <- tie_points[-c(1), ]

    for (i in 1:nrow(proxy_signal_max)) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      next_step <- "NO"
      while (next_step == "NO") {
        next_step <- "NO"
        dev.new(width = 20,
                height = 10,
                noRStudioGD = TRUE)
        par(mfrow = c(2, 1))
        plot(
          all_data_pre_max,
          type = "l",
          xlab = "depth/time",
          ylab = "proxy",
          xlim = rev(c(
            min(all_data_pre_max[, 1]), max(all_data_pre_max[, 1])
          ))
        )
        points(
          proxy_signal_max,
          type = "p",
          pch = 1,
          col = "black",
          lwd = "1",
          xlab = "depth/time",
          ylab = "proxy"
        )
        par(new = TRUE)
        plot(
          x = proxy_signal_max[i, 1],
          y = proxy_signal_max[i, 2],
          type = "p",
          pch = 19,
          cex = 1.5,
          col = "red",
          lwd = "1",
          xlim = rev(c(
            min(all_data_pre_max[, 1]), max(all_data_pre_max[, 1])
          )),
          ylim = c(min(all_data_pre_max[, 2]), max(all_data_pre_max[, 2])),
          xlab = "depth/time",
          ylab = "proxy"
        )
        par(new = TRUE)
        plot(
          x = tie_points[, 1],
          y = tie_points[, 2],
          type = "p",
          pch = 19,
          cex = 1.5,
          col = "green",
          lwd = "1",
          xlim = rev(c(
            min(all_data_pre_max[, 1]), max(all_data_pre_max[, 1])
          )),
          ylim = c(min(all_data_pre_max[, 2]), max(all_data_pre_max[, 2])),
          xlab = "depth/time",
          ylab = "proxy"
        )
        plot(
          astro_solution,
          type = "l",
          xlim = c(min(astro_solution[, 1]), max(astro_solution[, 1])),
          ylim = c(min(astro_solution[, 2]), max(astro_solution[, 2])),
          xlab = "depth/time",
          ylab = "tie_point value"
        )
        par(new = TRUE)
        plot(
          astro_maxdetect_error_corr[, c(1, 2)],
          type = "p",
          xlim = c(min(astro_solution[, 1]), max(astro_solution[, 1])),
          ylim = c(min(astro_solution[, 2]), max(astro_solution[, 2])),
          xlab = "depth/time",
          ylab = "tie_point value"
        )
        par(new = TRUE)
        plot(
          tie_points[, c(3, 4)],
          type = "p",
          pch = 19,
          cex = 1.5,
          col = "green",
          lwd = "1",
          xlim = c(min(astro_solution[, 1]), max(astro_solution[, 1])),
          ylim = c(min(astro_solution[, 2]), max(astro_solution[, 2])),
          xlab = "depth/time",
          ylab = "tie_point value"
        )


        x = astro_maxdetect_error_corr[, 1]
        y = astro_maxdetect_error_corr[, 2]

        defaultW <- getOption("warn")
        options(warn = -1)
        xy <- xy.coords(x, y)
        x <- xy$x
        y <- xy$y
        sel <- rep(FALSE, length(x))
        res <- integer(0)
        ans <- identify(x[!sel], y[!sel], n = 1, plot = F)
        ans <- which(!sel)[ans]
        points(
          x[ans],
          y[ans],
          pch = 19,
          col = "red",
          xlim = rev(c(
            min(all_data_pre_max[, 1]), max(all_data_pre_max[, 1])
          )),
          ylim = c(min(all_data_pre_max[, 2]), max(all_data_pre_max[, 2]))
        )
        sel[ans] <- TRUE
        res <- c(res, ans)
        pts <- res
        options(warn = defaultW)

        if (verbose == TRUE) {
          var = readline(prompt <-
                           "did you select the right point Y/N : ")
        }
        else
          (var = "Y")

        if (var == "Y" |
            var == "Yes" | var == "YES" | var == "yes" | var == "y")
        {
          if (verbose == TRUE) {
            cat("okay next point")
          }
          if (identical(pts, integer(0))) {
            if (verbose == TRUE) {
              cat("no points selected on to the next point")
            }
            sel_astro_maxdetect_error_corr <-
              matrix(data = NA, ncol = 2)
            sel_proxy_signal_max <-
              (as.data.frame(proxy_signal_max[i,]))
            tie_row <-
              cbind(sel_proxy_signal_max,
                    sel_astro_maxdetect_error_corr)
            colnames(tie_row) <-
              c("data_x",
                "data_y",
                "insolation_x",
                "insolation_y")
            tie_points <- rbind(tie_points, tie_row)
            graphics.off()
            next_step <- "YES"
          }
          else
          {
            sel_astro_maxdetect_error_corr <-
              as.data.frame(astro_maxdetect_error_corr[pts, c(1, 2)])
            sel_proxy_signal_max <-
              (as.data.frame(proxy_signal_max[i,]))
            tie_row <-
              cbind(sel_proxy_signal_max,
                    sel_astro_maxdetect_error_corr)
            colnames(tie_row) <-
              c("data_x",
                "data_y",
                "insolation_x",
                "insolation_y")
            tie_points <- rbind(tie_points, tie_row)
            graphics.off()
            next_step <- "YES"
          }
        }
        else if (var == "N" |
                 var == "NO" |
                 var == "No" | var == "no" | var == "n")
        {
          if (verbose == TRUE) {
            cat("okay selection of point will be redone")
          }
          graphics.off()
          next_step <- "NO"

        }
        else
        {
          if (verbose == TRUE) {
            cat("You did not type yes or no quiting of selection process")
          }
          graphics.off()
          next_step <- "Finished"

        }

      }
    }
  }

  tie_points <- na.omit(tie_points)
  tie_points <- tie_points[, c(1, 3, 2, 4)]

  if(genplot==TRUE){

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    proxy_signal <- all_data_pre_max
    anchor_pts <-  tie_points


    if (time_dir == TRUE){
      dev.new(width = 20,
              height = 10,
              noRStudioGD = TRUE)
      par(mfrow = c(2, 1),mai = c(0.5/2.54, 2/2.54, 2/2.54, 1/2.54),mgp=c(2,1,0))
      plot(x=proxy_signal[,1],
           y=proxy_signal[,2],
           type = "l",
           xlab = "",
           ylab = "Proxy value",xaxt = "n", xaxs = "i",lwd=2,
           ylim= c(min(proxy_signal[,2]),max(proxy_signal[,2])),
           xlim = (c(
             min(proxy_signal[, 1]), max(proxy_signal[, 1])
           )))
      mtext(text = "Depth (metres)",
            side = 3, #side 2 = left
            line = 2)
      box(lwd=2)
      axis(3)
      segments(x0=anchor_pts[,1],
               y0=rep(min(proxy_signal[,2])*0.2,times=(nrow(anchor_pts))),
               x1 = anchor_pts[,1],
               y1 =anchor_pts[,3],
               ylim= c(min(proxy_signal[,2]),max(proxy_signal[,2])),
               xlim = (c(
                 min(proxy_signal[, 1]), max(proxy_signal[, 1])
               )),col = "black", lwd = 2)

      points(x=anchor_pts[,1],y=anchor_pts[,3],col="green",pch=19,
             ylim= c(min(proxy_signal[,2]),max(proxy_signal[,2])),
             xlim = (c(
               min(proxy_signal[, 1]), max(proxy_signal[, 1])
             )),cex=2)




      par(new = FALSE,mai = c(2/2.54, 2/2.54, 0.5/2.54 , 1/2.54),mgp=c(2,1,0))


      plot(
        astro_solution,
        type = "l",
        xlab = "Time",
        ylab = "tie point value",xaxs = "i",yaxs = "i",lwd=2
      )
      box(lwd=2)

      segments(x0=anchor_pts[,2], y0=anchor_pts[,4], x1 = anchor_pts[,2],
               y1 = rep(max(astro_solution[,2])*1.2,times=(nrow(anchor_pts))),
               col = "black", lwd = 2)
      points(x=anchor_pts[,2],y=anchor_pts[,4],col="red",pch=19,cex=2)



      par(new = TRUE,mfrow = c(1, 1),mai = c(2/2.54, 2/2.54, 2/2.54, 1/2.54))
      plot.new()
      settings_dev <- dev.size("in")

      plot.window(xlim = c(min(astro_solution[,1]) , max(astro_solution[,1])),
                  ylim = c(2/2.54, settings_dev[2]-(3/2.54)), yaxs = "i", xaxs = "i")

      y0_val <- ((settings_dev[2]-((5)/2.54))/2)+1.5/2.54
      y1_val <-  ((settings_dev[2]-((5)/2.54))/2)+(2.5)/2.54


      fact <- (anchor_pts[,1]-min(proxy_signal[, 1]))/
        (max(proxy_signal[, 1])-min(proxy_signal[, 1]))
      x1_vals <- ((max(astro_solution[, 1])-min(astro_solution[, 1]))*(fact))+min(astro_solution[, 1])


      segments(x0=anchor_pts[,2],
               y0=rep(y0_val,times=(nrow(anchor_pts)))
               ,x1 = x1_vals,
               y1 = rep(y1_val,times=(nrow(anchor_pts))),
               col = "black", lwd = 2,lend=1)

    }else{
      dev.new(width = 20,
              height = 10,
              noRStudioGD = TRUE)
      par(mfrow = c(2, 1),mai = c(0.5/2.54, 2/2.54, 2/2.54, 1/2.54),mgp=c(2,1,0))
      plot(x=proxy_signal[,1],
           y=proxy_signal[,2],
           type = "l",
           xlab = "",
           ylab = "Proxy value",xaxt = "n", xaxs = "i",lwd=2,
           ylim= c(min(proxy_signal[,2]),max(proxy_signal[,2])),
           xlim = rev(c(
             min(proxy_signal[, 1]), max(proxy_signal[, 1])
           )))
      mtext(text = "Depth (metres)",
            side = 3, #side 2 = left
            line = 2)
      box(lwd=2)
      axis(3)


      segments(x0=anchor_pts[,1],
               y0=rep(min(proxy_signal[,2])*0.2,times=(nrow(anchor_pts))),
               x1 = anchor_pts[,1],
               y1 =anchor_pts[,3],
               ylim= c(min(proxy_signal[,2]),max(proxy_signal[,2])),
               xlim = rev(c(
                 min(proxy_signal[, 1]), max(proxy_signal[, 1])
               )),col = "black", lwd = 2)

      points(x=anchor_pts[,1],y=anchor_pts[,3],col="green",pch=19,
             ylim= c(min(proxy_signal[,2]),max(proxy_signal[,2])),
             xlim = rev(c(
               min(proxy_signal[, 1]), max(proxy_signal[, 1])
             )),cex=2)




      par(new = FALSE,mai = c(2/2.54, 2/2.54, 0.5/2.54 , 1/2.54),mgp=c(2,1,0))


      plot(
        astro_solution,
        type = "l",
        xlab = "Time",
        ylab = "tie point value",xaxs = "i",yaxs = "i",lwd=2
      )
      box(lwd=2)

      segments(x0=anchor_pts[,2], y0=anchor_pts[,4], x1 = anchor_pts[,2],
               y1 = rep(max(astro_solution[,2])*1.2,times=(nrow(anchor_pts))),
               col = "black", lwd = 2)
      points(x=anchor_pts[,2],y=anchor_pts[,4],col="red",pch=19,cex=2)

      par(new = TRUE,mfrow = c(1, 1),mai = c(2/2.54, 2/2.54, 2/2.54, 1/2.54))
      plot.new()
      settings_dev <- dev.size("in")

      plot.window(xlim = c(min(astro_solution[,1]) , max(astro_solution[,1])),
                  ylim = c(2/2.54, settings_dev[2]-(3/2.54)), yaxs = "i", xaxs = "i")

      y0_val <- ((settings_dev[2]-((5)/2.54))/2)+1.5/2.54
      y1_val <-  ((settings_dev[2]-((5)/2.54))/2)+(2.5)/2.54


      fact <- (anchor_pts[,1]-min(proxy_signal[, 1]))/
        (max(proxy_signal[, 1])-min(proxy_signal[, 1]))
      x1_vals <- ((max(astro_solution[, 1])-min(astro_solution[, 1]))*(1-fact))+min(astro_solution[, 1])


      segments(x0=anchor_pts[,2],
               y0=rep(y0_val,times=(nrow(anchor_pts)))
               ,x1 = x1_vals,
               y1 = rep(y1_val,times=(nrow(anchor_pts))),
               col = "black", lwd = 2,lend=1)
    }


    }




  return(tie_points)


}
