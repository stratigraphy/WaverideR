#' @title Calculate the uncertainty associated with the wavelet analysis based on the Gabor uncertainty principle
#'
#' @description The \code{\link{wavelet_uncertainty}} function is used to calculate uncertainties associated
#' with the wavelet analysis based on the Gabor uncertainty principle applied to the
#' continuous wavelet transform using a Morlet wavelet. The calculated uncertainty is the underlying
#' analytical uncertainty which is the result of applying the Gabor uncertainty principle to the
#' continuous wavelet transform using a Morlet wavelet.
#'
#'
#' @param tracked_cycle Curve of the cycle tracked using the \code{\link{track_period_wavelet}} function
#' Any input (matrix or data frame) in which the first column is depth or
#' time and the second column is period should work
#' @param period_of_tracked_cycle period of the tracked curve (in kyr).
#' @param wavelet wavelet object created using the \code{\link{analyze_wavelet}} function.
#' @param multi multiple of the standard deviation to be used for defining uncertainty \code{Default=1}.
#' @param verbose Print text \code{Default=FALSE}.
#' @param genplot_time plot time curves with a upper and lower uncertainty based on Gabor uncertainty principle applied to the
#' continuous wavelet transform using a Morlet wavelet, which uses which uses the omega number (number
#' of cycles in the wavelet) at one standard deviation to define the analytical uncertainty \code{Default=TRUE}
#' @param genplot_uncertainty Plot period curves with upper and lower uncertainty based on Gabor uncertainty principle applied to the
#' continuous wavelet transform using a Morlet wavelet, which uses which uses the omega number
#' (number of cycles in the wavelet) to define uncertainty at one standard deviation \code{Default=TRUE}
#' @param genplot_uncertainty_wt generate a wavelet plot with the uncertainty based on Gabor uncertainty
#' principle applied to the continuous wavelet transform using a Morlet wavelet superimposed on top of
#'original wavelet plot. The red curve is period of the tracked curve plus the analytical uncertainty.
#'The blue curve is period of the tracked curve min the analytical uncertainty.
#'The  black curve is the curve tracked using the '\code{Default=tracked_cycle_curve} function \code{Default=TRUE}
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#'
#'@return Results pertaining to the uncertainty calculated based on the Gabor uncertainty principle.\cr
#'If the genplot_time is TRUE then a depth time plot will be plotted with 3 lines, the mean age,age plus
#' x times the standard deviation and age minus x times the standard deviation . \cr
#'If the genplot_uncertainty is TRUE then a curve will be plotted with the mean period, the tracked period plus
#'x times the standard deviation and the tracked period minus x times the standard deviation. \cr
#'If the genplot_uncertainty_wt is TRUE a wavelet spectra will be plotted with the tracked period, the tracked period plus
#'x times the standard deviation,the tracked period minus x times the standard deviation and the area in between will be shaded in grey.\cr
#'
#'Returns a matrix with 8 columns.\cr
#'The first column is called "depth" eg. depth \cr
#'The second column is "period" of the originally tracked period. \cr
#'The third column is "frequency" of the originally tracked period. \cr
#'The fourth column "uncertainty in frequency FWHM" is the uncertainty in frequency based on the Gabor uncertainty principle defined as
#' (FWHM) full width at half maximum. \cr
#'The fifth column "uncertainty in frequency x_times SD" is the uncertainty in frequency based on the Gabor uncertainty principle defined as
#' times x standard deviations. \cr
#'The sixth column "time mean" is the mean time based on the tracked period. \cr
#'The seventh column "time plus x_times sd" is the time based on the tracked period plus x times the standard deviation. \cr
#'The eight column "time min x_times sd" is the time based on the tracked period min x times the standard deviation. \cr
#'
#' @author
#' Code based on the \link[WaveletComp]{analyze.wavelet} function of the 'WaveletComp' R package
#' and \link[biwavelet]{wt} function of the 'biwavelet' R package which are based on the
#' wavelet 'MATLAB' code written by Christopher Torrence and Gibert P. Compo (1998).
#' The assignment of the standard deviation of the uncertainty of the wavelet
#' is based on the work of Gabor (1946) and Russell et al., (2016)
#'
#' @references
#'Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
#'Wavelet Analysis. R package version 1.1.
#'\url{https://CRAN.R-project.org/package=WaveletComp}
#'
#'Gouhier TC, Grinsted A, Simko V (2021). R package biwavelet: Conduct Univariate and Bivariate Wavelet Analyses. (Version 0.20.21),
#'\url{https://github.com/tgouhier/biwavelet}
#'
#'Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#'Bulletin of the American Meteorological Society 79:61-78.
#'\url{https://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf}
#'
#'Gabor, Dennis. "Theory of communication. Part 1: The analysis of information."
#' Journal of the Institution of Electrical Engineers-part III: radio and
#' communication engineering 93, no. 26 (1946): 429-441.\url{http://genesis.eecg.toronto.edu/gabor1946.pdf}
#'
#'Russell, Brian, and Jiajun Han. "Jean Morlet and the continuous wavelet transform.
#'" CREWES Res. Rep 28 (2016): 115. \url{https://www.crewes.org/Documents/ResearchReports/2016/CRR201668.pdf}
#'
#'
#'Morlet, Jean, Georges Arens, Eliane Fourgeau, and Dominique Glard.
#'"Wave propagation and sampling theory—Part I: Complex signal and scattering in multilayered media.
#'" Geophysics 47, no. 2 (1982): 203-221.
#' \doi{<doi:10.1190/1.1441328>}
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#'  Geophysics 1982 47 (2): 222–236. \doi{<doi:10.1190/1.1441329>}
#'
#'
#'@examples
#'\donttest{
#'#calculate the Gabor uncertainty derived mathematical uncertainty of the
#'#magnetic susceptibility record of the Sullivan core of Pas et al., (2018)
#'
#'mag_wt <- analyze_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'#Track the 405 kyr eccentricity cycle in a wavelet spectra
#'
#'#mag_track <- track_period_wavelet(astro_cycle = 405,
#'#                                   wavelet=mag_wt,
#'#                                   n.levels = 100,
#'#                                   periodlab = "Period (metres)",
#'#                                   x_lab = "depth (metres)")
#'
#'#Instead of tracking, the tracked solution data set mag_track_solution is used
#'mag_track <- mag_track_solution
#'
#'mag_track_complete <- completed_series(
#'   wavelet = mag_wt,
#'   tracked_curve = mag_track,
#'   period_up = 1.2,
#'   period_down = 0.8,
#'   extrapolate = FALSE,
#'   genplot = FALSE,
#'   keep_editable=FALSE
#' )
#'
#' mag_track_complete <- loess_auto(time_series = mag_track_complete,
#'  genplot = FALSE, print_span = FALSE,keep_editable=FALSE)
#'
#'uncertainty <- wavelet_uncertainty(
#'   tracked_cycle = mag_track_complete,
#'   period_of_tracked_cycle = 405,
#'   wavelet = mag_wt,
#'   multi=1,
#'   verbose = FALSE,
#'   genplot_time = FALSE,
#'   genplot_uncertainty = FALSE,
#'   genplot_uncertainty_wt = FALSE,
#'   keep_editable=FALSE
#' )
#'}
#' @export
#' @importFrom grDevices dev.new
#' @importFrom WaveletComp analyze.wavelet
#' @importFrom biwavelet wt




wavelet_uncertainty <- function(tracked_cycle = NULL,
                                period_of_tracked_cycle = NULL,
                                wavelet = NULL,
                                multi = 1,
                                verbose = FALSE,
                                genplot_time = FALSE,
                                genplot_uncertainty = FALSE,
                                genplot_uncertainty_wt = FALSE,
                                keep_editable = FALSE) {
  data <- tracked_cycle[, c(1, 2)]
  ncycles <- wavelet$omega_nr
  b <- (2 * sqrt(2 * log(2)))
  a <- ((8 * log(2) / (2 * pi)))
  k <- (ncycles / (8 * log(2))) * 2
  data$f0 <- (1 / data[, 2])
  data$df <- (a * data$f0) / k
  data$sd_morlet <- data$df / b

  sedrate_min <- 1 / (data$f0 - data$sd_morlet * multi)
  sedrate_plus <-  1 / (data$f0 + data$sd_morlet * multi)


  sedrate_min <- cbind(data[, 1], sedrate_min)
  sedrate_plus <- cbind(data[, 1], sedrate_plus)
  sedrate <-  cbind(data[, 1], 1 / (data$f0))

  time_min <- curve2time(
    tracked_cycle_curve = sedrate_min,
    tracked_cycle_period = period_of_tracked_cycle,
    genplot = FALSE
  )

  time_plus <- curve2time(
    tracked_cycle_curve = sedrate_plus,
    tracked_cycle_period = period_of_tracked_cycle,
    genplot = FALSE
  )

  time <- curve2time(
    tracked_cycle_curve = sedrate,
    tracked_cycle_period = period_of_tracked_cycle,
    genplot = FALSE
  )


  if (verbose == TRUE) {
    cat("total duration is: ", time[nrow(time), 2])
    cat("uncertainty is: ", time_plus[nrow(time_plus), 2] - time[nrow(time), 2])
    cat("uncertainty is defined as ", multi, " standard deviations ")
    cat("uncertainty defined as standard deviation is: ",
        ((a / k) / b),
        " cycle per cycle")

  }


  if (genplot_time == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    dev.new(width = 7,
            height = 7,
            noRStudioGD = TRUE)
    plot(
      time_min,
      type = "l",
      col = "red",
      ylim = c(0, time_plus[nrow(time_plus), 2])
    )
    lines(time_plus, col = "blue")
    lines(time, col = "black", lwd = 3)
  }



  if (genplot_uncertainty == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    dev.new(width = 7,
            height = 7,
            noRStudioGD = TRUE)
    max_y <- max(1 / (data$f0 - data$sd_morlet * multi))
    min_y <- min(1 / (data$f0 + data$sd_morlet * multi))

    plot(
      data[, 1],
      1 / data$f0,
      type = "l",
      ylim = c(min_y, max_y),
      lwd = 2,
      xlab = "depth",
      ylab = "period (m)"
    )
    lines(data[, 1],
          1 / (data$f0 - data$sd_morlet * multi),
          col = "red",
          lwd = 2)
    lines(data[, 1],
          1 / (data$f0 + data$sd_morlet * multi),
          col = "blue",
          lwd = 2)
  }


  if (genplot_uncertainty_wt == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    plot_wavelet(
      wavelet = wavelet,
      plot.COI = TRUE,
      n.levels = 100,
      color.palette = "rainbow(n.levels, start = 0, end = 0.7)",
      useRaster = TRUE,
      periodlab = "Period (metres)",
      x_lab = "depth (metres)"
    )

    lines(data[, 1], log2(1 / data$f0), lwd = 2)
    lines(data[, 1], log2(1 / (data$f0 - data$sd_morlet * multi)), col =
            "red", lwd = 2)
          lines(data[, 1], log2(1 / (data$f0 + data$sd_morlet * multi)), col =
                  "blue", lwd = 2)

          combined_sedrate <-
            cbind(data[, 1], data$f0, data$sd_morlet * multi)

          xcords <-
            c(combined_sedrate[, 1],
              sort(combined_sedrate[, 1], decreasing = TRUE))
          xcords
          data_sort1 <-
            combined_sedrate[order(combined_sedrate[, 1], decreasing = TRUE), ]
          ycords <-
            c(1 / (combined_sedrate[, 2] + combined_sedrate[, 3]),
              1 / (data_sort1[, 2] - data_sort1[, 3]))

          polygon(
            x = xcords,
            y = log2(ycords),
            col = rgb(0.5, 0.5, 0.5, 0.5),
            border = "black"
          )

  }

  data$time_mean <- time[, 2]
  data$time_min_sd <- time_min[, 2]
  data$time_plus_sd <- time_plus[, 2]
  colnames(data) <- c(
    "depth",
    "period",
    "frequency",
    "uncertainty in frequency FWHM",
    "uncertainty in frequency x_times SD",
    "time mean",
    "time plus x_times sd",
    "time min x_times sd"
  )
return(data)
}
