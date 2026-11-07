#' @title Use a sedimentation curve to convert data to the time domain
#'
#' @description Convert a proxy record from the depth to time domain using
#' a sedimentation rate curve
#'
#' @param data Proxy record in the depth domain (required): matrix or data frame with depth in metres in column 1 and
#' the proxy value in column 2. Its depths need not be those of \code{sed_curve}: the sedimentation rate is
#' interpolated to every depth of the record. \code{Default=NULL}.
#' @param sed_curve Sedimentation rate curve (required): matrix or data frame with depth in metres in column 1 and the
#' sedimentation rate in cm/kyr in column 2, e.g. the output of \code{\link{minimal_tuning}}
#' with \code{output = 1}. It is sorted, cleaned of \code{NA} rows and linearly interpolated to the depths of
#' \code{data} (constant beyond its ends); the inverse rate is then integrated over depth to time (in kyr, 0 at the
#' first sample of \code{data}).
#' \code{Default=NULL}.
#' @param genplot Plot the proxy record against time (kyr) (logical). \code{Default=FALSE}.
#' @param keep_editable Keep the graphical parameters of the plot after the function returns, so that extra elements
#' can be added (logical); only used with \code{genplot = TRUE}. \code{Default=FALSE}.
#'
#'@author
#' Part of the code is based on the \link[astrochron]{sedrate2time}
#' function of the 'astrochron' R package
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#'
#' @examples
#' \donttest{
#'# Extract the 405kyr eccentricity cycle from the wavelet scalogram
#'# from the magnetic susceptibility record of the Sullivan core
#'# of Pas et al., (2018) and then create a age model using minimal tuning
#'# (e.g.) set the distance between peaks to 405 kyr. The age model
#'# (sedimentation rate curve) is then used to convert the data
#'# from the depth to the time domain
#'
#'mag_wt <- analyse_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'
#'mag_405 <- extract_signal_stable_V2(
#'  wavelet = mag_wt,
#'  period_max = 4,
#'  period_min = 2,
#'  add_mean = TRUE,
#'  plot_residual = FALSE,
#'  keep_editable = FALSE
#')
#'
#'mag_405_min_tuning <- minimal_tuning(data = mag_405,
#'pts = 5,
#'cycle = 405,
#'tune_opt = "max",
#'output = 1,
#'genplot = FALSE,
#'keep_editable = FALSE)
#'
#'mag_time <- sedrate2tune(
#'data=mag,
#'sed_curve=mag_405_min_tuning,
#'genplot=FALSE,
#'keep_editable=FALSE)
#'
#'}
#'@return
#'The output is a matrix with 2 columns.
#'The first column is time
#'The second column is the proxy value
#'If \code{genplot=TRUE} then a time vs proxy value plot will be plotted.
#'
#' @export
#' @importFrom astrochron sedrate2time

sedrate2tune <- function(data = NULL,
    sed_curve = NULL,
    genplot = FALSE,
    keep_editable = FALSE) {
  # the sedimentation rate (cm/kyr) at the depths of the record, and time as the integral of its inverse over depth
  # (trapezoid rule; time 0 at the first sample). Formerly the time of a regular grid of the curve was paired row by
  # row with the record, which was only right when both had the same depths.
  data <- as.data.frame(data)[, 1:2]
  data <- data[order(data[, 1]), ]
  sc <- stats::na.omit(as.matrix(sed_curve)[, 1:2])
  sc <- sc[order(sc[, 1]), , drop = FALSE]
  if (nrow(sc) == 1) {
    rate <- rep(sc[1, 2], nrow(data))
  } else {
    rate <- stats::approx(sc[, 1], sc[, 2], xout = data[, 1], rule = 2, ties = mean)$y
  }
  zz <- data[, 1] * 100
  inv <- 1 / rate
  hsum <- c(0, cumsum(diff(zz) * (inv[-1] + inv[-length(inv)]) / 2))
  out <- cbind(hsum, data[, 2])


  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    par(mfrow = c(1, 1))
    plot(
        x = out[, 1],
        y = out[, 2],
        cex = 1,
        xlab = "time (ka)",
        ylab = "proxy"
    )
    lines(out)
  }
  return(out)
}
