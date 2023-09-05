#' @title Convert the tracked curve to a depth time space
#'
#' @description Converts the tracked curve to a depth time space.
#'
#' @param tracked_cycle_curve  Curve of the cycle tracked using the
#' \code{\link{track_period_wavelet}} function \cr
#' Any input (matrix or data frame) in which the first column is depth in
#'  meters and the second column is period in meters can be used.
#' @param tracked_cycle_period Period of the tracked curve in kyr.
#' @param genplot Generates a plot with depth vs time \code{Default=FALSE}.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#' @author
#'Based on the \link[astrochron]{sedrate2time}
#'function of the 'astrochron' R package
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#' @examples
#' \donttest{
#' #Convert a tracked curve to a depth time space. The examples uses the
#' #magnetic susceptibility data set of Pas et al., (2018).
#'
#'#'# perform the CWT
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
#' mag_track_complete <- completed_series(
#'   wavelet = mag_wt,
#'   tracked_curve = mag_track,
#'   period_up = 1.2,
#'   period_down = 0.8,
#'   extrapolate = TRUE,
#'   genplot = FALSE
#' )
#'
#'# smooth the tracking of the 405 kyr eccentricity cycle
#' mag_track_complete <- loess_auto(time_series = mag_track_complete,
#' genplot = FALSE, print_span = FALSE)
#'
#'#convert period in meters to sedrate depth vs time
#'mag_track_time<- curve2time(tracked_cycle_curve=mag_track_complete,
#'tracked_cycle_period=405,
#'genplot=FALSE,
#'keep_editable=FALSE)
#'
#'}
#'@return
#'The output is a matrix with 2 columns.
#'The first column is depth.
#'The second column sedimentation rate in cm/kyr.
#'If \code{genplot=TRUE} then a depth vs time plot will be plotted.
#'
#' @export
#' @importFrom astrochron sedrate2time

curve2time <- function(tracked_cycle_curve = NULL,
                       tracked_cycle_period = NULL,
                       genplot = FALSE,
                       keep_editable = FALSE) {
  tracked_cycle_curve[, 2] <-
    tracked_cycle_curve[, 2] / (tracked_cycle_period / 100)
  sedrates <- data.frame(tracked_cycle_curve)
  dat <- as.matrix(tracked_cycle_curve)
  dat <- na.omit(dat)
  dat <- dat[order(dat[, 1], na.last = NA, decreasing = F),]
  npts <- length(dat[, 1])
  start <- dat[1, 1]
  end <- dat[length(dat[, 1]), 1]
  x1 <- dat[1:(npts - 1), 1]
  x2 <- dat[2:(npts), 1]
  dx = x2 - x1
  dt = median(dx)
  xout <- seq(start, end, by = dt)
  npts <- length(xout)
  interp <- approx(dat[, 1], dat[, 2], xout, method = "linear",
                   n = npts)
  sedrates <- as.data.frame(interp)

  npts <- length(sedrates[, 1])
  sedrates[1] = sedrates[1] * 100
  sedrates[2] = 1 / sedrates[2]
  dx = sedrates[2, 1] - sedrates[1, 1]
  midptx = (sedrates[2:npts, 1] + sedrates[1:(npts - 1), 1]) / 2
  slope = (sedrates[2:npts, 2] - sedrates[1:(npts - 1), 2]) / dx
  yint = sedrates[2:npts, 2] - (slope * sedrates[2:npts, 1])
  midpty = (slope * midptx) + yint
  hsum = cumsum(midpty * dx)
  hsum = append(0, hsum)
  out = data.frame(cbind(sedrates[, 1] / 100, hsum))
  colnames(out) <- c("meters", "ka")

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
      xlab = "meters",
      ylab = "Time (ka)",
      main = "Depth time"
    )
    lines(out)
  }
  return(out)
}
