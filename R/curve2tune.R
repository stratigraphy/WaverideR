#' @title Convert data from the depth to the time domain
#'
#' @description Converts a data set from the depth to the time domain
#' using a tracked curve/cycle to depth domain an assigning a duration (in kyr) set
#' tracked curve/cycle.
#'
#' @param data Record to convert to the time domain (matrix or data frame, two columns): depth in metres in column 1
#' and the proxy in column 2, normally the data set that was analysed with \code{\link{analyse_wavelet}} to track
#' the cycle. Depths outside the tracked curve get the age of its nearest end. Required (\code{Default=NULL} gives
#' an error).
#' @param tracked_cycle_curve Tracked period curve of one astronomical cycle (matrix or data frame, two columns):
#' depth in metres in column 1 and the period of the tracked cycle in metres in column 2, e.g. the output of
#' \code{\link{track_period_wavelet}} (preferably completed with \code{\link{completed_series}}). Rows with
#' \code{NA} are dropped. Required (\code{Default=NULL} gives an error).
#' @param tracked_cycle_period Duration of the tracked cycle in kyr (single numeric value, > 0, e.g. 405 for the
#' long eccentricity cycle), used to convert the tracked period to a sedimentation rate that is integrated to time.
#' Required (\code{Default=NULL} gives an error).
#' @param genplot Plot three panels stacked on top of each other (logical, TRUE or FALSE): the data in the depth
#' domain, the depth-time curve and the data in the time domain. \code{Default=FALSE}.
#' @param keep_editable Keep the graphical parameters of the plot so that extra elements can be added afterwards
#' (logical, TRUE or FALSE); \code{FALSE} restores the previous \code{par} settings. \code{Default=FALSE}.
#'
#'@author
#' Part of the code is based on the \link[astrochron]{sedrate2time}
#' function of the 'astrochron' R package
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#' @examples
#' \donttest{
#'#The example uses the magnetic susceptibility data set of Pas et al., (2018).
#'# perform the CWT
#'mag_wt <- analyse_wavelet(data = mag,
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
#'#                                   periodlab = "Period (meters)",
#'#                                   x_lab = "depth (meters)")
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
#convert period in meters to sedrate depth vs time
#'mag_track_time<- curve2tune(data=mag,
#'                            tracked_cycle_curve=mag_track_complete,
#'                            tracked_cycle_period=405,
#'                            genplot = FALSE,
#'                            keep_editable=FALSE)
#'}
#'@return
#'The output is a matrix with 2 columns.
#'The first column is time.
#'The second column sedimentation proxy value.
#'
#'If \code{genplot=TRUE} then 3 plots stacked on top of each other will be plotted.
#'Plot 1: the original data set.
#'Plot 2: the depth time plot.
#'Plot 3: the data set in the time domain.
#'
#' @export
#' @importFrom stats approx
#' @importFrom astrochron sedrate2time



curve2tune <- function(data = NULL,
    tracked_cycle_curve = NULL,
    tracked_cycle_period = NULL,
    genplot = FALSE,
    keep_editable = FALSE) {
  my.data <- data
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

  completed_series <- na.omit(out)
  yleft_comp <- completed_series[1, 2]
  yright_com <- completed_series[nrow(completed_series), 2]
  app <- approx(
      x = out[, 1],
      y = out[, 2],
      xout = my.data[, 1],
      method = "linear",
      yleft = yleft_comp,
      yright = yright_com
  )
  completed_series <- as.data.frame(cbind(app$y, my.data[, 2]))

  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    layout.matrix <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
    graphics::layout(
        mat = layout.matrix,
        heights = c(1, 1, 1),
        # Heights of the two rows
        widths = c(1, 1, 1)
    ) # Widths of the two columns
    par(mar = c(4, 4, 2, 2))
    plot(
        x = data[, 1],
        y = data[, 2],
        type = "l",
        main = "Data depth domain",
        xlab = "meters",ylab = "proxy value "
    )

    plot(
        x = out[, 1],
        y = out[, 2],
        type = "l",
        xlab = "meters",
        ylab = "Time (ka)",
        main = "Depth-time plot"
    )
    points(x = out[, 1], y = out[, 2], cex = 1)

    plot(
        completed_series[, 1],
        completed_series[, 2],
        type = "l",
        xlab = "Time (ka)",
        ylab = "proxy value",
        main = "Data time domain"
    )
  }

  return(completed_series)

}
