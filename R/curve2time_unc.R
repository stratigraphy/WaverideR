#' @title Convert the tracked curve to a depth time space with uncertainty
#'
#' @description Converts the tracked curve to a depth time space while also
#' taking into account the uncertainty of the tracked astronomical cycle
#'
#' @param tracked_cycle_curve  Curve of the cycle tracked using the
#' \code{\link{track_period_wavelet}} function \cr
#' Any input (matrix or data frame) in which the first column is depth in
#'  meters and the second column is period in meters can be used.
#'@param tracked_cycle_period Period of the tracked curve in kyr.
#'@param tracked_cycle_period_unc uncertainty in the period of the tracked cycle
#'@param tracked_cycle_period_unc_dist distribution of the uncertainty of the
#'tracked cycle value need to be either "u" for uniform distribution or
#'"n" for normal distribution  \code{Default="n"}
#'@param n_simulations number of time series to be modeled
#'@param output If output = 1 a matrix with the predicted ages given the input for each run
#'is given. If output = 2 a matrix with 6 columns is generated,
#'the first column is depth/height, the other columns are the quantile
#'(0.025,0.373,0.5,0.6827,0.975) age values of the runs
#'
#' @author
#'Based on the \link[astrochron]{sedrate2time}
#'function of the 'astrochron' R package
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <\doi{doi:10.1016/j.earscirev.2018.11.015}>
#'
#' @examples
#' \donttest{
#' #Convert a tracked curve to a depth time space. The examples uses the
#' #magnetic susceptibility data set of Pas et al., (2018).
#'
#'# perform the CWT
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
#'mag_track_time<- curve2time_unc(tracked_cycle_curve = mag_track_complete,
#'tracked_cycle_period = 405,
#'tracked_cycle_period_unc = 2.4,
#'tracked_cycle_period_unc_dist = "n",
#'n_simulations = 25,
#'output=1 )
#'
#'}
#'@return
#'If output = 1 a matrix with the predicted ages given the input for each run
#'is given
#'If output = 2 a matrix with 6 columns is generated, the first column is
#'depth/height, the other columns are the quantile
#'(0.025,0.373,0.5,0.6827,0.975) age values of the runs
#'
#' @export
#' @importFrom astrochron sedrate2time
#' @importFrom stats quantile
#' @importFrom stats runif
#' @importFrom stats rnorm

curve2time_unc <- function(tracked_cycle_curve = NULL,
                           tracked_cycle_period = NULL,
                           tracked_cycle_period_unc = NULL,
                           tracked_cycle_period_unc_dist = "n",
                           n_simulations = NULL,
                           output=1) {

  dat <- as.matrix(tracked_cycle_curve[,c(1,2)])
  dat <- na.omit(dat)
  dat <- dat[order(dat[, 1], na.last = NA, decreasing = F), ]
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

  interp_x  <- interp[[1]]

  out <- matrix(data=NA,nrow=length(interp_x),ncol=n_simulations)


  for (i in 1:n_simulations){

  if (tracked_cycle_period_unc_dist == "u") {
      tracked_cycle_period_new <-
        runif(
          1,
          min = tracked_cycle_period - tracked_cycle_period_unc,
          max = tracked_cycle_period + tracked_cycle_period_unc
        )
    }

    if (tracked_cycle_period_unc_dist == "n") {
      tracked_cycle_period_new <-
        rnorm(1, mean = tracked_cycle_period, sd = tracked_cycle_period_unc)
    }



    tracked_cycle_curve_2 <- tracked_cycle_curve
    tracked_cycle_curve_2[, 2] <-  tracked_cycle_curve[, 2] / (tracked_cycle_period_new / 100)

    dat <- as.matrix(tracked_cycle_curve_2[,c(1,2)])
    dat <- na.omit(dat)
    dat <- dat[order(dat[, 1], na.last = NA, decreasing = F), ]
    interp <- approx(dat[, 1], dat[, 2], interp_x, method = "linear",
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
    out[,i] <- hsum

    }

  if (output == 1) {
    res <- out

  }

  if (output == 2) {
    res <- interp_x
    res <- cbind(res, quantile(out,probs=c(0.025,0.373,0.5,0.6827,0.975),na.rm=TRUE))
  }


  return(res)

}
