#' @title Convert a tracked tracked to a sedimentation rate curve
#'
#' @description Converts the period of a tracked cycle to a sedimentation rate curve by
#' assigning a duration (in kyr) to the period of a tracked cycle
#'
#' @param tracked_cycle_curve Tracked period curve of one astronomical cycle (matrix or data frame, two columns):
#' depth in metres in column 1 and the period of the tracked cycle in metres in column 2, e.g. the output of
#' \code{\link{track_period_wavelet}} (preferably completed with \code{\link{completed_series}}). Required
#' (\code{Default=NULL} gives an error).
#' @param tracked_cycle_period Duration of the tracked cycle in kyr (single numeric value, > 0, e.g. 405 for the
#' long eccentricity cycle); the sedimentation rate in cm/kyr is the period in metres divided by
#' \code{tracked_cycle_period / 100}. Required (\code{Default=NULL} gives an error).
#'
#' @examples
#' \donttest{
#'#Conversion of the period (in meters) of a 405 kyr eccentricity cycle tracked
#'#in a wavelet spectra by assigning a duration of  405 kyr to the tracked cycle.
#'# the example uses the magnetic susceptibility data set of Pas et al., (2018)
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
#'#                                   periodlab = "Period (metres)",
#'#                                   x_lab = "depth (metres)")
#'
#'#Instead of tracking, the tracked solution data set \code{\link{mag_track_solution}} is used \cr
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
#'#convert period in meters to sedrate in cm/kyr
#'mag_track_sedrate <- curve2sedrate(tracked_cycle_curve=mag_track_complete,
#'tracked_cycle_period=405)
#'}
#'@return
#'The output is a matrix with 2 columns
#'The first column is depth
#'The second column sedimentation rate in cm/kyr
#' @export
#' @importFrom stats approx


curve2sedrate <-  function(tracked_cycle_curve = NULL,
    tracked_cycle_period = NULL) {
  tracked_cycle_curve[, 2] <-
      (tracked_cycle_curve[, 2] / (tracked_cycle_period / 100))
  return(tracked_cycle_curve)
}
