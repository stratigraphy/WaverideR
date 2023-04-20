#' @title Converts a tracked tracked to a sedimentation rate curve
#'
#' @description Converts the period of a tracked cycle to a sedimentation rate curve by assiging a duration
#' (in kyr) to the period of a tracked cycle
#'
#' @param tracked_cycle_curve  A cycle tracked using the \code{\link{track_period_wavelet}} function \cr
#' Any input (matrix or dataframe) in which the first column is depth in meters and the second column is period in meters
#' @param tracked_cycle_period Period of the tracked cycle (in kyr).
#'
#' @examples
#' \donttest{
#'#Conversion of the period (in meters) of a 405 kyr eccentricity cycle tracked \cr
#'#in a wavelet spectra by assigning a duration of  405 kyr to the tracked cycle.
#'# perform the CWT
#'mag_wt <- analyze_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = TRUE,
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
#'   genplot = TRUE
#' )
#'
#'# smooth the tracking of the 405 kyr eccentricity cycle
#' mag_track_complete <- loess_auto(time_series = mag_track_complete,
#' genplot = TRUE, print_span = TRUE)
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


curve2sedrate <-  function(
    tracked_cycle_curve= NULL,
    tracked_cycle_period = NULL
){
  tracked_cycle_curve[,2] <- tracked_cycle_curve[,2]/(tracked_cycle_period/100)}




