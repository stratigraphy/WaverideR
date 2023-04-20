#' @title Plot the average spectral power of a wavelet spectra
#'
#' @description Plot the average spectral power of a wavelet spectra using the results of
#' the \code{\link{analyze_wavelet}} function.
#'
#' @param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#' @param y_lab Label for the y-axis \code{Default="Power"}.
#' @param x_lab Label for the x-axis \code{Default="depth (metres)"}.
#'
#' @examples
#'\donttest{
#'#Example 1. Plot the average spectral power of the wavelet spectra of \cr
#'# the Total Solar Irradiance data set of Steinhilver et al., (2012)
#'TSI_wt <-
#'  analyze_wavelet(
#'    data = TSI,
#'    dj = 1/200,
#'    lowerPeriod = 16,
#'    upperPeriod = 8192,
#'    verbose = TRUE,
#'    omega_nr = 6
#'  )
#'
#' plot_avg_wavelet(wavelet=TSI_wt,
#'                  y_lab= "power",
#'                  x_lab="period (years)")
#'
#'
#'#Example 2. Plot the average spectral power of the wavelet spectra of \cr
#'# the magnetic susceptibility data set of De pas et al., (2018)
#'mag_wt <-
#'analyze_wavelet(
#'data = mag,
#'dj = 1/100,
#'lowerPeriod = 0.1,
#'upperPeriod = 254,
#'verbose = TRUE,
#'omega_nr = 10
#')
#' plot_avg_wavelet(wavelet=mag_wt,
#'                  y_lab= "power",
#'                  x_lab="period (metres)")
#'
#'
#'
#'#Example 3. Plot the average spectral power of the wavelet spectra of \cr
#'#the greyscale data set of Zeeden et al., (2013)
#'grey_wt <-
#'  analyze_wavelet(
#'    data = grey,
#'    dj = 1/200,
#'    lowerPeriod = 0.02,
#'    upperPeriod = 256,
#'    verbose = TRUE,
#'    omega_nr = 8
#'  )
#'
#' plot_avg_wavelet(wavelet=grey_wt,
#'                  y_lab= "power",
#'                  x_lab="period (metres)")
#'
#'}
#'@return
#'The output is a plot of the average spectral power of a wavelet spectra
#' @export


plot_avg_wavelet <- function(wavelet = NULL,
                         y_lab = "Power",
                         x_lab = "period (metres)") {

  plot(wavelet$Period,wavelet$Power.avg, log = "x",type="l",
       ylab = y_lab, xlab=x_lab)}
