#' @title Calculate average spectral power from red noise curves for a given percentile
#'
#' @description The \code{\link{percentile_from_red_noise}} function is
#' used to generate and average spectral power curve based on
#' a set percentile based. To generate the percentile curve the results of
#'  the \code{\link{model_red_noise_wt}} function are used.
#'
#'@param red_noise Red noise curves generated using the \code{\link{model_red_noise_wt}} function.
#'@param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#'@param percentile Percentile value (0-1).
#'
#'
#'@examples
#'\dontrun{
#'#'#generate red noise curves based on the magnetic susceptibility record of \cr
#'#the Sullivan core of De pas et al., (2018)
#'
#'mag_wt <- analyze_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'mag_wt_red_noise <- model_red_noise_wt(data=NULL,
#'n_simulations =1000,
#'verbose=FALSE)
#'
#'prob_curve <- percentile_from_red_noise(
#'red_noise = mag_wt_red_noise,
#'wavelet = mag_wt,
#'percentile = 0.9
#')}
#'
#' @return
#'Returns a matrix with 2 columns.\cr
#'The first column is the period (m). \cr
#'The second column is the spectral power at percentile x based on \cr
#'the red noise modelling runs. \cr
#'
#' @export
#' @importFrom matrixStats rowSds
#' @importFrom stats qnorm


percentile_from_red_noise <- function(red_noise=NULL,
                                       wavelet=NULL,
                                      percentile=NULL) {
noise_period <- cbind(wavelet$Period,red_noise)
noise_period2 <- as.data.frame(noise_period)
noise_period2$mean <- rowMeans(noise_period2[,2:(ncol(red_noise)+1)])
noise_period2$sd <- rowSds(as.matrix(noise_period2[,2:(ncol(red_noise)+1)]))
prob <- qnorm(p=percentile, mean=noise_period2$mean, sd=noise_period2$sd)
prob <- cbind(noise_period[,1],prob)
prob
}



