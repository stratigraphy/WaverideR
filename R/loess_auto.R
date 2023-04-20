#' @title Perform an automatically loess based smoothing of a timeseries
#'
#' @description Perform an automatically loess based smoothing of a timeseries.
#' The local polynomial regression with automatic smoothing parameter selection is based on an
#' optimization using the aicc  bias-corrected AIC criterion and the gcv generalized cross-validation criterion.
#'
#'
#'@param time_series Input is a time series with the first column being depth or time and the second column being a proxy
#'@param genplot Option to generate plot \code{Default=TRUE}. \cr
#'The plot will consist of the original signal in blue, the smoothed plot is displayed
#'in black and the + and - 1 sd bounds of the smoothing are displayed in red.
#'@param print_span Print span length as a fraction of the total length of the record.
#'
#'@author
#'Based on the the \code{\link{loess.as}} function of the fANCOVA package.
#'
#'@references
#'Cleveland, W. S. (1979) Robust locally weighted regression and smoothing scatterplots. Journal of the American Statistical Association. 74, 829–836.
#'Hurvich, C.M., Simonoff, J.S., and Tsai, C.L. (1998), Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion. Journal of the Royal Statistical Society B. 60, 271–293
#'Golub, G., Heath, M. and Wahba, G. (1979). Generalized cross validation as a method for choosing a good ridge parameter. Technometrics. 21, 215–224.
#'
#'
#'@examples
#'\donttest{
#'#'smooth the period curve of the 405 kyr eccentricity cycle extracted from \cr
#'# the magnetic susceptibility data set of De pas et al., (2018) \cr
#'#perform the CWT on the magnetic susceptibility data set of De pas et al., (2018)
#'
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
#'#Smooth the completed tracking of the 405 kyr eccentricity cycle as tracked in the wavelet spectra
#'mag_track_complete <- loess_auto(time_series = mag_track_complete,
#'genplot = TRUE, print_span = TRUE)
#'}
#'@return
#'A matrix with 3 columns.
#'The first column is depth/time.
#'The second column is the smoothed curve.
#'The third column is difference between the original curve and the smoothed curve.
#'
#' @export
#' @importFrom fANCOVA loess.as
#' @importFrom stats loess
#' @importFrom stats predict

loess_auto <- function(time_series=NULL,genplot=TRUE,print_span=TRUE){
  completed_series_test <- na.omit(time_series)
  parameters_to_use <- loess.as(y=completed_series_test[,2],x=completed_series_test[,1],
                                degree = 2, criterion = c("aicc", "gcv")[2],
                                family = c( "symmetric"), plot = FALSE)
  span <- parameters_to_use$pars$span
  time_series <- as.data.frame(time_series)

  colnames(time_series) <- c("V1","V2")
  loessMod <- loess(V2 ~ V1, data=time_series, span=span,model=TRUE,
                    degree = 1,family = c("gaussian", "symmetric"))


  smoothed <- as.data.frame(predict(loessMod,time_series,se = FALSE,
                                    interval="predict",level=0.68))

  smoothed <- cbind(time_series[,1],smoothed)
  colnames(smoothed) <- c("depth","fit")
  smoothed$SD <- (abs(smoothed$fit-time_series[,2]))
  smoothed <- na.omit(smoothed)


  if (genplot==TRUE){
    max_y <-  (max(smoothed[,2]+smoothed[,3]))
    min_y <-  (min(smoothed[,2]-smoothed[,3]))
    plot(smoothed[,1],smoothed[,2],type="l",lwd=2,col="black",ylim=c(min_y,max_y))
    lines(completed_series_test[,1],completed_series_test[,2],col="blue")
    lines(x=smoothed[,1],y=(smoothed[,2]+smoothed[,3]),col="red",lty=2)
    lines(x=smoothed[,1],y=(smoothed[,2]-smoothed[,3]),col="red",lty=2)
  }

  if (print_span==TRUE){
    print(paste0("span is: ", parameters_to_use$pars$span))
  }
  smoothed
}


