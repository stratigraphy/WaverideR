#' @title Perform an automatically loess based smoothing of a time series
#'
#' @description Perform an automatically loess based smoothing of a time series.
#' The local polynomial regression with automatic smoothing parameter selection is based on an
#' optimization using the 'aicc' bias-corrected 'AIC' criterion and the 'gcv' generalized cross-validation criterion.
#'
#'
#'@param time_series Input is a time series with the first column being depth or time and the second column being a proxy
#'@param genplot Option to generate plot \code{Default=TRUE}. \cr
#'The plot will consist of the original signal in blue, the smoothed plot is displayed
#'in black and the + and - 1 sd bounds of the smoothing are displayed in red.
#'@param print_span Print span length as a fraction of the total length of the record.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#'@author
#'Based on the the loess.as function of the 'fANCOVA' R package.
#'
#'@references
#'Cleveland, W. S. (1979) Robust locally weighted regression and smoothing scatter plots. Journal of the American Statistical Association. 74, 829–836. <doi:10.1080/01621459.1979.10481038>
#'Hurvich, C.M., Simonoff, J.S., and Tsai, C.L. (1998), Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion. Journal of the Royal Statistical Society B. 60, 271–293 <doi:10.1111/1467-9868.00125>
#'Golub, G., Heath, M. and Wahba, G. (1979). Generalized cross validation as a method for choosing a good ridge parameter. Technometrics. 21, 215–224. <doi:10.2307/1268518>
#'
#'
#'@examples
#'\donttest{
#'#'smooth the period curve of the 405 kyr eccentricity cycle extracted from
#'# the magnetic susceptibility data set of Pas et al., (2018)
#'#perform the CWT on the magnetic susceptibility data set of Pas et al., (2018)
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
#' mag_track_complete <- completed_series(
#'   wavelet = mag_wt,
#'   tracked_curve = mag_track,
#'   period_up = 1.2,
#'   period_down = 0.8,
#'   extrapolate = TRUE,
#'   genplot = FALSE,
#'   keep_editable=FALSE
#' )
#'
#'#Smooth the completed tracking of the 405 kyr eccentricity cycle as tracked in the wavelet spectra
#'mag_track_complete <- loess_auto(time_series = mag_track_complete,
#'genplot = FALSE, print_span = FALSE,keep_editable=FALSE)
#'}
#'@return
#'A matrix with 3 columns.
#'The first column is depth/time.
#'The second column is the smoothed curve.
#'The third column is difference between the original curve and the smoothed curve.
#'
#' @export
#' @importFrom stats loess
#' @importFrom stats predict
#' @importFrom fANCOVA loess.as


loess_auto <-
  function(time_series = NULL,
           genplot = FALSE,
           print_span = FALSE,
           keep_editable = FALSE) {
    completed_series_test <- na.omit(time_series)
    parameters_to_use <-
      loess.as(
        y = completed_series_test[, 2],
        x = completed_series_test[, 1],
        degree = 2,
        criterion = c("aicc", "gcv")[2],
        family = c("symmetric"),
        plot = FALSE
      )
    span <- parameters_to_use$pars$span
    time_series <- as.data.frame(time_series)

    colnames(time_series) <- c("V1", "V2")
    loessMod <- loess(
      V2 ~ V1,
      data = time_series,
      span = span,
      model = TRUE,
      degree = 1,
      family = c("gaussian", "symmetric")
    )


    smoothed <- as.data.frame(predict(
      loessMod,
      time_series,
      se = FALSE,
      interval = "predict",
      level = 0.68
    ))

    smoothed <- cbind(time_series[, 1], smoothed)
    colnames(smoothed) <- c("depth", "fit")
    smoothed$SD <- (abs(smoothed$fit - time_series[, 2]))
    smoothed <- na.omit(smoothed)


    if (genplot == TRUE) {
      if (keep_editable == FALSE) {
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
      }
      max_y <-  (max(smoothed[, 2] + smoothed[, 3]))
      min_y <-  (min(smoothed[, 2] - smoothed[, 3]))
      plot(
        smoothed[, 1],
        smoothed[, 2],
        type = "l",
        lwd = 2,
        col = "black",
        ylim = c(min_y, max_y)
      )
      lines(completed_series_test[, 1],
            completed_series_test[, 2],
            col = "blue")
      lines(
        x = smoothed[, 1],
        y = (smoothed[, 2] + smoothed[, 3]),
        col = "red",
        lty = 2
      )
      lines(
        x = smoothed[, 1],
        y = (smoothed[, 2] - smoothed[, 3]),
        col = "red",
        lty = 2
      )
    }

    if (print_span == TRUE) {
      cat("span is: ", parameters_to_use$pars$span)
    }
    return(smoothed)
  }
