#' @title Extract signal from a wavelet spectra using a traced period curve
#'
#' @description Extract signal power from the wavelet in the depth domain using the traced period.
#'
#' @param tracked_cycle_curve Tracked period of the cycle along the record (matrix or data frame: depth/time in
#' column 1, period in column 2, in the depth/time unit of the wavelet), e.g. the output of
#' \code{\link{track_period_wavelet}} completed with \code{\link{completed_series}}, optionally smoothed with
#' \code{\link{loess_auto}}. Rows with \code{NA} are dropped and the curve is linearly interpolated (and
#' extrapolated) to the depth/time axis of the wavelet. Required (\code{Default=NULL}).
#' @param wavelet Wavelet object (list) created with \code{\link{analyse_wavelet}}; the signal is reconstructed from
#' its coefficients. Required (\code{Default=NULL}).
#' @param period_up Upper boundary of the band that is reconstructed, as a factor of the (rescaled) tracked period
#' (numeric, > \code{period_down}). \code{Default=1.2}.
#' @param period_down Lower boundary of the band that is reconstructed, as a factor of the (rescaled) tracked period
#' (numeric, > 0). \code{Default=0.8}.
#' @param add_mean Add the mean of the original data series to the extracted cycle (logical, TRUE or FALSE);
#' \code{FALSE} returns it centred on zero. \code{Default=TRUE}.
#' @param tracked_cycle_period Duration in time of the tracked cycle, e.g. 405 (kyr) for the long eccentricity cycle
#' (numeric, > 0), in the same unit as \code{extract_cycle}. The tracked curve is multiplied by
#' \code{extract_cycle / tracked_cycle_period}; also used for the tuning when \code{tune=TRUE}. Required
#' (\code{Default=NULL}).
#' @param extract_cycle Duration of the cycle to extract (numeric, > 0), in the same unit as
#' \code{tracked_cycle_period}; use the same value to extract the tracked cycle itself. Required
#' (\code{Default=NULL}).
#' @param tune Convert the extracted cycle from the depth to the time domain with \code{\link{curve2tune}}, using
#' \code{tracked_cycle_curve} and \code{tracked_cycle_period} (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param plot_residual Plot the residual (extracted cycle minus data) against depth and its histogram (logical,
#' TRUE or FALSE). \code{Default=FALSE}.
#'
#' @author
#' Code based on the reconstruct function of the 'WaveletComp' R package
#' which is based on the wavelet 'MATLAB' code written by Christopher Torrence and Gibert P. Compo (1998).
#'
#' @references
#'Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
#'Wavelet Analysis. R package version 1.1.
#'\url{https://CRAN.R-project.org/package=WaveletComp}
#'
#'Gouhier TC, Grinsted A, Simko V (2021). R package biwavelet: Conduct Univariate and Bivariate Wavelet Analyses. (Version 0.20.21),
#'\url{https://github.com/tgouhier/biwavelet}
#'
#'Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#'Bulletin of the American Meteorological Society 79:61-78.
#'\url{https://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf}
#'
#'@examples
#'\donttest{
#'#Extract the 405 kyr eccentricity cycle from the the magnetic susceptibility \cr
#'#record of the Sullivan core and use the Gabor uncertainty principle to define \cr
#'#the mathematical uncertainty of the analysis and use a factor of that standard \cr
#'#deviation to define boundaries.
#'
#'#Perform the CWT
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
#'# extract the 405 kyr eccentricity cycle from the wavelet spectrum and use the \cr
#'# tracked cycle curve and set factors of the extracted cycle as boundaries
#'
#'mag_405_ecc  <- extract_signal(
#'tracked_cycle_curve = mag_track_complete,
#'wavelet = mag_wt,
#'period_up = 1.2,
#'period_down = 0.8,
#'add_mean = TRUE,
#'tracked_cycle_period = 405,
#'extract_cycle = 405,
#'tune = FALSE,
#'plot_residual = FALSE
#')
#'}
#'@return
#'Returns a matrix with 2 columns
#'The first column is depth/time
#'The second column is extracted signal
#'
#'
#' @export
#' @importFrom stats na.omit


extract_signal <- function(tracked_cycle_curve = NULL,
    wavelet=NULL,
    period_up =1.2,
    period_down = 0.8,
    add_mean=TRUE,
    tracked_cycle_period=NULL,
    extract_cycle=NULL,
    tune=FALSE,
    plot_residual=FALSE){




  my.w <- wavelet
  my.data <- cbind(wavelet$x,wavelet$y)
  filtered_cycle <- my.data[,1]
  filtered_cycle <- as.data.frame(filtered_cycle)
  filtered_cycle$value <- NA

  completed_series <- na.omit(tracked_cycle_curve)
  completed_series[,2] <- completed_series[,2]*(extract_cycle/tracked_cycle_period)
  app <- wr_approx_extrap(x=completed_series[,1],y=completed_series[,2],xout=my.data[,1],
      method="linear")
  interpolated <- cbind(app$x,app$y)

  Wave = my.w$Wave
  Power = my.w$Power

  nc = my.w$nc
  nr = my.w$nr
  dt = my.w$dt
  dj = my.w$dj

  Scale = my.w$Scale
  Period = my.w$Period
  loess.span = my.w$loess.span
  rec.waves = matrix(0, nrow = nr, ncol = nc)


  for (s.ind in seq_len(nr)) {
    rec.waves[s.ind, ] = (Re(Wave[s.ind, ])/sqrt(Scale[s.ind])) *
        dj * sqrt(dt)/.wr_psi0(my.w)
  }


  interpolated <- as.data.frame(interpolated)
  interpolated$high <- interpolated[,2]*(period_up)
  interpolated$low <- interpolated[,2]*(period_down)

  for (i in 1:nrow(filtered_cycle)){
    row_nr_high <- wr_closest(Period[],interpolated[i,3],which=TRUE)
    row_nr_low <- wr_closest(Period[],interpolated[i,4],which=TRUE)

    row_nr_high <-row_nr_high[1]
    row_nr_low <- row_nr_low[1]

    value <- rec.waves[c(row_nr_low:row_nr_high),i]
    value  <- sum(value, na.rm = T)
    value <- as.numeric(value)
    filtered_cycle[i,2] <- value
  }

  rec_value  <- colSums(rec.waves, na.rm = T)

  filtered_cycle[,2] <- (filtered_cycle[,2]) * sd(my.data[,2])/sd(rec_value)

  if(add_mean==TRUE){
    filtered_cycle[,2] <- filtered_cycle[,2] + mean(my.data[,2])
  }

  if(plot_residual==TRUE){
    residual <- filtered_cycle[,2]-my.data[,2]
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    layout.matrix <- matrix(c(1,2), nrow = 2, ncol = 1)
    graphics::layout(mat = layout.matrix,
        heights = c(1, 1), # Heights of the two rows
        widths = c(1,1))
    par(mar = c(4, 4, 1, 1))
    plot(x=filtered_cycle[,1],y=residual,xlab="depth m")
    hist(residual)
  }



  if(tune==TRUE){
    filtered_cycle <- curve2tune(
        data = filtered_cycle,
        tracked_cycle_curve = tracked_cycle_curve,
        tracked_cycle_period = tracked_cycle_period,
        genplot = FALSE)

  }


  filtered_cycle
}
