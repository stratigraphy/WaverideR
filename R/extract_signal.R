#' @title Extract signal from a wavelet spectra using a traced period curve
#'
#' @description Extract signal power from the wavelet in the depth domain using the traced period.
#'
#'@param tracked_cycle_curve Traced period result from the \code{track_period_wavelet}
#'function completed using the \code{completed_series}.
#'The input can be pre-smoothed using the the \code{loess_auto} function.
#'@param wavelet wavelet object created using the \code{\link{analyze_wavelet}} function.
#'@param period_up Upper period as a factor of the to be extracted cycle \code{Default=1.2}.
#'@param period_down Lower period as a factor of the to be extracted cycle \code{Default=0.8}.
#'@param add_mean Add mean to the extracted cycle \code{Default=TRUE}.
#'@param tracked_cycle_period Period in time of the traced cycle.
#'@param extract_cycle Period of the to be extracted cycle.
#'@param tune Convert record from the depth to the time domain using the traced period \code{Default=FALSE}.
#'@param plot_residual Plot the residual signal after extraction of set cycle \code{Default=FALSE}.
#'
#'@examples
#'\donttest{
#'#Extract the 405 kyr eccentricity cycle from the the magnetic susceptibility \cr
#'#record of the Sullivan core and use the Gabor uncertainty principle to define \cr
#'#the mathematical uncertainty of the analysis and use a factor of that standard \cr
#'#deviation to define boundaries.
#'
#'#Perform the CWT
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
#'# extract the 405 kyr eccentricty cycle from the wavelet spectrum and use the \cr
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
#' @importFrom Hmisc approxExtrap
#' @importFrom stats na.omit
#' @importFrom DescTools Closest


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
  app <- approxExtrap(x=completed_series[,1],y=completed_series[,2],xout=my.data[,1],
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
      dj * sqrt(dt)/(pi^(-1/4))}


  interpolated <- as.data.frame(interpolated)
  interpolated$high <- interpolated[,2]*(period_up)
  interpolated$low <- interpolated[,2]*(period_down)

  for (i in 1:nrow(filtered_cycle)){
    row_nr_high <- Closest(Period[],interpolated[i,3],which=TRUE)
    row_nr_low <- Closest(Period[],interpolated[i,4],which=TRUE)

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
