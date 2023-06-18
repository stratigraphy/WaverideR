#' @title Complete the tracking of cycle in a wavelet spectra
#'
#'@description Use the traced series and the existing wavelet
#' spectra to complete the tracking of a cycle of the wavelet spectra.
#' The selected points using the \code{\link{track_period_wavelet}} function form a incomplete line
#' unless every point is tracked. However clicking every individual point along a wavelet ridge is time
#' intensive and error prone.
#' To avoid errors and save time the \code{\link{completed_series}} function can be used to
#' complete the tracing of a cycle in a wavelet spectra.The \code{\link{completed_series}} function interpolates the
#' data points selected using the \code{\link{track_period_wavelet}}. A a search a algorithm then looks up and replaces
#' the interpolated curve values with the values of the nearest spectral peak in the wavelet spectra.
#'
#'@param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#'@param tracked_curve Traced period result from the \code{track_period_wavelet} function.
#'@param period_up The period_up parameter is the factor with which the linear interpolated tracked_curve
#'curve is multiplied by. This linear interpolated tracked_curve multiplied by the period_up factor is
#'the upper boundary which is used  for detecting the spectral peak nearest to the linear interpolated tracked_curve
#'curve. If no spectral peak is detected within the specified boundary the interpolated
#'value is used instead. between spectral peaks \code{Default=1.2},
#'@param period_down  The period_down parameter is the factor with which the linear interpolated tracked_curve
#'curve is multiplied by. This linear interpolated tracked_curve multiplied by the period_down factor is
#'the lower boundary which is used  for detecting the spectral peak nearest to the linear interpolated tracked_curve
#'curve. If no spectral peak is detected within the specified boundary the interpolated
#'value is used instead. between spectral peaks \code{Default=0.8},
#'@param extrapolate Extrapolate the completed curve when through parts where no spectral peaks could be traced
#' \code{Default=TRUE}.
#'@param genplot Generate a plot \code{Default=TRUE}. The red curve is the completed curve,
#'the black curve is the original curve.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#'@examples
#'\donttest{
#'#Use the grey_track example points to complete the tracking of the
#'# precession cycle in the wavelet spectra of the grey scale data set
#'# of Zeeden et al., (2013).
#'
#'
#'
#'grey_wt <-
#'  analyze_wavelet(
#'    data = grey,
#'    dj = 1/200,
#'    lowerPeriod = 0.02,
#'    upperPeriod = 256,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'#The ~22kyr precession cycle is between 0.25 and 1m The grey_track data
#'#set is a pre-loaded uncompleted tracking of the precession cycle
#'
#'#grey_track <- track_period_wavelet(
#'#astro_cycle = 22,
#'#wavelet = NULL,
#'#n.levels = 100,
#'#periodlab = "Period (meters)",
#'#x_lab = "depth (meters)"
#'#)
#'
#'
#'
#'grey_track <- completed_series(
#'  wavelet = grey_wt,
#'  tracked_curve  = grey_track,
#'  period_up  = 1.25,
#'  period_down  = 0.75,
#'  extrapolate = TRUE,
#' genplot = FALSE,
#' keep_editable=FALSE
#')
#'
#'}
#' @return
#' Returns a matrix with 2 columns
#' The first column is the depth axis
#' The second column is the completed tracking of the period a cycle of the wavelet spectra
#'
#' @export
#' @importFrom stats approx


completed_series <- function(wavelet = NULL,
                             tracked_curve = NULL,
                             period_up = 1.2,
                             period_down = 0.8,
                             extrapolate = TRUE,
                             genplot = FALSE,
                             keep_editable = FALSE) {
  my.w <- wavelet
  my.data <- cbind(wavelet$x, wavelet$y)
  Pwert <- my.w$Power
  maxdetect <- matrix(nrow = (nrow(Pwert)), ncol = ncol(Pwert), 0)

  for (j in 1:ncol(Pwert)) {
    for (i in 2:(nrow(maxdetect) - 1)) {
      if ((Pwert[i, j] - Pwert[(i + 1), j] > 0) &
          (Pwert[i, j] - Pwert[(i - 1), j]  > 0))
      {
        maxdetect[i, j] <- 1
      }
    }
  }

  out <-  tracked_curve
  out <- na.omit(out)
  yleft_out <- out[1, 2]
  yright_out <- out[nrow(out), 2]
  seq <-
    seq(
      from = min(my.data[, 1]),
      to = max(my.data[, 1]),
      by = my.data[, 1][2] - my.data[, 1][1]
    )
  app <- approx(
    x = out[, 1],
    y = out[, 2],
    xout = seq,
    method = "linear",
    yleft = yleft_out,
    yright = yright_out
  )
  app <- as.data.frame(cbind(app$x, app$y))
  periods <- as.data.frame(my.w$Period)

  completed_series <- matrix(data = NA,
                             nrow = nrow(app[,]),
                             ncol = 2)
  completed_series[, 1] <- app[, 1]

  for (i  in 1:nrow(completed_series)) {
    row_nr <- Closest(periods[, 1], app[i, 2], which = TRUE)
    row_nr <- row_nr[1]
    if (maxdetect[row_nr, i] == 1) {
      completed_series[i, 2] <- periods[row_nr,]
    } else {
      sel_row <- as.data.frame(maxdetect[, i])
      sel_row$period <- periods[, 1]
      sel_row <- sel_row[sel_row[, 1] > 0, ]
      row_nr_sel_row <-
        Closest(sel_row[, 2], app[i, 2], which = TRUE)
      row_nr_closest <-
        Closest(periods[, 1], sel_row[row_nr_sel_row, 2], which = TRUE)
      closest_period <- periods[row_nr_closest, 1]
      if ((closest_period < (app[i, 2] * period_up)) &
          (closest_period > (app[i, 2] * period_down))) {
        completed_series[i, 2] <- closest_period
      }
    }
  }

  if (extrapolate == TRUE) {
    completed_series <- na.omit(completed_series)
    yleft_comp <- completed_series[1, 2]
    yright_com <- completed_series[nrow(completed_series), 2]
    seq <-
      seq(
        from = min(my.data[, 1]),
        to = max(my.data[, 1]),
        by = my.data[, 1][2] - my.data[, 1][1]
      )
    app <-
      approx(
        x = completed_series[, 1],
        y = completed_series[, 2],
        xout = seq,
        method = "linear",
        yleft = yleft_comp,
        yright = yright_com
      )
    completed_series <- as.data.frame(cbind(app$x, app$y))
  }


  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    plot(completed_series, type = "l", col = "red")
    lines(tracked_curve, col = "black")
  }

  return(completed_series)
}
