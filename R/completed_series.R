#' @title Complete the tracking of the period a cycle in a scalogram
#'
#'@description Use the traced series and the existing wavelet or superlet scalogram
#'to complete the tracking of a cycle in scalogram.
#' The selected points using the \code{\link{track_period_wavelet}} or \code{\link{track_period}}
#' function form a incomplete line unless every point is tracked. However
#' clicking every individual point along a wavelet ridge is time intensive and error prone.
#' To avoid errors and save time the \code{\link{completed_series}} function can be used to
#' complete the tracing of a cycle in a wavelet spectra.
#' The \code{\link{completed_series}} function interpolates the
#' data points selected using the \code{\link{track_period_wavelet}}
#' or \code{\link{track_period}}. A a search a algorithm then looks up and replaces
#' the interpolated curve values with the values of the nearest spectral peak in the wavelet spectra.
#'
#'@param scalogram Scalogram in which the cycle was tracked: the output of \code{\link{analyse_wavelet}} (class
#'\code{analyze.wavelet}) or of \code{\link{analyse_superlet}} (class \code{analyze.superlet}). Its depth axis sets
#'the positions of the completed curve. Required (\code{Default=NULL} gives an error).
#'@param tracked_curve Incomplete tracked period curve (matrix or data frame, two columns: depth in column 1 and the
#'tracked period in column 2, in the depth units of the scalogram), e.g. the output of
#'\code{\link{track_period_wavelet}} or \code{\link{track_period}}; rows with \code{NA} are dropped. It is linearly
#'interpolated onto the depth axis of the scalogram before the search for spectral peaks. Required
#'(\code{Default=NULL} gives an error).
#'@param period_up Factor (numeric, > 1) that sets the upper limit of the search window: at every depth the
#'nearest spectral peak (local power maximum along the period axis) is used only if its period is below the
#'interpolated period times \code{period_up}; otherwise no value is assigned there (see \code{extrapolate}).
#'\code{Default=1.2}.
#'@param period_down Factor (numeric, between 0 and 1) that sets the lower limit of the search window: the nearest
#'spectral peak is used only if its period is above the interpolated period times \code{period_down}.
#'\code{Default=0.8}.
#'@param extrapolate Fill the depths where no spectral peak was found within the window by linear interpolation
#'between the found values (with the first and last found values held constant towards the ends) (logical, TRUE
#'or FALSE); \code{FALSE} leaves \code{NA} at those depths. \code{Default=TRUE}.
#'@param genplot Plot the completed curve (red) together with the original tracked curve (black) (logical, TRUE or
#'FALSE). \code{Default=FALSE}.
#'@param keep_editable Keep the graphical parameters of the plot so that extra elements can be added afterwards
#'(logical, TRUE or FALSE); \code{FALSE} restores the previous \code{par} settings. \code{Default=FALSE}.
#'@param ... Arguments of older versions: \code{wavelet} is accepted as an alias of \code{scalogram}; other
#'arguments are ignored.
#'
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
#'  analyse_wavelet(
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
#'  scalogram = grey_wt,
#'  tracked_curve  = grey_track,
#'  period_up  = 1.25,
#'  period_down  = 0.75,
#'  extrapolate = TRUE,
#' genplot = FALSE,
#' keep_editable=FALSE
#')
#'}
#' @return
#' Returns a matrix with 2 columns
#' The first column is the depth axis
#' The second column is the completed tracking of the period a cycle of the wavelet spectra
#'
#' @export
#' @importFrom stats approx

completed_series <- function (scalogram = NULL, tracked_curve = NULL, period_up = 1.2,
    period_down = 0.8, extrapolate = TRUE, genplot = FALSE,
    keep_editable = FALSE,...)
{


  dots <- list(...)

  # Allow alias: wavelet = scalogram
  if (is.null(scalogram) && !is.null(dots$wavelet)) {
    scalogram <- dots$wavelet
  }
  my.w <- scalogram
  my.data <- cbind(scalogram$x, scalogram$y)

  if(inherits(scalogram,"analyze.superlet") == TRUE){
    Pwert <- t(my.w$Power)
  }
  if(inherits(scalogram,"analyze.wavelet") == TRUE ){
    Pwert <- my.w$Power
  }

  maxdetect <- matrix(nrow = (nrow(Pwert)), ncol = ncol(Pwert), 0)
  for (j in 1:ncol(Pwert)) {
    for (i in 2:(nrow(maxdetect) - 1)) {
      if ((Pwert[i, j] - Pwert[(i + 1), j] > 0) &
          (Pwert[i, j] - Pwert[(i - 1), j] > 0)) {
        maxdetect[i, j] <- 1
      }
    }
  }

  out <- tracked_curve
  out <- na.omit(out)
  yleft_out <- out[1, 2]
  yright_out <- out[nrow(out), 2]
  seq <- seq(
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


  if(inherits(scalogram,"analyze.superlet") == TRUE){
    periods <- as.data.frame(rev(my.w$Period))
  }
  if(inherits(scalogram,"analyze.wavelet") == TRUE ){
    periods <- as.data.frame(my.w$Period)
  }
  completed_series <- matrix(data = NA,
      nrow = nrow(app[, ]),
      ncol = 2)
  completed_series[, 1] <- app[, 1]

  for (i in 1:nrow(completed_series)) {
    row_nr <- wr_closest(periods[, 1], app[i, 2], which = TRUE)
    row_nr
    row_nr <- row_nr[1]
    if (maxdetect[row_nr, i] == 1) {
      completed_series[i, 2] <- periods[row_nr, ]
    } else {
      sel_row <- as.data.frame(maxdetect[, i])
      sel_row$period <- periods[, 1]
      sel_row <- sel_row[sel_row[, 1] > 0, ]
      row_nr_sel_row <- wr_closest(sel_row[, 2], app[i, 2], which = TRUE)
      row_nr_closest <- wr_closest(periods[, 1], sel_row[row_nr_sel_row, 2], which = TRUE)
      closest_period <- periods[unlist(row_nr_closest[1]), 1]
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
    seq <- seq(
        from = min(my.data[, 1]),
        to = max(my.data[, 1]),
        by = my.data[, 1][2] - my.data[, 1][1]
    )
    app <- approx(
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
    plot(completed_series,
        type = "l",
        col = "red")
    lines(tracked_curve, col = "black")
  }
  return(completed_series)
}
