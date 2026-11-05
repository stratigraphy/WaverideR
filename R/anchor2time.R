#'@title Convert a proxy record to the time domain using anchor points
#'
#' @description
#' Convert a proxy record to the time domain using anchor points made using the \code{\link{astro_anchor}} function.
#'
#'
#'@param anchor_points Anchor points: the output of \code{\link{astro_anchor}}, or a matrix or data frame with
#'depth in column 1 and time in column 2 (at least two rows, sorted by increasing depth). Time is linearly
#'interpolated between the anchor points; above the first and below the last anchor point it is extrapolated
#'with the sedimentation rate of the first/last anchored interval. Required (\code{Default=NULL} gives an error).
#'@param data Proxy record to convert from the depth to the time domain: a two-column matrix or data frame with
#'depth in column 1 (same unit as column 1 of \code{anchor_points}, sorted by increasing depth) and the proxy
#'value in column 2. Required (\code{Default=NULL} gives an error).
#'@param genplot Plot the result (logical, TRUE or FALSE). If \code{TRUE}, 3 plots are stacked on top of each
#'other: 1) the original record in the depth domain, 2) the depth-time plot of the (extended) anchor points,
#'3) the record in the time domain. \code{Default=FALSE}.
#'@param keep_editable Keep the plot layout and graphical parameters after plotting so that extra features can
#'be added to the plot (logical, TRUE or FALSE); if \code{FALSE} the previous \code{par()} settings are restored.
#'Only used when \code{genplot=TRUE}. \code{Default=FALSE}.
#'
#'@examples
#'# Use the age_model_zeeden example anchor points of Zeeden et al., (2013)
#'#to anchor the grey data set of Zeeden et al., (2013) in the time domain.
#'
#'grey_time <- anchor2time(anchor_points=age_model_zeeden,
#' data=grey,
#' genplot=FALSE,
#' keep_editable=FALSE)
#'
#' @return
#'The output is a matrix with 2 columns.
#'The first column is time.
#'The second column sedimentation proxy value.
#'
#'If \code{genplot=TRUE} then 3 plots stacked on top of each other will be plotted.
#'Plot 1: the original data set.
#'Plot 2: the depth time plot.
#'Plot 3: the data set in the time domain.
#'
#' @export
#' @importFrom stats approx
anchor2time <- function(anchor_points = NULL,
    data = NULL,
    genplot = FALSE,
    keep_editable = FALSE) {

  my.data <-  cbind(data)
  completed_series <- cbind(anchor_points[, 1], anchor_points[, 2])
  out <- completed_series

  if (my.data[1,1] < completed_series[1,1] ){
    m_begin <- completed_series[1,1]-my.data[1,1]
    sedrate_begin <- (completed_series[2,2]-completed_series[1,2])/(completed_series[2,1]-completed_series[1,1])
    topvals <- c(my.data[1,1],completed_series[1,2]-(sedrate_begin*m_begin))
  }


  if (my.data[nrow(my.data),1] > completed_series[nrow(completed_series),1] ){
    m_end <- my.data[nrow(my.data),1]-completed_series[nrow(completed_series),1]
    sedrate_end <- (completed_series[nrow(completed_series),2]-completed_series[nrow(completed_series)-1,2])/(completed_series[nrow(completed_series),1]-completed_series[nrow(completed_series)-1,1])
    bottomvals <- c(my.data[nrow(my.data),1],completed_series[nrow(completed_series),2]+(sedrate_end*m_end))
  }


  if ((my.data[1,1] < completed_series[1,1] )&(my.data[nrow(my.data),1] > completed_series[nrow(completed_series),1]) ){
    out <- rbind(topvals,out,bottomvals)
  }


  if ((my.data[1,1] > completed_series[1,1] )&(my.data[nrow(my.data),1] > completed_series[nrow(completed_series),1]) ){
    out <- rbind(out,bottomvals)
  }

  if ((my.data[1,1] < completed_series[1,1] )&(my.data[nrow(my.data),1] < completed_series[nrow(completed_series),1]) ){
    out <- rbind(topvals,out)
  }


  app <- approx(
      x = out[, 1],
      y = out[, 2],
      xout = my.data[, 1])


  completed_series <- as.data.frame(cbind(app$y, my.data[, 2]))


  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    layout.matrix <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
    graphics::layout(
        mat = layout.matrix,
        heights = c(1, 1, 1),
        widths = c(1, 1, 1)
    )
    par(mar = c(4, 2, 1, 1))
    plot(
        x = data[, 1],
        y = data[, 2],
        type = "l",
        main = "Data depth domain",
        xlab = "meters",
        ylab = "prxoy"
    )

    plot(
        x = out[, 1],
        y = out[, 2],
        type = "l",
        xlab = "meters",
        ylab = "Time (ka)",
        main = "Depth-time plot"
    )
    points(x = out[, 1], y = out[, 2], cex = 1)

    plot(
        completed_series[, 1],
        completed_series[, 2],
        type = "l",
        xlab = "Time (ka)",
        ylab = "proxy",
        main = "Data time domain"
    )
  }
  return(completed_series)

}
