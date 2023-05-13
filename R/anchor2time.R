#'@title Convert a proxy record to the time domain using anchor points
#'
#' @description
#' Convert a proxy record to the time domain using anchor points made using the \code{\link{astro_anchor}} function.
#'
#'
#'@param anchor_points Anchor points made using the \code{\link{astro_anchor}} function or a matrix in which the first column is depth
#'and the second column is time.
#'@param data Data set which needs to be converted from the depth to time domain using set anchor points.
#'The data set should consist of a matrix with 2 column the first column should be depth
#'and the second column should be a proxy value.
#' @param genplot
#'If  \code{genplot=FALSE} then 3 plots stacked on top of each other will be plotted.
#'Plot 1: the original data set
#'Plot 2: the depth time plot
#'Plot 3: the data set in the time domain
#'set to TRUE to allow for anchoring using the GUI
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#'@examples
#'# Use the age_model_zeeden example anchor points of Zeeden et al., (2013)
#'#to anchor the grey data set of Zeeden et al., (2013) in the time domain.
#'
#'anchored2time <- anchor2time(anchor_points=age_model_zeeden,
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
  my.data <- cbind(data)
  completed_series <- cbind(anchor_points[, 1], anchor_points[, 2])
  yleft_comp <- completed_series[1, 2]
  yright_com <- completed_series[nrow(completed_series), 2]
  out <- completed_series
  app <-
    approx(
      x = out[, 1],
      y = out[, 2],
      xout = my.data[, 1],
      method = "linear",
      yleft = yleft_comp,
      yright = yright_com
    )
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
      ylab = "prxoy",
      main = "Data time domain"
    )
  }
return(completed_series)

}
