#' @title Plot DYNOT results
#'
#' @description The \code{\link{plot_dynot}} function plots the input data
#' and the DYNOT noise model (median and confidence envelopes, reversed
#' y-axis) as in the Acycle DYNOT window.
#' Translated from the plotting part of Acycle v3.1 DYNOS.m (Mingsong Li).
#'
#' @param dynot_res DYNOT result to plot: the list returned by \code{\link{dynot}} (its \code{data},
#' \code{percentile}, \code{interval_count} and \code{median_index} elements are used). Required
#' (\code{Default=NULL} gives an error).
#' @param plot_data Also plot the input data in an upper panel above the noise model (logical, TRUE or FALSE); with
#' \code{FALSE} only the noise panel (median and up to five confidence envelopes, reversed y-axis) is drawn.
#' \code{Default=TRUE}.
#' @param x_lab Label of the time axis of both panels (character). \code{Default="time (kyr)"}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' (Acycle).
#'
#' @references
#' Li, M., Hinnov, L.A., Huang, C., Ogg, J.G., 2018. Sedimentary noise and
#' sea levels linked to land-ocean water exchange and obliquity forcing.
#' Nature Communications 9, 1004. \doi{10.1038/s41467-018-03454-y}
#'
#' @examples
#' \donttest{
#' set.seed(2)
#' t <- cumsum(runif(900, 1, 3))
#' x <- sin(2 * pi * t / 405) + sin(2 * pi * t / 41) + rnorm(900)
#' res <- dynot(data = cbind(t, x), n_sim = 20, n_out = 200, seed = 1)
#' plot_dynot(res)
#' }
#' @return No return value; a plot is drawn.
#'
#' @export
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics polygon
#' @importFrom graphics lines
#' @importFrom grDevices rgb

plot_dynot <- function(dynot_res = NULL,
    plot_data = TRUE,
    x_lab = "time (kyr)"){
  # Internal helper functions (local to this function)

  # @keywords internal
  # @noRd
  dynot_band_plot <- function(percentile, n_interval, median_index, x_range,
      x_lab = "depth/time", y_lab = "", reverse = FALSE){
    colorcode <- grDevices::rgb(c(221, 201, 176, 126, 67),
        c(234, 227, 219, 201, 180),
        c(224, 209, 188, 146, 100), maxColorValue = 255)
    x <- percentile[, 1]
    pm <- percentile[, -1, drop = FALSE]
    np <- ncol(pm)
    y_lim <- range(pm, na.rm = TRUE)
    if (reverse == TRUE){
      y_lim <- rev(y_lim)
    }
    graphics::plot(x, pm[, median_index], type = "n", xlim = x_range,
        ylim = y_lim, xlab = x_lab, ylab = y_lab)
    for (i in seq_len(min(n_interval, 5))){
      graphics::polygon(c(x, rev(x)), c(pm[, np + 1 - i], rev(pm[, i])),
          col = colorcode[i], border = NA)
    }
    graphics::lines(x, pm[, median_index], col = grDevices::rgb(0, 120, 0,
        maxColorValue = 255), lwd = 1.5, lty = 2)
    invisible(NULL)
  }

  # Main code
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  if (plot_data == TRUE){
    graphics::par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
    graphics::plot(dynot_res$data[, 1], dynot_res$data[, 2], type = "l",
        col = grDevices::rgb(0.2, 0.6, 1), xlab = x_lab,
        ylab = "value")
  }
  dynot_band_plot(dynot_res$percentile, dynot_res$interval_count,
      dynot_res$median_index, range(dynot_res$data[, 1]),
      x_lab = x_lab, y_lab = "Noise", reverse = TRUE)
  invisible(NULL)
}
