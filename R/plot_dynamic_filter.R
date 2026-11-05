#' @title Plot the result of a dynamic filter
#'
#' @description The \code{\link{plot_dynamic_filter}} function plots the
#' output of \code{\link{dynamic_filter}}: the (demeaned) data with the
#' filtered signal, and the depth-varying pass band (window lower, upper and
#' centre frequencies or periods, with the control points), reproducing the
#' overlay Acycle draws on its evolutionary spectrum
#' (dynamic_filter_lang.m). Written for the R translation of Acycle v3.1
#' dynamic filtering (Mingsong Li).
#'
#' @param data Record that was filtered, the same as given to \code{\link{dynamic_filter}}: a two-column matrix or
#' data frame with depth/time in column 1 and the proxy value in column 2. It is drawn demeaned in the top panel.
#' Required (\code{Default=NULL} gives an error).
#' @param dyn_filter Dynamic filter result to plot: the list returned by \code{\link{dynamic_filter}} (its
#' \code{filtered}, \code{windows}, \code{lower_controls} and \code{upper_controls} elements are used). Required
#' (\code{Default=NULL} gives an error).
#' @param period Unit of the pass band in the bottom panel (logical): \code{TRUE} draws the band limits, centre and
#' control points as periods (1/frequency, depth/time units), \code{FALSE} as frequencies (cycles per depth/time
#' unit). \code{Default=FALSE}.
#' @param x_lab Label of the depth/time axis of both panels (character). \code{Default="depth"}.
#' @param y_lab Label of the y-axis of the data panel (character). \code{Default="proxy"}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Nicolas
#' Thibault, Giovanni Rizzi and Mingsong Li (Acycle).
#'
#' @references
#' Li, M., Hinnov, L., Kump, L., 2019. Acycle: Time-series analysis software
#' for paleoclimate research and education. Computers & Geosciences 127,
#' 12-22. <doi:10.1016/j.cageo.2019.02.011>
#'
#' @examples
#' \donttest{
#' depth <- seq(0, 300, by = 0.25)
#' period <- 5 + 10 * depth / 300
#' proxy <- sin(2 * pi * cumsum(1 / period) * 0.25)
#' ctrl <- dynamic_filter_controls(tracked = cbind(depth, period),
#'   rel_bandwidth = 0.3, n_points = 30)
#' res <- dynamic_filter(data = cbind(depth, proxy), controls = ctrl,
#'   window_length = 40, step_length = 1)
#' plot_dynamic_filter(data = cbind(depth, proxy), dyn_filter = res)
#' }
#'
#' @return No return value; a plot is produced.
#'
#' @export
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom graphics legend

plot_dynamic_filter <- function(data = NULL,
    dyn_filter = NULL,
    period = FALSE,
    x_lab = "depth",
    y_lab = "proxy"){
  if (is.null(data) || is.null(dyn_filter)){
    stop("data and dyn_filter are required")
  }
  data <- as.matrix(data)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
  y <- data[, 2] - mean(data[, 2])
  plot(data[, 1], y, type = "l", col = "grey50", xlab = x_lab, ylab = y_lab)
  lines(dyn_filter$filtered[, 1], dyn_filter$filtered[, 2], col = "red",
      lwd = 1.5)
  legend("topright", legend = c("data (demeaned)", "dynamic filter"),
      col = c("grey50", "red"), lty = 1, bty = "n", cex = 0.8)

  win <- dyn_filter$windows
  lo <- dyn_filter$lower_controls
  up <- dyn_filter$upper_controls
  if (period == TRUE){
    band <- cbind(1 / win[, "f_high"], 1 / win[, "f_low"],
        1 / win[, "f_centre"])
    lo_y <- 1 / lo[, 2]
    up_y <- 1 / up[, 2]
    band_lab <- "period"
  } else {
    band <- cbind(win[, "f_low"], win[, "f_high"], win[, "f_centre"])
    lo_y <- lo[, 2]
    up_y <- up[, 2]
    band_lab <- "frequency"
  }
  plot(win[, "centre"], band[, 3], type = "l", col = "red",
      xlim = range(data[, 1]), ylim = range(c(band, lo_y, up_y)),
      xlab = x_lab, ylab = band_lab)
  lines(win[, "centre"], band[, 1], col = "darkgreen")
  lines(win[, "centre"], band[, 2], col = "darkgreen")
  points(lo[, 1], lo_y, pch = 1, cex = 0.7)
  points(up[, 1], up_y, pch = 1, cex = 0.7)
  legend("top", legend = c("pass band", "centre", "control points"), horiz = TRUE,
      col = c("darkgreen", "red", "black"), lty = c(1, 1, NA),
      pch = c(NA, NA, 1), bty = "n", cex = 0.8)
  invisible(NULL)
}
