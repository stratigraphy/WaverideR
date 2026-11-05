#' @title Plot a maximum entropy (MESA) spectrum
#'
#' @description Plots the power spectrum of \code{\link{mesa_spectrum}} on a logarithmic power
#' axis, against period or frequency, with the strongest peaks labelled.
#'
#' @param x Maximum entropy spectrum to plot: the list returned by \code{\link{mesa_spectrum}}. Required, no default.
#' @param xaxis Quantity on the x axis (character): \code{"period"} (in the depth/time unit of the input record) or
#' \code{"frequency"} (cycles per unit). \code{Default="period"} (the first element of
#' \code{c("period", "frequency")}).
#' @param peaks Mark the spectral peaks of \code{x$peaks} that fall inside \code{xlim} with red triangles and label
#' the strongest ones with their period (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param n_label Number of strongest peaks (inside \code{xlim}) that are labelled with their period (integer, >= 0).
#' \code{Default=5}.
#' @param xlim Range of the x axis, in the unit of \code{xaxis} (numeric vector of length 2); \code{NULL} uses the
#' full range of positive frequencies (or the corresponding periods). \code{Default=NULL}.
#' @param log_x Draw the x axis on a logarithmic scale (logical, TRUE or FALSE); the power axis is always
#' logarithmic. \code{Default=FALSE}.
#' @param col Colour of the spectrum line (a colour name, hexadecimal string or number). \code{Default="black"}.
#' @param main Title of the plot (character); \code{NULL} gives a title with the AR order and the order-selection
#' criterion. \code{Default=NULL}.
#' @param ... Further graphical arguments passed to \code{\link[graphics]{plot}} (e.g. \code{lwd}, \code{ylim}).
#'
#' @return Called for its plot; invisibly returns \code{x}.
#' @examples
#' set.seed(1)
#' t <- seq(0, 600, by = 2)
#' y <- sin(2 * pi * t / 100) + rnorm(length(t), sd = 0.5)
#' plot_mesa_spectrum(mesa_spectrum(data.frame(t, y)))
#' @importFrom graphics plot lines points text
#' @export
plot_mesa_spectrum <- function(x,
    xaxis = c("period", "frequency"),
    peaks = TRUE,
    n_label = 5,
    xlim = NULL,
    log_x = FALSE,
    col = "black",
    main = NULL,
    ...) {
  xaxis <- match.arg(xaxis)
  s <- x$spectrum
  s <- s[s$frequency > 0, ]
  xv <- if (xaxis == "period") { s$period } else { s$frequency }
  xl <- if (xaxis == "period") { "Period" } else { "Frequency" }
  if (is.null(main)) {
    main <- sprintf("MESA (Burg), AR order %d (%s)", x$order, x$criterion)
  }
  if (is.null(xlim)) {
    xlim <- range(xv)
    if (xaxis == "period") {
      xlim <- c(min(xv), max(xv))
    }
  }
  plot(xv, s$power, type = "l", log = paste0(if (log_x) { "x" } else { "" }, "y"), xlim = xlim,
      xlab = xl, ylab = "Power spectral density", main = main, col = col, ...)
  if (peaks && nrow(x$peaks)) {
    pxv <- if (xaxis == "period") { x$peaks$period } else { x$peaks$frequency }
    ok <- pxv >= min(xlim) & pxv <= max(xlim)
    points(pxv[ok], x$peaks$power[ok], pch = 25, bg = "red3", col = "red3", cex = 0.8)
    lab <- utils::head(which(ok), n_label)
    if (length(lab)) {
      text(pxv[lab], x$peaks$power[lab], labels = format(signif(x$peaks$period[lab], 3)),
          pos = 3, cex = 0.7, col = "red3", xpd = NA)
    }
  }
  invisible(x)
}
