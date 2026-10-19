#' @title Plot a maximum entropy (MESA) spectrum
#'
#' @description Plots the power spectrum of \code{\link{mesa_spectrum}} on a logarithmic power
#' axis, against period or frequency, with the strongest peaks labelled.
#'
#' @param x output of \code{\link{mesa_spectrum}}.
#' @param xaxis \code{"period"} (default) or \code{"frequency"}.
#' @param peaks If TRUE (default) the peaks of \code{x$peaks} are marked and labelled with their period.
#' @param n_label Number of strongest peaks to label (default 5).
#' @param xlim Optional limits of the x axis (in the units of \code{xaxis}).
#' @param log_x If TRUE the x axis is logarithmic (default FALSE).
#' @param col Line colour.
#' @param main Title. The default describes the order and criterion.
#' @param ... further arguments passed to \code{plot}.
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
  xv <- if (xaxis == "period") s$period else s$frequency
  xl <- if (xaxis == "period") "Period" else "Frequency"
  if (is.null(main)) {
    main <- sprintf("MESA (Burg), AR order %d (%s)", x$order, x$criterion)
  }
  if (is.null(xlim)) {
    xlim <- range(xv)
    if (xaxis == "period") {
      xlim <- c(min(xv), max(xv))
    }
  }
  plot(xv, s$power, type = "l", log = paste0(if (log_x) "x" else "", "y"), xlim = xlim,
       xlab = xl, ylab = "Power spectral density", main = main, col = col, ...)
  if (peaks && nrow(x$peaks)) {
    pxv <- if (xaxis == "period") x$peaks$period else x$peaks$frequency
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
