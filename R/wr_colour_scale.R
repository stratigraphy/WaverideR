#' Colour breaks for image plots (one rule for all plots of the package)
#'
#' Every image plot of WaverideR (scalograms, superlet and synchrosqueezed maps, de-chirped maps, evolutive spectra,
#' evidence maps, rate-grams) takes the same argument \code{scale}, handled here:
#' \describe{
#' \item{\code{"auto"}}{(default) \code{"log"} for a map that is mostly empty (synchrosqueezed, superlet-SST,
#' chirplet-SST, reassigned: more than 30 \% of the cells are zero) or that spans more than three decades between its
#' 5th and 99th percentile (e.g. maximum-entropy spectra, red-noise spectra), \code{"quantile"} otherwise.}
#' \item{\code{"log"}}{logarithmic from the 2nd percentile of the positive cells (at most eight decades below the
#' maximum) to the maximum; empty cells get the lowest colour.}
#' \item{\code{"quantile"}}{histogram-equalised: every colour covers the same number of cells.}
#' \item{\code{"sqrt"}}{linear in amplitude (square root of power).}
#' \item{\code{"linear"}}{linear in the value.}
#' }
#' @param M Matrix of values (power, amplitude, score).
#' @param scale One of \code{"auto"}, \code{"log"}, \code{"quantile"}, \code{"sqrt"}, \code{"linear"}.
#' @param n Number of colours.
#' @param sparse Force the sparse-map rule of \code{"auto"} (\code{TRUE}/\code{FALSE}); \code{NULL} decides from the data.
#' @return A strictly increasing vector of \code{n + 1} breaks in the units of \code{M} (for \code{image(breaks = )}),
#' with attribute \code{scale} (the scale used).
#' @keywords internal
#' @noRd
.wr_colour_breaks <- function(M, scale = "auto", n = 64, sparse = NULL) {
  scale <- match.arg(scale, c("auto", "log", "quantile", "sqrt", "linear"))
  v <- as.numeric(M); v <- v[is.finite(v)]
  if (!length(v)) { br <- seq(0, 1, length.out = n + 1); attr(br, "scale") <- scale; return(br) }
  mn <- min(v); mx <- max(v)
  if (mx <= mn) { br <- seq(mn - 0.5, mx + 0.5, length.out = n + 1); attr(br, "scale") <- scale; return(br) }
  pos <- v[v > 0]
  if (scale == "auto") {
    if (is.null(sparse)) sparse <- mean(v <= 0) > 0.3 && length(pos) > 10
    wide <- length(pos) > 10 && mn >= 0 &&
      stats::quantile(pos, 0.99, names = FALSE) > 1e3 * stats::quantile(pos, 0.05, names = FALSE)
    scale <- if (isTRUE(sparse) || wide) "log" else "quantile"
  }
  if (scale == "log" && (length(pos) < 10 || mn < -1e-12 * mx)) scale <- "quantile"   # signed values: no log
  br <- switch(scale,
    log = {
      lo <- max(1e-8 * mx, stats::quantile(pos, 0.02, names = FALSE))
      c(min(mn, 0), exp(seq(log(lo), log(mx), length.out = n)))
    },
    quantile = stats::quantile(v, seq(0, 1, length.out = n + 1), names = FALSE),
    sqrt = mn + (mx - mn) * seq(0, 1, length.out = n + 1)^2,
    linear = seq(mn, mx, length.out = n + 1))
  # strictly increasing, covering the whole range
  eps <- 1e-9 * (mx - mn)
  br <- cummax(br) + seq_along(br) * eps
  br[1] <- mn - eps; br[length(br)] <- max(br[length(br)], mx + eps)
  attr(br, "scale") <- scale
  br
}

#' Replace the missing cells of a map by its lowest value, so that empty space is drawn in the lowest colour
#' @noRd
.wr_fill_na <- function(M) {
  ok <- is.finite(M)
  if (any(ok) && !all(ok)) M[!ok] <- min(M[ok])
  M
}

#' Draw a matrix with the common colour scale: image() with breaks from .wr_colour_breaks
#' @param x,y Axes; \code{M} has one row per \code{y} and one column per \code{x} (as the scalograms).
#' @noRd
.wr_image <- function(x, y, M, scale = "auto", palette = "viridis", n = 64, fill_na = TRUE, sparse = NULL, ...) {
  if (fill_na) M <- .wr_fill_na(M)
  br <- .wr_colour_breaks(M, scale, n, sparse)
  graphics::image(x, y, t(M), col = grDevices::hcl.colors(length(br) - 1, palette), breaks = br, useRaster = TRUE, ...)
  invisible(br)
}

#' Draw a log10 significance ratio (0 = threshold) with a diverging palette whose middle is the threshold
#' @noRd
.wr_image_signif <- function(x, y, M, ...) {
  M <- .wr_fill_na(M)
  up <- max(0.5, stats::quantile(M[M > 0], 0.995, na.rm = TRUE, names = FALSE), na.rm = TRUE)
  M[M < -1] <- -1; M[M > up] <- up
  br <- c(seq(-1, 0, length.out = 33), seq(0, up, length.out = 33)[-1])
  graphics::image(x, y, t(M), col = grDevices::hcl.colors(64, "Blue-Red 3"), breaks = br, useRaster = TRUE, ...)
  invisible(br)
}
