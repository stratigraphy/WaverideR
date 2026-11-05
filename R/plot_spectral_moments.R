#' @title Plot the result of spectral_moments
#'
#' @description Plots the mean frequency and bandwidth (and, if available,
#' the bandwidth trend and the sedimentation-rate estimate) returned by
#' \code{\link{spectral_moments}}, as in the Acycle spectral moments GUI.
#' Translated from the plotting part of Acycle v3.1 SpectralMomentsGUI.m
#' (Mingsong Li).
#'
#' @param sm Spectral moments result to plot: the list returned by \code{\link{spectral_moments}} (its
#' \code{result} table, and \code{spectra} and \code{anchor} when present, are used). Required (\code{Default=NULL}
#' gives an error).
#' @param x_lab Label of the depth/time axis (character). \code{Default="depth"}.
#' @param y_lab Label of the frequency axis of the moments panel (and of the image when \code{yaxis = "frequency"})
#' (character). \code{Default="frequency (cycles/unit)"}.
#' @param image Draw the periodograms of all windows as an image (frequency or period against depth, the evolutive
#' spectrum the moments are computed from) above the moments, with the mean frequency and bandwidth overlain
#' (logical, TRUE or FALSE); only drawn when \code{sm} contains the window spectra. \code{Default=TRUE}.
#' @param yaxis Spectral axis of the image (character): \code{"frequency"} (cycles per unit, linear axis) or
#' \code{"period"} (log2 axis). \code{Default="frequency"} (the first element of \code{c("frequency", "period")}).
#' @param scale Colour scale of the image (character): \code{"log"} (log10 power, floored at the 2nd percentile of
#' the positive power), \code{"quantile"} (histogram-equalised, 64 classes) or \code{"linear"}.
#' \code{Default="log"} (the first element of \code{c("log", "quantile", "linear")}).
#' @param palette Name of a \code{\link[grDevices]{hcl.colors}} palette for the colour scale, a character string, for
#' example \code{"viridis"}, \code{"Plasma"}, \code{"Inferno"}, \code{"YlOrRd"}, \code{"Spectral"}, \code{"Blue-Red"}, \code{"Zissou 1"}; run \code{grDevices::hcl.pals()} to list all
#' 115 names. Low values get the first colour of the palette. Names of the viridis family (\code{"viridis"},
#' \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"cividis"}, \code{"rocket"}, \code{"mako"}, \code{"turbo"}) give the colours of the
#' viridis package, as in the wavelet plots. \code{Default="viridis"}.
#' @param freq_max Highest frequency shown in the image, in cycles per unit of depth/time (numeric, > 0); \code{NULL}
#' shows all frequencies. \code{Default=NULL}.
#' @param sr_lab Label of the sedimentation-rate axis of the sedimentation-rate panel, drawn when \code{sm}
#' contains a sedimentation-rate estimate (character). \code{Default="sed. rate"}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li.
#'
#' @references
#' Sinnesael, M., Zivanovic, M., De Vleeschouwer, D., Claeys, P. (2018).
#' Spectral moments in cyclostratigraphy: advantages and disadvantages
#' compared to more classic approaches. Paleoceanography and
#' Paleoclimatology, 33, 493-510. <doi:10.1029/2017PA003293>
#'
#' @examples
#' \donttest{
#' depth <- seq(0, 20, by = 0.01)
#' proxy <- sin(2 * pi * depth / (0.4 + 0.02 * depth))
#' sm <- spectral_moments(data = cbind(depth, proxy), window = 4, step = 0.1)
#' plot_spectral_moments(sm)
#' }
#'
#' @return No return value, called for its side effect (a base graphics plot).
#'
#' @export
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics legend
plot_spectral_moments <- function(sm = NULL,
    x_lab = "depth",
    y_lab = "frequency (cycles/unit)",
    sr_lab = "sed. rate",
    image = TRUE,
    yaxis = c("frequency", "period"),
    scale = c("log", "quantile", "linear"),
    palette = "viridis",
    freq_max = NULL){
  yaxis <- match.arg(yaxis)
  scale <- match.arg(scale)
  res <- sm$result
  has_sr <- "sedrate" %in% colnames(res)
  if (has_sr == TRUE){
    if (all(!is.finite(res[, "sedrate"]))){
      has_sr <- FALSE
    }
  }
  has_img <- isTRUE(image) && !is.null(sm$spectra) && !is.null(sm$spectra$power)
  has_an <- "sedrate_anchored" %in% colnames(res)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  n_panel <- 1 + has_img + has_sr + has_an
  if (n_panel > 1) {
    graphics::layout(matrix(seq_len(n_panel), ncol = 1), heights = c(if (has_img) { 2.2 }, 1, if (has_sr) { 1 }, if (has_an) { 1 }))
  }
  graphics::par(mar = c(4, 4.5, 1, 1))
  if (has_img){
    sp <- sm$spectra
    f <- sp$frequency
    P <- sp$power
    z <- sp$depth
    keep <- f > 0 & (if (is.null(freq_max)) { TRUE } else { f <= freq_max })
    f <- f[keep]
    P <- P[keep, , drop = FALSE]
    P[!is.finite(P)] <- 0
    if (yaxis == "period") {
      ya <- log2(1 / f)
      o <- order(ya)
      ya <- ya[o]
      P <- P[o, , drop = FALSE]
    } else {
      ya <- f
    }
    L <- switch(scale,
        log = {
          v <- P[P > 0]
          lo <- stats::quantile(v, 0.02)
          log10(pmax(P, lo))
        },
        quantile = {
          br <- unique(stats::quantile(P, seq(0, 1, length.out = 65)))
          matrix(findInterval(P, br, all.inside = TRUE), nrow(P), ncol(P))
        },
        linear = P)
    graphics::image(z, ya, t(L), col = .wr_hcl_cols(64, palette), xlab = "", xaxs = "i",
        ylab = if (yaxis == "period") { "period" } else { y_lab }, yaxt = if (yaxis == "period") { "n" } else { "s" },
        main = "Window periodograms; mean frequency (red) and bandwidth (dashed)")
    if (yaxis == "period") {
      tk <- 2^(floor(min(ya)):ceiling(max(ya)))
      graphics::axis(2, at = log2(tk), labels = signif(tk, 3))
    }
    conv <- function(fr) if (yaxis == "period") { log2(1 / pmax(fr, 1e-12)) } else { fr }
    graphics::lines(res[, "depth"], conv(res[, "mean_frequency"]), col = "red", lwd = 2)
    graphics::lines(res[, "depth"], conv(res[, "mean_frequency"] + res[, "bandwidth"]), col = "red", lty = 2)
    graphics::lines(res[, "depth"], conv(pmax(res[, "mean_frequency"] - res[, "bandwidth"], 1e-12)), col = "red", lty = 2)
  }
  yr <- range(c(res[, "mean_frequency"], res[, "bandwidth"]), finite = TRUE)
  graphics::plot(res[, "depth"], res[, "mean_frequency"], type = "l", xaxs = if (has_img) { "i" } else { "r" },
      col = "red", ylim = yr, xlab = x_lab, ylab = y_lab)
  graphics::lines(res[, "depth"], res[, "bandwidth"], col = "blue", lty = 4)
  if (has_an) {
    graphics::lines(res[, "depth"], res[, "moment_trend"], col = "black", lwd = 2)
    graphics::legend("topright", legend = c("mean frequency", "bandwidth", paste("trend of the", sub("_", " ", sm$anchor$moment))),
        col = c("red", "blue", "black"), lty = c(1, 4, 1), lwd = c(1, 1, 2), bty = "n")
    graphics::plot(res[, "depth"], res[, "sedrate_anchored"], type = "l", lwd = 2, xaxs = if (has_img) { "i" } else { "r" },
        xlab = x_lab, ylab = sr_lab, main = "Sedimentation rate from the anchored trend (Sinnesael et al. 2018)")
    graphics::points(sm$anchor$depth, sm$anchor$rate, pch = 23, bg = "red", cex = 1.4)
    graphics::abline(v = sm$anchor$depth, lty = 3)
    if (has_sr == TRUE) {
      graphics::lines(res[, "depth"], res[, "sedrate"], col = "grey50", lty = 2)
    }
  } else
  if (has_sr == TRUE){
    graphics::lines(res[, "depth"], res[, "bandwidth_trend"], col = "green3")
    graphics::legend("topright", legend = c("mean frequency", "bandwidth",
        "bandwidth trend"),
        col = c("red", "blue", "green3"), lty = c(1, 4, 1),
        bty = "n")
    graphics::plot(res[, "depth"], res[, "sedrate"], type = "l", xaxs = if (has_img) { "i" } else { "r" },
        xlab = x_lab, ylab = sr_lab)
  } else {
    graphics::legend("topright", legend = c("mean frequency", "bandwidth"),
        col = c("red", "blue"), lty = c(1, 4), bty = "n")
  }
  invisible(NULL)
}
