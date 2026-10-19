#' Plot a chirp-matched EHA map, with or without the surrogate significance
#'
#' Two styles of "enhanced spectrum" for a \code{\link{analyse_chirp_eha}} result: \code{"significance"} shows the
#' line-to-continuum ratio divided by the 95th percentile of the warped surrogates (\code{Power_sig}, from
#' \code{\link{chirp_eha_significance}}; above 1, i.e. above 0 on the log scale, is significant), \code{"outline"}
#' shows the de-chirped map itself (\code{Power_norm}) with the significant cells outlined, and \code{"plain"} the
#' de-chirped map alone.
#'
#' @param ce Result of \code{\link{analyse_chirp_eha}} (after \code{\link{chirp_eha_significance}} for the first
#' two styles).
#' @param style \code{"outline"} (default), \code{"significance"} or \code{"plain"}.
#' @param track Optional track whose periods are drawn (\code{show} cycles in kyr at the tracked rate).
#' @param show Astronomical periods (kyr) to draw with \code{track}.
#' @param level Significance level for the outline \code{Default=0.95}.
#' @param palette \code{hcl.colors} palette \code{Default="viridis"}.
#' @param scale Colour scale as in all image plots of the package (see \code{\link{plot_wavelet}}) \code{Default="auto"}.
#' @param main Title.
#' @param xlab,ylab Axis labels.
#' @return Invisibly, the matrix shown.
#' @export
plot_chirp_eha <- function(ce, style = c("outline", "significance", "plain"), track = NULL, show = NULL, level = 0.95,
                           palette = "viridis", main = NULL, xlab = "Depth", ylab = "Period", scale = "auto") {
  style <- match.arg(style)
  if (inherits(track, "wr_track_auto")) track <- track$track
  if (style != "plain" && is.null(ce$Power_sig)) stop("run chirp_eha_significance() first, or use style = 'plain'")
  M <- if (style == "significance") log10(pmax(ce$Power_sig, 1e-6)) else ce$Power_norm
  M <- .wr_fill_na(M)
  if (is.null(main)) main <- switch(style,
    significance = "Lines against warped AR1 surrogates (log10 of the ratio to the 95th percentile; > 0 significant)",
    outline = paste0("De-chirped EHA, cells significant at ", round(100 * level), " % against warped surrogates outlined"),
    plain = "De-chirped EHA (lines against the continuum of each window)")
  if (style == "significance") {
    .wr_image_signif(ce$axis.1, ce$axis.2, M, yaxt = "n", xlab = xlab, ylab = ylab, main = main, xaxs = "i")
  } else {
    .wr_image(ce$axis.1, ce$axis.2, M, scale = scale, palette = palette, fill_na = FALSE, yaxt = "n", xlab = xlab,
              ylab = ylab, main = main, xaxs = "i")
  }
  ticks <- 2^(floor(min(ce$axis.2)):ceiling(max(ce$axis.2)))
  graphics::axis(2, at = log2(ticks), labels = ticks)
  if (style == "outline") {
    sig <- (ce$Power_pval <= 1 - level) * 1
    graphics::contour(ce$axis.1, ce$axis.2, t(sig), levels = 0.5, drawlabels = FALSE, add = TRUE, col = "white", lwd = 0.6)
  }
  if (!is.null(track) && !is.null(show)) {
    for (p in show) graphics::lines(track$depth, log2(p * track$sedrate), col = if (style == "significance") "black" else "white", lwd = 0.9, lty = 3)
  }
  invisible(M)
}

#' Rate-gram: the astronomical template as a spectrum over sedimentation rate, per window
#'
#' A purely spectral estimate of the sedimentation rate along a record, without any tracking: for every window of
#' \code{\link{analyse_chirp_eha}} (\code{criterion = "template"}) the template score (the summed log line-to-continuum
#' ratio at the periods of all astronomical cycles, the cyclostratigraphic counterpart of the harmonic product
#' spectrum used for pitch detection in speech and music) is taken at its best gradient and shown against the rate;
#' its maximum per window is the spectral rate estimate, optionally median-filtered over a few windows. On CIP Case 3
#' the maximum per window correlates 0.56 with the true rate, the median over 11 windows 0.73 (the tracked rate of
#' \code{\link{track_chirp_hmm}}, which adds continuity, 0.97).
#'
#' @param ce Result of \code{\link{analyse_chirp_eha}} with \code{criterion = "template"}.
#' @param k Width (windows) of the running median of the maximum \code{Default=11}.
#' @param plot Draw the rate-gram \code{Default=TRUE}.
#' @param rate_unit Factor and label of the rate axis, e.g. \code{c(100, "cm/kyr")}.
#' @param reference Optional data frame (depth, rate) drawn for comparison.
#' @param main Title.
#' @return A list with \code{depth}, \code{sedrate} (rate grid), \code{map} (rates x windows), \code{best} (the
#' maximum per window) and \code{smooth} (its running median).
#' @export
rate_gram <- function(ce, k = 11, plot = TRUE, rate_unit = c(1, "depth units / kyr"), reference = NULL, main = NULL) {
  TS <- ce$template_score_all
  if (is.null(TS)) stop("ce must come from analyse_chirp_eha(criterion = 'template')")
  RG <- apply(TS, c(1, 3), max, na.rm = TRUE)
  r <- ce$sedrate_grid; z <- if (!is.null(ce$positions)) ce$positions else ce$axis.1
  best <- r[apply(RG, 2, which.max)]
  k <- max(1, 2 * floor(k / 2) + 1)
  sm <- if (k > 1 && length(best) > k) exp(stats::runmed(log(best), k)) else best
  if (plot) {
    f <- as.numeric(rate_unit[1])
    if (is.null(main)) main <- "Rate-gram: astronomical template score against sedimentation rate per window"
    graphics::image(z, log2(r), t(RG), col = grDevices::hcl.colors(64, "viridis"), yaxt = "n", xlab = "Depth",
                    ylab = paste0("Sedimentation rate (", rate_unit[2], ")"), main = main, xaxs = "i")
    rng <- range(r) * f
    ticks <- c(1, 2, 3, 5) * rep(10^((floor(log10(rng[1])) - 1):(ceiling(log10(rng[2])) + 1)), each = 4)
    ticks <- ticks[ticks >= rng[1] & ticks <= rng[2]]
    graphics::axis(2, at = log2(ticks / f), labels = ticks)
    graphics::points(z, log2(best), pch = 16, cex = 0.3, col = "white")
    graphics::lines(z, log2(sm), col = "white", lwd = 2)
    if (!is.null(reference)) graphics::lines(reference[, 1], log2(reference[, 2]), col = "red", lwd = 1.3, lty = 2)
  }
  invisible(list(depth = z, sedrate = r, map = RG, best = best, smooth = sm))
}
