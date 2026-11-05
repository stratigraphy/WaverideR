#' Plot a chirp-matched EHA map, with or without the surrogate significance
#'
#' Two styles of "enhanced spectrum" for a \code{\link{analyse_chirp_eha}} result: \code{"significance"} shows the
#' line-to-continuum ratio divided by the 95th percentile of the warped surrogates (\code{Power_sig}, from
#' \code{\link{chirp_eha_significance}}; above 1, i.e. above 0 on the log scale, is significant), \code{"outline"}
#' shows the de-chirped map itself (\code{Power_norm}) with the significant cells outlined, and \code{"plain"} the
#' de-chirped map alone.
#'
#' @param ce Chirp-matched EHA result to plot: the list returned by \code{\link{analyse_chirp_eha}}; for
#' \code{style = "outline"} or \code{"significance"} it must first be passed through
#' \code{\link{chirp_eha_significance}} (which adds \code{Power_sig} and \code{Power_pval}). Required, no default.
#' @param style Plot style (character), one of:
#' \itemize{
#' \item \code{"outline"}: the de-chirped map (\code{Power_norm}) with the cells significant at \code{level} outlined
#' in white;
#' \item \code{"significance"}: log10 of the ratio of the observed map to the 95th percentile of the warped
#' surrogates (\code{Power_sig}), with a diverging blue-red scale centred on 0 (above 0 is significant);
#' \item \code{"plain"}: the de-chirped map alone (no significance needed).
#' }
#' \code{Default="outline"} (the first element of \code{c("outline", "significance", "plain")}).
#' @param track Sedimentation-rate track used to draw the expected periods of the astronomical cycles in \code{show}:
#' a \code{wr_track} object or a \code{\link{track_astro_auto}} result (its \code{track} element is used), i.e. a list
#' or data frame with \code{depth} and \code{sedrate} (depth units per kyr). Nothing is drawn when \code{track} or
#' \code{show} is \code{NULL}. \code{Default=NULL}.
#' @param show Astronomical periods in kyr (numeric vector, e.g. \code{c(405, 100, 40, 20)}) drawn as dotted lines at
#' period = \code{show} times \code{track$sedrate} along the record; only used with \code{track}. \code{Default=NULL}.
#' @param level Confidence level of the significance outline of \code{style = "outline"} (numeric, between 0 and 1):
#' cells with surrogate p-value \code{Power_pval <= 1 - level} are outlined. \code{Default=0.95}.
#' @param palette Name of a \code{\link[grDevices]{hcl.colors}} palette for the colour scale, a character string, for
#' example \code{"viridis"}, \code{"Plasma"}, \code{"Inferno"}, \code{"YlOrRd"}, \code{"Spectral"}, \code{"Blue-Red"}, \code{"Zissou 1"}; run \code{grDevices::hcl.pals()} to list all
#' 115 names. Low values get the first colour of the palette. Names of the viridis family (\code{"viridis"},
#' \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"cividis"}, \code{"rocket"}, \code{"mako"}, \code{"turbo"}) give the colours of the
#' viridis package, as in the wavelet plots. \code{Default="viridis"}.
#' @param scale Spacing of the colour breaks of the colour scale, as in all image plots of the package (character):
#' \code{"quantile"} (histogram-equalised), \code{"log"}, \code{"sqrt"}, \code{"linear"} or \code{"auto"}
#' (logarithmic for mostly empty maps or values spanning more than three orders of magnitude, quantiles otherwise);
#' see \code{\link{plot_wavelet}}. Not used for \code{style = "significance"}, which has a fixed diverging scale.
#' \code{Default="auto"}.
#' @param main Title of the plot (character); \code{NULL} gives a title that describes the chosen \code{style}.
#' \code{Default=NULL}.
#' @param xlab Label of the depth axis (character). \code{Default="Depth"}.
#' @param ylab Label of the period axis (character). \code{Default="Period"}.
#' @return Invisibly, the matrix shown.
#' @export
plot_chirp_eha <- function(ce, style = c("outline", "significance", "plain"), track = NULL, show = NULL, level = 0.95,
    palette = "viridis", main = NULL, xlab = "Depth", ylab = "Period", scale = "auto") {
  style <- match.arg(style)
  if (inherits(track, "wr_track_auto")) {
    track <- track$track
  }
  if (style != "plain" && is.null(ce$Power_sig)) {
    stop("run chirp_eha_significance() first, or use style = 'plain'")
  }
  M <- if (style == "significance") { log10(pmax(ce$Power_sig, 1e-6)) } else { ce$Power_norm }
  M <- .wr_fill_na(M)
  if (is.null(main)) {
    main <- switch(style,
        significance = "Lines against warped AR1 surrogates (log10 of the ratio to the 95th percentile; > 0 significant)",
        outline = paste0("De-chirped EHA, cells significant at ", round(100 * level), " % against warped surrogates outlined"),
        plain = "De-chirped EHA (lines against the continuum of each window)")
  }
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
    for (p in show) {
      graphics::lines(track$depth, log2(p * track$sedrate), col = if (style == "significance") { "black" } else { "white" }, lwd = 0.9, lty = 3)
    }
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
#' @param ce Chirp-matched EHA result: the list returned by \code{\link{analyse_chirp_eha}} run with
#' \code{criterion = "template"} (it must contain \code{template_score_all} and \code{sedrate_grid}). Required, no
#' default.
#' @param k Width, in windows, of the running median applied to the best rate per window (integer, >= 1; even values
#' are rounded up to the next odd number, 1 means no smoothing). \code{Default=11}.
#' @param plot Draw the rate-gram (map of the template score against sedimentation rate and depth, with the best rate
#' per window as white points and its running median as a white line) (logical, TRUE or FALSE); with \code{FALSE}
#' only the list is returned. \code{Default=TRUE}.
#' @param rate_unit Conversion factor and label of the rate axis: a vector \code{c(factor, "label")} whose first
#' element multiplies the rates (in depth units per kyr) for the axis labels, e.g. \code{c(100, "cm/kyr")} for a
#' depth axis in metres. Only changes the labelling. \code{Default=c(1, "depth units / kyr")}.
#' @param reference Reference sedimentation-rate curve drawn as a red dashed line for comparison: a two-column matrix
#' or data frame with depth in column 1 and rate (depth units per kyr) in column 2; \code{NULL} draws none.
#' \code{Default=NULL}.
#' @param main Title of the plot (character); \code{NULL} gives a default title. \code{Default=NULL}.
#' @return A list with \code{depth}, \code{sedrate} (rate grid), \code{map} (rates x windows), \code{best} (the
#' maximum per window) and \code{smooth} (its running median).
#' @export
rate_gram <- function(ce, k = 11, plot = TRUE, rate_unit = c(1, "depth units / kyr"), reference = NULL, main = NULL) {
  TS <- ce$template_score_all
  if (is.null(TS)) {
    stop("ce must come from analyse_chirp_eha(criterion = 'template')")
  }
  RG <- matrix(NA_real_, dim(TS)[1], dim(TS)[3])
  for (i in seq_len(dim(TS)[1])) {
    for (j in seq_len(dim(TS)[3])) {
      RG[i, j] <- max(TS[i, , j], na.rm = TRUE)
    }
  }
  r <- ce$sedrate_grid
  z <- if (!is.null(ce$positions)) { ce$positions } else { ce$axis.1 }
  best_i <- integer(ncol(RG))
  for (j in seq_len(ncol(RG))) {
    best_i[j] <- which.max(RG[, j])
  }
  best <- r[best_i]
  k <- max(1, 2 * floor(k / 2) + 1)
  sm <- if (k > 1 && length(best) > k) { exp(stats::runmed(log(best), k)) } else { best }
  if (plot) {
    f <- as.numeric(rate_unit[1])
    if (is.null(main)) {
      main <- "Rate-gram: astronomical template score against sedimentation rate per window"
    }
    graphics::image(z, log2(r), t(RG), col = .wr_hcl_cols(64, "viridis"), yaxt = "n", xlab = "Depth",
        ylab = paste0("Sedimentation rate (", rate_unit[2], ")"), main = main, xaxs = "i")
    rng <- range(r) * f
    ticks <- c(1, 2, 3, 5) * rep(10^((floor(log10(rng[1])) - 1):(ceiling(log10(rng[2])) + 1)), each = 4)
    ticks <- ticks[ticks >= rng[1] & ticks <= rng[2]]
    graphics::axis(2, at = log2(ticks / f), labels = ticks)
    graphics::points(z, log2(best), pch = 16, cex = 0.3, col = "white")
    graphics::lines(z, log2(sm), col = "white", lwd = 2)
    if (!is.null(reference)) {
      graphics::lines(reference[, 1], log2(reference[, 2]), col = "red", lwd = 1.3, lty = 2)
    }
  }
  invisible(list(depth = z, sedrate = r, map = RG, best = best, smooth = sm))
}
