#' @title Plot an evolutive MESA
#'
#' @description Plots the spectral image of an evolutive maximum entropy spectrum
#' (\code{\link{mesa_evolutive}}): one panel for every window length, side by side, with the
#' position (depth or time) on one axis and the period or frequency on the other, as in the
#' spectral image columns of log-analysis software. The spectrum of every window position is
#' divided by its maximum, so that the ridges of the spectral peaks show along the whole record.
#' The peaks of every position, the tracked ridge (\code{\link{track_mesa}}) and the average
#' amplitude spectrum (\code{\link{mesa_average}}) can be added.
#'
#' @param x Evolutive maximum entropy spectrum to plot: the list returned by \code{\link{mesa_evolutive}} (one or
#' more window lengths). Required, no default.
#' @param vertical Orientation of the plot (logical). \code{TRUE}: the position (depth/time) is the vertical axis,
#' increasing downwards as for a depth log, and the window lengths are panels side by side; \code{FALSE}: the
#' position is the horizontal axis and the window lengths are stacked panels. \code{Default=TRUE}.
#' @param xaxis Spectral axis (character): \code{"period"} or \code{"frequency"} (cycles per unit of position).
#' \code{Default="period"} (the first element of \code{c("period", "frequency")}).
#' @param period_lim Range of the period axis, in the unit of the position axis (numeric vector of length 2,
#' \code{c(min, max)}); used when \code{xaxis = "period"}. \code{NULL} runs from twice the sampling interval to the
#' longest window length. \code{Default=NULL}.
#' @param frequency_lim Range of the frequency axis, in cycles per unit of position (numeric vector of length 2);
#' used when \code{xaxis = "frequency"}. \code{NULL} uses the full frequency range. \code{Default=NULL}.
#' @param log_scale Kept for compatibility (logical): \code{TRUE} is the same as \code{scale = "log"}; \code{FALSE}
#' uses \code{scale}. \code{Default=FALSE}.
#' @param peaks Draw the spectral peaks of every window position as white points (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param peak_min Smallest relative power of a peak that is drawn with \code{peaks = TRUE} (numeric, 0 to 1, power
#' relative to the maximum of the spectrum at the same position). \code{Default=0.3}.
#' @param track Ridge track(s) drawn as red lines on the panel of the matching window length: a \code{mesa_track}
#' object returned by \code{\link{track_mesa}} or a list of them; \code{NULL} draws none. \code{Default=NULL}.
#' @param average Draw the average amplitude spectrum of every window length (\code{\link{mesa_average}}) beside or
#' below its image, with the periods (or frequencies) of its three strongest peaks labelled (logical).
#' \code{Default=FALSE}.
#' @param colors Colours to interpolate for the colour scale (character vector of colours, e.g.
#' \code{c("white", "orange", "red")}); if given it overrides \code{palette}. \code{NULL} uses \code{palette}.
#' \code{Default=NULL}.
#' @param position_lim Range of the position (depth/time) axis (numeric vector of length 2); \code{NULL} uses the
#' full range of window positions. \code{Default=NULL}.
#' @param main Overall title, drawn in the outer top margin (character); \code{NULL} draws none. \code{Default=NULL}.
#' @param log_axis Draw the period or frequency axis on a log2 scale (logical): the spectrum of a
#' cyclostratigraphic record spans several octaves, and on a linear axis the long periods take most of the plot.
#' \code{Default=TRUE}.
#' @param scale Spacing of the colour breaks of the colour scale, as in all image plots of the package (character):
#' \code{"quantile"} (histogram-equalised), \code{"log"}, \code{"sqrt"}, \code{"linear"} or \code{"auto"}
#' (logarithmic for mostly empty maps or values spanning more than three orders of magnitude, quantiles otherwise;
#' see \code{\link{plot_wavelet}}). A maximum-entropy spectrum spans many decades, so \code{"auto"} is usually
#' logarithmic for it. \code{Default="auto"}.
#' @param palette Name of a \code{\link[grDevices]{hcl.colors}} palette for the colour scale, a character string, for
#' example \code{"viridis"}, \code{"Plasma"}, \code{"Inferno"}, \code{"YlOrRd"}, \code{"Spectral"}, \code{"Blue-Red"}, \code{"Zissou 1"}; run \code{grDevices::hcl.pals()} to list all
#' 115 names. Low values get the first colour of the palette. Names of the viridis family (\code{"viridis"},
#' \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"cividis"}, \code{"rocket"}, \code{"mako"}, \code{"turbo"}) give the colours of the
#' viridis package, as in the wavelet plots. \code{Default="viridis"}.
#' @param normalise How the spectrum is scaled before plotting (character), one of:
#' \code{"window"} (the spectrum of each window position divided by its maximum, so the ridges show along the whole
#' record), \code{"period"} (each period divided by its mean along the record, which removes the red-noise slope and
#' shows where a cycle is stronger or weaker) or \code{"none"} (raw power).
#' \code{Default="window"} (the first element of \code{c("window", "period", "none")}).
#' @param ... Not used; accepted (and ignored) so that extra arguments do not give an error.
#'
#' @return Called for its plot; invisibly returns \code{x}.
#' @examples
#' set.seed(2)
#' t <- seq(0, 800, by = 2)
#' y <- sin(2 * pi * cumsum(1 / ifelse(t < 400, 40, 80)) * 2) + rnorm(length(t), sd = 0.4)
#' ev <- mesa_evolutive(data.frame(t, y), window = c(150, 250), order = 20)
#' plot_mesa_evolutive(ev, average = TRUE)
#' @importFrom graphics image layout par axis box abline lines points mtext polygon
#' @importFrom grDevices colorRampPalette
#' @export
plot_mesa_evolutive <- function(x,
    vertical = TRUE,
    xaxis = c("period", "frequency"),
    period_lim = NULL,
    frequency_lim = NULL,
    log_scale = FALSE,
    peaks = FALSE,
    peak_min = 0.3,
    track = NULL,
    average = FALSE,
    colors = NULL,
    position_lim = NULL,
    main = NULL,
    log_axis = TRUE,
    scale = "auto",
    palette = "viridis",
    normalise = c("window", "period", "none"),
    ...) {
  xaxis <- match.arg(xaxis)
  normalise <- match.arg(normalise)
  if (isTRUE(log_scale)) {
    scale <- "log"
  }
  nw <- length(x$windows)
  freq <- x$frequency
  pos <- x$position
  if (is.null(position_lim)){
    position_lim <- range(pos)
  }
  lens <- numeric(nw)
  for (i in seq_len(nw)) {
    lens[i] <- x$windows[[i]]$window
  }
  longest <- max(lens)
  if (is.null(period_lim)){
    period_lim <- c(2 * x$dt, longest)
  }
  if (is.null(frequency_lim)){
    frequency_lim <- range(freq)
  }
  if (!is.null(track) && inherits(track, "mesa_track")){
    track <- list(track)
  }
  keep <- if (xaxis == "period") {
    freq > 0 & 1 / freq >= period_lim[1] & 1 / freq <= period_lim[2]
  } else {
    freq >= frequency_lim[1] & freq <= frequency_lim[2]
  }
  if (log_axis) {
    keep <- keep & freq > 0
  }
  if (sum(keep) < 2){
    stop("the spectral axis limits contain fewer than 2 frequencies")
  }
  fk <- freq[keep]
  sv <- if (xaxis == "period") { rev(1 / fk) } else { fk }
  slab <- if (xaxis == "period") { "Period" } else { "Frequency" }
  if (log_axis) {
    sv_lab <- sv
    sv <- log2(sv)
    slab <- paste(slab, "(log2 axis)")
  }
  log_ticks <- function(side) {
    if (!log_axis) {
      return(graphics::axis(side, las = if (side == 2) { 1 } else { 0 }))
    }
    tk <- 2^(floor(min(sv)):ceiling(max(sv)))
    graphics::axis(side, at = log2(tk), labels = signif(tk, 3), las = if (side == 2) { 1 } else { 0 })
  }
  tx <- function(v) if (log_axis) { log2(v) } else { v }
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  if (average){
    if (vertical){
      graphics::layout(matrix(seq_len(2 * nw), 2, nw), heights = c(3, 1))
    } else {
      graphics::layout(matrix(seq_len(2 * nw), nw, 2, byrow = TRUE), widths = c(3, 1))
    }
  } else if (vertical){
    graphics::layout(matrix(seq_len(nw), 1, nw))
  } else {
    graphics::layout(matrix(seq_len(nw), nw, 1))
  }
  for (wi in seq_len(nw)){
    w <- x$windows[[wi]]
    pw <- w$power[keep, , drop = FALSE]
    rel <- switch(normalise,
        window = {                                                              # each window relative to its maximum
          cmax <- numeric(ncol(pw))
          for (j in seq_len(ncol(pw))) {
            cmax[j] <- max(pw[, j], na.rm = TRUE)
          }
          sweep(pw, 2, cmax, "/")
        },
        period = pw / rowMeans(pw, na.rm = TRUE),                               # each period relative to its mean
        none = pw)
    rel[!is.finite(rel)] <- NA
    if (xaxis == "period"){
      rel <- rel[nrow(rel):1, , drop = FALSE]
    }
    rel <- .wr_fill_na(rel)
    br <- .wr_colour_breaks(rel, scale)
    cols <- if (is.null(colors)) { .wr_hcl_cols(length(br) - 1, palette) } else { grDevices::colorRampPalette(colors)(length(br) - 1) }
    ttl <- sprintf("Window = %s", format(w$window))
    if (vertical){
      graphics::par(mar = c(if (average) { 0.5 } else { 3.5 }, 3.8, 2.5, 0.5), mgp = c(2.2, 0.6, 0), las = 0)
      graphics::image(sv, pos, rel, col = cols, breaks = br, xlab = "", ylab = if (wi == 1) { "Position" } else { "" },
          ylim = rev(position_lim), axes = FALSE, useRaster = FALSE)
      graphics::axis(2, las = 1)
      if (!average){
        log_ticks(1)
        graphics::mtext(slab, 1, 2.2, cex = 0.8)
      }
    } else {
      graphics::par(mar = c(if (average) { 3.5 } else { 3.5 }, 4, 1.8, 0.5), mgp = c(2.2, 0.6, 0), las = 0)
      graphics::image(pos, sv, t(rel), col = cols, breaks = br, xlab = "Position", ylab = slab,
          xlim = position_lim, axes = FALSE, useRaster = FALSE)
      graphics::axis(1)
      log_ticks(2)
    }
    graphics::box()
    graphics::mtext(ttl, 3, 0.4, font = 2, cex = 0.8)
    if (peaks && nrow(w$peaks)){
      pk <- w$peaks
      sx <- tx(if (xaxis == "period") { pk$period } else { pk$frequency })
      ok <- sx >= min(sv) & sx <= max(sv) & pk$rel_power >= peak_min
      if (vertical){
        graphics::points(sx[ok], pk$position[ok], pch = 16, cex = 0.35, col = "white")
      } else {
        graphics::points(pk$position[ok], sx[ok], pch = 16, cex = 0.35, col = "white")
      }
    }
    if (!is.null(track)){
      for (tr in track){
        if (abs(attr(tr, "window") - w$window) > 1e-9 * w$window){
          next
        }
        sx <- tx(if (xaxis == "period") { tr$period } else { tr$frequency })
        if (vertical){
          graphics::lines(sx, tr$position, col = "red3", lwd = 2)
          graphics::lines(sx, tr$position, col = "white", lwd = 0.6, lty = 3)
        } else {
          graphics::lines(tr$position, sx, col = "red3", lwd = 2)
          graphics::lines(tr$position, sx, col = "white", lwd = 0.6, lty = 3)
        }
      }
    }
    if (average){
      av <- mesa_average(x, wi)
      av <- av[keep, ]
      sa <- tx(if (xaxis == "period") { av$period } else { av$frequency })
      o <- order(sa)
      sa <- sa[o]
      am <- av$amplitude[o]
      if (vertical){
        graphics::par(mar = c(3.5, 3.8, 0.3, 0.5))
        plot(sa, am, type = "n", xlim = range(sv), ylim = c(0, max(am) * 1.15), xlab = "", ylab = if (wi == 1) { "Amplitude" } else { "" },
            las = 1, xaxs = "i", xaxt = "n")
        log_ticks(1)
        graphics::polygon(c(sa[1], sa, sa[length(sa)]), c(0, am, 0), col = "grey25", border = NA)
        graphics::mtext(paste("Average amplitude spectrum -", slab), 1, 2.2, cex = 0.7)
      } else {
        graphics::par(mar = c(3.5, 0.5, 1.8, 3.5))
        plot(am, sa, type = "n", ylim = range(sv), xlim = c(0, max(am) * 1.15), xlab = "", ylab = "",
            yaxt = "n", yaxs = "i")
        graphics::polygon(c(0, am, 0), c(sa[1], sa, sa[length(sa)]), col = "grey25", border = NA)
        graphics::mtext("Average amplitude", 1, 2.2, cex = 0.7)
      }
      pks <- .mesa_peaks(sa, am, 3)
      if (nrow(pks)){
        if (vertical){
          graphics::text(pks$frequency, pks$power, labels = format(signif(if (log_axis) { 2^pks$frequency } else { pks$frequency }, 3)),
              pos = 3, cex = 0.65, col = "red3", xpd = NA)
        } else {
          graphics::text(pks$power, pks$frequency, labels = format(signif(if (log_axis) { 2^pks$frequency } else { pks$frequency }, 3)),
              pos = 4, cex = 0.65, col = "red3", xpd = NA)
        }
      }
    }
  }
  if (!is.null(main)){
    graphics::mtext(main, 3, outer = TRUE, font = 2)
  }
  invisible(x)
}
