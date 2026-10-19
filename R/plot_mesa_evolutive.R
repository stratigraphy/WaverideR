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
#' @param x output of \code{\link{mesa_evolutive}}.
#' @param vertical If TRUE (default) the position is the vertical axis, increasing downwards (as
#' for a depth log), and the windows are panels side by side. If FALSE the position is horizontal
#' and the windows are stacked.
#' @param xaxis \code{"period"} (default) or \code{"frequency"}: the spectral axis.
#' @param period_lim Limits of the period axis (default: from twice the sampling interval to the
#' longest window).
#' @param frequency_lim Limits of the frequency axis, used when \code{xaxis = "frequency"}.
#' @param log_scale Kept for compatibility: \code{TRUE} is \code{scale = "log"}.
#' @param peaks If TRUE the peaks of every position are drawn as white points.
#' @param peak_min Smallest relative power (0 to 1) of a peak that is drawn (default 0.3).
#' @param track A \code{mesa_track} from \code{\link{track_mesa}} or a list of them, drawn as lines.
#' @param average If TRUE the average amplitude spectrum of every window is drawn next to its image,
#' with the period of its strongest peaks marked.
#' @param colors Optional vector of colours to interpolate (overrides \code{palette}).
#' @param log_axis Period or frequency axis on a log2 scale \code{Default=TRUE}: the spectrum of a cyclostratigraphic
#' record spans several octaves, and on a linear axis the long periods take most of the plot.
#' @param scale Colour scale as in all image plots of the package (see \code{\link{plot_wavelet}}) \code{Default="auto"};
#' a maximum-entropy spectrum spans many decades, so \code{"auto"} is logarithmic for it.
#' @param palette \code{hcl.colors} palette \code{Default="viridis"}.
#' @param normalise \code{"window"} (each analysis window relative to its maximum, default), \code{"period"} (each
#' period relative to its mean along the record, which removes the red-noise slope and shows where a cycle is
#' stronger or weaker) or \code{"none"}.
#' @param position_lim Limits of the position axis.
#' @param main Overall title.
#' @param ... currently unused.
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
  xaxis <- match.arg(xaxis); normalise <- match.arg(normalise)
  if (isTRUE(log_scale)) scale <- "log"
  nw <- length(x$windows)
  freq <- x$frequency
  pos <- x$position
  if (is.null(position_lim)){
    position_lim <- range(pos)
  }
  longest <- max(vapply(x$windows, function(w) w$window, numeric(1)))
  if (is.null(period_lim)){
    period_lim <- c(2 * x$dt, longest)
  }
  if (is.null(frequency_lim)){
    frequency_lim <- range(freq)
  }
  if (!is.null(track) && inherits(track, "mesa_track")){
    track <- list(track)
  }
  keep <- if (xaxis == "period") freq > 0 & 1 / freq >= period_lim[1] & 1 / freq <= period_lim[2] else
    freq >= frequency_lim[1] & freq <= frequency_lim[2]
  if (log_axis) keep <- keep & freq > 0
  if (sum(keep) < 2){
    stop("the spectral axis limits contain fewer than 2 frequencies")
  }
  fk <- freq[keep]
  sv <- if (xaxis == "period") rev(1 / fk) else fk
  slab <- if (xaxis == "period") "Period" else "Frequency"
  if (log_axis) {
    sv_lab <- sv; sv <- log2(sv)
    slab <- paste(slab, "(log2 axis)")
  }
  log_ticks <- function(side) {
    if (!log_axis) return(graphics::axis(side, las = if (side == 2) 1 else 0))
    tk <- 2^(floor(min(sv)):ceiling(max(sv)))
    graphics::axis(side, at = log2(tk), labels = signif(tk, 3), las = if (side == 2) 1 else 0)
  }
  tx <- function(v) if (log_axis) log2(v) else v
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
      window = sweep(pw, 2, apply(pw, 2, max, na.rm = TRUE), "/"),           # each window relative to its maximum
      period = pw / rowMeans(pw, na.rm = TRUE),                               # each period relative to its mean
      none = pw)
    rel[!is.finite(rel)] <- NA
    if (xaxis == "period"){
      rel <- rel[nrow(rel):1, , drop = FALSE]
    }
    rel <- .wr_fill_na(rel)
    br <- .wr_colour_breaks(rel, scale)
    cols <- if (is.null(colors)) grDevices::hcl.colors(length(br) - 1, palette) else grDevices::colorRampPalette(colors)(length(br) - 1)
    ttl <- sprintf("Window = %s", format(w$window))
    if (vertical){
      graphics::par(mar = c(if (average) 0.5 else 3.5, 3.8, 2.5, 0.5), mgp = c(2.2, 0.6, 0), las = 0)
      graphics::image(sv, pos, rel, col = cols, breaks = br, xlab = "", ylab = if (wi == 1) "Position" else "",
                      ylim = rev(position_lim), axes = FALSE, useRaster = FALSE)
      graphics::axis(2, las = 1)
      if (!average){
        log_ticks(1)
        graphics::mtext(slab, 1, 2.2, cex = 0.8)
      }
    } else {
      graphics::par(mar = c(if (average) 3.5 else 3.5, 4, 1.8, 0.5), mgp = c(2.2, 0.6, 0), las = 0)
      graphics::image(pos, sv, t(rel), col = cols, breaks = br, xlab = "Position", ylab = slab,
                      xlim = position_lim, axes = FALSE, useRaster = FALSE)
      graphics::axis(1)
      log_ticks(2)
    }
    graphics::box()
    graphics::mtext(ttl, 3, 0.4, font = 2, cex = 0.8)
    if (peaks && nrow(w$peaks)){
      pk <- w$peaks
      sx <- tx(if (xaxis == "period") pk$period else pk$frequency)
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
        sx <- tx(if (xaxis == "period") tr$period else tr$frequency)
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
      sa <- tx(if (xaxis == "period") av$period else av$frequency)
      o <- order(sa)
      sa <- sa[o]
      am <- av$amplitude[o]
      if (vertical){
        graphics::par(mar = c(3.5, 3.8, 0.3, 0.5))
        plot(sa, am, type = "n", xlim = range(sv), ylim = c(0, max(am) * 1.15), xlab = "", ylab = if (wi == 1) "Amplitude" else "",
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
          graphics::text(pks$frequency, pks$power, labels = format(signif(if (log_axis) 2^pks$frequency else pks$frequency, 3)),
                         pos = 3, cex = 0.65, col = "red3", xpd = NA)
        } else {
          graphics::text(pks$power, pks$frequency, labels = format(signif(if (log_axis) 2^pks$frequency else pks$frequency, 3)),
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
