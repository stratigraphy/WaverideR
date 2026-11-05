# Internal building blocks of plot_wavelet(). The 16 layouts of plot_wavelet() (data panel, global spectrum and
# MTM panels switched on/off, vertical or horizontal) are built from the same few panels. Every helper reproduces
# the base graphics calls of the original, long-hand implementation exactly.
# `p` is the list with the prepared plot state that plot_wavelet() builds.

.pw_layout <- function(m, nrow, ncol, heights, widths) {
  layout.matrix <- matrix(m, nrow = nrow, ncol = ncol, byrow = TRUE)
  graphics::layout(mat = layout.matrix, heights = heights, widths = widths)
}

.pw_levels <- function(p) {
  # default colour scale: quantiles for a standard CWT, "auto" (log for the mostly empty map) for a synchrosqueezed one
  sq <- isTRUE(p$wavelet$synchrosqueezed)
  sc <- if (!is.null(p$scale)) { p$scale } else if (sq) { "auto" } else { "quantile" }
  sparse <- if (sq) { TRUE } else { NULL }
  br <- .wr_colour_breaks(p$pmax_avg_sel, sc, p$n.levels, sparse = sparse)
  as.numeric(br)
}

# Colour keys

# colour key below the scalogram (vertical plots). power_lab: where the "Power" label goes
.pw_key_bottom <- function(p, power_lab = c("side2", "below3", "below4")) {
  power_lab <- match.arg(power_lab)
  image(x = seq(from = 0, to = p$n.levels), y = 1, z = t(matrix(p$levels, nrow = 1)),
      col = p$key.cols, breaks = p$levels, useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  axis(1, lwd = p$lwd.axis, at = p$key.marks, labels = NA, tck = 0.02, tcl = 1.24)
  mtext(p$key.labels, side = 1, at = p$key.marks, line = 0.5, las = 2, font = par()$font.axis, cex = par()$cex.axis)
  if (power_lab == "side2") {
    mtext(c("Power"), side = 2, at = mean(p$key.marks), line = 0.5, las = 2, font = par()$font.axis, cex = par()$cex.axis)
  } else if (power_lab == "below3") {
    mtext(c("Power"), side = 1, at = stats::median(p$key.marks), line = 3, las = 2, font = par()$font.axis,
        cex = par()$cex.axis, xpd = NA)
  } else {
    mtext(c("Power"), side = 1, at = mean(p$key.marks), line = 4, font = par()$font.axis, cex = par()$cex.axis, las = 1)
  }
}

# colour key to the right of the scalogram (horizontal plots)
.pw_key_right <- function(p) {
  image(y = seq(from = 0, to = p$n.levels), x = 1, z = (matrix(p$levels, nrow = 1)), col = p$key.cols,
      breaks = p$levels, useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  axis(4, lwd = p$lwd.axis, at = p$key.marks, labels = NA, tck = 0.02, tcl = 1.24)
  mtext(p$key.labels, side = 4, at = p$key.marks, line = 0.5, las = 2, font = par()$font.axis, cex = par()$cex.axis)
  mtext(c("Power"), side = 4, at = mean(p$key.marks), line = 0.5, las = 2, font = par()$font.axis, cex = par()$cex.axis)
}

# colour key to the left of the global spectrum (vertical plots with data, avg and no MTM)
.pw_key_left <- function(p) {
  image(y = seq(from = 0, to = p$n.levels), x = 1, z = (matrix(p$levels, nrow = 1)), col = p$key.cols,
      breaks = p$levels, useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  axis(2, lwd = p$lwd.axis, at = p$key.marks, labels = NA, tck = 0.02, tcl = 1.24)
  mtext(p$key.labels, side = 2, at = p$key.marks, line = 0.5, las = 2, font = par()$font.axis, cex = par()$cex.axis)
  mtext(c("Power"), side = 1, at = 1, line = 0.5, font = par()$font.axis, cex = par()$cex.axis, las = 1, xpd = NA)
}

# colour key above the scalogram (horizontal plots with a data panel on top)
.pw_key_top <- function(p) {
  image(x = seq(from = 0, to = p$n.levels), y = 1, z = t(matrix(p$levels, nrow = 1)), col = p$key.cols,
      breaks = p$levels, useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "Power", ylab = "")
  axis(1, lwd = p$lwd.axis, at = p$key.marks, labels = NA, tck = 0.02, tcl = 1.2)
  mtext(p$key.labels, side = 1, at = p$key.marks, line = 0.1, las = 2, cex = 0.75)
}

# Data panel

.pw_data_v <- function(p, abline_h = FALSE) {
  if (!is.null(p$data2)) {                        # cross spectra: two records, the second in red on top
    plot(x = p$data_y, y = p$data_x, type = "l", col = "black", yaxs = "i", ylim = p$xlim_vals, xlab = "proxy_1",
        ylab = p$x_lab, xaxt = "s", yaxt = "s")
    par(new = TRUE)
    plot(x = p$data2$y, y = p$data2$x, type = "l", col = "red", yaxs = "i", ylim = p$xlim_vals, xlab = "", ylab = "",
        xaxt = "n", yaxt = "n")
    axis(side = 3, col = "red", col.axis = "red")
    mtext("proxy_2", side = 3, line = 2, col = "red")
  } else {
    plot(x = p$data_y, y = p$data_x, type = "l", yaxs = "i", xlab = "proxy value", ylab = p$x_lab,
        ylim = p$xlim_vals)
  }
  if (abline_h && is.null(p$add_abline_h) != TRUE) {
    abline(h = p$add_abline_h)
  }
}

.pw_data_h <- function(p) {
  if (!is.null(p$data2)) {
    plot(y = p$data_y, x = p$data_x, type = "l", col = "black", xaxs = "i", xlim = p$xlim_vals, ylab = "proxy_1",
        xlab = p$x_lab, yaxt = "s", xaxt = "s")
    par(new = TRUE)
    plot(y = p$data2$y, x = p$data2$x, type = "l", col = "red", xaxs = "i", xlim = p$xlim_vals, ylab = "", xlab = "",
        yaxt = "n", xaxt = "n")
    axis(side = 4, col = "red", col.axis = "red")
    mtext("proxy_2", side = 4, line = 2, col = "red")
  } else {
    plot(y = p$data_y, x = p$data_x, type = "l", xaxs = "i", xlab = "", ylab = "proxy value", xaxt = "n",
        xlim = p$xlim_vals)
  }
  if (is.null(p$add_abline_v) != TRUE) {
    abline(v = p$add_abline_v)
  }
}

# Scalogram

# vertical: depth on the y axis. ylab = NULL leaves the argument out; yaxt_n suppresses the depth axis
.pw_env <- new.env(parent = emptyenv())   # geometry of the last scalogram (sizes the phase arrows)

.pw_scalogram_v <- function(p, ylab = NULL, yaxt_n = TRUE) {
  args <- list(y = p$depth, x = p$ax2, z = (p$z), col = p$key.cols, breaks = p$levels, useRaster = TRUE,
      xlab = p$periodlab)
  if (!is.null(ylab)) {
    args$ylab <- ylab
  }
  if (yaxt_n) {
    args$yaxt <- "n"
  }
  args <- c(args, list(xaxt = "n", main = p$main, ylim = p$xlim_vals, xlim = log2(p$ylim_vals)))
  do.call(image, args)
  .pw_env$geom <- list(usr = par("usr"), pin = par("pin"), horizontal = FALSE)
  .wr_add_signif(p$add_signif, p$wavelet, horizontal = FALSE)
}

.pw_scalogram_h <- function(p) {
  image(x = p$depth, y = p$ax2, z = t(p$z), col = p$key.cols, breaks = p$levels, useRaster = TRUE,
      ylab = p$periodlab, xlab = p$x_lab, axes = TRUE, yaxt = "n", main = p$main, xlim = p$xlim_vals,
      ylim = log2(p$ylim_vals))
  .pw_env$geom <- list(usr = par("usr"), pin = par("pin"), horizontal = TRUE)
  .wr_add_signif(p$add_signif, p$wavelet, horizontal = TRUE)
}

# p-value overlay (white veil on the cells that are not significant) or, without p-values, the cone of influence
.pw_overlay <- function(p, horizontal) {
  w <- p$wavelet
  if (!isTRUE(p$overlay)) {
    return(invisible(NULL))
  }
  veil <- c(rgb(red = 255, green = 255, blue = 255, alpha = 125, maxColorValue = 255),
      rgb(red = 255, green = 255, blue = 255, alpha = 0, maxColorValue = 255))
  if (p$add_pval == TRUE & is.null(w$Power.pval) != TRUE) {
    par(new = TRUE)
    if (horizontal) {
      image(x = p$depth, y = p$ax2, z = t(w$Power.pval < p$pval_cutoff) * 1, col = veil, useRaster = TRUE,
          xlab = "", ylab = "", yaxt = "n", xaxt = "n", main = "", xlim = p$xlim_vals, ylim = log2(p$ylim_vals))
    } else {
      image(y = p$depth, x = p$ax2, z = (w$Power.pval < p$pval_cutoff) * 1, col = veil, useRaster = TRUE,
          xlab = "", ylab = "", yaxt = "n", xaxt = "n", main = "", ylim = p$xlim_vals, xlim = log2(p$ylim_vals))
    }
  } else if (is.null(w$coi.1)) {
    invisible(NULL)
  } else if (horizontal) {
    polygon(x = w$coi.1, y = w$coi.2, border = NA, col = rgb(1, 1, 1, 0.5), xlim = p$xlim_vals,
        ylim = log2(p$ylim_vals))
  } else {
    polygon(y = w$coi.1, x = w$coi.2, border = NA, col = rgb(1, 1, 1, 0.5), ylim = p$xlim_vals,
        xlim = log2(p$ylim_vals))
  }
}

# frame, period axis and the user supplied lines/points/ablines of the scalogram
.pw_tail <- function(p, horizontal) {
  box(lwd = p$lwd.axis)
  period.tick <- unique(trunc(p$ax2))
  period.tick[period.tick < p$tick_min] <- NA
  if (!is.null(p$tick_max)) {
    period.tick[period.tick > p$tick_max] <- NA
  }
  period.tick <- na.omit(period.tick)
  period.tick.label <- 2^(period.tick)
  if (horizontal) {
    axis(2, lwd = p$lwd.axis, at = period.tick, labels = NA, tck = p$periodtck, tcl = p$periodtcl)
    axis(4, lwd = p$lwd.axis, at = period.tick, labels = NA, tck = p$periodtck, tcl = p$periodtcl)
    mtext(period.tick.label, side = 2, at = period.tick, las = 2, line = par()$mgp[2] - 0.5,
        font = par()$font.axis, cex = par()$cex.axis)
    if (is.null(p$add_lines) != TRUE) {
      for (i in 2:ncol(p$add_lines)) {
        lines(p$add_lines[, 1], log2(p$add_lines[, i]))
      }
    }
    if (is.null(p$add_points) != TRUE) {
      for (i in 2:ncol(p$add_points)) {
        points(p$add_points[, 1], log2(p$add_points[, i]))
      }
    }
    if (p$add_MTM_peaks == TRUE & is.na(p$MTM_res_2[1, 2]) == FALSE) {
      abline(h = log2(1 / p$MTM_res_2[, 1]), col = "black", lty = 3)
    }
    if (is.null(p$add_abline_h) != TRUE) {
      abline(h = log2(p$add_abline_h))
    }
    if (is.null(p$add_abline_v) != TRUE) {
      abline(v = p$add_abline_v)
    }
  } else {
    axis(1, lwd = p$lwd.axis, at = period.tick, labels = NA, tck = p$periodtck, tcl = p$periodtcl)
    mtext(period.tick.label, side = 1, at = period.tick, las = 2, line = par()$mgp[2] - 0.5,
        font = par()$font.axis, cex = par()$cex.axis)
    if (is.null(p$add_lines) != TRUE) {
      for (i in 2:ncol(p$add_lines)) {
        lines(y = p$add_lines[, 1], x = log2(p$add_lines[, i]))
      }
    }
    if (is.null(p$add_points) != TRUE) {
      for (i in 2:ncol(p$add_points)) {
        points(y = p$add_points[, 1], x = log2(p$add_points[, i]))
      }
    }
    if (p$add_MTM_peaks == TRUE & is.na(p$MTM_res_2[1, 2]) == FALSE) {
      abline(v = log2(1 / p$MTM_res_2[, 1]), col = "black", lty = 3)
    }
    if (is.null(p$add_abline_h) != TRUE) {
      abline(h = (p$add_abline_h))
    }
    if (is.null(p$add_abline_v) != TRUE) {
      abline(v = log2(p$add_abline_v))
    }
  }
}

# Global wavelet spectrum

# global spectrum beside a vertical scalogram (period on the x axis)
# xlab_empty: xlab = "" given; title_ylab: ylab drawn with title(); pval_xpd: p-value lines may leave the panel;
# abline_v: user abline at the end
.pw_gws_v <- function(p, xlab_empty = FALSE, title_ylab = FALSE, pval_xpd = FALSE, abline_v = FALSE) {
  w <- p$wavelet
  args <- list(y = p$avg, x = p$per, log = "x", type = "l", xaxs = "i", xaxt = "n", ylab = "Wt. power")
  if (xlab_empty) {
    args$xlab <- ""
  }
  args <- c(args, list(xaxs = "i", xlim = p$ylim_vals, ylim = .wr_gws_lim(p$add_signif, p$avg)))
  do.call(plot, args)
  .wr_add_signif_gws(p$add_signif, horizontal_panel = FALSE)
  if (title_ylab) {
    title(ylab = "Wt. power", xpd = NA)
  }
  if (p$add_pval == TRUE) {
    par(new = TRUE)
    plot(x = p$per, y = w$Power.avg.pval, col = "red", xlim = p$ylim_vals, ylim = c(0, 1), axes = FALSE, type = "l",
        xlab = "", ylab = "", log = "x", xaxs = "i", xaxt = "n", yaxs = "i")
    axis(side = 4, at = pretty(c(0, 1)))
    if (pval_xpd) {
      abline(h = p$pval_abline, col = "grey", lwd = 2, lty = 2, xpd = TRUE)
    } else {
      abline(h = p$pval_abline, col = "grey", lwd = 2, lty = 2)
    }
  }
  if ((abline_v || isTRUE(p$gws_abline_all)) && is.null(p$add_abline_v) != TRUE) {
    abline(v = (p$add_abline_v), xpd = FALSE)
  }
}

# global spectrum beside a horizontal scalogram (period on the y axis)
.pw_gws_h <- function(p) {
  w <- p$wavelet
  plot(x = p$avg, y = p$per, log = "y", type = "l", yaxs = "i", yaxt = "n", xlab = "Wt. power", xaxs = "i",
      ylim = p$ylim_vals, xlim = .wr_gws_lim(p$add_signif, p$avg))
  .wr_add_signif_gws(p$add_signif, horizontal_panel = TRUE)
  if (p$add_pval == TRUE) {
    par(new = TRUE)
    plot(y = p$per, x = w$Power.avg.pval, col = "red", ylim = p$ylim_vals, xlim = c(0, 1), axes = FALSE, type = "l",
        xlab = "", ylab = "", log = "y", yaxs = "i", yaxt = "n", xaxs = "i")
    # p-value axis drawn inside the panel (red), so that it does not collide with the colour key above
    axis(side = 3, at = pretty(c(0, 1)), line = -1.1, lwd = 0, lwd.ticks = 1, col.axis = "red", col.ticks = "red",
        cex.axis = 0.7, mgp = c(3, 0.2, 0))
    abline(v = p$pval_abline, col = "grey", lwd = 2, lty = 2)
  }
  if (p$add_MTM_peaks == TRUE & is.na(p$MTM_res_2[1, 2]) == FALSE) {
    abline(h = 1 / p$MTM_res_2[, 1], col = "black", lty = 3)
  }
  if (is.null(p$add_abline_h) != TRUE) {
    abline(h = p$add_abline_h)
  }
}

# MTM panels

.pw_mtm_peaks_red <- function(p, horizontal) {
  if (p$add_MTM_peaks == TRUE & is.na(p$MTM_res_2[1, 2]) == FALSE) {
    if (horizontal) {
      abline(h = 1 / p$MTM_res_2[, 1], col = "red", lty = 3)
    } else {
      abline(v = 1 / p$MTM_res_2[, 1], col = "red", lty = 3)
    }
  }
}

# vertical plots: MTM power, AR(1) and harmonic confidence levels as stacked panels (period on the x axis)
# kind: "power", "ar", "har"; ylab_in_title: label drawn with title() instead of the plot call
.pw_mtm_v <- function(p, kind = c("power", "ar", "har"), ylab_in_title = FALSE) {
  kind <- match.arg(kind)
  r <- p$MTM_res_1
  per <- 1 / r[, 1]
  if (kind == "power") {
    plot(x = per, y = r[, 2], type = "l", xaxs = "i", xaxt = "n", ylab = "MTM power", xlab = "", log = "x",
        xlim = p$ylim_vals)
    for (i in 5:8) {
      lines(x = per, y = r[, i], xlim = p$ylim_vals, lty = 3, col = "grey", lwd = 2)
    }
  } else if (kind == "ar") {
    plot(x = per, y = r[, 4], type = "l", xaxs = "i", xaxt = "n", ylab = if (ylab_in_title) { "" } else { "Ar. conf lvl" },
        xlab = "", log = "x", ylim = c(80, 101), xaxs = "i", xlim = p$ylim_vals)
    if (ylab_in_title) {
      title(ylab = "Ar. conf lvl", xpd = NA)
    }
    abline(h = c(90, 95, 99), lty = 3, col = "grey", lwd = 2)
  } else {
    plot(x = per, y = r[, 3], type = "l", xaxs = "i", xaxt = "n", ylab = if (ylab_in_title) { "" } else { "Har. conf lvl" },
        xlab = "", log = "xy", ylim = c(80, 101), yaxs = "i", xlim = p$ylim_vals)
    if (ylab_in_title) {
      title(ylab = "Har. conf lvl", xpd = NA)
    }
    abline(h = c(90, 95, 99), lty = 3, col = "grey", lwd = 2)
  }
  .pw_mtm_peaks_red(p, horizontal = FALSE)
}

# horizontal plots: the same panels side by side (period on the y axis)
.pw_mtm_h <- function(p, kind = c("power", "ar", "har")) {
  kind <- match.arg(kind)
  r <- p$MTM_res_1
  per <- 1 / r[, 1]
  if (kind == "power") {
    plot(y = per, x = r[, 2], type = "l", yaxs = "i", yaxt = "n", xlab = "MTM power", ylab = "", log = "y",
        ylim = p$ylim_vals)
    for (i in 5:8) {
      lines(y = per, x = r[, i], ylim = p$ylim_vals, lty = 3, col = "grey", lwd = 2)
    }
  } else if (kind == "ar") {
    plot(y = per, x = r[, 4], type = "l", yaxs = "i", yaxt = "n", xlab = "Ar. conf lvl", ylab = "", log = "y",
        xlim = c(80, 101), xaxs = "i", ylim = p$ylim_vals)
    abline(v = c(90, 95, 99), lty = 3, col = "grey", lwd = 2)
  } else {
    plot(y = per, x = r[, 3], type = "l", yaxs = "i", yaxt = "n", xlab = "Har. conf lvl", ylab = "", log = "xy",
        xlim = c(80, 101), xaxs = "i", ylim = p$ylim_vals)
    abline(v = c(90, 95, 99), lty = 3, col = "grey", lwd = 2)
  }
  .pw_mtm_peaks_red(p, horizontal = TRUE)
}

# The 16 layouts

# draws the complete figure; layouts differ in which panels exist and in their margins
.wr_fix_key_labels <- function(labels, levels, key.marks) {
  # fixed decimals hide small values (0.001, 0.001, 0.000): use significant digits when labels are not distinct
  lev <- as.numeric(levels)[key.marks + 1]
  lev[lev <= 0] <- 0
  if (anyDuplicated(labels) > 0 || any(as.numeric(labels) == 0 & lev != 0)) {
    return(ifelse(lev == 0, "0", formatC(lev, digits = 2, format = "g")))
  }
  labels
}

.pw_draw <- function(p) {
  hz <- p$plot_horizontal
  d <- p$add_data
  a <- p$add_avg
  m <- p$add_MTM
  if (is.null(p$levels)) {
    p$levels <- .pw_levels(p)
  }
  if (!is.null(p$key.labels) && !is.null(p$key.marks)) {
    p$key.labels <- .wr_fix_key_labels(p$key.labels, p$levels, p$key.marks)
  }
  V <- function(...) par(...)
  if (!hz) {
    if (d == FALSE & a == FALSE & m == FALSE) {
      .pw_layout(c(1, 2), 2, 1, c(1, 4), c(1))
      V(mar = c(2, 4, 2, 3))
      .pw_key_bottom(p, "side2")
      box(lwd = p$lwd.axis)
      V(new = FALSE, mar = c(4, 4, 2, 3))
      .pw_scalogram_v(p, ylab = p$x_lab, yaxt_n = FALSE)
    } else if (d == FALSE & a == FALSE & m == TRUE) {
      .pw_layout(c(2, 1, 3, 0, 4, 0, 5, 0), 4, 2, c(1, 1, 1, 4), c(4, 1))
      V(mar = c(2, 2, 2, 3))
      .pw_key_bottom(p, "side2")
      box(lwd = p$lwd.axis)
      V(mar = c(0, 4, 2, 2))
      .pw_mtm_v(p, "power")
      V(mar = c(0, 4, 0, 2))
      .pw_mtm_v(p, "ar")
      .pw_mtm_v(p, "har")
      V(new = FALSE, mar = c(4, 4, 0, 2))
      .pw_scalogram_v(p, ylab = p$x_lab, yaxt_n = FALSE)
    } else if (d == FALSE & a == TRUE & m == FALSE) {
      .pw_layout(c(1, 2, 0, 3), 2, 2, c(1, 4), c(1, 4))
      V(mar = c(2, 2, 2, 3))
      .pw_key_bottom(p, "below3")
      box(lwd = p$lwd.axis)
      V(mar = c(0, 4, 2, 2))
      .pw_gws_v(p, abline_v = TRUE)
      V(new = FALSE, mar = c(4, 4, 0, 2))
      .pw_scalogram_v(p, ylab = p$x_lab, yaxt_n = FALSE)
    } else if (d == FALSE & a == TRUE & m == TRUE) {
      .pw_layout(c(2, 1, 3, 0, 4, 0, 5, 0, 6, 0), 5, 2, c(1, 1, 1, 1, 4), c(4, 1))
      V(mar = c(2, 2, 2, 3))
      .pw_key_bottom(p, "side2")
      box(lwd = p$lwd.axis)
      V(mar = c(0, 4, 2, 2))
      .pw_mtm_v(p, "power")
      V(mar = c(0, 4, 0, 2))
      .pw_mtm_v(p, "ar")
      .pw_mtm_v(p, "har")
      .pw_gws_v(p)
      V(new = FALSE, mar = c(4, 4, 0, 2))
      .pw_scalogram_v(p, ylab = p$x_lab, yaxt_n = FALSE)
    } else if (d == TRUE & a == FALSE & m == FALSE) {
      .pw_layout(c(1, 3, 0, 2), 2, 2, c(1, 0.25), c(1, 4))
      V(mar = c(4, 4, 2, 0))
      .pw_data_v(p, abline_h = TRUE)
      V(new = FALSE, mar = c(3, 0, 2, 2))
      .pw_key_bottom(p, "side2")
      box(lwd = p$lwd.axis)
      V(new = FALSE, mar = c(4, 0, 2, 2))
      .pw_scalogram_v(p, ylab = "")
    } else if (d == TRUE & a == FALSE & m == TRUE) {
      .pw_layout(c(1, 2, 0, 3, 0, 4, 5, 6), 4, 2, c(1, 1, 1, 4), c(1, 4))
      V(mar = c(0, 1, 2, 6), xpd = NA)
      .pw_key_bottom(p, "below4")
      box(lwd = p$lwd.axis)
      V(mar = c(0, 0, 2, 2))
      .pw_mtm_v(p, "power")
      V(mar = c(0, 0, 0, 2), xpd = FALSE)
      .pw_mtm_v(p, "ar", TRUE)
      .pw_mtm_v(p, "har", TRUE)
      V(mar = c(4, 4, 0, 0))
      .pw_data_v(p)
      V(new = FALSE, mar = c(4, 0, 0, 2))
      .pw_scalogram_v(p, ylab = NULL)
    } else if (d == TRUE & a == TRUE & m == FALSE) {
      .pw_layout(c(1, 2, 3, 4), 2, 2, c(1, 4), c(1, 4))
      V(mar = c(2, 0.5, 2, 6), xpd = NA)
      .pw_key_left(p)
      box(lwd = p$lwd.axis)
      V(mar = c(0, 0, 2, 2))
      .pw_gws_v(p, xlab_empty = TRUE, title_ylab = TRUE, pval_xpd = TRUE)
      V(mar = c(4, 4, 0, 0), xpd = TRUE)
      .pw_data_v(p)
      V(new = FALSE, mar = c(4, 0, 0, 2), xpd = FALSE)
      .pw_scalogram_v(p, ylab = "")
    } else {
      .pw_layout(c(1, 2, 0, 3, 0, 4, 0, 5, 6, 7), 5, 2, c(1, 1, 1, 1, 4), c(1, 4))
      V(mar = c(0, 1, 2, 6), xpd = FALSE)
      .pw_key_bottom(p, "below4")
      box(lwd = p$lwd.axis)
      V(mar = c(0, 0, 2, 2), xpd = FALSE)
      .pw_mtm_v(p, "power")
      V(mar = c(0, 0, 0, 2), xpd = FALSE)
      .pw_mtm_v(p, "ar", TRUE)
      .pw_mtm_v(p, "har", TRUE)
      .pw_gws_v(p, title_ylab = TRUE)
      V(mar = c(4, 4, 0, 0), xpd = FALSE)
      .pw_data_v(p)
      V(new = FALSE, mar = c(4, 0, 0, 2))
      .pw_scalogram_v(p, ylab = NULL)
    }
    .pw_overlay(p, FALSE)
    .pw_tail(p, FALSE)
  } else {
    if (d == FALSE & a == FALSE & m == FALSE) {
      .pw_layout(c(2, 1), 1, 2, c(1), c(10, 2.25))
      V(mar = c(4, 0, 2, 5))
      .pw_key_right(p)
      box(lwd = p$lwd.axis)
      V(new = FALSE, mar = c(4, 4, 2, 0.5))
      .pw_scalogram_h(p)
      .pw_overlay(p, TRUE)
      .pw_tail(p, TRUE)
    } else if (d == FALSE & a == FALSE & m == TRUE) {
      .pw_layout(c(5, 4, 3, 2, 1), 1, 5, c(1), c(7, 1, 1, 1, 1))
      V(mar = c(4, 0, 2, 5))
      .pw_key_right(p)
      box(lwd = p$lwd.axis)
      V(new = FALSE, mar = c(4, 0, 2, 0.5))
      .pw_mtm_h(p, "har")
      .pw_mtm_h(p, "ar")
      .pw_mtm_h(p, "power")
      V(new = FALSE, mar = c(4, 4, 2, 0), mgp = c(2, 1, 0))
      .pw_scalogram_h(p)
      .pw_overlay(p, TRUE)
      .pw_tail(p, TRUE)
    } else if (d == FALSE & a == TRUE & m == FALSE) {
      .pw_layout(c(3, 2, 1), 1, 3, c(1), c(6, 1, 1))
      V(mar = c(4, 0, 2, 5))
      .pw_key_right(p)
      box(lwd = p$lwd.axis)
      V(new = FALSE, mar = c(4, 0, 2, 0.5))
      .pw_gws_h(p)
      V(new = FALSE, mar = c(4, 4, 2, 0), mgp = c(2, 1, 0))
      .pw_scalogram_h(p)
      .pw_overlay(p, TRUE)
      .pw_tail(p, TRUE)
    } else if (d == FALSE & a == TRUE & m == TRUE) {
      .pw_layout(c(6, 5, 4, 3, 2, 1), 1, 6, c(1), c(6, 1, 1, 1, 1, 1))
      V(mar = c(4, 0, 2, 5))
      .pw_key_right(p)
      box(lwd = p$lwd.axis)
      V(new = FALSE, mar = c(4, 0, 2, 0.5))
      .pw_mtm_h(p, "har")
      .pw_mtm_h(p, "ar")
      .pw_mtm_h(p, "power")
      .pw_gws_h(p)
      V(new = FALSE, mar = c(4, 4, 2, 0), mgp = c(2, 1, 0))
      .pw_scalogram_h(p)
      .pw_overlay(p, TRUE)
      .pw_tail(p, TRUE)
    } else if (d == TRUE & a == FALSE & m == FALSE) {
      .pw_layout(c(1, 0, 3, 2), 2, 2, c(0.25, 1), c(8, 2))
      V(mar = c(0, 4, 2, 2))
      .pw_data_h(p)
      V(new = FALSE, mar = c(4, 0, 2, 5))
      .pw_key_right(p)
      box(lwd = p$lwd.axis)
      V(new = FALSE, mar = c(4, 4, 0, 2))
      .pw_scalogram_h(p)
      .pw_overlay(p, TRUE)
      .pw_tail(p, TRUE)
    } else if (d == TRUE & a == FALSE & m == TRUE) {
      .pw_layout(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 5, 6), 2, 9, c(0.25, 1), c(1))
      V(mar = c(0, 4, 2, 0))
      .pw_data_h(p)
      V(new = FALSE, mar = c(3, 2, 2, 2), mgp = c(2, 1, 0))
      .pw_key_top(p)
      box(lwd = p$lwd.axis)
      V(new = FALSE, mar = c(4, 4, 0, 0), mgp = c(2, 1, 0))
      .pw_scalogram_h(p)
      .pw_overlay(p, TRUE)
      .pw_tail(p, TRUE)
      V(new = FALSE, mar = c(4, 0, 0, 0.5))
      .pw_mtm_h(p, "power")
      .pw_mtm_h(p, "ar")
      .pw_mtm_h(p, "har")
    } else if (d == TRUE & a == TRUE & m == FALSE) {
      .pw_layout(c(1, 2, 4, 3), 2, 2, c(0.25, 1), c(8, 2))
      V(mar = c(0, 4, 2, 0))
      .pw_data_h(p)
      V(new = FALSE, mar = c(3, 2, 2, 2), mgp = c(2, 1, 0))
      .pw_key_top(p)
      box(lwd = p$lwd.axis)
      V(new = FALSE, mar = c(4, 0, 0, 0.5))
      .pw_gws_h(p)
      V(new = FALSE, mar = c(4, 4, 0, 0), xpd = FALSE)
      .pw_scalogram_h(p)
      .pw_overlay(p, TRUE)
      .pw_tail(p, TRUE)
    } else {
      .pw_layout(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7), 2, 10, c(0.25, 1), c(1))
      V(mar = c(0, 4, 2, 0))
      .pw_data_h(p)
      V(new = FALSE, mar = c(3, 2, 2, 2), mgp = c(2, 1, 0))
      .pw_key_top(p)
      box(lwd = p$lwd.axis)
      V(new = FALSE, mar = c(4, 4, 0, 0), mgp = c(2, 1, 0))
      .pw_scalogram_h(p)
      .pw_overlay(p, TRUE)
      .pw_tail(p, TRUE)
      V(new = FALSE, mar = c(4, 0, 0, 0.5))
      .pw_gws_h(p)
      .pw_mtm_h(p, "power")
      .pw_mtm_h(p, "ar")
      .pw_mtm_h(p, "har")
    }
  }
  .pw_arrows(p)
  invisible(NULL)
}

# phase arrows of the cross spectra on top of the scalogram (maxima of power above a quantile)
.pw_arrows <- function(p) {
  a <- p$arrows
  if (is.null(a)) {
    return(invisible(NULL))
  }
  Power <- a$Power
  Phase <- a$Phase
  maxdetect <- matrix(0, nrow = a$nr, ncol = a$nc)
  for (j in 1:a$nc) {
    for (i in 2:(a$nr - 1)) {
      if ((Power[i, j] > Power[i + 1, j]) & (Power[i, j] > Power[i - 1, j])) {
        maxdetect[i, j] <- 1
      }
    }
  }
  # arrow size: T_unit / P_unit (depth and log2-period units) when given; otherwise a fixed visual length
  # (0.3 inch) whose slope on the screen equals the phase angle, whatever the record length
  T_unit <- a$T_unit
  P_unit <- a$P_unit
  if (is.null(T_unit) || is.null(P_unit)) {
    g <- .pw_env$geom
    len <- 0.3
    if (g$horizontal) {
      T_unit <- len * diff(g$usr[1:2]) / g$pin[1]
      P_unit <- len * diff(g$usr[3:4]) / g$pin[2]
    } else {
      T_unit <- len * diff(g$usr[3:4]) / g$pin[2]
      P_unit <- len * diff(g$usr[1:2]) / g$pin[1]
    }
  }
  t_idx <- seq(1, a$nc, by = max(1, round(a$nc / a$n_arrows)))
  s_idx <- 1:a$nr
  for (i in t_idx) {
    for (j in s_idx) {
      if (Power[j, i] >= a$p_thresh && maxdetect[j, i] == 1) {
        x0 <- a$depth[i]
        y0 <- a$ax2[j]
        angle <- Phase[j, i]
        x1 <- x0 + cos(angle) * T_unit
        y1 <- y0 + sin(angle) * P_unit
        if (p$plot_horizontal == TRUE) {
          arrows(x0, y0, x1, y1, length = a$length, lwd = a$lwd, angle = a$angle, col = "black")
        } else {
          arrows(y0, x0, y1, x1, length = a$length, lwd = a$lwd, angle = a$angle, col = "black")
        }
      }
    }
  }
}
