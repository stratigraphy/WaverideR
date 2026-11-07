#' @title plot the windowed timeOpt sedimentation rate estimation
#'
#' @description The \code{\link{plot_win_timeOpt}} function plots a windowed
#'  timeOpt sedimentation rate estimation of \code{\link{win_timeOpt}}.
#'@param win_timeOpt_result Result to plot (required): the list returned by \code{\link{win_timeOpt}}. \code{Default=NULL}.
#'@param proxy_name Label of the y-axis of the proxy panel on top (character, e.g. \code{"magnetic susceptibility"});
#'\code{NULL} leaves it empty. \code{Default=NULL}.
#'@param abline_h Sedimentation rates at which dashed horizontal lines are drawn on the map, in the units of
#'\code{win_timeOpt_result$sedrate} (numeric vector, e.g. \code{c(2, 3, 5)}, or \code{NULL} for none). \code{Default=NULL}.
#'@param abline_v Depths at which dashed vertical lines are drawn on the map, in the depth unit of the record (numeric
#'vector or \code{NULL} for none). \code{Default=NULL}.
#'@param add_lines Lines to draw on the map, e.g. a sedimentation rate curve (matrix or data frame, or \code{NULL} for
#'none): depth in column 1 and one or more columns of sedimentation rates after it; each column is drawn as a separate
#'line. \code{Default=NULL}.
#'@param fig_lts Text (e.g. a figure letter such as \code{"A"}) written in large type in the upper left corner of the
#'proxy panel (character, or \code{NULL} for none). \code{Default=NULL}.
#'@param xlab Label of the depth axis (character). \code{Default="depth (m)"}.
#'@param ylab Label of the sedimentation rate axis (character). \code{Default="sedrate (cm/kyr)"}.
#'@param sel_parameter Which timeOpt fit to plot (integer): \code{1} envelope fit, \code{2} spectral power fit,
#'\code{3} combined fit (opt). \code{Default=3}.
#'@param n.levels Number of colour levels of the map (integer, > 1); the breaks are quantiles of the plotted values.
#'\code{Default=100}.
#'@param corrected Plot the bias-corrected z-scores (\code{z_envelope}, \code{z_power}, \code{z_opt}) instead of the
#'raw r-squared values, when the result contains them (\code{bias_correction = "null"} in \code{\link{win_timeOpt}})
#'(logical). The median z threshold is then drawn as a dotted line on the average panel. \code{Default=TRUE}.
#'@param show_untestable Grey out the cells that cannot be tested (too few samples per cycle, too few cycles in the
#'window or window outside the record), as flagged in \code{testable} of the result (logical). \code{Default=TRUE}.
#'@param show_best Mark the best sedimentation rate of every window (logical); filled symbols are significant at the
#'family-wise level of \code{\link{win_timeOpt}}, open symbols are not. \code{Default=TRUE}.
#'
#' @author
#'Based on the \code{\link[astrochron:eTimeOpt]{eTimeOpt}}
#'function of the 'astrochron' R package.
#'
#'@references
#' Meyers, S.R. (2019). Cyclostratigraphy and the problem of astrochronologic testing.
#' Earth-Science Reviews 190, 190-223 <doi:10.1016/j.earscirev.2018.11.015>
#'
#'@examples
#'\donttest{
#'#plot the windowed timeOpt of the magnetic susceptibility record
#'#of the Sullivan core of Pas et al., (2018).
#'set.seed(1)
#'mag_win_timeOpt <-win_timeOpt(
#'data = mag[1:250, ],
#'window_size = 6,
#'sedmin = 0.1,
#'sedmax = 1,
#'numsed = 5,
#'limit = FALSE,
#'fit = 2,
#'fitModPwr = TRUE,
#'flow = NULL,
#'fhigh = NULL,
#'roll = 10 ^ 6,
#'targetE = c(405.7, 130.7, 123.8, 98.9, 94.9),
#'targetP = c(20.9, 19.9, 17.1, 17.2),
#'detrend = TRUE,
#'normalize =TRUE,
#'linLog = 1,
#'bias_correction = "null",
#'n_null = 6,
#'noise_model = "ar1",
#'step = 10,
#'run_multicore = FALSE,
#'verbose=FALSE)
#'
#'plot_win_timeOpt(win_timeOpt_result = mag_win_timeOpt,
#'proxy_name= "mag",
#'abline_h=NULL,
#'abline_v = NULL,
#'add_lines=NULL,
#'fig_lts = NULL,
#'xlab="depth (m)",
#'ylab= "sedrate (cm/kyr)",
#'sel_parameter=3,
#'n.levels=100,
#'corrected = TRUE)
#'}
#'
#' @return
#'The output is a plot of a windowed timeOpt with the average curve over all windows.
#'
#' @export
#' @importFrom stats quantile
#' @importFrom graphics par image axis mtext text box layout abline lines points plot.new legend
#' @importFrom grDevices rainbow adjustcolor

plot_win_timeOpt <- function(win_timeOpt_result = NULL,
    proxy_name = NULL,
    abline_h = NULL,
    abline_v = NULL,
    add_lines = NULL,
    fig_lts = NULL,
    xlab = "depth (m)",
    ylab = "sedrate (cm/kyr)",
    sel_parameter = 3,
    n.levels = 100,
    corrected = TRUE,
    show_untestable = TRUE,
    show_best = TRUE) {
  res <- win_timeOpt_result
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  wr_plot_win_timeopt_panels(res, proxy_name = proxy_name, abline_h = abline_h,
      abline_v = abline_v, add_lines = add_lines, fig_lts = fig_lts, xlab = xlab, ylab = ylab,
      sel_parameter = sel_parameter, n.levels = n.levels, corrected = corrected,
      show_untestable = show_untestable, show_best = show_best)
  invisible(NULL)
}

# matrix, average curve and labels of one result
wr_win_timeopt_sel <- function(res, sel_parameter, corrected) {
  sel_name <- c("envelope", "power", "opt")[sel_parameter]
  use_z <- isTRUE(corrected) && !is.null(res[[paste0("z_", sel_name)]])
  # no window could be tested against noise (all z values NA): show the uncorrected R^2 instead of failing
  untested <- use_z && !any(is.finite(as.matrix(res[[paste0("z_", sel_name)]])))
  if (untested) {
    return(list(M = res[[paste0("r_2_", sel_name)]], avg = res[[paste0("r_2_", sel_name, "_avg")]],
        key_lab = "R^2 (no window testable vs noise)", avg_lab = "avg. R^2", use_z = FALSE))
  }
  if (use_z) {
    list(M = res[[paste0("z_", sel_name)]], avg = res[[paste0("z_", sel_name, "_avg")]],
        key_lab = "z (vs noise)", avg_lab = "avg. z", use_z = TRUE)
  } else {
    list(M = res[[paste0("r_2_", sel_name)]], avg = res[[paste0("r_2_", sel_name, "_avg")]],
        key_lab = "R^2", avg_lab = "avg. R^2", use_z = FALSE)
  }
}

# strictly increasing colour breaks from the finite values of M
wr_breaks <- function(M, n.levels) {
  v <- M[is.finite(M)]
  if (length(v) == 0L) {
    v <- c(0, 1)
  }
  b <- stats::quantile(v, probs = seq(0, 1, length.out = n.levels + 1), names = FALSE)
  eps <- max(1e-12, 1e-9 * diff(range(b)))
  for (i in 2:length(b)) {
    if (b[i] <= b[i - 1]) {
      b[i] <- b[i - 1] + eps
    }
  }
  b
}

wr_plot_win_timeopt_panels <- function(res, proxy_name, abline_h, abline_v, add_lines, fig_lts,
    xlab, ylab, sel_parameter, n.levels, corrected, show_untestable, show_best) {
  sel <- wr_win_timeopt_sel(res, sel_parameter, corrected)
  M <- as.matrix(sel$M)
  sedrate <- res$sedrate
  depth <- res$x
  if (!is.null(res$data)) {
    proxy <- res$data
  } else {
    proxy <- cbind(res$x, res$y)
  }
  log_y <- isTRUE(res$linLog == 1)
  yv <- if (log_y) { log10(sedrate) } else { sedrate }
  brk <- wr_breaks(M, n.levels)
  key.cols <- .wr_hcl_cols(n.levels, "viridis")
  n.ticks <- 6
  key.marks <- round(seq(from = 0, to = 1, length.out = n.ticks) * n.levels)
  key.labels <- formatC(brk, digits = 2, format = "f")[key.marks + 1]
  xlim_vals <- range(proxy[, 1])

  graphics::plot.new()
  graphics::layout(mat = matrix(c(1, 2, 4, 3), nrow = 2, ncol = 2, byrow = TRUE),
      heights = c(0.25, 1), widths = c(8, 2))
  graphics::par(mar = c(0, 4, 2, 0))
  plot(y = proxy[, 2], x = proxy[, 1], type = "l", xaxs = "i", xlab = "", ylab = proxy_name,
      xaxt = "n", xlim = xlim_vals)
  if (!is.null(fig_lts)) {
    graphics::text(x = min(proxy[, 1]) + 0.05 * diff(xlim_vals), y = max(proxy[, 2]) * 0.9, fig_lts,
        col = "black", cex = 3)
  }
  graphics::par(mar = c(3, 2, 2, 2), mgp = c(2, 1, 0))
  graphics::image(x = seq(from = 0, to = n.levels), y = 1, z = t(matrix(brk, nrow = 1)),
      col = key.cols, breaks = brk, useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = sel$key_lab,
      ylab = "")
  graphics::axis(1, at = key.marks, labels = NA, tck = 0.02, tcl = 1.2)
  graphics::mtext(key.labels, side = 1, at = key.marks, line = 0.1, las = 2, cex = 0.75)
  graphics::box()
  graphics::par(mar = c(4, 0, 0, 0.5))
  avg <- sel$avg
  avg[!is.finite(avg)] <- NA
  plot(x = avg, y = yv, type = "l", yaxs = "i", yaxt = "n", xlab = sel$avg_lab, xaxs = "i",
      ylim = range(yv), xlim = if (any(is.finite(avg))) { range(avg, na.rm = TRUE) } else { c(0, 1) })
  if (sel$use_z && !is.null(res$z_threshold) && any(is.finite(res$z_threshold))) {
    graphics::abline(v = stats::median(res$z_threshold, na.rm = TRUE), lty = 3, col = "grey40")
  }
  graphics::par(mar = c(4, 4, 0, 0), xpd = FALSE)
  regular_x <- length(depth) > 2 && all(abs(diff(diff(depth))) < 1e-6 * abs(diff(depth)[1]))
  graphics::image(x = depth, y = yv, z = t(M), col = key.cols, breaks = brk,
      useRaster = regular_x, yaxt = "n", xlab = xlab, ylab = ylab, xlim = xlim_vals)
  if (show_untestable && !is.null(res$testable)) {
    U <- ifelse(res$testable, NA, 1)
    if (any(is.finite(U))) {
      graphics::image(x = depth, y = yv, z = t(U), col = grDevices::adjustcolor("grey30", 0.55),
          useRaster = regular_x, add = TRUE)
    }
  }
  tick <- unique(signif(sedrate[round(seq(1, length(sedrate), length.out = 8))], 2))
  if (log_y) {
    graphics::axis(2, at = log10(tick), labels = tick, las = 2)
  } else {
    graphics::axis(2, at = tick, labels = tick, las = 2)
  }
  if (show_best && !is.null(res$best_sedrate)) {
    by <- if (log_y) { log10(res$best_sedrate) } else { res$best_sedrate }
    sig <- res$significant
    if (is.null(sig)) {
      sig <- rep(FALSE, length(by))
    }
    graphics::points(depth[!sig], by[!sig], pch = 1, cex = 0.5, col = "black")
    graphics::points(depth[sig], by[sig], pch = 19, cex = 0.5, col = "black")
  }
  if (!is.null(abline_h)) {
    graphics::abline(h = if (log_y) { log10(abline_h) } else { abline_h }, lty = 2)
  }
  if (!is.null(abline_v)) {
    graphics::abline(v = abline_v, lty = 2)
  }
  if (!is.null(add_lines)) {
    for (i in 2:ncol(add_lines)) {
      graphics::lines(x = add_lines[, 1], y = if (log_y) { log10(add_lines[, i]) } else { add_lines[, i] })
    }
  }
  invisible(NULL)
}
