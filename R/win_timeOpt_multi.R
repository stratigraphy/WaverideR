#' @title Windowed timeOpt with precession, short-eccentricity and obliquity interpretations
#'
#' @description Runs \code{\link{win_timeOpt}} with fit 1 (precession amplitude modulation),
#' fit 2 (short-eccentricity amplitude modulation) and fit 3 (obliquity amplitude modulation)
#' on the same record and compares their best sedimentation rates window by window.
#'
#' @details A record that holds only one cycle family can be read in two ways: obliquity cycles
#' read as precession give a peak at about twice the true rate (fit 1), precession cycles read
#' as obliquity give a peak at about half the true rate (fit 3). For every window the function
#' flags
#' \describe{
#'   \item{agree_12, agree_13, agree_23}{both fits are significant and their best rates agree within
#'   \code{ratio_tol}: independent support for the same rate}
#'   \item{ambiguous_13}{both fit 1 and fit 3 are significant and the fit-1 rate is about twice the
#'   fit-3 rate: precession/obliquity 2:1 ambiguity, which needs other evidence (eccentricity
#'   cycles, an age tie point) to resolve}
#'   \item{ambiguous_23}{both fit 2 and fit 3 are significant and the fit-2 rate is about 0.41 times
#'   the fit-3 rate: obliquity (41 kyr with its 173-kyr modulation) read as short eccentricity
#'   (~100 kyr with the 405-kyr modulation)}
#' }
#'
#' @param data Record to analyse (matrix or data frame with two columns: depth in metres in column 1, proxy value in
#' column 2), evenly spaced; passed to \code{\link{win_timeOpt}}. Must be given (the \code{NULL} default fails).
#' \code{Default=NULL}.
#' @param fits The timeOpt fits to run (numeric vector with any of 1, 2, 3): 1 precession amplitude modulation, 2
#' short-eccentricity amplitude modulation, 3 obliquity amplitude modulation (the \code{fit} argument of
#' \code{\link{win_timeOpt}}). The agreement and ambiguity flags need the corresponding pairs. \code{Default=c(1, 2, 3)}.
#' @param ratio_tol Relative tolerance used to compare the best rates of two fits (numeric, > 0; 0.15 counts rates
#' within a factor 1.15 as equal); used for the \code{agree_*} and \code{ambiguous_*} flags. \code{Default=0.15}.
#' @param ... Further arguments passed to \code{\link{win_timeOpt}} for every fit (e.g. \code{window_size},
#' \code{sedmin}, \code{sedmax}, \code{numsed}, \code{n_null}); \code{fit} is not allowed (use \code{fits}).
#'
#' @return A list of class \code{win_timeOpt_multi.result} with the \code{\link{win_timeOpt}}
#' result of every fit (\code{fit1}, \code{fit2}, \code{fit3}) and \code{summary}, a data frame
#' with per window the best rate and significance of each fit and the flags above.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' mag_multi <- win_timeOpt_multi(data = mag[1:250, ], fits = c(1, 2), window_size = 6,
#'     sedmin = 0.1, sedmax = 1, numsed = 5, n_null = 6, noise_model = "ar1", step = 25)
#' head(mag_multi$summary)
#' }
#' @export
win_timeOpt_multi <- function(data = NULL, fits = c(1, 2, 3), ratio_tol = 0.15, ...) {
  dots <- list(...)
  if (!is.null(dots$fit)) {
    stop("win_timeOpt_multi: use 'fits', not 'fit'", call. = FALSE)
  }
  out <- list()
  for (f in fits) {
    args <- c(list(data = data, fit = f), dots)
    out[[paste0("fit", f)]] <- do.call(win_timeOpt, args)
  }
  first <- out[[1]]
  sm <- data.frame(x = first$x)
  for (f in fits) {
    r <- out[[paste0("fit", f)]]
    sm[[paste0("best_fit", f)]] <- r$best_sedrate
    sm[[paste0("z_fit", f)]] <- r$best_z
    if (is.null(r$significant)) {
      sm[[paste0("sig_fit", f)]] <- rep(NA, nrow(sm))
    } else {
      sm[[paste0("sig_fit", f)]] <- r$significant
    }
  }
  near <- function(a, b) {
    is.finite(a) & is.finite(b) & abs(log(a / b)) < log(1 + ratio_tol)
  }
  both <- function(a, b) {
    sa <- sm[[paste0("sig_fit", a)]]
    sb <- sm[[paste0("sig_fit", b)]]
    !is.na(sa) & !is.na(sb) & sa & sb
  }
  pairs <- list(c(1, 2), c(1, 3), c(2, 3))
  for (p in pairs) {
    if (all(p %in% fits)) {
      sm[[paste0("agree_", p[1], p[2])]] <- both(p[1], p[2]) &
          near(sm[[paste0("best_fit", p[1])]], sm[[paste0("best_fit", p[2])]])
    }
  }
  if (all(c(1, 3) %in% fits)) {
    sm$ambiguous_13 <- both(1, 3) & near(sm$best_fit1, 2 * sm$best_fit3)
  }
  if (all(c(2, 3) %in% fits)) {
    # obliquity (41 kyr, 173-kyr modulation) read as short eccentricity (~100 kyr, 405-kyr
    # modulation) at about 0.41 times the obliquity rate
    sm$ambiguous_23 <- both(2, 3) & near(sm$best_fit2, 0.41 * sm$best_fit3)
  }
  out$summary <- sm
  out$fits <- fits
  class(out) <- "win_timeOpt_multi.result"
  invisible(out)
}

#' @title Plot the fits of win_timeOpt_multi side by side
#'
#' @description Plots the (bias-corrected) z_opt of every fit of
#' \code{\link{win_timeOpt_multi}} as depth x sedimentation rate panels (untestable cells
#' greyed out, best rates marked), the average curves of all fits on one rate axis, and the
#' best rates per window of all fits with the 2:1 precession/obliquity ambiguity marked.
#'
#' @param multi_result Result of \code{\link{win_timeOpt_multi}} (list of class \code{win_timeOpt_multi.result}).
#' Must be given (the \code{NULL} default fails). \code{Default=NULL}.
#' @param corrected What is plotted (logical, TRUE or FALSE): \code{TRUE} the bias-corrected z-scores
#' (\code{z_opt}; without them, i.e. with \code{bias_correction = "none"}, the raw values are used),
#' \code{FALSE} the raw r-squared values (\code{r_2_opt}). \code{Default=TRUE}.
#' @param n.levels Number of colour levels of the colour scale (integer, > 1). \code{Default=100}.
#' @param xlab Label of the depth axis (character). \code{Default="depth (m)"}.
#' @param ylab Label of the sedimentation rate axis (character). \code{Default="sedrate (cm/kyr)"}.
#'
#' @return A plot.
#' @examples
#' \donttest{
#' set.seed(1)
#' mag_multi <- win_timeOpt_multi(data = mag[1:250, ], fits = c(1, 2), window_size = 6,
#'     sedmin = 0.1, sedmax = 1, numsed = 5, n_null = 6, noise_model = "ar1", step = 25)
#' plot_win_timeOpt_multi(mag_multi)
#' }
#' @export
#' @importFrom graphics par image axis mtext box layout abline lines points plot.new legend title rect
#' @importFrom grDevices rainbow adjustcolor
plot_win_timeOpt_multi <- function(multi_result = NULL, corrected = TRUE, n.levels = 100,
    xlab = "depth (m)", ylab = "sedrate (cm/kyr)") {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  fits <- multi_result$fits
  nf <- length(fits)
  fit_lab <- c("fit 1: precession AM", "fit 2: short-ecc. AM", "fit 3: obliquity AM")
  fit_col <- c("#2a78d6", "#eb6834", "#1baf7a")
  first <- multi_result[[paste0("fit", fits[1])]]
  sedrate <- first$sedrate
  log_y <- isTRUE(first$linLog == 1)
  yv <- if (log_y) { log10(sedrate) } else { sedrate }
  depth <- first$x
  regular_x <- length(depth) > 2 && all(abs(diff(diff(depth))) < 1e-6 * abs(diff(depth)[1]))
  tick <- unique(signif(sedrate[round(seq(1, length(sedrate), length.out = 8))], 2))
  at_tick <- if (log_y) { log10(tick) } else { tick }
  # common colour scale
  allv <- NULL
  for (f in fits) {
    allv <- c(allv, as.vector(wr_win_timeopt_sel(multi_result[[paste0("fit", f)]], 3, corrected)$M))
  }
  brk <- wr_breaks(matrix(allv, ncol = 1), n.levels)
  key.cols <- .wr_hcl_cols(n.levels, "viridis")

  lay <- matrix(0, nf + 1, 2)
  for (i in seq_len(nf)) {
    lay[i, 1] <- i
    lay[i, 2] <- nf + 2
  }
  lay[nf + 1, 1] <- nf + 1
  lay[nf + 1, 2] <- nf + 3
  graphics::plot.new()
  graphics::layout(lay, widths = c(7, 2.4), heights = c(rep(1, nf), 0.9))
  for (i in seq_len(nf)) {
    r <- multi_result[[paste0("fit", fits[i])]]
    sel <- wr_win_timeopt_sel(r, 3, corrected)
    graphics::par(mar = c(1.5, 4.5, 1.8, 0.5))
    graphics::image(x = depth, y = yv, z = t(sel$M), col = key.cols, breaks = brk,
        useRaster = regular_x, yaxt = "n", xlab = "", ylab = ylab)
    if (!is.null(r$testable)) {
      U <- ifelse(r$testable, NA, 1)
      if (any(is.finite(U))) {
        graphics::image(x = depth, y = yv, z = t(U), col = grDevices::adjustcolor("grey30", 0.55),
            useRaster = regular_x, add = TRUE)
      }
    }
    graphics::axis(2, at = at_tick, labels = tick, las = 2)
    sig <- r$significant
    if (is.null(sig)) {
      sig <- rep(FALSE, length(depth))
    }
    by <- if (log_y) { log10(r$best_sedrate) } else { r$best_sedrate }
    graphics::points(depth[sig], by[sig], pch = 19, cex = 0.45)
    graphics::points(depth[!sig], by[!sig], pch = 1, cex = 0.45)
    graphics::title(main = fit_lab[fits[i]], font.main = 1, cex.main = 1, line = 0.4, adj = 0)
  }
  # best rates of all fits
  graphics::par(mar = c(4, 4.5, 1.5, 0.5))
  plot(NA, xlim = range(depth), ylim = range(yv), yaxt = "n", xlab = xlab, ylab = ylab, xaxs = "i")
  graphics::axis(2, at = at_tick, labels = tick, las = 2)
  sm <- multi_result$summary
  amb <- rep(FALSE, nrow(sm))
  if (!is.null(sm$ambiguous_13)) {
    amb <- amb | sm$ambiguous_13
  }
  if (!is.null(sm$ambiguous_23)) {
    amb <- amb | sm$ambiguous_23
  }
  if (any(amb)) {
    dd <- if (length(depth) > 1) { stats::median(diff(depth)) } else { 1 }
    for (k in which(amb)) {
      graphics::rect(depth[k] - dd / 2, min(yv), depth[k] + dd / 2, max(yv),
          col = grDevices::adjustcolor("#e34948", 0.18), border = NA)
    }
  }
  for (i in seq_len(nf)) {
    f <- fits[i]
    b <- sm[[paste0("best_fit", f)]]
    s <- sm[[paste0("sig_fit", f)]]
    s[is.na(s)] <- FALSE
    by <- if (log_y) { log10(b) } else { b }
    graphics::points(depth[s], by[s], pch = 19, cex = 0.6, col = fit_col[f])
  }
  graphics::title(main = "significant best rates (shaded: precession/obliquity or obliquity/eccentricity ambiguity)",
      font.main = 1, cex.main = 0.9, line = 0.3, adj = 0)
  # average curves
  graphics::par(mar = c(1.5, 0.5, 1.8, 4))
  xr <- NULL
  for (f in fits) {
    xr <- c(xr, wr_win_timeopt_sel(multi_result[[paste0("fit", f)]], 3, corrected)$avg)
  }
  xr <- range(xr[is.finite(xr)])
  plot(NA, xlim = xr, ylim = range(yv), yaxt = "n", xlab = "", ylab = "", yaxs = "i")
  graphics::axis(4, at = at_tick, labels = tick, las = 2)
  graphics::mtext(ylab, side = 4, line = 2.8, cex = 0.7)
  for (f in fits) {
    a <- wr_win_timeopt_sel(multi_result[[paste0("fit", f)]], 3, corrected)$avg
    a[!is.finite(a)] <- NA
    graphics::lines(a, yv, col = fit_col[f], lwd = 2)
  }
  graphics::legend("topright", legend = fit_lab[fits], col = fit_col[fits], lwd = 2, bty = "n",
      cex = 0.75)
  graphics::title(main = if (corrected) { "avg. z" } else { "avg. R^2" }, font.main = 1, cex.main = 0.9,
      line = 0.4)
  graphics::par(mar = c(4, 0.5, 1.5, 1))
  graphics::image(x = seq(0, n.levels), y = 1, z = t(matrix(brk, nrow = 1)), col = key.cols,
      breaks = brk, useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = if (corrected) { "z" } else { "R^2" },
      ylab = "")
  km <- round(seq(0, 1, length.out = 5) * n.levels)
  graphics::axis(1, at = km, labels = formatC(brk[km + 1], digits = 1, format = "f"), cex.axis = 0.8)
  graphics::box()
  invisible(NULL)
}
