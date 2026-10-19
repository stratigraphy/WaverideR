#' @title Plot a cvCOCO result
#'
#' @description Plots the diagnostics of a \code{\link{cvcoco}} result in the
#' style of Acycle's plotcvcoco.m, using base graphics. Four pages can be
#' drawn: (1) the two regularized folds and the held-out validation spectra
#' with their frozen targets, (2) the stacked result panels: directional and
#' consensus correlation, global p, local p and the number of resolvable
#' astronomical periods, (3) the consensus pCOCO curve and (4) the Monte
#' Carlo audit (null distributions of the consensus and directional maxima).
#' Translated from Acycle v3.1 plotcvcoco.m (Mingsong Li).
#'
#' @param result Output of \code{\link{cvcoco}}.
#' @param pages Pages to draw, any of "spectra", "result", "pcoco",
#' "audit". \code{Default=c("spectra", "result", "pcoco", "audit")}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li.
#'
#' @references
#' Li, M., Kump, L.R., Hinnov, L.A., Mann, M.E., 2018. Tracking variable
#' sedimentation rates and astronomical forcing in Phanerozoic
#' paleoclimate proxy series with evolutionary correlation coefficients and
#' hypothesis testing. Earth and Planetary Science Letters 501, 165-179.
#'
#' @examples
#' \donttest{
#' set.seed(3)
#' t_kyr <- seq(0, 1600, by = 2)
#' etp <- sin(2 * pi * t_kyr / 405) + 0.6 * sin(2 * pi * t_kyr / 125) +
#'   0.8 * sin(2 * pi * t_kyr / 41) + 0.5 * sin(2 * pi * t_kyr / 23.7)
#' dat <- cbind(t_kyr * 5 / 100, etp + 0.5 * stats::rnorm(length(t_kyr)))
#' res <- cvcoco(data = dat, sr1 = 2, sr2 = 10, srstep = 0.5, n_sim = 20,
#'   seed = 1)
#' plot_cvcoco(res, pages = "result")
#' }
#'
#' @return Invisibly returns \code{NULL}; called for its plots.
#'
#' @export
#' @importFrom graphics par plot lines points abline legend axis mtext
#' @importFrom graphics hist text rect box title
#' @importFrom graphics plot.new
#' @importFrom grDevices adjustcolor
plot_cvcoco <- function(result = NULL,
                        pages = c("spectra", "result", "pcoco", "audit")){
  # ---- internal helper functions (local to this function) ----

  # @keywords internal
  # @noRd
  cvcoco_curve_max <- function(curve, sr_grid){
    curve <- as.numeric(curve)
    finite_index <- which(is.finite(curve))
    if (length(finite_index) == 0){
      return(list(score = NA_real_, rate = NA_real_, index = NA_integer_))
    }
    relative <- which.max(curve[finite_index])
    index <- finite_index[relative]
    list(score = curve[index], rate = sr_grid[index], index = index)
  }

  # @keywords internal
  # @noRd
  cvcoco_mark <- function(rate, score, col){
    if (length(rate) == 1 && length(score) == 1 && is.finite(rate) &&
        is.finite(score)){
      graphics::points(rate, score, pch = 19, col = col, cex = 0.7)
    }
    invisible(NULL)
  }

  # @keywords internal
  # @noRd
  cvcoco_format_p <- function(value){
    if (length(value) != 1 || !is.finite(value)){
      return("NaN")
    }
    if (value == 0){
      return("0.000")
    }
    exponent <- floor(log10(abs(value)))
    digits <- max(0, 4 - exponent - 1)
    sprintf(paste0("%0.", digits, "f"), value)
  }

  # @keywords internal
  # @noRd
  cvcoco_plot_spectrum <- function(diagnostic, validation, train_label,
                                   validation_label, fmax, target_legend,
                                   show_legend){
    main <- sprintf("%s train -> %s held out at %.4g cm/kyr\nrho = %.4g",
                    train_label, validation_label, validation$best_rate,
                    validation$score)
    if (is.finite(validation$participating_period_count)){
      main <- paste0(main, sprintf("; %d/9 periods",
                                   round(validation$participating_period_count)))
    }
    if (length(diagnostic$frequency) == 0){
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Spectrum unavailable for the selected rate grid")
      graphics::title(main, cex.main = 0.9)
      return(invisible(NULL))
    }
    f <- diagnostic$frequency
    data_power <- cvcoco_normalize_visible(f, diagnostic$data_power, fmax)
    target_power <- cvcoco_normalize_visible(f, diagnostic$target_power, fmax)
    graphics::plot(f, data_power, type = "l", xlim = c(0, fmax),
                   ylim = c(0, 1.05), xlab = "Frequency (cycle/kyr)",
                   ylab = "Normalized power", main = main, cex.main = 0.9)
    graphics::lines(f, target_power, col = "red")
    if (show_legend == TRUE){
      graphics::legend("topright", legend = c("Validation data", target_legend),
                       col = c("black", "red"), lty = 1, bty = "n", cex = 0.75)
    }
    invisible(NULL)
  }

  # @keywords internal
  # @noRd
  cvcoco_normalize_visible <- function(x, y, fmax){
    ok <- is.finite(x) & is.finite(y) & x >= 0 & x <= fmax
    scale <- NA
    if (any(ok)){
      scale <- max(y[ok])
    }
    if (!is.finite(scale) || scale <= 0){
      scale <- max(c(y[is.finite(y)], -Inf))
    }
    if (!is.finite(scale) || scale <= 0){
      scale <- 1
    }
    y / scale
  }

  # @keywords internal
  # @noRd
  cvcoco_plot_p_panel <- function(sr, xl, p_BtoA, p_AtoB, p_cons, n_mc, title,
                                  best_BtoA, best_AtoB, best_cons, dir_labels,
                                  thresholds){
    if (n_mc <= 0 || (!any(is.finite(p_AtoB)) && !any(is.finite(p_BtoA)) &&
                      !any(is.finite(p_cons)))){
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Monte Carlo was not run")
      graphics::title(title)
      return(invisible(NULL))
    }
    p_floor <- 1 / (n_mc + 1)
    s_BtoA <- cvcoco_p_to_score(p_BtoA, p_floor)
    s_AtoB <- cvcoco_p_to_score(p_AtoB, p_floor)
    s_cons <- cvcoco_p_to_score(p_cons, p_floor)
    all_s <- c(s_BtoA, s_AtoB, s_cons)
    top <- max(c(all_s[is.finite(all_s)], -log10(p_floor), 0.5))
    graphics::plot(sr, s_BtoA, type = "l", col = "blue", xlim = xl,
                   ylim = c(-log10(1.01), top * 1.12), yaxt = "n", xlab = "",
                   ylab = title, main = title)
    graphics::lines(sr, s_AtoB, col = "red")
    for (th in thresholds){
      if (th >= p_floor){
        lty <- 3
        if (th == 0.05 || (title == "Local p" && th == 0.01)){
          lty <- 2
        }
        graphics::abline(h = -log10(th), lty = lty, col = "grey30")
      }
    }
    graphics::lines(sr, s_cons, col = "black", lwd = 2)
    cvcoco_mark(best_BtoA[1], -log10(max(best_BtoA[2], p_floor)), "blue")
    cvcoco_mark(best_AtoB[1], -log10(max(best_AtoB[2], p_floor)), "red")
    cvcoco_mark(best_cons[1], -log10(max(best_cons[2], p_floor)), "black")
    ticks <- c(1, 0.1, 0.05, 0.02, 0.01, 0.002, 0.001)
    ticks <- ticks[ticks >= p_floor & -log10(ticks) <= top * 1.12]
    graphics::axis(2, at = -log10(ticks), labels = format(ticks), las = 1,
                   cex.axis = 0.8)
    graphics::legend("topright",
                     legend = c(paste0("p_", dir_labels[1], "=",
                                       cvcoco_format_p(best_BtoA[2])),
                                paste0("p_", dir_labels[2], "=",
                                       cvcoco_format_p(best_AtoB[2])),
                                paste0("p_cons=", cvcoco_format_p(best_cons[2]))),
                     col = c("blue", "red", "black"), lwd = c(1, 1, 2),
                     bty = "n", cex = 0.8)
    invisible(NULL)
  }

  # @keywords internal
  # @noRd
  cvcoco_p_to_score <- function(p, p_floor){
    s <- rep(NA_real_, length(p))
    ok <- is.finite(p) & p > 0
    s[ok] <- -log10(pmax(p[ok], p_floor))
    s
  }

  # @keywords internal
  # @noRd
  cvcoco_plot_null <- function(null_scores, observed, p_value, title, xlab){
    null_scores <- null_scores[is.finite(null_scores)]
    if (length(null_scores) == 0){
      graphics::plot.new()
      graphics::text(0.5, 0.5, sprintf("Monte Carlo was not run\nObserved statistic = %.4g",
                                       observed))
      graphics::title(title)
      return(invisible(NULL))
    }
    h <- graphics::hist(null_scores, plot = FALSE)
    prop <- h$counts / sum(h$counts)
    xr <- range(c(h$breaks, observed[is.finite(observed)]))
    graphics::plot(NA, xlim = xr, ylim = c(0, max(prop) * 1.1), xlab = xlab,
                   ylab = "Proportion",
                   main = sprintf("%s\np=%.4g, N=%d", title, p_value,
                                  length(null_scores)))
    graphics::rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], prop,
                   col = grDevices::adjustcolor("#A6B3C7"), border = NA)
    if (is.finite(observed)){
      graphics::abline(v = observed, col = "red", lwd = 2)
    }
    graphics::legend("topright", legend = c("Stationary AR(1) null",
                                            "Observed statistic"),
                     fill = c("#A6B3C7", NA), border = c(NA, NA),
                     col = c(NA, "red"), lwd = c(NA, 2), bty = "n", cex = 0.8)
    invisible(NULL)
  }

  # ---- main code ----
  if (is.null(result) || is.null(result$sr_grid)){
    stop("result must be the output of cvcoco()")
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  sr <- result$sr_grid
  name <- result$name
  if (isTRUE(result$degraded_mode)){
    name <- paste0(name, " (partial-orbit exploratory)")
  }
  lab_A <- result$fold_labels[1]
  lab_B <- result$fold_labels[2]
  dir_A <- "A"
  dir_B <- "B"
  if (lab_A != "Segment A"){
    dir_A <- lab_A
    dir_B <- lab_B
  }
  target_legend <- "Frozen coherent 9-period target"
  if (result$target_model == "four-group"){
    target_legend <- "Frozen phase-averaged four-group target"
  }
  fmax <- result$config$max_frequency
  n_mc <- max(c(result$n_sim_valid_AtoB, result$n_sim_valid_BtoA,
                result$n_sim_valid_consensus, 0))
  xl <- range(sr)
  if (xl[1] == xl[2]){
    xl <- xl + c(-0.5, 0.5)
  }

  if ("spectra" %in% pages){
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3.5, 1), oma = c(0, 0, 2, 0))
    graphics::plot(result$data_A[, 1], result$data_A[, 2], type = "l",
                   xlab = "Depth (m)", ylab = "Proxy value",
                   main = paste(lab_A, "depth series"))
    graphics::plot(result$data_B[, 1], result$data_B[, 2], type = "l",
                   xlab = "Depth (m)", ylab = "Proxy value",
                   main = paste(lab_B, "depth series"))
    cvcoco_plot_spectrum(result$spectra$validate_BtoA, result$validate_BtoA,
                         dir_B, dir_A, fmax, target_legend, TRUE)
    cvcoco_plot_spectrum(result$spectra$validate_AtoB, result$validate_AtoB,
                         dir_A, dir_B, fmax, target_legend, FALSE)
    graphics::mtext(paste(name, "data and held-out validation spectra"),
                    outer = TRUE, font = 2)
  }

  if ("result" %in% pages){
    graphics::par(mfrow = c(4, 1), mar = c(2.5, 4.5, 2, 1), oma = c(2, 0, 2, 0))
    curves <- c(result$validate_BtoA$curve, result$validate_AtoB$curve,
                result$consensus$curve)
    yl <- range(curves[is.finite(curves)])
    if (!all(is.finite(yl))){
      yl <- c(-1, 1)
    }
    graphics::plot(sr, result$validate_BtoA$curve, type = "l", col = "blue",
                   xlim = xl, ylim = yl, xlab = "", ylab = expression(rho),
                   main = "Correlation coefficient")
    graphics::lines(sr, result$validate_AtoB$curve, col = "red")
    graphics::lines(sr, result$consensus$curve, col = "black", lwd = 2)
    cvcoco_mark(result$validate_BtoA$best_rate, result$validate_BtoA$score,
                "blue")
    cvcoco_mark(result$validate_AtoB$best_rate, result$validate_AtoB$score,
                "red")
    cvcoco_mark(result$consensus$best_rate, result$consensus$best_correlation,
                "black")
    graphics::legend("topright", legend = c(dir_A, dir_B, "Consensus"),
                     col = c("blue", "red", "black"), lwd = c(1, 1, 2),
                     bty = "n", cex = 0.8)
    cvcoco_plot_p_panel(sr, xl, result$p_curve_BtoA, result$p_curve_AtoB,
                        result$p_curve_consensus, n_mc, "Global p",
                        c(result$validate_BtoA$best_rate, result$p_A),
                        c(result$validate_AtoB$best_rate, result$p_B),
                        c(result$consensus$best_rate, result$p_consensus),
                        c(dir_A, dir_B), c(0.5, 0.1, 0.05, 0.02, 0.01, 0.001))
    cvcoco_plot_p_panel(sr, xl, result$p_local_curve_BtoA,
                        result$p_local_curve_AtoB,
                        result$p_local_curve_consensus, n_mc, "Local p",
                        c(result$validate_BtoA$best_rate,
                          result$validate_BtoA$p_local_at_best),
                        c(result$validate_AtoB$best_rate,
                          result$validate_AtoB$p_local_at_best),
                        c(result$consensus$best_rate,
                          result$consensus$p_local_at_best),
                        c(dir_A, dir_B), c(0.1, 0.05, 0.02, 0.01, 0.001))
    graphics::plot(sr, result$orbit_count_A, type = "n", xlim = xl,
                   ylim = c(0, 9.5), xlab = "", ylab = "#",
                   main = "Number of contributing astronomical parameters")
    rng9 <- result$all_nine_rate_range_shared
    if (all(is.finite(rng9)) && rng9[1] < rng9[2]){
      usr <- graphics::par("usr")
      if (rng9[1] > usr[1]){
        graphics::rect(usr[1], usr[3], rng9[1], usr[4], col = "grey90",
                       border = NA)
      }
      if (rng9[2] < usr[2]){
        graphics::rect(rng9[2], usr[3], usr[2], usr[4], col = "grey90",
                       border = NA)
      }
    }
    graphics::lines(sr, result$orbit_count_A, col = "blue")
    graphics::lines(sr, result$orbit_count_B, col = "red", lty = 2)
    graphics::legend("bottomright", legend = c(lab_A, lab_B),
                     col = c("blue", "red"), lty = c(1, 2), bty = "n",
                     cex = 0.8)
    graphics::box()
    graphics::mtext("Sedimentation rate (cm/kyr)", side = 1, line = 0.5,
                    outer = TRUE)
    graphics::mtext(paste(name, "bidirectional held-out result"),
                    outer = TRUE, font = 2)
  }

  if ("pcoco" %in% pages){
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 0, 0))
    pc <- result$pcoco
    if (!any(is.finite(pc))){
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Monte Carlo was not run\nConsensus pCOCO is unavailable")
      graphics::title(paste(name, ": pCOCO"))
    } else {
      graphics::plot(sr, pc, type = "l", col = "red", lwd = 2, xlim = xl,
                     xlab = "Sedimentation rate (cm/kyr)", ylab = "pCOCO",
                     main = paste(name, ": pCOCO"))
      best <- cvcoco_curve_max(pc, sr)
      if (is.finite(best$rate)){
        graphics::abline(v = best$rate, col = "red", lty = 2)
        graphics::points(best$rate, best$score, pch = 19, col = "red",
                         cex = 0.8)
        graphics::mtext(sprintf("%.4g cm/kyr", best$rate), side = 3,
                        line = -1.5, at = best$rate, col = "red", cex = 0.8)
      }
    }
  }

  if ("audit" %in% pages){
    graphics::par(mfrow = c(3, 1), mar = c(4.5, 4.5, 3.5, 1), oma = c(0, 0, 2, 0))
    cvcoco_plot_null(result$null_consensus, result$score_consensus,
                     result$p_consensus, "Consensus global p",
                     "Null consensus maximum")
    cvcoco_plot_null(result$null_BtoA, result$validate_BtoA$score, result$p_A,
                     sprintf("%s held out; global p; best %.4g cm/kyr", dir_A,
                             result$validate_BtoA$best_rate),
                     "Null directional maximum")
    cvcoco_plot_null(result$null_AtoB, result$validate_AtoB$score, result$p_B,
                     sprintf("%s held out; global p; best %.4g cm/kyr", dir_B,
                             result$validate_AtoB$best_rate),
                     "Null directional maximum")
    graphics::mtext(paste(name, "Monte Carlo audit"), outer = TRUE, font = 2)
  }
  invisible(NULL)
}
