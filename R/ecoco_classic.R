# 4. eCOCO (sliding window; 2019 Acycle driver)

#' @title Evolutionary COCO (Li et al. 2018)
#' @description COCO in a sliding window (2019 Acycle ecoco.m driver) with
#'   best-rate summary and ridge tracking per window.
#' @param data Record to analyse (matrix or data frame; the first two columns are used): depth in metres in column 1
#' and the proxy in column 2, evenly spaced (it is not interpolated; incomplete rows are dropped and the rows sorted
#' by depth). Required, no default.
#' @param target Astronomical target periodogram (two-column matrix or data frame): frequency in cycles/kyr in
#' column 1 and power in column 2, e.g. the output of \code{\link{coco_target_from_periods}} or
#' \code{\link{coco_target_from_series}}; passed to \code{\link{coco_classic}} for every window. Required, no
#' default.
#' @param orbit Astronomical periods in kyr (numeric vector), used to count how many astronomical parameters
#' contribute at each rate (and, with \code{adjust=1}, to define the target bands; then seven periods are needed).
#' \code{NULL} uses the periods stored with the target (attribute \code{"orbit"}), or else the seven periods of Li
#' et al. (2018) (400, 120.5, 93.1, 38.6, 23.0, 21.8 and 18.5 kyr). \code{Default=NULL}.
#' @param sr Sedimentation rates to test in every window, in cm/kyr (numeric vector, > 0).
#' \code{Default=seq(0.2, 20, by = 0.2)}.
#' @param dt Sampling interval of \code{data}, in metres (numeric, > 0); \code{NULL} uses the median spacing of the
#' depth column. It also sets the number of samples per window, \code{trunc(window / dt)}. \code{Default=NULL}.
#' @param pad Length (number of points) to which each window is zero-padded for the periodogram (integer, at least
#' the number of samples per window); a larger value gives a finer frequency grid. \code{Default=5000}.
#' @param red Red-noise treatment of the window periodograms, passed to \code{\link{coco_classic}} (integer):
#' \code{0} none; \code{1} the AR(1) background is subtracted (as in Li et al. 2018); \code{2} the periodogram is
#' divided by the AR(1) background minus 1; \code{3} a robust AR(1) background (median-smoothed spectrum) is
#' subtracted; a value from \code{50} to \code{99} divides by that AR(1) confidence level (in \%) minus 1. Negative
#' values are set to zero. \code{Default=1}.
#' @param nsim Number of Monte Carlo simulations of the null hypothesis of no astronomical forcing, used for the
#' p-value \code{p_h0} and pCOCO (integer, >= 0; 2000 in Li et al. 2018); \code{0} skips them and the score becomes
#' \code{norbit * rho}. See \code{mc} for where the null is simulated. \code{Default=0}.
#' @param slices Number of equal-depth slices of each window whose periodograms are averaged (integer, >= 1);
#' \code{1} uses the whole window. \code{Default=1}.
#' @param adjust Adjust the power of the target in four bands (long eccentricity, short eccentricity, obliquity,
#' precession) to the power of the data in those bands at each rate (integer): \code{0} no, \code{1} yes (only for a
#' seven-period target). \code{Default=0}.
#' @param window Length of the sliding window, in metres (numeric, > 0; it must contain at least 3 samples and at
#' most the whole record). Required, no default.
#' @param step Step between successive windows, in number of samples (integer, >= 1 and smaller than half the number
#' of samples). \code{Default=1}.
#' @param delinear Remove the mean of each window before the analysis, as in Acycle ecoco.m (numeric): \code{1} yes,
#' any other value no. \code{Default=1}.
#' @param method Correlation coefficient between the data and target spectra (character): \code{"Pearson"} (as in
#' Li et al. 2018) or \code{"Spearman"}. \code{Default=c("Pearson", "Spearman")}, i.e. \code{"Pearson"}.
#' @param mc Where the Monte Carlo null distribution is simulated (character): \code{"first-window"} reproduces
#' Acycle 2018-2019, the null is simulated once in the first window and reused for all windows;
#' \code{"each-window"} simulates a separate null for every window (slower). Only used when \code{nsim > 0}.
#' \code{Default=c("first-window", "each-window")}, i.e. \code{"first-window"}.
#' @param seed Seed of the random number generator (integer), so that the Monte Carlo results can be reproduced;
#' \code{NULL} leaves the current random state unchanged. \code{Default=NULL}.
#' @param track_ridge Track a smooth sedimentation-rate ridge through the score map with
#' \code{\link{ecoco_track_ridge}} (logical, TRUE or FALSE); \code{FALSE} returns the best rate per window
#' (\code{\link{ecoco_best_per_window}}) as \code{track}. \code{Default=TRUE}.
#' @param verbose Print progress messages for every window (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @return list of matrices (rows = sr, columns = windows): rho, p_h0, norbit,
#'   pcoco, score (= norbit/norb * pcoco) plus `track` (best rate per window)
#' @export
ecoco_classic <- function(data, target, orbit = NULL, window,
    step = 1, sr = seq(0.2, 20, by = 0.2), dt = NULL,
    pad = 5000, red = 1, nsim = 0, slices = 1, adjust = 0,
    delinear = 1, method = c("Pearson", "Spearman"),
    mc = c("first-window", "each-window"), seed = NULL,
    track_ridge = TRUE, verbose = TRUE){
  # Internal helper functions (local to this function)

  `%||%` <- function(a, b){
    if (is.null(a)){
      b
    } else {
      a
    }
  }

  # Monte Carlo p-value and pCOCO
  pcoco_value <- function(rho, p){
    ps <- p
    ps[!is.finite(ps) | ps <= 0] <- NA
    ps[ps > 1] <- 1
    rho * abs(log10(ps))
  }

  # Seven periods used in Li et al. (2018), 55-57 Ma La2004 target.
  orbit7_default <- function(){
    c(400, 1 / 0.0083, 1 / 0.01074, 1 / 0.02588, 1 / 0.04346, 1 / 0.0459, 1 / 0.0542)
  }

  # Main code
  method <- match.arg(method)
  mc <- match.arg(mc)
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (is.null(orbit)){
    orbit <- attr(target, "orbit") %||% orbit7_default()
  }
  dat <- as.matrix(data)[, 1:2]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  if (is.null(dt)){
    dt <- stats::median(diff(dat[, 1]))
  }
  nrow_d <- nrow(dat)
  npts <- trunc(window / dt)
  if (npts < 3 || npts > nrow_d){
    stop("Window must contain between 3 and N samples")
  }
  if (step >= nrow_d / 2){
    stop("Sliding step is too large")
  }
  m <- trunc((nrow_d - npts) / step) + 1
  starts <- seq(1, by = step, length.out = m)
  starts <- starts[starts + npts - 1 <= nrow_d]
  m <- length(starts)
  nsr <- length(sr)
  norb <- length(orbit)
  win <- function(i){
    ii <- starts[i]:(starts[i] + npts - 1)
    w <- dat[ii, , drop = FALSE]
    if (delinear == 1){
      w[, 2] <- w[, 2] - mean(w[, 2])
    }
    w
  }
  rho <- p_param <- p_h0 <- norbit <- matrix(NA_real_, nsr, m)
  null <- NULL
  if (nsim > 0 && mc == "first-window"){
    first <- coco_classic(win(1), target, orbit, sr, dt, pad, red, nsim, slices,
        adjust, method, keep_null = TRUE, verbose = verbose)
    null <- first$null
  }
  for (i in seq_len(m)){
    each <- nsim > 0 && mc == "each-window"
    if (each){
      nsim_i <- nsim
    } else {
      nsim_i <- 0
    }
    r <- coco_classic(win(i), target, orbit, sr, dt, pad, red,
        nsim_i, slices, adjust, method,
        verbose = FALSE)$result
    rho[, i] <- r$rho
    p_param[, i] <- r$p_param
    norbit[, i] <- r$norbit
    if (each){
      p_h0[, i] <- r$p_h0
    } else if (!is.null(null)){
      for (k in seq_len(nsr)){
        s <- null[k, ]
        s <- s[is.finite(s)]
        if (is.finite(rho[k, i]) && length(s)){
          pv <- (length(s) - sum(s < rho[k, i])) / length(s)   # findperct.m
          if (pv == 0){
            p_h0[k, i] <- 1 / (nsim + 1)
          } else {
            p_h0[k, i] <- pv
          }
        }
      }
    }
    if (verbose){
      message(sprintf("eCOCO window %d of %d", i, m))
    }
  }
  depth <- seq(dat[1, 1] + window / 2, dat[nrow_d, 1] - window / 2, length.out = m)
  pc <- pcoco_value(rho, p_h0)
  if (nsim > 0){
    score <- norbit / norb * pc
  } else {
    score <- norbit * rho
  }
  best <- ecoco_best_per_window(sr, depth, rho, p_h0, norbit, score)
  if (track_ridge){
    track <- ecoco_track_ridge(score, sr, depth, rho, p_h0, norbit, best)
  } else {
    track <- best
  }
  out <- list(sr = sr, depth = depth, rho = rho, p_param = p_param, p_h0 = p_h0,
      norbit = norbit, pcoco = pc, score = score, best = best, track = track,
      window = window, step = step, nsim = nsim, mc = mc, orbit = orbit)
  class(out) <- "ecoco_classic"
  out
}

#' @title Plot an eCOCO result (Li et al. 2018 method)
#' @description Maps of correlation, H0 significance, contributing parameters
#'   and window-normalised score with the tracked sedimentation-rate ridge.
#' @param x eCOCO result to plot: the output of \code{\link{ecoco_classic}} (an object of class
#' \code{ecoco_classic}). Required, no default.
#' @param ridge Draw the tracked sedimentation-rate ridge (red line) on the last map, the window-normalised score
#' (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param reverse_depth Draw depth increasing downwards, i.e. reverse the depth (y) axis of the maps (logical, TRUE or
#' FALSE). \code{Default=FALSE}.
#' @param ... Not used; present for compatibility with the generic \code{plot} function.
#' @return invisible NULL
#' @export
plot.ecoco_classic <- function(x, ridge = TRUE, reverse_depth = FALSE, ...){
  # Internal helper functions (local to this function)

  .ecoco_map <- function(sr, depth, z, main, zlim = NULL, col = .wr_hcl_cols(64, "viridis"),
      reverse_depth = FALSE, legend_lab = ""){
    if (is.null(zlim)){
      fz <- z[is.finite(z)]
      if (length(fz)){
        zlim <- range(fz)
      } else {
        zlim <- c(0, 1)
      }
    }
    if (zlim[1] == zlim[2]){
      zlim <- zlim + c(-0.5, 0.5)
    }
    if (reverse_depth){
      yl <- rev(range(depth))
    } else {
      yl <- range(depth)
    }
    if (length(depth) == 1){
      graphics::plot(sr, z[, 1], type = "l", main = main, xlab = "Sedimentation rate (cm/kyr)", ylab = legend_lab)
      return(invisible())
    }
    graphics::image(sr, depth, z, zlim = zlim, col = col, ylim = yl, useRaster = TRUE,
        xlab = "Sedimentation rate (cm/kyr)", ylab = "Depth (m)", main = main)
    graphics::box()
  }

  .ecoco_plot <- function(x, maps, ridge = TRUE, reverse_depth = FALSE){
    n <- length(maps)
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mfrow = c(1, n), mar = c(3.4, 3.6, 2, 0.6), mgp = c(2.2, 0.7, 0))
    for (i in seq_len(n)){
      mp <- maps[[i]]
      .ecoco_map(x$sr, x$depth, mp$z, mp$title, mp$zlim, reverse_depth = reverse_depth)
      if (!is.null(mp$contour) && length(x$depth) > 1 && any(is.finite(mp$z))){
        try(graphics::contour(x$sr, x$depth, mp$z, levels = mp$contour, add = TRUE,
            drawlabels = FALSE, col = "white", lwd = 0.6), silent = TRUE)
      }
      if (ridge && i == n && !is.null(x$track)){
        graphics::lines(x$track$sr, x$track$depth, col = "red", lwd = 1.5)
        graphics::points(x$track$sr, x$track$depth, col = "red", pch = 16, cex = 0.4)
      }
    }
  }

  .norm_by_window <- function(s){
    # min-max normalisation of every column (window) over its finite values
    out <- matrix(NA_real_, nrow(s), ncol(s))
    for (j in seq_len(ncol(s))){
      v <- s[, j]
      ok <- is.finite(v)
      oj <- rep(NA_real_, length(v))
      if (any(ok)){
        r <- range(v[ok])
        if (diff(r) > 0){
          oj[ok] <- (v[ok] - r[1]) / diff(r)
        } else {
          oj[ok] <- 0
        }
      }
      out[, j] <- oj
    }
    if (!is.null(colnames(s))){
      dimnames(out) <- list(NULL, colnames(s))
    }
    if (nrow(s) == 1){
      # a single rate gives one value per window, returned as a plain vector
      out <- out[1, ]
    }
    out
  }

  # Main code
  maps <- list(list(z = x$rho, title = "Correlation coefficient"))
  if (x$nsim > 0){
    maps[[length(maps) + 1]] <- list(z = -log10(x$p_h0), title = "H0 significance (-log10 p)", contour = 2)
  }
  maps[[length(maps) + 1]] <- list(z = x$norbit, title = "No. astronomical parameters", zlim = c(0, length(x$orbit)))
  maps[[length(maps) + 1]] <- list(z = .norm_by_window(x$score), title = "Score (window-normalised)", zlim = c(0, 1))
  .ecoco_plot(x, maps, ridge, reverse_depth)
}

#' @export
#' @noRd
print.ecoco_classic <- function(x, ...){
  cat(sprintf("eCOCO (Li et al. 2018): %d windows x %d rates, window %.4g m, %d MC sims (%s null)\n",
      length(x$depth), length(x$sr), x$window, x$nsim, x$mc))
  print(utils::head(x$track, 10))
  invisible(x)
}
