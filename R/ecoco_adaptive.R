# Adaptive eCOCO (ecoco.m adaptive branch + ecocoAdaptiveCore.m)

#' @title Adaptive eCOCO
#' @description Adaptive COCO in physical-depth sliding windows with per-window
#'   AR(1) Monte Carlo nulls and ridge tracking (Acycle ecoco.m adaptive branch,
#'   ecocoAdaptiveCore.m; Mingsong Li).
#' @param data Record to analyse (matrix or data frame, two columns): depth in metres in column 1 and the proxy in
#' column 2. Non-finite rows are dropped, duplicated depths averaged, and unevenly spaced data are linearly
#' interpolated to the median sample spacing (see \code{\link{coco_prepare_regular_data}}). Required, no default.
#' @param orbit9 The nine astronomical target periods in kyr (numeric vector of nine distinct values: long
#' eccentricity, four short eccentricity terms, obliquity, three precession terms), e.g. from
#' \code{\link{orbit9_at_age}} or \code{\link{calculate_orbit9}} for the age of the record. The default uses the
#' present-day periods (405.7, 130.7, 123.9, 98.9, 94.9, 41.0, 23.7, 22.4 and 19.0 kyr).
#' \code{Default=orbit9_default()}.
#' @param window Length of each sliding window, in metres (numeric, > 0, shorter than the record; every window must
#' contain at least 4 data points). Each window is resampled to an odd number of points around its centre and
#' analysed with Adaptive COCO. Required, no default.
#' @param sr Sedimentation rates to test in every window, in cm/kyr (numeric vector, > 0).
#' \code{Default=seq(0.5, 20, by = 0.1)}.
#' @param step Step between successive window centres, in samples of the (regularised) record (numeric, > 0); only
#' used when \code{step_depth=NULL}, giving \code{step_depth = step * sample interval}. \code{Default=1}.
#' @param step_depth Exact spacing of the window centres, in metres (numeric, > 0, or \code{NULL}); overrides
#' \code{step}. At most 10000 windows are allowed. \code{Default=NULL}.
#' @param pad Length (number of points) to which every window is zero-padded for the periodogram (integer, at least
#' the number of points per window); a larger value gives a finer frequency grid. \code{Default=5000}.
#' @param red Background (red-noise) removal applied to every window spectrum before the correlation (integer):
#' \code{0} no removal, \code{1} classical AR(1) background subtracted, \code{2} robust AR(1) background
#' (median-smoothed spectrum, Mann and Lees 1996 style), \code{3} smoothed-window-average (SWA) background (needs at
#' least 33 spectral bins). Negative power after subtraction is set to zero. \code{Default=1}.
#' @param nsim Number of AR(1) Monte Carlo simulations per window used for the global and local p-values (integer,
#' >= 0); the AR(1) coefficient and standard deviation are estimated in each window. \code{0} skips the
#' simulations, so the p-values are \code{NA}. \code{Default=500}.
#' @param method Correlation coefficient between the window spectrum and the target spectrum (character):
#' \code{"Pearson"} or \code{"Spearman"}. \code{Default=c("Pearson", "Spearman")}, i.e. \code{"Pearson"}.
#' @param max_frequency Highest frequency of the target spectrum that enters the correlation, in cycles/kyr
#' (numeric, must be at least the highest orbital frequency \code{max(1/orbit9)}); \code{NULL} uses
#' \code{1.2 * max(1/orbit9)}. Passed to \code{\link{coco_adaptive_evaluate}}. \code{Default=NULL}.
#' @param amplitude_mode How the amplitudes of the nine target sinusoids are set (character):
#' \code{"four-group-area"} fits one amplitude per group (long eccentricity, short eccentricity, obliquity,
#' precession) to the band areas of the data spectrum (Acycle default), \code{"per-orbit"} fits each of the nine
#' amplitudes to the spectral peak at its frequency, \code{"fixed"} uses fixed weights (1.0 eccentricity, 0.8
#' obliquity, 0.6 precession; see \code{\link{coco_fixed_target_weights}}). The target is always synthesised as
#' \code{"coherent-nine"}. \code{Default=c("four-group-area", "per-orbit", "fixed")}, i.e. \code{"four-group-area"}.
#' @param seed Seed of the random number generator (integer), so that the Monte Carlo results can be reproduced;
#' \code{NULL} leaves the current random state unchanged. \code{Default=1}.
#' @param center_limits Depths of the first and last window centre, in metres (numeric vector \code{c(first, last)},
#' or \code{NULL}); the windows must still fit inside the record. \code{NULL} uses
#' \code{c(first depth + window/2, last depth - window/2)}. \code{Default=NULL}.
#' @param batch Number of Monte Carlo simulations processed together in one batch (integer, >= 1); larger batches are
#' faster but use more memory. \code{NULL} chooses \code{min(50, nsim, 3e6 / (frequencies * windows))}.
#' \code{Default=NULL}.
#' @param verbose Print messages about the windows and the progress of the Monte Carlo simulations (logical, TRUE or
#' FALSE). \code{Default=TRUE}.
#' @return list with matrices (rate-by-window): rho, p_global, p_local, norbit,
#'   pcoco, score; `track` = ridge-tracked sedimentation-rate path
#' @param centre_limits British spelling of \code{center_limits}; if given (not \code{NULL}) it overrides
#' \code{center_limits}. \code{Default=NULL}.
#' @export
ecoco_adaptive <- function(data, orbit9 = orbit9_default(), window,
    sr = seq(0.5, 20, by = 0.1), step = 1, step_depth = NULL,
    pad = 5000, red = 1, nsim = 500,
    method = c("Pearson", "Spearman"), max_frequency = NULL,
    amplitude_mode = c("four-group-area", "per-orbit", "fixed"),
    seed = 1, center_limits = NULL, batch = NULL,
    verbose = TRUE,
    centre_limits = NULL){
  if (!is.null(centre_limits)) {
    center_limits <- centre_limits
  }  # British spelling synonym
  # Internal helper functions (local to this function)

  detrend_linear <- function(x){
    # Least-squares removal of a straight line from every column (MATLAB detrend(x,1)).
    xm <- as.matrix(x)
    n <- nrow(xm)
    if (n < 2){
      if (is.matrix(x)){
        return(xm * 0)
      } else {
        return(as.vector(xm) * 0)
      }
    }
    t <- seq_len(n)
    X <- cbind(1, t)
    q <- qr(X)
    res <- qr.resid(q, xm)
    if (is.matrix(x)){
      res
    } else {
      as.vector(res)
    }
  }

  # MATLAB periodogram (rectangular window, one-sided PSD)
  ml_periodogram <- function(x, nfft, fs){
    # x: vector or samples-by-series matrix. Every column is zero padded to nfft.
    xm <- as.matrix(x)
    n <- nrow(xm)
    if (nfft < n){
      stop("pad (nfft) must be >= number of samples")
    }
    if (nfft > n){
      xm <- rbind(xm, matrix(0, nfft - n, ncol(xm)))
    }
    X <- mvfft(xm)
    if (nfft %% 2 == 0){
      nf <- nfft / 2 + 1
    } else {
      nf <- (nfft + 1) / 2
    }
    P <- (Mod(X[seq_len(nf), , drop = FALSE])^2) / (fs * n)
    if (nf > 2){
      if (nfft %% 2 == 0){
        P[2:(nf - 1), ] <- 2 * P[2:(nf - 1), ]
      } else {
        P[2:nf, ] <- 2 * P[2:nf, ]
      }
    } else if (nf == 2 && nfft %% 2 == 1){
      P[2, ] <- 2 * P[2, ]
    }
    f <- (0:(nf - 1)) * fs / nfft
    if (is.matrix(x)){
      p_out <- P
    } else {
      p_out <- as.vector(P)
    }
    list(p = p_out, f = f)
  }

  # Interpolation
  interp1_na <- function(x, y, xi){
    # Linear interpolation, NA outside range (MATLAB interp1 default).
    o <- order(x)
    x <- x[o]
    y <- y[o]
    keep <- !duplicated(x)
    stats::approx(x[keep], y[keep], xout = xi, method = "linear", rule = 1,
        ties = "ordered")$y
  }

  interp1_extrap <- function(x, y, xi){
    # Linear interpolation with linear extrapolation (MATLAB 'linear','extrap').
    out <- interp1_na(x, y, xi)
    n <- length(x)
    if (n >= 2){
      lo <- xi < x[1]
      hi <- xi > x[n]
      if (any(lo)){
        out[lo] <- y[1] + (xi[lo] - x[1]) * (y[2] - y[1]) / (x[2] - x[1])
      }
      if (any(hi)){
        out[hi] <- y[n] + (xi[hi] - x[n]) * (y[n] - y[n - 1]) / (x[n] - x[n - 1])
      }
    }
    out
  }

  # Moving median with shrinking endpoints
  move_median <- function(x, w){
    x <- as.vector(x)
    n <- length(x)
    halfw <- floor(w / 2)
    if (w %% 2 == 1){
      back <- halfw
      ahead <- halfw
    } else {
      back <- max(0, halfw - 1)
      ahead <- halfw
    }
    out <- numeric(n)
    for (i in seq_len(n)){
      out[i] <- stats::median(x[max(1, i - back):min(n, i + ahead)])
    }
    out
  }

  # AR(1) red-noise backgrounds
  rho_ar1_ml <- function(x){
    # Mudelsee (2002) lag-1 estimator on demeaned data (Husson 2012 / Li 2017).
    x <- as.vector(x) - mean(x)
    n <- length(x)
    sum(x[2:n] * x[1:(n - 1)]) / sum(x[1:(n - 1)]^2)
  }

  theored_ar1_ml <- function(datax, f, meanp, dt){
    rho <- rho_ar1_ml(datax)
    fnyq <- 1 / (2 * dt)
    theored <- (1 - rho^2) / (1 - 2 * rho * cos(pi * f / fnyq) + rho^2)
    tmp <- theored
    tmp[1] <- 0
    art <- mean(tmp)
    theored * meanp / art
  }

  minirhos0 <- function(s0, fn, ft, pxxsmooth, linlog = 2){
    cospara <- cos(pi * as.vector(ft) / fn)
    power <- as.vector(pxxsmooth)
    dist_grid <- function(rho_g, s0_g){
      # AR(1) shape for every (frequency, rho) pair: rows = frequencies
      cmat <- matrix(cospara, length(cospara), length(rho_g))
      rmat <- matrix(rho_g, length(cospara), length(rho_g), byrow = TRUE)
      shape <- (1 - rmat^2) / (1 - 2 * cmat * rmat + rmat^2)
      if (linlog == 1){
        a <- colSums(shape^2)
        b <- as.vector(crossprod(shape, power))
        cc <- sum(power^2)
        outer(a, s0_g^2) - 2 * outer(b, s0_g) + cc
      } else {
        base <- log(shape) - log(power)
        rs <- colSums(base)
        rss <- colSums(base^2)
        ls <- log(s0_g)
        matrix(rss, length(rho_g), length(s0_g)) + 2 * outer(rs, ls) +
            length(power) * matrix(ls^2, length(rho_g), length(s0_g), byrow = TRUE)
      }
    }
    nn <- 50
    rhoi <- seq(0.001, 0.999, length.out = nn)
    s0i <- seq(0.2 * s0, 5 * s0, length.out = nn)
    d <- dist_grid(rhoi, s0i)
    k <- which.min(d)
    x <- (k - 1) %% nn + 1
    y <- (k - 1) %/% nn + 1
    mm <- nn / 2
    for (it in 1:3){
      ds <- s0i[2] - s0i[1]
      dr <- rhoi[2] - rhoi[1]
      s0max <- s0i[y] + ds
      s0min <- s0i[y] - ds
      rmax <- min(rhoi[x] + dr, 0.9999)
      rmin <- rhoi[x] - dr
      if (rmin <= 0){
        rmin <- 1e-4
      }
      rhoi <- seq(rmin, rmax, length.out = mm)
      s0i <- seq(s0min, s0max, length.out = mm)
      d <- dist_grid(rhoi, s0i)
      k <- which.min(d)
      x <- (k - 1) %% mm + 1
      y <- (k - 1) %/% mm + 1
    }
    list(rho = rhoi[x], s0 = s0i[y])
  }

  redconf_any <- function(f_angular, pxx, dt, smoothwin = 0.25, linlog = 2){
    # Robust AR(1) background (Mann & Lees 1996 style): median-smoothed spectrum,
    # grid-search fit of rho and S0. f_angular is radians/sample (2*pi*f*dt).
    fn <- 1 / (2 * dt)
    ft <- f_angular / pi * fn
    smoothn <- max(1, round(smoothwin * length(pxx)))
    pxxsmooth <- move_median(pxx, smoothn)
    pos <- pxxsmooth[pxxsmooth > 0]
    if (length(pos)){
      floor_val <- max(.Machine$double.xmin, max(pos) * .Machine$double.eps)
    } else {
      floor_val <- 1e-300
    }
    pxxsmooth <- pmax(pxxsmooth, floor_val)
    s0 <- mean(pxxsmooth)
    fit <- minirhos0(s0, fn, ft, pxxsmooth, linlog)
    fit$s0 * (1 - fit$rho^2) / (1 - 2 * fit$rho * cos(pi * ft / fn) + fit$rho^2)
  }

  # Monte Carlo p-value and pCOCO
  pcoco_value <- function(rho, p){
    ps <- p
    ps[!is.finite(ps) | ps <= 0] <- NA
    ps[ps > 1] <- 1
    rho * abs(log10(ps))
  }

  # Present-day periods; falls back to constants if the data files are not found.
  orbit9_default <- function(){
    tryCatch(orbit9_at_age(0),
        error = function(e){
          c(405.6912, 130.6979, 123.8532, 98.8517, 94.8856, 40.9897, 23.6820, 22.3758, 18.9519)
        })
  }

  # Batched spectra with background removal (ecocoWindowSpectra.m)
  coco_window_spectra <- function(values, dt, pad, red = 1, already_detrended = FALSE){
    values <- as.matrix(values)
    if (!already_detrended){
      values <- detrend_linear(values)
    }
    sp <- ml_periodogram(values, pad, 1 / dt)
    P <- as.matrix(sp$p)
    f <- sp$f
    P[P < 0] <- 0
    if (red == 0){
      return(list(power = P, f = f))
    }
    if (red == 3){
      return(coco_swa_subtract(P, f, nrow(values)))
    }
    for (j in seq_len(ncol(P))){
      p <- P[, j]
      if (red == 1){
        bg <- theored_ar1_ml(values[, j], f, mean(p), dt)
      } else {
        bg <- redconf_any(2 * pi * f * dt, p, dt, 0.25, 2)
      }
      p <- p - bg
      p[p < 0] <- 0
      P[, j] <- p
    }
    list(power = P, f = f)
  }

  # @title Remove the smoothed-window-average background (red = 3)
  # @description Helper of \code{coco_window_spectra}. Acycle red = 3: every
  # one-sided periodogram (at least 33 bins) is floored at
  # \code{max(realmin, max(positive power) * eps)}, converted to log10, passed
  # to specswa (\code{swa_background}), and the SWA background is subtracted
  # with negative values set to zero. Translated from Acycle v3.1
  # corrcoefslices_rankNew.m, redNoisePeriodogramMC.m and ecocoWindowSpectra.m
  # (Mingsong Li); verified in verification/ecoco_ext (red = 3 cases).
  # @keywords internal
  # @noRd
  coco_swa_subtract <- function(P, f, n_data){
    if (nrow(P) < 33){
      stop(paste0("red = 3 (smoothed-window average) requires at least 33 ",
          "one-sided spectral bins; increase pad"))
    }
    for (j in seq_len(ncol(P))){
      p <- P[, j]
      pos <- p[is.finite(p) & p > 0]
      if (length(pos) < 1){
        sc <- 1
      } else {
        sc <- max(pos)
      }
      fl <- max(2.2250738585072014e-308, sc * .Machine$double.eps)
      bg <- as.numeric(swa_background(f, log10(pmax(p, fl)), n_data,
          verbose = FALSE)$swa)
      if (length(bg) != length(p) || any(!is.finite(bg) | bg <= 0)){
        stop("The smoothed-window background is nonpositive or nonfinite.")
      }
      p <- p - bg
      p[!is.finite(p) | p < 0] <- 0
      P[, j] <- p
    }
    list(power = P, f = f)
  }

  # Physical-depth windows (ecocoPhysicalCenters.m / ecocoPhysicalWindowValues.m)
  ecoco_physical_centers <- function(depth, window, step_depth, center_limits = NULL){
    if (is.null(center_limits)){
      lim <- c(depth[1] + window / 2, depth[length(depth)] - window / 2)
    } else {
      lim <- center_limits
    }
    if (lim[2] < lim[1]){
      stop("Data span is shorter than the window")
    }
    ratio <- (lim[2] - lim[1]) / step_depth
    nint <- max(0, floor(ratio + 64 * .eps * max(1, abs(ratio))))
    if (nint + 1 > 10000){
      stop("More than 10000 windows; increase the step")
    }
    cen <- lim[1] + (0:nint) * step_depth
    cen[cen <= lim[2] + 64 * .eps * max(1, abs(lim), abs(cen))]
  }

  ecoco_physical_window_values <- function(data, window, centers, nominal_dt, min_source = 4){
    depth <- data[, 1]
    proxy <- data[, 2]
    h <- round(window / (2 * nominal_dt))
    if (h < 1){
      stop("Window too short for the sampling interval")
    }
    npt <- 2 * h + 1
    adt <- window / (2 * h)
    rel <- (-h:h) * adt
    rel[1] <- -window / 2
    rel[h + 1] <- 0
    rel[npt] <- window / 2
    V <- matrix(NA_real_, npt, length(centers))
    tol0 <- function(lo, hi){
      64 * .eps * max(1, abs(depth[1]), abs(depth[length(depth)]), abs(lo), abs(hi), window) *
          max(1, length(depth))
    }
    for (w in seq_along(centers)){
      lo <- centers[w] - window / 2
      hi <- centers[w] + window / 2
      tol <- tol0(lo, hi)
      if (lo < depth[1] - tol || hi > depth[length(depth)] + tol){
        stop(sprintf("Window %d exceeds the data support", w))
      }
      ii <- which(depth >= lo - tol & depth <= hi + tol)
      if (length(ii) < min_source){
        stop(sprintf("Window %d has fewer than %d observations", w, min_source))
      }
      q <- centers[w] + rel
      q[1] <- lo
      q[h + 1] <- centers[w]
      q[npt] <- hi
      V[, w] <- interp1_extrap(depth[ii], proxy[ii], q)
    }
    list(values = V, n = npt, dt = adt)
  }

  # Small helpers
  .eps <- .Machine$double.eps

  # Main code
  method <- match.arg(method)
  amplitude_mode <- match.arg(amplitude_mode)
  if (!is.null(seed)){
    set.seed(seed)
  }
  prep <- coco_prepare_regular_data(data, verbose)
  dat <- prep$data
  dt <- prep$dt
  if (is.null(step_depth)){
    step_depth <- step * dt
  }
  centers <- ecoco_physical_centers(dat[, 1], window, step_depth, center_limits)
  wv <- ecoco_physical_window_values(dat, window, centers, dt)
  V <- wv$values
  nW <- wv$n
  wdt <- wv$dt
  if (pad < nW){
    stop(sprintf("pad (%d) must be >= %d samples per window", pad, nW))
  }
  nWin <- length(centers)
  nRate <- length(sr)
  f <- (0:(floor(pad / 2))) / (pad * wdt)
  if (verbose){
    message(sprintf("Adaptive eCOCO: %d windows, %d points per window, %d rates", nWin, nW, nRate))
  }
  ev <- function(Pm){
    coco_adaptive_evaluate(Pm, f, nW, wdt, pad, orbit9, sr, method,
        max_frequency, "coherent-nine", amplitude_mode)
  }
  Pobs <- coco_window_spectra(V, wdt, pad, red)$power
  obs <- ev(Pobs)
  rho <- obs$rho
  # per-window AR(1) (ecocoEstimateAr1.m)
  X <- detrend_linear(V)
  rhoM <- sdv <- rep(NA_real_, nWin)
  validW <- logical(nWin)
  for (w in seq_len(nWin)){
    x <- X[, w]
    s <- sqrt(sum(x^2) / (nW - 1))
    rel <- sqrt(mean(x^2)) / max(abs(V[, w]))
    if (!is.finite(s) || s <= 0 || !(rel > 256 * .eps)){
      next
    }
    a <- x[-nW]
    b <- x[-1]
    rhoM[w] <- min(max(sum(a * b) / sum(a^2), -0.999), 0.999)
    sdv[w] <- s
    validW[w] <- TRUE
  }
  p_global <- p_local <- matrix(NA_real_, nRate, nWin)
  if (nsim > 0){
    nf <- floor(pad / 2) + 1
    if (is.null(batch)){
      batch <- max(1, min(50, nsim, floor(3e6 / (nf * nWin))))
    }
    gEx <- lEx <- matrix(0, nRate, nWin)
    done <- 0
    while (done < nsim){
      nb <- min(batch, nsim - done)
      rr <- rep(rhoM, nb)
      ss <- rep(sdv, nb)
      vv <- rep(validW, nb)
      Z <- matrix(stats::rnorm(nW * nWin * nb), nW)
      Y <- matrix(0, nW, nWin * nb)
      Y[1, vv] <- Z[1, vv]
      is <- sqrt(pmax(0, 1 - rr[vv]^2))
      for (r in 2:nW){
        Y[r, vv] <- rr[vv] * Y[r - 1, vv] + is * Z[r, vv]
      }
      Y[, vv] <- sweep(Y[, vv, drop = FALSE], 2, ss[vv], "*")
      Pn <- coco_window_spectra(Y, wdt, pad, red)$power
      nr <- ev(Pn)$rho
      dim(nr) <- c(nRate, nWin, nb)
      for (bi in seq_len(nb)){
        cur <- nr[, , bi, drop = FALSE]
        dim(cur) <- c(nRate, nWin)
        # maximum over rates of the finite null correlations, per window
        cmax <- rep(NA_real_, nWin)
        for (j in seq_len(nWin)){
          v <- cur[, j]
          if (any(is.finite(v))){
            cmax[j] <- max(v[is.finite(v)])
          }
        }
        gEx <- gEx + (sweep(rho, 2, cmax, "<=") %in% TRUE)
        lEx <- lEx + ((cur >= rho) %in% TRUE)
      }
      done <- done + nb
      if (verbose){
        message(sprintf("Adaptive eCOCO Monte Carlo %d of %d", done, nsim))
      }
    }
    p_global <- (gEx + 1) / (nsim + 1)
    p_local <- (lEx + 1) / (nsim + 1)
    p_global[!is.finite(rho)] <- NA
    p_local[!is.finite(rho)] <- NA
    p_global[, !validW] <- NA
    p_local[, !validW] <- NA
  }
  norbit <- matrix(length(orbit9) - obs$nMissing, nRate, nWin)
  pc <- pcoco_value(rho, p_global)
  score <- pc * norbit / length(orbit9)
  best <- ecoco_best_per_window(sr, centers, rho, p_global, norbit, score)
  track <- ecoco_track_ridge(score, sr, centers, rho, p_global, norbit, best)
  out <- list(sr = sr, depth = centers, rho = rho, p_global = p_global, p_local = p_local,
      norbit = norbit, pcoco = pc, score = score, best = best, track = track,
      window = window, step_depth = step_depth, window_points = nW,
      analysis_dt = wdt, rhoM = rhoM, nsim = nsim, orbit9 = orbit9,
      amplitude_mode = amplitude_mode, method = method)
  class(out) <- "ecoco_adaptive"
  out
}

#' @title Plot an adaptive eCOCO result
#' @description Maps of correlation, global and local p, contributing periods and
#'   normalised ridge score with the tracked sedimentation-rate ridge.
#' @param x Adaptive eCOCO result to plot: the output of \code{\link{ecoco_adaptive}} (an object of class
#' \code{ecoco_adaptive}). Required, no default.
#' @param ridge Draw the ridge-tracked sedimentation-rate path (red line) on the last map, the normalised ridge score
#' (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param reverse_depth Draw depth increasing downwards, i.e. reverse the depth (y) axis of the maps (logical, TRUE or
#' FALSE). \code{Default=FALSE}.
#' @param ... Not used; present for compatibility with the generic \code{plot} function.
#' @return invisible NULL
#' @export
plot.ecoco_adaptive <- function(x, ridge = TRUE, reverse_depth = FALSE, ...){
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
    maps[[length(maps) + 1]] <- list(z = -log10(x$p_global), title = "Global p (-log10)", contour = -log10(0.05))
    maps[[length(maps) + 1]] <- list(z = -log10(x$p_local), title = "Local p (-log10)", contour = 2)
  }
  maps[[length(maps) + 1]] <- list(z = x$norbit, title = "Contributing orbital periods", zlim = c(0, 9))
  maps[[length(maps) + 1]] <- list(z = .norm_by_window(x$score), title = "Ridge score (normalised)", zlim = c(0, 1))
  .ecoco_plot(x, maps, ridge, reverse_depth)
}

#' @export
#' @noRd
print.ecoco_adaptive <- function(x, ...){
  cat(sprintf("Adaptive eCOCO: %d windows x %d rates, window %.4g m (%d pts), %d MC sims\n",
      length(x$depth), length(x$sr), x$window, x$window_points, x$nsim))
  print(utils::head(x$track, 10))
  invisible(x)
}
