# 3. COCO (one record)

#' @title COCO analysis (Li et al. 2018)
#' @description Correlation coefficient between the (AR(1)-corrected) data
#'   periodogram and an astronomical target over a grid of sedimentation rates,
#'   with the H0 (no orbital forcing) Monte Carlo test of Li et al. (2018).
#'   Translated from Acycle corrcoefslices.m / cyclecorr.m (Mingsong Li).
#' @param data Record to analyse (matrix or data frame; the first two columns are used): depth in metres in column 1
#' and the proxy in column 2, evenly spaced (it is not interpolated; incomplete rows are dropped and the rows sorted
#' by depth). Required, no default.
#' @param target Astronomical target periodogram (two-column matrix or data frame): frequency in cycles/kyr in
#' column 1 and power in column 2, e.g. the output of \code{\link{coco_target_from_periods}},
#' \code{\link{coco_target_from_series}}. Its peaks higher than half the mean
#' power are rebuilt as Gaussian peaks with the Rayleigh width of the data at every tested rate. Required, no
#' default.
#' @param orbit Astronomical periods in kyr (numeric vector), used to count how many astronomical parameters
#' contribute at each rate (and, with \code{adjust=1}, to define the target bands; then seven periods are needed).
#' \code{NULL} uses the periods stored with the target (attribute \code{"orbit"}, set by
#' \code{\link{coco_target_from_periods}}), or else the seven periods of Li et al. (2018) (400, 120.5, 93.1, 38.6,
#' 23.0, 21.8 and 18.5 kyr). \code{Default=NULL}.
#' @param sr Sedimentation rates to test, in cm/kyr (numeric vector, > 0). \code{Default=seq(0.2, 20, by = 0.2)}.
#' @param dt Sampling interval of \code{data}, in metres (numeric, > 0); \code{NULL} uses the median spacing of the
#' depth column. \code{Default=NULL}.
#' @param pad Length (number of points) to which the data (or each slice) is zero-padded for the periodogram
#' (integer, at least the number of samples); a larger value gives a finer frequency grid. \code{Default=5000}.
#' @param red Red-noise treatment of the data periodogram (integer): \code{0} none; \code{1} the AR(1) background
#' is subtracted (as in Li et al. 2018); \code{2} the periodogram is divided by the AR(1) background minus 1;
#' \code{3} a robust AR(1) background (median-smoothed spectrum) is subtracted; a value from \code{50} to
#' \code{99} divides by that AR(1) confidence level (in \%) minus 1, so only power above it remains. Negative values
#' are set to zero. \code{Default=1}.
#' @param nsim Number of Monte Carlo simulations of the null hypothesis of no astronomical forcing (random
#' spectra with as many peaks as the data spectrum), used for the p-value \code{p_h0} and pCOCO (integer, >= 0;
#' 2000 in Li et al. 2018); \code{0} skips them. \code{Default=0}.
#' @param slices Number of equal-depth slices whose periodograms are averaged (integer, >= 1); \code{1} uses the
#' whole record. \code{Default=1}.
#' @param adjust Adjust the power of the target in four bands (long eccentricity, short eccentricity, obliquity,
#' precession) to the power of the data in those bands at each rate (integer): \code{0} no, \code{1} yes (only for a
#' seven-period target). \code{Default=0}.
#' @param method Correlation coefficient between the data and target spectra (character): \code{"Pearson"} (as in
#' Li et al. 2018) or \code{"Spearman"}. \code{Default=c("Pearson", "Spearman")}, i.e. \code{"Pearson"}.
#' @param seed Seed of the random number generator (integer), so that the Monte Carlo results can be reproduced;
#' \code{NULL} leaves the current random state unchanged. \code{Default=NULL}.
#' @param keep_null Return the Monte Carlo null correlations (matrix, rates in rows, simulations in columns) as
#' element \code{null} of the output (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param verbose Print progress messages of the Monte Carlo simulations (logical, TRUE or FALSE).
#' \code{Default=TRUE}.
#' @return list with data.frame `result` (sr, rho, p_param, p_h0, norbit, pcoco)
#' @export
coco_classic <- function(data, target, orbit = NULL,
    sr = seq(0.2, 20, by = 0.2), dt = NULL, pad = 5000,
    red = 1, nsim = 0, slices = 1, adjust = 0,
    method = c("Pearson", "Spearman"), seed = NULL,
    keep_null = FALSE, verbose = TRUE){
  # Internal helper functions (local to this function)

  `%||%` <- function(a, b){
    if (is.null(a)){
      b
    } else {
      a
    }
  }

  # Detrending
  detrend_mean <- function(x){
    if (is.matrix(x)){
      sweep(x, 2, colMeans(x))
    } else {
      x - mean(x)
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

  # Peak finder (strict local maxima, plateaus counted once)
  ml_findpeaks <- function(y){
    y <- as.vector(y)
    n <- length(y)
    if (n < 3){
      return(list(pks = numeric(0), locs = integer(0)))
    }
    # collapse plateaus
    r <- rle(y)
    vals <- r$values
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    k <- length(vals)
    if (k < 3){
      return(list(pks = numeric(0), locs = integer(0)))
    }
    idx <- which(vals[2:(k - 1)] > vals[1:(k - 2)] & vals[2:(k - 1)] > vals[3:k]) + 1
    idx <- idx[is.finite(vals[idx])]
    list(pks = vals[idx], locs = starts[idx])
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

  # Correlation helpers
  cor_with_p <- function(x, y, method = c("Pearson", "Spearman")){
    method <- match.arg(method)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    n <- length(x)
    if (n < 3){
      return(c(rho = NA_real_, p = NA_real_))
    }
    if (stats::sd(x) == 0 || stats::sd(y) == 0){
      return(c(rho = NA_real_, p = NA_real_))
    }
    if (method == "Pearson"){
      cm <- "pearson"
    } else {
      cm <- "spearman"
    }
    r <- stats::cor(x, y, method = cm)
    tt <- r * sqrt((n - 2) / max(1e-300, 1 - r^2))
    p <- 2 * stats::pt(-abs(tt), n - 2)
    c(rho = r, p = p)
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

  # Peaks of the target higher than half of the mean power (targetrebuilt.m)
  target_rebuilt <- function(target){
    psh <- 0.5 * mean(target[, 2])
    pk <- ml_findpeaks(target[, 2])
    keep <- pk$pks >= psh
    cbind(power = pk$pks[keep], freq = target[pk$locs[keep], 1])
  }

  # Gaussian-peak reconstruction of the target at a given Rayleigh width (freq2target.m)
  freq2target <- function(f, p, fx, fray){
    sigma <- 0.5 * fray
    px <- numeric(length(fx))
    for (i in seq_along(f)){
      g <- stats::dnorm(fx, f[i], sigma)
      if (max(g) > 0){
        px <- px + p[i] * g / max(g)
      }
    }
    px
  }

  # Number of orbital terms NOT involved at sed. rate sr (norbits.m)
  norbits_missing <- function(target, f, power, orbit, rayleigh, sr){
    aimarx <- f * sr / 100
    P_nyquist <- 1 / aimarx[length(aimarx)]
    P_dfr <- 1 / (rayleigh * sr / 100)
    no <- length(orbit)
    orbit8 <- numeric(no + 1)
    if (no > 1){
      for (i in 1:(no - 1)){
        orbit8[i + 1] <- 2 / (orbit[i + 1] + orbit[i])
      }
    }
    orbit8[no + 1] <- target[nrow(target), 1]
    nm <- 0
    for (i in 1:no){
      pow <- power[aimarx > orbit8[i] & aimarx <= orbit8[i + 1]]
      if (any(pow > 0)){
        po <- mean(pow[pow > 0])
      } else {
        po <- NaN
      }
      if (is.nan(po) || P_dfr < orbit[i] || P_nyquist >= orbit[i]){
        nm <- nm + 1
      }
    }
    nm
  }

  # Adjust target band powers to data (targetadj_real.m; 7-period target only)
  target_adjust_real <- function(target, f, power, orbit7, rayleigh, sr){
    orbit5 <- c(0, 2 / (orbit7[1] + orbit7[2]), 2 / (orbit7[3] + orbit7[4]),
        2 / (orbit7[4] + orbit7[6]), 0.06)
    aimarx <- f * sr / 100
    p4 <- numeric(4)
    for (i in 1:4){
      pow <- power[aimarx > orbit5[i] & aimarx <= orbit5[i + 1]]
      if (any(pow >= 1)){
        v <- mean(pow[pow >= 1])
      } else {
        v <- NaN
      }
      if (is.nan(v)){
        p4[i] <- 0
      } else {
        p4[i] <- v
      }
    }
    tx <- target[, 1]
    tp <- target[, 2]
    out <- numeric(length(tx))
    bands <- list(tx <= orbit5[2], tx > orbit5[2] & tx <= orbit5[3],
        tx > orbit5[3] & tx <= orbit5[4], tx > orbit5[4] & tx <= orbit5[5])
    for (b in 1:4){
      idx <- which(bands[[b]])
      if (length(idx) && is.finite(p4[b])){
        d <- max(tp[idx])
        if (is.finite(d) && d > 0){
          out[idx] <- p4[b] * tp[idx] / d
        }
      }
    }
    out
  }

  # Random half-sine-peak spectrum for the null hypothesis (randspec_sin.m)
  randspec_sin <- function(f, n, dat_ray){
    nf <- length(f)
    out <- numeric(nf)
    if (n == 0){
      return(out)
    }
    df <- f[2] - f[1]
    f1 <- seq(0, 2 * dat_ray, by = df)
    nf1 <- length(f1)
    p1 <- sin(2 * pi * (0.25 / dat_ray) * f1)
    if (nf - nf1 < 1){
      stop("Frequency grid too short for random spectrum")
    }
    for (i in 1:n){
      loc <- sample.int(nf - nf1, 1)
      out[loc:(loc + nf1 - 1)] <- out[loc:(loc + nf1 - 1)] + stats::runif(1) * p1
    }
    out
  }

  # 2. Data spectrum with AR(1) background removed
  coco_data_spectrum <- function(dat, dt, pad, red = 1, slices = 1){
    x <- dat[, 1]
    v <- dat[, 2]
    if (slices > 1){
      b <- seq(x[1], x[length(x)], length.out = slices + 1)
      # select_interval() in Acycle is inclusive at both ends
      idx <- vector("list", slices)
      for (j in 1:slices){
        idx[[j]] <- which(x >= b[j] & x <= b[j + 1])
      }
    } else {
      idx <- list(seq_along(x))
    }
    nslice <- min(lengths(idx))
    ps <- NULL
    for (j in seq_along(idx)){
      vj <- detrend_mean(v[idx[[j]]])
      sp <- ml_periodogram(vj, pad, 1 / dt)
      p <- sp$p
      f <- sp$f
      if (red == 1){
        p <- p - theored_ar1_ml(v, f, mean(p), dt)
        p[p < 0] <- 0
      } else if (red == 2){
        p <- p / theored_ar1_ml(v, f, mean(p), dt) - 1
        p[p < 0] <- 0
      } else if (red == 3){
        p <- p - redconf_any(2 * pi * f * dt, p, dt, 0.25, 2)
        p[p < 0] <- 0
      } else if (red >= 50 && red < 100){
        th <- theored_ar1_ml(v, f, mean(p), dt) * stats::qgamma(red / 100, 2) / 2
        p <- p / th - 1
        p[p < 0] <- 0
      }
      ps <- cbind(ps, p)
    }
    list(f = f, p = rowMeans(ps), p_slices = ps,
        rayleigh = 1 / (nslice * dt), nyquist = 1 / (2 * dt))
  }

  # Main code
  method <- match.arg(method)
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (is.null(orbit)){
    orbit <- attr(target, "orbit") %||% orbit7_default()
  }
  orbit <- sort(orbit, decreasing = TRUE)
  dat <- as.matrix(data)[, 1:2]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  if (is.null(dt)){
    dt <- stats::median(diff(dat[, 1]))
  }
  target <- as.matrix(target)
  spec <- coco_data_spectrum(dat, dt, pad, red, slices)
  xx <- spec$f
  yy <- spec$p
  dat_ray <- spec$rayleigh
  dat_nyq <- spec$nyquist
  tn <- target_rebuilt(target)
  lax <- target[, 1]
  sr0 <- lax[length(lax)] * 100 / dat_nyq
  npks <- length(ml_findpeaks(yy)$pks)
  if (method == "Pearson"){
    cm <- "pearson"
  } else {
    cm <- "spearman"
  }

  nsr <- length(sr)
  # Pre-build per-rate target on its evaluation grid (independent of data spectrum)
  prep <- vector("list", nsr)
  for (k in seq_len(nsr)){
    s <- sr[k]
    la <- freq2target(tn[, "freq"], tn[, "power"], lax, s * dat_ray / 100)
    prep[[k]] <- list(la = la, y = s * xx / 100, low = s < sr0)
  }
  eval_rate <- function(k, power){
    pr <- prep[[k]]
    la <- pr$la
    if (adjust == 1){
      la <- target_adjust_real(cbind(lax, la), xx, power, orbit, dat_ray, sr[k])
    }
    if (pr$low){
      a <- interp1_na(lax, la, pr$y)
      b <- power
    } else {
      a <- la
      b <- interp1_na(pr$y, power, lax)
    }
    list(a = a, b = b)
  }

  rho <- p_param <- nmi <- numeric(nsr)
  for (k in seq_len(nsr)){
    ab <- eval_rate(k, yy)
    cp <- cor_with_p(ab$a, ab$b, method)
    rho[k] <- cp["rho"]
    p_param[k] <- cp["p"]
    nmi[k] <- norbits_missing(cbind(lax, prep[[k]]$la), xx, yy, orbit, dat_ray, sr[k])
  }

  p_h0 <- rep(NA_real_, nsr)
  null <- NULL
  if (nsim > 0){
    null <- matrix(NA_real_, nsr, nsim)
    for (i in seq_len(nsim)){
      rs <- randspec_sin(xx, npks, dat_ray)
      for (k in seq_len(nsr)){
        ab <- eval_rate(k, rs)
        ok <- is.finite(ab$a) & is.finite(ab$b)
        null[k, i] <- suppressWarnings(stats::cor(ab$a[ok], ab$b[ok], method = cm))
      }
      if (verbose && i %% max(1, round(nsim / 10)) == 0){
        message(sprintf("COCO Monte Carlo %d of %d", i, nsim))
      }
    }
    for (k in seq_len(nsr)){
      s <- null[k, ]
      s <- s[is.finite(s)]
      if (is.finite(rho[k])){
        p_h0[k] <- (sum(s >= rho[k]) + 1) / (length(s) + 1)
      }
    }
  }
  res <- data.frame(sr = sr, rho = rho, p_param = p_param, p_h0 = p_h0,
      norbit = length(orbit) - nmi)
  res$pcoco <- pcoco_value(res$rho, res$p_h0)
  out <- list(result = res, spectrum = cbind(f = xx, power = yy), target = target,
      target_peaks = tn, sr0 = sr0, rayleigh = dat_ray, npeaks = npks,
      nsim = nsim, orbit = orbit, method = method, dt = dt)
  if (keep_null){
    out$null <- null
  }
  class(out) <- "coco_classic"
  out
}

#' @title Plot a COCO result (Li et al. 2018 method)
#' @description Correlation, H0 significance, number of contributing parameters
#'   and pCOCO against sedimentation rate.
#' @param x COCO result to plot: the output of \code{\link{coco_classic}} (an object of class
#' \code{coco_classic}). Required, no default.
#' @param ... Not used; present for compatibility with the generic \code{plot} function.
#' @return invisible NULL
#' @export
plot.coco_classic <- function(x, ...){
  # Internal helper functions (local to this function)

  .coco_panels <- function(sr, rho, p_list, norbit, pcoco, norb, main = "COCO"){
    np <- 2 + length(p_list) + !is.null(pcoco)
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mfrow = c(np, 1), mar = c(3.2, 4.2, 1.8, 1), mgp = c(2.1, 0.7, 0))
    graphics::plot(sr, rho, type = "l", col = "red", lwd = 1.2,
        xlab = "Sedimentation rate (cm/kyr)", ylab = "Correlation coefficient", main = main)
    for (nm in names(p_list)){
      p <- p_list[[nm]]
      if (all(!is.finite(p))){
        next
      }
      yl <- c(1, max(1e-4, 0.5 * min(p, na.rm = TRUE)))
      graphics::plot(sr, p, type = "l", col = "red", lwd = 1.2, log = "y", ylim = yl,
          xlab = "Sedimentation rate (cm/kyr)", ylab = nm, main = nm)
      graphics::abline(h = c(0.1, 0.05, 0.01, 0.001), lty = c(3, 3, 2, 3), col = "grey30")
    }
    graphics::plot(sr, norbit, type = "s", col = "blue", lwd = 1.2, ylim = c(0, norb + 0.5),
        xlab = "Sedimentation rate (cm/kyr)", ylab = "#",
        main = "Number of contributing astronomical parameters")
    if (!is.null(pcoco) && any(is.finite(pcoco))){
      graphics::plot(sr, pcoco, type = "l", col = "red", lwd = 1.5,
          xlab = "Sedimentation rate (cm/kyr)", ylab = "pCOCO", main = "pCOCO")
      b <- which.max(pcoco)
      graphics::points(sr[b], pcoco[b], pch = 19, col = "darkred")
      graphics::text(sr[b], pcoco[b], sprintf("%.3g cm/kyr", sr[b]), pos = 4, cex = 0.9)
    }
    invisible(NULL)
  }

  # Main code
  r <- x$result
  if (x$nsim > 0){
    pc <- r$pcoco
  } else {
    pc <- NULL
  }
  .coco_panels(r$sr, r$rho, list("H0 significance level" = r$p_h0), r$norbit,
      pc, length(x$orbit), "COCO (Li et al. 2018)")
}

#' @export
#' @noRd
print.coco_classic <- function(x, ...){
  r <- x$result
  if (x$nsim > 0){
    b <- which.max(r$pcoco)
  } else {
    b <- which.max(r$rho)
  }
  cat(sprintf("COCO (Li et al. 2018): %d rates, %d MC sims, %s correlation\n", nrow(r), x$nsim, x$method))
  cat(sprintf("Best rate %.4g cm/kyr: rho = %.3f, H0 p = %.3g, %d of %d parameters\n",
      r$sr[b], r$rho[b], r$p_h0[b], r$norbit[b], length(x$orbit)))
  invisible(x)
}
