# Adaptive COCO (one record) -- corrcoefslices_rankNew.m adaptive path

#' @title Adaptive COCO
#' @description Modern Acycle adaptive COCO: nine-term orbital target
#'   synthesised on the finite record at every tested rate, amplitudes fitted to
#'   the spectrum, AR(1) Monte Carlo null with global (max-over-rate) and local
#'   p-values. Translated from Acycle corrcoefslices_rankNew.m (Mingsong Li).
#' @param data Record to analyse (matrix or data frame, two columns): depth in metres in column 1 and the proxy in
#' column 2. Non-finite rows are dropped, duplicated depths averaged, and unevenly spaced data are linearly
#' interpolated to the median sample spacing (see \code{\link{coco_prepare_regular_data}}). Required, no default.
#' @param orbit9 The nine astronomical target periods in kyr (numeric vector of nine distinct values: long
#' eccentricity, four short eccentricity terms, obliquity, three precession terms), e.g. from
#' \code{\link{orbit9_at_age}} or \code{\link{calculate_orbit9}} for the age of the record. The default uses the
#' present-day periods (405.7, 130.7, 123.9, 98.9, 94.9, 41.0, 23.7, 22.4 and 19.0 kyr).
#' \code{Default=orbit9_default()}.
#' @param sr Sedimentation rates to test, in cm/kyr (numeric vector, > 0); the correlation is computed for every
#' value. \code{Default=seq(0.5, 20, by = 0.1)}.
#' @param pad Length (number of points) to which every slice is zero-padded for the periodogram (integer, at least
#' the number of samples in the longest slice); a larger value gives a finer frequency grid. \code{Default=5000}.
#' @param red Background (red-noise) removal applied to every spectrum before the correlation (integer): \code{0} no
#' removal, \code{1} classical AR(1) background subtracted, \code{2} robust AR(1) background (median-smoothed
#' spectrum, Mann and Lees 1996 style), \code{3} smoothed-window-average (SWA) background (needs at least 33
#' spectral bins). Negative power after subtraction is set to zero. \code{Default=1}.
#' @param nsim Number of AR(1) Monte Carlo simulations used for the global and local p-values (integer, >= 0);
#' \code{0} skips the simulations, so no p-values or pCOCO are returned. \code{Default=1000}.
#' @param slices Number of equal-depth slices into which the record is cut; each slice is standardised and
#' detrended and their spectra are averaged (integer, >= 1; every slice needs at least 4 points).
#' \code{Default=1}.
#' @param method Correlation coefficient between the data spectrum and the target spectrum (character):
#' \code{"Pearson"} or \code{"Spearman"}. \code{Default=c("Pearson", "Spearman")}, i.e. \code{"Pearson"}.
#' @param max_frequency Highest frequency of the target spectrum that enters the correlation, in cycles/kyr
#' (numeric, must be at least the highest orbital frequency \code{max(1/orbit9)}); \code{NULL} uses
#' \code{1.2 * max(1/orbit9)}. \code{Default=NULL}.
#' @param amplitude_mode How the amplitudes of the nine target sinusoids are set (character):
#' \code{"four-group-area"} fits one amplitude per group (long eccentricity, short eccentricity, obliquity,
#' precession) to the band areas of the data spectrum (Acycle default; needs \code{target_model="coherent-nine"}),
#' \code{"per-orbit"} fits each of the nine amplitudes to the spectral peak at its frequency, \code{"fixed"} uses
#' fixed weights (1.0 eccentricity, 0.8 obliquity, 0.6 precession; see \code{\link{coco_fixed_target_weights}}).
#' \code{Default=c("four-group-area", "per-orbit", "fixed")}, i.e. \code{"four-group-area"}.
#' @param target_model How the target spectrum is synthesised on the finite record (character):
#' \code{"coherent-nine"} (the spectrum of the sum of the nine sinusoids, including their interference) or
#' \code{"phase-averaged"} (the sum of the individual spectra, averaged over sine and cosine phase).
#' \code{Default=c("coherent-nine", "phase-averaged")}, i.e. \code{"coherent-nine"}.
#' @param seed Seed of the random number generator (integer), so that the Monte Carlo results can be reproduced;
#' \code{NULL} leaves the current random state unchanged. \code{Default=1}.
#' @param batch Number of Monte Carlo simulations processed together in one batch (integer, >= 1); larger batches
#' are faster but use more memory. \code{NULL} chooses it from \code{pad} (about 4 million spectral values per
#' batch). \code{Default=NULL}.
#' @param verbose Print messages about the preprocessing and the progress of the Monte Carlo simulations
#' (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @return list with `result` data.frame: sr, rho, p_param, p_global, p_local,
#'   norbit, pcoco
#' @export
coco_adaptive <- function(data, orbit9 = orbit9_default(), sr = seq(0.5, 20, by = 0.1),
    pad = 5000, red = 1, nsim = 1000, slices = 1,
    method = c("Pearson", "Spearman"), max_frequency = NULL,
    amplitude_mode = c("four-group-area", "per-orbit", "fixed"),
    target_model = c("coherent-nine", "phase-averaged"),
    seed = 1, batch = NULL, verbose = TRUE){
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

  # Main code
  method <- match.arg(method)
  amplitude_mode <- match.arg(amplitude_mode)
  target_model <- match.arg(target_model)
  if (!is.null(seed)){
    set.seed(seed)
  }
  prep <- coco_prepare_regular_data(data, verbose)
  dat <- prep$data
  dt <- prep$dt
  n <- nrow(dat)
  # slices: standardise + linearly detrend each (data_slices in rankNew)
  b <- seq(dat[1, 1], dat[n, 1], length.out = slices + 1)
  sidx <- vector("list", slices)
  for (j in 1:slices){
    if (j < slices){
      upper_ok <- dat[, 1] < b[j + 1]
    } else {
      upper_ok <- dat[, 1] <= b[j + 1]
    }
    sidx[[j]] <- which(dat[, 1] >= b[j] & upper_ok)
  }
  if (any(lengths(sidx) < 4)){
    stop("Every slice needs at least 4 points")
  }
  if (pad < max(lengths(sidx))){
    stop("pad must be >= samples in the longest slice")
  }
  svals <- vector("list", length(sidx))
  for (j in seq_along(sidx)){
    v <- dat[sidx[[j]], 2]
    svals[[j]] <- detrend_linear((v - mean(v)) / stats::sd(v))
  }
  Ps <- NULL
  for (j in seq_along(svals)){
    pw <- coco_window_spectra(svals[[j]], dt, pad, red, already_detrended = TRUE)$power
    if (is.null(Ps)){
      Ps <- matrix(0, length(pw), length(svals))
    }
    Ps[, j] <- pw
  }
  Pobs <- matrix(rowMeans(as.matrix(Ps)), ncol = 1)
  f <- (0:(floor(pad / 2))) / (pad * dt)
  if (slices == 1){
    ndata <- n
  } else {
    ndata <- max(lengths(sidx))
  }
  ev <- function(Pm, obsp = FALSE){
    coco_adaptive_evaluate(Pm, f, ndata, dt, pad, orbit9, sr, method, max_frequency,
        target_model, amplitude_mode, observed_p = obsp)
  }
  obs <- ev(Pobs, TRUE)
  rho <- obs$rho[, 1]
  norbit <- length(orbit9) - obs$nMissing
  p_local <- p_global <- rep(NA_real_, length(sr))
  null <- NULL
  if (nsim > 0){
    # pooled conditional LS AR(1) within slices
    num <- 0
    den <- 0
    for (v in svals){
      x <- v - mean(v)
      num <- num + sum(x[-length(x)] * x[-1])
      den <- den + sum(x[-length(x)]^2)
    }
    rhoM <- min(max(num / den, -0.999), 0.999)
    sdv <- sqrt(sum(detrend_linear(dat[, 2])^2) / (n - 1))
    if (is.null(batch)){
      batch <- max(1, min(nsim, floor(4e6 / (floor(pad / 2) + 1))))
    }
    null <- matrix(NA_real_, length(sr), nsim)
    done <- 0
    while (done < nsim){
      nb <- min(batch, nsim - done)
      z <- matrix(stats::rnorm(n * nb), n, nb)
      z[-1, ] <- sqrt(1 - rhoM^2) * z[-1, ]
      rs <- sdv * stats::filter(z, rhoM, method = "recursive")
      rs <- matrix(rs, n, nb)
      Pnull <- 0
      for (ii in sidx){
        x <- rs[ii, , drop = FALSE]
        if (slices > 1){
          sds <- numeric(ncol(x))
          for (j in seq_len(ncol(x))){
            sds[j] <- stats::sd(x[, j])
          }
          x <- sweep(sweep(x, 2, colMeans(x)), 2, sds, "/")
        }
        Pnull <- Pnull + coco_window_spectra(x, dt, pad, red)$power
      }
      Pnull <- Pnull / slices
      null[, done + seq_len(nb)] <- ev(Pnull)$rho
      done <- done + nb
      if (verbose){
        message(sprintf("Adaptive COCO Monte Carlo %d of %d", done, nsim))
      }
    }
    okr <- is.finite(rho)
    null_ok <- null[okr, , drop = FALSE]
    nullMax <- numeric(ncol(null_ok))
    for (j in seq_len(ncol(null_ok))){
      nullMax[j] <- max(null_ok[, j], na.rm = TRUE)
    }
    for (k in which(okr)){
      p_local[k] <- (sum(null[k, ] >= rho[k], na.rm = TRUE) + 1) / (nsim + 1)
      p_global[k] <- (sum(nullMax >= rho[k]) + 1) / (nsim + 1)
    }
  }
  res <- data.frame(sr = sr, rho = rho, p_param = obs$p, p_global = p_global,
      p_local = p_local, norbit = norbit)
  res$pcoco <- pcoco_value(res$rho, res$p_global)
  if (nsim > 0){
    null_max <- nullMax
  } else {
    null_max <- NULL
  }
  out <- list(result = res, spectrum = cbind(f = f, power = Pobs[, 1]), dt = dt,
      orbit9 = orbit9, nsim = nsim, amplitude_mode = amplitude_mode,
      method = method, null_max = null_max)
  class(out) <- "coco_adaptive"
  out
}

#' @title Plot an adaptive COCO result
#' @description Correlation, global and local p, contributing periods and pCOCO
#'   against sedimentation rate.
#' @param x Adaptive COCO result to plot: the output of \code{\link{coco_adaptive}} (an object of class
#' \code{coco_adaptive}). Required, no default.
#' @param ... Not used; present for compatibility with the generic \code{plot} function.
#' @return invisible NULL
#' @export
plot.coco_adaptive <- function(x, ...){
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
  .coco_panels(r$sr, r$rho, list("Global p" = r$p_global, "Local p" = r$p_local),
      r$norbit, pc, length(x$orbit9), "Adaptive COCO")
}

#' @export
#' @noRd
print.coco_adaptive <- function(x, ...){
  r <- x$result
  if (x$nsim > 0){
    b <- which.max(r$pcoco)
  } else {
    b <- which.max(r$rho)
  }
  cat(sprintf("Adaptive COCO (%s): %d rates, %d MC sims\n", x$amplitude_mode, nrow(r), x$nsim))
  cat(sprintf("Best rate %.4g cm/kyr: rho = %.3f, global p = %.3g, local p = %.3g, %d of 9 periods\n",
      r$sr[b], r$rho[b], r$p_global[b], r$p_local[b], r$norbit[b]))
  invisible(x)
}
