# Internal helpers of flmw() and flmw_path() --------------------------------------------------------------------------

# Spectral peaks of one wavelet spectrum (column of the power matrix). Reproduces the peak selection that flmw() always
# used: local maxima/minima are found by comparing a point with the points two steps away, after which consecutive
# maxima (or minima) are reduced so that maxima and minima alternate. Returns the periods and powers of the maxima.
.flmw_peaks <- function(period, power) {
  n <- length(power)
  empty <- list(period = numeric(0), power = numeric(0))
  if (n < 5) return(empty)
  i <- 3:(n - 2)
  is_max <- (power[i] - power[i + 2] > 0) & (power[i] - power[i - 2] > 0)
  is_min <- (power[i] - power[i + 2] < 0) & (power[i] - power[i - 2] < 0)
  sel <- which(is_max | is_min)
  if (length(sel) == 0) return(empty)
  idx <- i[sel]
  per <- period[idx]
  pw <- power[idx]
  ty <- ifelse(is_max[sel], 1, -1)
  yes <- function(z) !is.na(z) && z
  i <- 1
  while (i < length(ty)) {
    if (ty[i] == ty[i + 1]) {
      dropped <- FALSE
      # the four checks of the original code, each one evaluated on the already reduced set
      drops <- list(
        list(t = 1,  cmp = `>`, which = 1),
        list(t = 1,  cmp = `<`, which = 0),
        list(t = -1, cmp = `<`, which = 1),
        list(t = -1, cmp = `>`, which = 0)
      )
      for (d in drops) {
        if (yes(ty[i] == d$t) && yes(ty[i + 1] == d$t) && yes(d$cmp(pw[i], pw[i + 1]))) {
          k <- i + d$which
          per <- per[-k]; pw <- pw[-k]; ty <- ty[-k]
          dropped <- TRUE
        }
      }
      if (!dropped) {                      # identical power: keep the first one
        per <- per[-(i + 1)]; pw <- pw[-(i + 1)]; ty <- ty[-(i + 1)]
      }
    }
    if (is.na(ty[i + 1]) || ty[i] != ty[i + 1]) i <- i + 1
  }
  keep <- ty > 0
  list(period = per[keep], power = pw[keep])
}

# Classic flmw metrics for one set of peak periods and all trial sedimentation rates (vectorised over the rates).
# Returns a matrix [rate, 4]: |1 - slope|, r squared, number of matched cycles, 1 - |intercept|
.flmw_metrics <- function(pk, cycles, srs) {
  nsr <- length(srs)
  out <- matrix(0, nsr, 4)
  nc <- length(cycles)
  if (length(pk) == 0 || nc == 0) return(out)
  tgt <- outer(cycles, srs) / 100                              # cycles in depth [cycle, rate]
  K <- matrix(max.col(-abs(outer(as.vector(tgt), pk, "-")), ties.method = "first"), nc, nsr)
  X <- matrix(pk[K], nc, nsr)                                  # nearest peak
  D <- abs(X - tgt)
  keep <- matrix(TRUE, nc, nsr)
  if (nc > 1) {
    for (k in seq_len(nc)) {
      for (k2 in seq_len(nc)) {
        if (k2 == k) next
        keep[k, ] <- keep[k, ] & !((K[k, ] == K[k2, ]) & (D[k2, ] < D[k, ] | (D[k2, ] == D[k, ] & k2 > k)))
      }
    }
  }
  n <- colSums(keep)
  out[, 3] <- n
  ok <- n >= 2
  if (any(ok)) {
    w <- keep * 1
    xbar <- colSums(w * X) / pmax(n, 1)
    ybar <- colSums(w * tgt) / pmax(n, 1)
    dx <- (X - rep(xbar, each = nc)) * w
    dy <- (tgt - rep(ybar, each = nc)) * w
    sxx <- colSums(dx^2); syy <- colSums(dy^2); sxy <- colSums(dx * dy)
    slope <- sxy / sxx
    icpt <- ybar - slope * xbar
    out[ok, 1] <- abs(1 - slope[ok])
    out[ok, 2] <- (sxy^2 / (sxx * syy))[ok]
    out[ok, 4] <- 1 - abs(icpt[ok])
  }
  out
}

# All eight classic flmw columns for one depth level
.flmw_depth_fit <- function(period, power, cycles, srs, run_random, rand_simulations, randomize_mat) {
  pk <- .flmw_peaks(period, power)
  fits <- matrix(NA_real_, length(srs), 8)
  fits[, 1:4] <- .flmw_metrics(pk$period, cycles, srs)
  if (run_random) {
    np <- length(pk$period)
    rnd <- array(NA_real_, c(length(srs), 4, rand_simulations))
    for (rt in seq_len(rand_simulations)) {
      rp <- if (np > 0) sort(pk$period * randomize_mat[seq_len(np), rt]) else numeric(0)
      rnd[, , rt] <- .flmw_metrics(rp, cycles, srs)
    }
    pm <- function(j) rowMeans(matrix(rnd[, j, ], nrow = length(srs)) <= fits[, j])
    fits[, 5] <- 1 - pm(1)
    fits[, 6] <- pm(2)
    fits[, 7] <- pm(3)
    fits[, 8] <- pm(4)
  }
  fits
}

# Weighted one-to-one matching of the predicted cycles to the spectral peaks of one depth level.
# Every predicted cycle and every peak can be used once. The strength of a pair is the (relative) power of the peak times a
# Gaussian in the relative (log) period misfit. Returns the score per trial sedimentation rate (0-1).
.flmw_weighted_depth <- function(period, power, cycles, srs, tolerance) {
  pk <- .flmw_peaks(period, power)
  nsr <- length(srs)
  nc <- length(cycles)
  score <- numeric(nsr)
  np <- length(pk$period)
  if (np == 0 || !is.finite(max(power)) || max(power) <= 0) return(score)
  a <- pk$power / max(power)
  lp <- log(pk$period)
  for (s in seq_len(nsr)) {
    tg <- log(cycles * srs[s] / 100)
    M <- exp(-0.5 * (outer(tg, lp, "-") / tolerance)^2) * rep(a, each = nc)
    tot <- 0
    for (r in seq_len(min(nc, np))) {
      mx <- which.max(M)
      if (M[mx] <= 0) break
      tot <- tot + M[mx]
      rr <- (mx - 1) %% nc + 1
      cc <- (mx - 1) %/% nc + 1
      M[rr, ] <- 0
      M[, cc] <- 0
    }
    score[s] <- tot / nc
  }
  score
}

# Template correlation score for all depths and rates at once (no peak picking).
.flmw_template_scores <- function(period, Power, cycles, srs, tolerance) {
  lp <- log(period)
  nsr <- length(srs)
  nc <- length(cycles)
  Tm <- matrix(0, length(lp), nsr)
  frac <- numeric(nsr)
  for (s in seq_len(nsr)) {
    tg <- log(cycles * srs[s] / 100)
    inside <- tg >= min(lp) & tg <= max(lp)
    frac[s] <- sum(inside) / nc
    for (k in which(inside)) Tm[, s] <- Tm[, s] + exp(-0.5 * ((lp - tg[k]) / tolerance)^2)
  }
  cs <- function(M) {
    M <- sweep(M, 2, colMeans(M), "-")
    nrm <- sqrt(colSums(M^2))
    nrm[nrm == 0] <- NA
    sweep(M, 2, nrm, "/")
  }
  R <- crossprod(cs(Tm), cs(Power))
  R[is.na(R)] <- 0
  R[R < 0] <- 0
  R * frac
}
