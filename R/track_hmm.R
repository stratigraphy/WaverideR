## ---------------------------------------------------------------------------------------------
## Hidden Markov model ridge tracking of (synchrosqueezed) scalograms
## ---------------------------------------------------------------------------------------------

# Engine. logE: matrix states x positions of log emission scores. Transition between neighbouring
# positions: a discrete Gaussian random walk over the states with sd sig (in states) mixed with a
# uniform jump with probability rho (abrupt changes, e.g. at a hiatus). Returns the Viterbi path,
# the posterior state probabilities and n_samples paths sampled from the posterior (forward filtering,
# backward sampling).
.wr_hmm <- function(logE, sig, rho = 1e-3, n_samples = 0, start = NULL) {
  ns <- nrow(logE)
  nt <- ncol(logE)
  # rho may vary along the record: rho[t] is the jump probability between positions t - 1 and t
  rho <- rep_len(rho, nt)
  k <- max(1, ceiling(4 * sig))
  offs <- -k:k
  K <- exp(-offs^2 / (2 * sig^2))
  K <- K / sum(K)
  shift <- function(v, d, fill = 0) {        # out[j] = v[j - d]
    if (d == 0) return(v)
    if (d > 0) c(rep(fill, d), v[seq_len(ns - d)]) else c(v[(1 - d):ns], rep(fill, -d))
  }
  conv <- function(v, r) {                   # sum_i v[i] T(i, j)
    out <- numeric(ns)
    for (m in seq_along(offs)) out <- out + (1 - r) * K[m] * shift(v, offs[m])
    out + r * sum(v) / ns
  }
  E <- exp(logE - matrix(apply(logE, 2, max), ns, nt, byrow = TRUE))
  prior <- if (is.null(start)) rep(1 / ns, ns) else start / sum(start)
  alpha <- matrix(0, ns, nt)
  a <- prior * E[, 1]
  alpha[, 1] <- a / sum(a)
  for (t in 2:nt) {
    a <- conv(alpha[, t - 1], rho[t]) * E[, t]
    s <- sum(a)
    alpha[, t] <- if (s > 0) a / s else rep(1 / ns, ns)
  }
  beta <- matrix(0, ns, nt)
  beta[, nt] <- 1
  for (t in (nt - 1):1) {
    b <- conv(E[, t + 1] * beta[, t + 1], rho[t + 1])     # T is symmetric apart from the edges
    s <- sum(b)
    beta[, t] <- if (s > 0) b / s else rep(1, ns)
  }
  post <- alpha * beta
  post <- post / matrix(colSums(post), ns, nt, byrow = TRUE)
  # Viterbi
  delta <- log(prior + 1e-300) + logE[, 1]
  from <- matrix(0L, ns, nt)
  for (t in 2:nt) {
    lK <- log((1 - rho[t]) * K)
    best <- rep(-Inf, ns)
    arg <- rep(NA_integer_, ns)
    for (m in seq_along(offs)) {
      cand <- shift(delta, offs[m], -Inf) + lK[m]
      better <- cand > best
      best[better] <- cand[better]
      arg[better] <- (seq_len(ns) - offs[m])[better]
    }
    jm <- which.max(delta)
    jump <- delta[jm] + log(rho[t] / ns)
    better <- jump > best
    best[better] <- jump
    arg[better] <- jm
    from[, t] <- arg
    delta <- best + logE[, t]
  }
  path <- integer(nt)
  path[nt] <- which.max(delta)
  for (t in nt:2) path[t - 1] <- from[path[t], t]
  samples <- NULL
  if (n_samples > 0) {
    samples <- matrix(0L, nt, n_samples)
    g <- function(p) {                       # sample one state per column of a probability matrix
      lp <- log(p + 1e-300) - log(-log(matrix(stats::runif(length(p)), nrow(p))))
      max.col(t(lp), ties.method = "first")
    }
    samples[nt, ] <- g(matrix(alpha[, nt], ns, n_samples))
    idx <- seq_len(ns)
    for (t in (nt - 1):1) {
      d <- outer(idx, samples[t + 1, ], function(i, z) z - i)
      tr <- ifelse(abs(d) <= k, (1 - rho[t + 1]) * K[pmin(pmax(d + k + 1, 1), length(K))], 0) + rho[t + 1] / ns
      samples[t, ] <- g(alpha[, t] * tr)
    }
  }
  list(path = path, post = post, samples = samples)
}

# Viterbi path only (used for the search of the Lagrange multipliers of duration constraints)
.wr_viterbi <- function(logE, sig, rho) {
  ns <- nrow(logE)
  nt <- ncol(logE)
  rho <- rep_len(rho, nt)
  k <- max(1, ceiling(4 * sig))
  offs <- -k:k
  K <- exp(-offs^2 / (2 * sig^2))
  K <- K / sum(K)
  shift <- function(v, d) {
    if (d == 0) return(v)
    if (d > 0) c(rep(-Inf, d), v[seq_len(ns - d)]) else c(v[(1 - d):ns], rep(-Inf, -d))
  }
  delta <- logE[, 1]
  from <- matrix(0L, ns, nt)
  for (t in 2:nt) {
    lK <- log((1 - rho[t]) * K)
    best <- rep(-Inf, ns)
    arg <- rep(NA_integer_, ns)
    for (m in seq_along(offs)) {
      cand <- shift(delta, offs[m]) + lK[m]
      better <- cand > best
      best[better] <- cand[better]
      arg[better] <- (seq_len(ns) - offs[m])[better]
    }
    jm <- which.max(delta)
    jump <- delta[jm] + log(rho[t] / ns)
    better <- jump > best
    best[better] <- jump
    arg[better] <- jm
    from[, t] <- arg
    delta <- best + logE[, t]
  }
  path <- integer(nt)
  path[nt] <- which.max(delta)
  for (t in nt:2) path[t - 1] <- from[path[t], t]
  path
}

# Exact most probable path under duration constraints: inside every constrained interval the state is extended by
# the elapsed time (in bins of bin_frac * duration), so that the path that maximises the evidence among all paths
# whose interval duration lies within the tolerance is found by dynamic programming. Intervals must not overlap.
.wr_viterbi_constrained <- function(logE, sig, rho, Tm, zc, cons, bin_frac = 0.002) {
  ns <- nrow(logE)
  nt <- ncol(logE)
  rho <- rep_len(rho, nt)
  k <- max(1, ceiling(4 * sig))
  offs <- -k:k
  K <- exp(-offs^2 / (2 * sig^2))
  K <- K / sum(K)
  cons <- cons[order(cons$from), , drop = FALSE]
  iv <- lapply(seq_len(nrow(cons)), function(i) which(zc >= cons$from[i] & zc <= cons$to[i]))
  start_of <- rep(0L, nt); in_iv <- rep(0L, nt)
  for (i in seq_along(iv)) if (length(iv[[i]])) { start_of[iv[[i]][1]] <- i; in_iv[iv[[i]]] <- i }
  # elapsed-time increments in bins: differences of the floored cumulative time of each state, so that a
  # path that stays near a state accumulates the right total (rounding each step separately would count
  # steps much shorter than a bin as zero time)
  Cum <- t(apply(Tm, 1, cumsum))
  if (nt == 1) Cum <- matrix(Cum, ncol = 1)
  inc_at <- function(t, bw) {
    prev <- if (t > 1) floor(Cum[, t - 1] / bw) else 0
    as.integer(floor(Cum[, t] / bw) - prev)
  }
  shift_s <- function(v, d) {                      # v[j - d] along the state axis (vector or matrix rows)
    if (is.matrix(v)) {
      if (d == 0) return(v)
      fill <- matrix(-Inf, abs(d), ncol(v))
      if (d > 0) rbind(fill, v[seq_len(ns - d), , drop = FALSE]) else rbind(v[(1 - d):ns, , drop = FALSE], fill)
    } else {
      if (d == 0) return(v)
      if (d > 0) c(rep(-Inf, d), v[seq_len(ns - d)]) else c(v[(1 - d):ns], rep(-Inf, -d))
    }
  }
  # best predecessor (state) for every target state, for a vector or a matrix (columns = time bins)
  best_pred <- function(V, t) {
    lK <- log((1 - rho[t]) * K)
    best <- if (is.matrix(V)) matrix(-Inf, ns, ncol(V)) else rep(-Inf, ns)
    arg <- if (is.matrix(V)) matrix(NA_integer_, ns, ncol(V)) else rep(NA_integer_, ns)
    for (m in seq_along(offs)) {
      cand <- shift_s(V, offs[m]) + lK[m]
      better <- cand > best
      best[better] <- cand[better]
      src <- seq_len(ns) - offs[m]
      if (is.matrix(V)) arg[better] <- matrix(src, ns, ncol(V))[better] else arg[better] <- src[better]
    }
    lj <- log(rho[t] / ns)
    if (is.matrix(V)) {
      jm <- max.col(t(V), ties.method = "first")
      jv <- V[cbind(jm, seq_len(ncol(V)))] + lj
      better <- matrix(jv, ns, ncol(V), byrow = TRUE) > best
      best[better] <- matrix(jv, ns, ncol(V), byrow = TRUE)[better]
      arg[better] <- matrix(jm, ns, ncol(V), byrow = TRUE)[better]
    } else {
      jm <- which.max(V); jv <- V[jm] + lj
      better <- jv > best; best[better] <- jv; arg[better] <- jm
    }
    list(best = best, arg = arg)
  }
  ptr <- vector("list", nt)                      # plain positions: integer vector; constrained: integer matrix
  binptr <- vector("list", nt)                   # at the closing of an interval: chosen bin per state
  delta <- logE[, 1]
  cur <- 0L; D <- NULL; bw <- NULL; nb <- NULL; inc <- NULL
  close_iv <- function(D, i) {
    lo <- floor(cons$duration[i] * (1 - cons$tolerance[i]) / bw); hi <- ceiling(cons$duration[i] * (1 + cons$tolerance[i]) / bw)
    ok <- seq_len(nb) - 1L
    valid <- ok >= lo & ok <= hi
    Dv <- D; Dv[, !valid] <- -Inf
    b <- max.col(Dv, ties.method = "first")
    list(delta = Dv[cbind(seq_len(ns), b)], bin = b)
  }
  t_first <- if (start_of[1] > 0) 1L else NA
  if (!is.na(t_first)) {
    i <- start_of[1]; bw <- cons$duration[i] * bin_frac
    nb <- ceiling(cons$duration[i] * (1 + cons$tolerance[i]) / bw) + 2L
    D <- matrix(-Inf, ns, nb)
    b0 <- pmin(inc_at(1, bw), nb - 1L) + 1L
    D[cbind(seq_len(ns), b0)] <- logE[, 1]
    cur <- i
  }
  for (t in 2:nt) {
    if (cur > 0 && in_iv[t] == cur) {
      # inside the interval: extend the elapsed time
      bp <- best_pred(D, t)
      inc <- inc_at(t, bw)
      Dn <- matrix(-Inf, ns, nb); An <- matrix(NA_integer_, ns, nb)
      for (s in seq_len(ns)) {
        sh <- inc[s]
        if (sh >= nb) next
        idx <- (1 + sh):nb
        Dn[s, idx] <- bp$best[s, idx - sh]
        An[s, idx] <- bp$arg[s, idx - sh]
      }
      D <- Dn + logE[, t]
      ptr[[t]] <- An
      next
    }
    if (cur > 0) {
      # the interval ended at t - 1
      cl <- close_iv(D, cur); binptr[[t - 1]] <- cl$bin; delta <- cl$delta; D <- NULL; cur <- 0L
    }
    if (start_of[t] > 0) {
      i <- start_of[t]; bw <- cons$duration[i] * bin_frac
      nb <- ceiling(cons$duration[i] * (1 + cons$tolerance[i]) / bw) + 2L
      bp <- best_pred(delta, t)
      D <- matrix(-Inf, ns, nb)
      b0 <- pmin(inc_at(t, bw), nb - 1L) + 1L
      D[cbind(seq_len(ns), b0)] <- bp$best + logE[, t]
      ptr[[t]] <- bp$arg
      cur <- i
      next
    }
    bp <- best_pred(delta, t)
    delta <- bp$best + logE[, t]
    ptr[[t]] <- bp$arg
  }
  if (cur > 0) { cl <- close_iv(D, cur); binptr[[nt]] <- cl$bin; delta <- cl$delta }
  # backtracking
  path <- integer(nt)
  path[nt] <- which.max(delta)
  b <- if (!is.null(binptr[[nt]])) binptr[[nt]][path[nt]] else NA
  for (t in nt:2) {
    s <- path[t]
    if (is.matrix(ptr[[t]])) {
      # inside an interval (not its first position): pointer depends on the bin
      if (is.na(b)) b <- binptr[[t]][s]
      prev <- ptr[[t]][s, b]
      b <- b - inc_at(t, cons$duration[in_iv[t]] * bin_frac)[s]
      path[t - 1] <- prev
    } else {
      path[t - 1] <- ptr[[t]][s]
      b <- NA
    }
    if (!is.null(binptr[[t - 1]])) b <- binptr[[t - 1]][path[t - 1]]
  }
  path
}

# Duration constraints (geological constraints such as a biostratigraphic bracket or radiometric ages): for every
# constraint the log evidence is tilted by lambda * (time spent at each position), and lambda is found by bisection
# such that the most probable path gives the required duration of the interval (Lagrangian relaxation of the
# constrained Viterbi problem). inv_rate: time per unit depth of every state.
.wr_hmm_constrained <- function(logE, sig, rho, n_samples, inv_rate, zc, constraints, tol = 0.025,
                               n_draw = 600, min_kept = 100, max_batches = 8, verbose = FALSE) {
  if (is.null(constraints) || !nrow(constraints)) {
    return(list(h = .wr_hmm(logE, sig, rho, n_samples = n_samples), constraints = NULL))
  }
  constraints <- as.data.frame(constraints)
  names(constraints)[1:3] <- c("from", "to", "duration")
  if (is.null(constraints$tolerance)) constraints$tolerance <- tol
  nt <- ncol(logE)
  ns <- nrow(logE)
  step <- c(diff(zc), utils::tail(diff(zc), 1))
  Tm <- outer(inv_rate, step)
  masks <- lapply(seq_len(nrow(constraints)), function(i) zc >= min(constraints$from[i], constraints$to[i]) &
                    zc <= max(constraints$from[i], constraints$to[i]))
  typ <- stats::median(Tm)
  tilt <- function(l) {
    out <- logE
    for (i in seq_along(l)) if (l[i] != 0) out[, masks[[i]]] <- out[, masks[[i]]] + l[i] * Tm[, masks[[i]]] / typ
    out
  }
  durs <- function(S, i) colSums(matrix(Tm[cbind(as.vector(S), rep(seq_len(nt), ncol(S)))], nt)[masks[[i]], , drop = FALSE])
  # The constraints condition the posterior: paths are drawn from the posterior tilted by Lagrange multipliers chosen so
  # that the drawn durations centre on the targets, and only draws within the tolerance are kept. Within that narrow
  # window the tilt is (nearly) constant, so the kept draws follow the untilted posterior conditioned on the constraints.
  lam <- rep(0, nrow(constraints))
  for (round in 1:4) {
    for (i in seq_len(nrow(constraints))) {
      f <- function(u) { l <- lam; l[i] <- u; stats::median(durs(.wr_hmm(tilt(l), sig, rho, n_samples = 60)$samples, i)) }
      target <- constraints$duration[i]
      lo <- -6; hi <- 6
      if (f(lo) >= target) { lam[i] <- lo; next }
      if (f(hi) <= target) { lam[i] <- hi; next }
      for (it in 1:18) {
        mid <- (lo + hi) / 2
        dm <- f(mid)
        if (abs(dm - target) <= 0.01 * target) break
        if (dm < target) lo <- mid else hi <- mid
      }
      lam[i] <- mid
    }
    if (nrow(constraints) == 1) break
  }
  # draw batches until enough paths satisfy the constraints (few do when the posterior is multimodal and the target
  # lies between the modes)
  Lt <- tilt(lam)
  S <- NULL; D <- NULL
  for (batch in seq_len(max_batches)) {
    Sb <- .wr_hmm(Lt, sig, rho, n_samples = n_draw)$samples
    Db <- matrix(sapply(seq_len(nrow(constraints)), function(i) durs(Sb, i)), ncol = nrow(constraints))
    S <- cbind(S, Sb); D <- rbind(D, Db)
    ok <- rep(TRUE, ncol(S))
    for (i in seq_len(nrow(constraints))) ok <- ok & abs(D[, i] - constraints$duration[i]) <= constraints$tolerance[i] * constraints$duration[i]
    if (sum(ok) >= min_kept) break
  }
  if (verbose) cat("draws:", ncol(S), " within tolerance:", sum(ok), "\n")
  if (sum(ok) < 10) {
    # fall back to the draws closest to the targets
    dist <- rowSums(abs(sweep(D, 2, constraints$duration)) / matrix(constraints$duration, nrow(D), ncol(D), byrow = TRUE))
    ok <- rank(dist, ties.method = "first") <= 10
    warning("few posterior paths satisfy the duration constraints; the 10 closest paths are used (achieved: ",
            paste(round(colMeans(D[ok, , drop = FALSE])), collapse = ", "), ")")
  }
  Sk <- S[, ok, drop = FALSE]
  # untilted log probability of every kept path: evidence plus the random-walk transition terms
  k <- max(1, ceiling(4 * sig))
  K <- exp(-((-k):k)^2 / (2 * sig^2)); K <- K / sum(K)
  rr <- rep_len(rho, nt)
  lp <- apply(Sk, 2, function(p) {
    dd <- diff(p)
    tr <- ifelse(abs(dd) <= k, (1 - rr[-1]) * K[pmin(pmax(dd + k + 1, 1), length(K))], 0) + rr[-1] / ns
    sum(logE[cbind(p, seq_len(nt))]) + sum(log(tr))
  })
  # exact constrained Viterbi on a quantised elapsed-time axis; the path may exploit the rounding of the
  # time per step, so the target is corrected with the rounding error of the previous solution
  achieved_of <- function(path) vapply(seq_len(nrow(constraints)), function(i) sum(Tm[cbind(path, seq_len(nt))][masks[[i]]]), numeric(1))
  cons_dp <- constraints
  cons_dp$tolerance <- constraints$tolerance / 2
  best <- NULL
  for (it in 1:6) {
    cand <- tryCatch(.wr_viterbi_constrained(logE, sig, rho, Tm, zc, cons_dp, bin_frac = getOption("WaverideR.dp_bin", 0.0005)), error = function(e) NULL)
    if (is.null(cand) || any(is.na(cand))) break
    a <- achieved_of(cand)
    rel <- abs(a - constraints$duration) / constraints$duration
    if (is.null(best) || max(rel) < best_rel) { best <- cand; best_rel <- max(rel) }
    if (isTRUE(getOption("WaverideR.debug_dp"))) message("dp iter ", it, " ", format(Sys.time()), ": achieved ", paste(round(a), collapse = ","))
    if (all(rel <= constraints$tolerance / 2)) break
    cons_dp$duration <- pmax(cons_dp$duration - (a - constraints$duration), 0.2 * constraints$duration)
  }
  if (is.null(best) || best_rel > max(constraints$tolerance)) {
    if (is.null(best) || best_rel > max(constraints$tolerance)) best <- Sk[, which.max(lp)]
  }
  post <- matrix(0, ns, nt)
  for (t in seq_len(nt)) post[, t] <- tabulate(Sk[t, ], ns) / ncol(Sk)
  samp <- Sk[, sample.int(ncol(Sk), n_samples, replace = TRUE), drop = FALSE]
  constraints$achieved <- vapply(seq_len(nrow(constraints)), function(i) sum(Tm[cbind(best, seq_len(nt))][masks[[i]]]), numeric(1))
  constraints$lambda <- lam
  constraints$n_paths <- sum(ok)
  if (verbose) cat("lambda:", signif(lam, 3), " kept paths:", sum(ok), "\n")
  list(h = list(path = best, post = post, samples = if (n_samples > 0) samp else NULL), constraints = constraints)
}

# Per-position jump probabilities: the base value everywhere, raised to jump_prob_at near the given depths
# (e.g. hiatuses, erosion surfaces or abrupt facies changes known from the field)
.wr_jump_rho <- function(zc, jump_prob, jump_depths, jump_prob_at) {
  nt <- length(zc)
  rho <- rep(jump_prob, nt)
  if (length(jump_depths)) {
    step <- if (nt > 1) stats::median(diff(zc)) else 1
    for (d in jump_depths) {
      near <- which(abs(zc - d) <= 1.5 * step)
      rho[near] <- pmax(rho[near], jump_prob_at)
    }
  }
  rho
}

# Long-window evidence (eCOCO-like): the wavelet power (Morlet, omega 12) averaged over long windows, with the
# red-noise slope removed (linear fit of log power against log period), correlated with a template of the
# astronomical periods for every sedimentation rate. Correlating spectral shapes over a long window is much less
# sensitive to noise than the local power, but has a poor resolution along the record; it is used as an additional
# soft term that keeps the track on the right cycles.
.wr_long_evidence <- function(wavelet, zc, long_window, ls, astro_periods, weights, sigma = 0.06) {
  x <- wavelet$x
  y <- wavelet$y
  dt <- wavelet$dt
  len <- diff(range(x))
  lo <- max(2 * dt, min(astro_periods) * 2^min(ls) * 0.8)
  hi <- min(len / 2, max(astro_periods) * 2^max(ls) * 1.25)
  cw <- analyse_wavelet(data.frame(x, y), dj = 1 / 50, lowerPeriod = lo, upperPeriod = hi, omega_nr = 12)
  lp <- log2(cw$Period)
  cols <- vapply(zc, function(v) which.min(abs(cw$axis.1 - v)), integer(1))
  Ps <- .wr_smooth_map(cw$Power, 0, max(1, round(long_window / dt)))[, cols, drop = FALSE]
  L <- log(Ps + 1e-12 * max(Ps))
  X <- cbind(1, lp)
  L <- L - X %*% solve(crossprod(X), crossprod(X, L))
  R <- scale(L)
  Tm <- vapply(ls, function(l) {
    v <- numeric(length(lp))
    # only cycles that repeat at least three times within the window are resolved by its averaged spectrum
    for (k in seq_along(astro_periods)) if (astro_periods[k] * 2^l <= long_window / 3)
      v <- v + weights[k] * exp(-(lp - (log2(astro_periods[k]) + l))^2 / (2 * sigma^2))
    v
  }, numeric(length(lp)))
  sdv <- apply(Tm, 2, stats::sd)
  Tm <- sweep(sweep(Tm, 2, colMeans(Tm)), 2, ifelse(sdv > 0, sdv, 1), "/")
  C <- t(Tm) %*% R / (length(lp) - 1)
  C[!is.finite(C)] <- 0
  C
}

#' @title Astronomical periods and weights for tracking, for any age
#'
#' @description Returns the main eccentricity, obliquity and precession periods with weights proportional to their
#' typical amplitude, for use in \code{\link{track_astro_hmm}}. The obliquity and precession periods follow from the
#' precession constant k at the given age (AstroGeo22, via \code{\link{astro_terms}}) or from a k given directly (e.g.
#' estimated with \code{\link{test_tuned}}), so that the periods fit deep-time records such as the Silurian and Devonian.
#' The eccentricity periods are the La2004 values; the 405 kyr period is stable through the Phanerozoic, the short
#' eccentricity periods are less certain before about 60 Ma.
#'
#' @param age_ma Age in Ma \code{Default=0}.
#' @param k_yr Optional precession constant (arcsec/yr) instead of the AstroGeo22 value for \code{age_ma}.
#' @param include Which families to include \code{Default=c("eccentricity", "obliquity", "precession")}.
#' When only \code{"eccentricity"} is included, the 405-kyr cycle gets weight 1 (instead of 0.5): in records in which
#' only the eccentricity is resolved, the stable 405-kyr cycle carries the age model.
#' @return A data frame with \code{term}, \code{family}, \code{period} (kyr) and \code{weight}.
#' @examples
#' astro_periods_tracking(0)
#' astro_periods_tracking(420)        # Silurian
#' @export
astro_periods_tracking <- function(age_ma = 0, k_yr = NULL,
                                   include = c("eccentricity", "obliquity", "precession")) {
  at <- astro_terms(age_ma = age_ma, k_yr = k_yr)
  pick <- data.frame(term = c("g2-g5", "g4-g2", "g3-g2", "g4-g5", "g3-g5", "k+s3", "k+s4", "k+s6",
                              "k+g5", "k+g2", "k+g4", "k+g3"),
                     family = c(rep("eccentricity", 5), rep("obliquity", 3), rep("precession", 4)),
                     weight = c(0.5, 0.3, 0.2, 0.3, 0.2, 1, 0.4, 0.3, 1, 0.8, 0.6, 0.4),
                     stringsAsFactors = FALSE)
  pick$period <- at$period_kyr[match(pick$term, at$term)]
  pick <- pick[pick$family %in% include & is.finite(pick$period), c("term", "family", "period", "weight")]
  # when only the eccentricity is tracked, the stable 405-kyr cycle carries the age model
  if (all(pick$family == "eccentricity")) pick$weight[pick$term == "g2-g5"] <- 1
  rownames(pick) <- NULL
  attr(pick, "age_ma") <- age_ma
  attr(pick, "k") <- attr(at, "k_kyr")
  pick
}

# Gaussian smoothing of the rows (period axis, sd in rows) and moving average of the columns
.wr_smooth_map <- function(P, sd_rows = 0, n_cols = 1) {
  if (sd_rows > 0) {
    kh <- max(1, ceiling(3 * sd_rows))
    kern <- exp(-((-kh):kh)^2 / (2 * sd_rows^2))
    kern <- kern / sum(kern)
    Pp <- rbind(matrix(0, kh, ncol(P)), P, matrix(0, kh, ncol(P)))
    out <- matrix(0, nrow(P), ncol(P))
    for (i in seq_along(kern)) out <- out + kern[i] * Pp[i:(i + nrow(P) - 1), , drop = FALSE]
    P <- out
  }
  if (n_cols > 1) {
    k <- rep(1 / n_cols, n_cols)
    P <- t(apply(P, 1, function(v) {
      f <- as.numeric(stats::filter(v, k, sides = 2))
      f[is.na(f)] <- v[is.na(f)]
      f
    }))
  }
  P
}

# Evidence map of a scalogram on a depth grid: the power, smoothed, normalised per position over the
# periods (so the average is one), optionally mixed with the normalised wavelet power as a fallback
.wr_evidence <- function(wavelet, cols, smooth_period, smooth_depth, fallback, fallback_weight, whiten = FALSE,
                         background = "none", normalise = "column", local_octaves = 0.5) {
  P <- wavelet$Power
  bgdiv <- function(M) M
  if (identical(background, "ar1")) {
    # divide by the theoretical AR1 (red noise) background of the record (Torrence and Compo 1998, eq. 16),
    # so that red noise has the same expected level at every period and cycles stand out as excess power
    yy <- wavelet$y
    r1 <- stats::acf(yy, lag.max = 1, plot = FALSE)$acf[2]
    r1 <- min(max(r1, 0), 0.99)
    bg <- (1 - r1^2) / (1 + r1^2 - 2 * r1 * cos(2 * pi * wavelet$dt / wavelet$Period))
    bgdiv <- function(M) M / bg
  }
  dj <- wavelet$dj
  dt <- wavelet$dt
  per <- wavelet$Period
  ncol_s <- max(1, round(smooth_depth / dt))
  norm_cols <- function(M) {
    if (identical(normalise, "local")) {
      # power relative to its local background along the period axis (Gaussian, sd local_octaves): removes the
      # red-noise slope and makes noise have an expected level of one at every period
      bgm <- .wr_smooth_map(M, local_octaves / dj, 1)
      tiny <- max(M) * 1e-8 + 1e-300
      return(M / (bgm + tiny))
    }
    if (identical(normalise, "median")) {
      # robust background: running median over +- local_octaves (not pulled up by the cycles themselves)
      h <- max(1, round(local_octaves / dj))
      bgm <- apply(M, 2, function(v) stats::runmed(v, 2 * h + 1, endrule = "constant"))
      bgm <- pmax(bgm, max(M) * 1e-6)
      return(M / bgm)
    }
    if (identical(normalise, "loglog")) {
      # smooth background: quadratic fit of log power against log period at every position (red-noise slope)
      lp <- log2(per)
      X <- cbind(1, lp, lp^2)
      L <- log(M + max(M) * 1e-8 + 1e-300)
      fit <- X %*% solve(crossprod(X), crossprod(X, L))
      return(exp(L - fit))
    }
    if (whiten) {
      # divide every period by its typical power along the record (removes the red-noise slope)
      rm <- apply(M, 1, function(v) mean(v))
      M <- M / pmax(rm, max(rm) * 1e-6)
    }
    M / matrix(pmax(colMeans(M), 1e-300), nrow(M), ncol(M), byrow = TRUE)
  }
  Ps <- bgdiv(.wr_smooth_map(P, smooth_period / dj, ncol_s)[, cols, drop = FALSE])
  ev <- norm_cols(Ps)
  if (fallback_weight > 0) {
    Fb <- NULL
    if (identical(fallback, "cwt")) {
      if (!is.null(wavelet$Power_cwt)) {
        Fb <- wavelet$Power_cwt
      } else {
        cw <- analyse_wavelet(data.frame(wavelet$x, wavelet$y), dj = dj, lowerPeriod = min(per),
                              upperPeriod = max(per) * 1.0001, omega_nr = 8)
        Fb <- cw$Power[seq_len(min(nrow(cw$Power), length(per))), , drop = FALSE]
        if (nrow(Fb) < length(per)) Fb <- rbind(Fb, matrix(0, length(per) - nrow(Fb), ncol(Fb)))
        if (ncol(Fb) != length(wavelet$axis.1)) {
          # the regridded record can differ by one sample: put the fallback on the positions of the scalogram
          Fb <- t(apply(Fb, 1, function(v) stats::approx(cw$axis.1, v, wavelet$axis.1, rule = 2)$y))
        }
      }
    } else if (is.list(fallback) && !is.null(fallback$Power)) {
      Fb <- fallback$Power
    }
    if (!is.null(Fb)) {
      Fb <- bgdiv(.wr_smooth_map(Fb, 0, ncol_s)[, cols, drop = FALSE])
      Fb[!is.finite(Fb)] <- 0
      ev <- (1 - fallback_weight) * ev + fallback_weight * norm_cols(Fb)
    }
  }
  ev
}

.wr_track_positions <- function(n, n_positions) {
  unique(round(seq(1, n, length.out = min(n, n_positions))))
}

#' @title Continuous tracking of a cycle in a (synchrosqueezed) scalogram with a hidden Markov model
#'
#' @description Tracks the period of one cycle continuously through a scalogram made by
#' \code{\link{analyse_wavelet}}, \code{\link{analyse_wavelet_sst}}, \code{\link{analyse_superlet_sst}},
#' \code{\link{analyse_conceft}} or any object of class \code{analyze.wavelet}. Sharp synchrosqueezed scalograms can have
#' gaps where the cycle is weak (low eccentricity, noise, a facies change), which makes manual or ridge-by-ridge tracking
#' discontinuous. Here the period at every position is a hidden state that follows a random walk (the sedimentation rate
#' changes gradually, with an optional small probability of an abrupt change, e.g. at a hiatus) and the scalogram is the
#' evidence. The forward-backward algorithm combines the evidence of the whole record: where the cycle is clear it fixes the
#' period, and through gaps the track is bridged by the neighbouring evidence and by the fallback scalogram (the ordinary
#' wavelet power, which is never empty). The result is a continuous curve with
#' \itemize{
#'   \item the most probable path (Viterbi),
#'   \item the posterior mean and a credible band of the period at every position,
#'   \item a confidence value (the posterior probability within 2 percent of the tracked period), which shows where the
#'   track is constrained by data and where it is interpolated,
#'   \item optionally (\code{period_of_tracked_cycle}) the sedimentation rate and an age model with uncertainty from paths
#'   sampled from the posterior.
#' }
#' The first two columns (depth, period) can be used directly by \code{\link{curve2time}} and the other functions that
#' take a tracked curve.
#'
#' @param wavelet Scalogram object (class \code{analyze.wavelet}).
#' @param period Band \code{c(lower, upper)} in which the cycle is searched (in depth units). Alternatively give
#' \code{guide}.
#' @param guide Optional rough track (data frame: depth, period), e.g. from \code{\link{track_period_wavelet}} or a
#' previous run; the search is restricted to \code{guide * (1 +- guide_tolerance)}.
#' @param guide_tolerance Relative width of the search band around the guide \code{Default=0.25}.
#' @param sd_change Standard deviation of the change of log2(period) over a distance of one period
#' \code{Default=0.1} (about 7 percent per cycle). Larger values follow faster changes of the sedimentation rate.
#' @param jump_prob Probability of an abrupt change of the period between neighbouring positions
#' \code{Default=1e-3}.
#' @param smooth_period Smoothing of the scalogram along the period axis before tracking (sd in octaves)
#' \code{Default=0.02}.
#' @param smooth_depth Moving average of the scalogram along the record before tracking (in depth units)
#' \code{Default=NULL}: the lower limit of the band (one shortest period).
#' @param fallback \code{"cwt"} (default) to mix in the normalised Morlet wavelet power (from the object or computed),
#' another scalogram object, or \code{"none"}.
#' @param fallback_weight Weight of the fallback evidence \code{Default=0.2}.
#' @param sharpness Weight of the evidence \code{Default=0.5}. Because neighbouring positions of a scalogram are not
#' independent, the log evidence is scaled so that one period of the cycle counts as one observation; larger values
#' trust the scalogram more relative to the smoothness of the track and give narrower credible bands; with 0.5 the 95
#' percent bands of the age model contained the true age 95 to 100 percent of the time in synthetic tests.
#' @param n_positions Number of positions at which the track is computed (it is interpolated to all positions)
#' \code{Default=1500}; fewer are used when the positions would be closer than a quarter of the (reference) period.
#' @param jump_depths Optional depths where an abrupt change of the sedimentation rate is plausible (a hiatus,
#' erosion surface, abrupt facies change or tectonic contact known from the field or from a change-point analysis such
#' as \code{\link{bayes_changepoint}}). Near these depths the jump probability is raised to \code{jump_prob_at}.
#' Abrupt changes are hard to detect from the scalogram alone in noisy records, so giving them helps the track to
#' follow the jump instead of smoothing it.
#' @param jump_prob_at Jump probability near \code{jump_depths} \code{Default=0.05}.
#' @param constraints Optional geological constraints: a data frame with the columns \code{from}, \code{to} (depths) and
#' \code{duration} (time, in the units of the astronomical periods, e.g. kyr), for instance a biostratigraphic bracket,
#' the time between two dated ash beds or the duration of a stage known from another section. The track is then the most
#' probable one that satisfies them: the evidence is tilted by a Lagrange multiplier times the time spent in each
#' interval, with the multiplier found by bisection so that the interval has the required duration. The achieved
#' durations and multipliers are returned in the attribute \code{constraints}.
#' @param period_of_tracked_cycle Optional duration of the tracked cycle (e.g. 405 kyr): adds the sedimentation rate
#' and an age model with uncertainty.
#' @param n_samples Number of posterior paths used for the uncertainty of the age model \code{Default=200}.
#' @param level Width of the credible bands \code{Default=0.95}.
#' @param seed Seed of the posterior sampling \code{Default=1}.
#'
#' @return A data frame of class \code{wr_track} with the columns \code{depth}, \code{period} (Viterbi path),
#' \code{period_mean}, \code{period_lower}, \code{period_upper}, \code{confidence} and, with
#' \code{period_of_tracked_cycle}, \code{sedrate} (depth units per time unit), \code{sedrate_lower}, \code{sedrate_upper},
#' \code{time}, \code{time_lower} and \code{time_upper} (time since the top/start of the record).
#'
#' @references
#' Rabiner, L.R., 1989. A tutorial on hidden Markov models and selected applications in speech recognition. Proceedings
#' of the IEEE 77, 257-286.
#'
#' Carmona, R.A., Hwang, W.L., Torresani, B., 1997. Characterization of signals by the ridges of their wavelet
#' transforms. IEEE Transactions on Signal Processing 45, 2586-2590.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' z <- 0:1999
#' Pt <- 40 + 10 * sin(2 * pi * z / 1500)
#' amp <- 0.3 + abs(sin(2 * pi * z / 600))
#' x <- amp * sin(2 * pi * cumsum(1 / Pt)) + rnorm(length(z), sd = 0.5)
#' s <- analyse_superlet_sst(data.frame(z, x), lowerPeriod = 15, upperPeriod = 120, dj = 1 / 50)
#' tr <- track_ridge_hmm(s, period = c(25, 70), period_of_tracked_cycle = 20)
#' plot_track(tr, s)
#' }
#' @export
track_ridge_hmm <- function(wavelet,
                            period = NULL,
                            guide = NULL,
                            guide_tolerance = 0.25,
                            sd_change = 0.1,
                            jump_prob = 1e-3,
                            smooth_period = 0.02,
                            smooth_depth = NULL,
                            fallback = "cwt",
                            fallback_weight = 0.2,
                            sharpness = 0.5,
                            n_positions = 1500,
                            period_of_tracked_cycle = NULL,
                            n_samples = 200,
                            level = 0.95,
                            jump_depths = NULL,
                            constraints = NULL,
                            jump_prob_at = 0.05,
                            seed = 1) {
  per <- wavelet$Period
  z <- wavelet$axis.1
  n <- length(z)
  dt <- wavelet$dt
  if (is.null(period) && is.null(guide)) stop("give a period band or a guide")
  if (!is.null(guide)) {
    gz <- stats::approx(guide[, 1], guide[, 2], z, rule = 2)$y
    lo_band <- min(gz) * (1 - guide_tolerance)
    hi_band <- max(gz) * (1 + guide_tolerance)
  } else {
    lo_band <- min(period)
    hi_band <- max(period)
  }
  # positions about a quarter of a period apart: neighbouring positions are then nearly independent and the
  # random walk can move by whole period bins between them
  pref <- exp(mean(log(c(lo_band, hi_band))))
  cols <- .wr_track_positions(n, min(n_positions, ceiling(diff(range(z)) / (pref / 4)) + 1))
  zc <- z[cols]
  if (!is.null(guide)) g <- gz[cols]
  rows <- which(per >= lo_band & per <= hi_band)
  if (length(rows) < 3) stop("the band holds fewer than 3 periods of the scalogram")
  if (is.null(smooth_depth)) smooth_depth <- lo_band
  ev <- .wr_evidence(wavelet, cols, smooth_period, smooth_depth, fallback, fallback_weight)
  ev <- ev[rows, , drop = FALSE]
  ev <- ev / matrix(pmax(colMeans(ev), 1e-300), nrow(ev), ncol(ev), byrow = TRUE)
  logE <- log(ev + 0.05)
  if (!is.null(guide)) {
    lg <- log2(per[rows])
    outside <- abs(outer(lg, log2(g), "-")) > log2(1 + guide_tolerance)
    logE[outside] <- logE[outside] - 50
  }
  # random walk: sd of the change of log2(period) per step of the track, scaled from one period
  step <- if (length(cols) > 1) stats::median(diff(zc)) else dt
  sig <- sd_change * sqrt(step / pref) / wavelet$dj
  # neighbouring positions see (nearly) the same part of the record: temper the evidence so that one
  # period (the correlation length of the scalogram) counts as one independent observation
  logE <- logE * sharpness * min(1, step / max(pref, smooth_depth))
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  rho <- .wr_jump_rho(zc, jump_prob, jump_depths, jump_prob_at)
  if (!is.null(constraints) && is.null(period_of_tracked_cycle)) stop("duration constraints need period_of_tracked_cycle")
  inv_rate <- if (is.null(period_of_tracked_cycle)) rep(1, length(rows)) else period_of_tracked_cycle / per[rows]
  hc <- .wr_hmm_constrained(logE, max(sig, 0.3), rho, if (is.null(period_of_tracked_cycle)) 0 else n_samples,
                            inv_rate, zc, constraints)
  h <- hc$h
  lp <- log2(per[rows])
  pm <- colSums(h$post * lp)
  cdf <- apply(h$post, 2, cumsum)
  qf <- function(p) apply(cdf, 2, function(cc) lp[min(which(cc >= p))])
  a <- (1 - level) / 2
  vit <- lp[h$path]
  conf <- sapply(seq_along(cols), function(t) sum(h$post[abs(lp - vit[t]) <= log2(1.02), t]))
  out <- data.frame(depth = z,
                    period = 2^stats::approx(zc, vit, z, rule = 2)$y,
                    period_mean = 2^stats::approx(zc, pm, z, rule = 2)$y,
                    period_lower = 2^stats::approx(zc, qf(a), z, rule = 2)$y,
                    period_upper = 2^stats::approx(zc, qf(1 - a), z, rule = 2)$y,
                    confidence = stats::approx(zc, conf, z, rule = 2)$y)
  if (!is.null(period_of_tracked_cycle)) {
    out <- .wr_track_age(out, out$period / period_of_tracked_cycle, zc,
                         2^lp[h$samples] / period_of_tracked_cycle, level)
  }
  attr(out, "band") <- c(lo_band, hi_band)
  attr(out, "jump_depths") <- jump_depths
  attr(out, "constraints") <- hc$constraints
  attr(out, "period_of_tracked_cycle") <- period_of_tracked_cycle
  class(out) <- c("wr_track", "data.frame")
  out
}

# Sedimentation rate and age model (with uncertainty from posterior paths) of a track
.wr_track_age <- function(out, sr, zc, sr_samples, level) {
  z <- out$depth
  a <- (1 - level) / 2
  integ <- function(zz, s) c(0, cumsum(diff(zz) * (1 / s[-1] + 1 / s[-length(s)]) / 2))
  out$sedrate <- sr
  sr_s <- matrix(sr_samples, nrow = length(zc))
  out$sedrate_lower <- stats::approx(zc, apply(sr_s, 1, stats::quantile, a), z, rule = 2)$y
  out$sedrate_upper <- stats::approx(zc, apply(sr_s, 1, stats::quantile, 1 - a), z, rule = 2)$y
  out$time <- integ(z, sr)
  tt <- apply(sr_s, 2, function(s) integ(zc, s))
  out$time_lower <- stats::approx(zc, apply(tt, 1, stats::quantile, a), z, rule = 2)$y
  out$time_upper <- stats::approx(zc, apply(tt, 1, stats::quantile, 1 - a), z, rule = 2)$y
  out$time_median <- stats::approx(zc, apply(tt, 1, stats::median), z, rule = 2)$y
  out
}

#' @title Astronomically constrained continuous tracking of the sedimentation rate in a scalogram
#'
#' @description Tracks the sedimentation rate (and with it the periods of all astronomical cycles) continuously through a
#' scalogram, using all cycles at once. At every position the hidden state is the sedimentation rate s; the evidence
#' for a value of s is the (normalised, smoothed) power of the scalogram at the depth periods \code{P_k * s} of all
#' astronomical cycles \code{P_k} together (e.g. eccentricity, obliquity and precession). Where one cycle is weak (a node
#' of the eccentricity modulation of precession, a gap in a synchrosqueezed ridge) the others keep the track on course,
#' and the fixed ratios between the periods prevent the confusion of one cycle with another (e.g. an obliquity ridge read
#' as precession) that single ridge tracking suffers from. The sedimentation rate follows a random walk with an optional
#' small probability of abrupt changes. The forward-backward algorithm gives a continuous most probable sedimentation rate
#' curve with a credible band, and paths sampled from the posterior give an age model with uncertainty.
#'
#' Performance in synthetic tests (3000 kyr; eccentricity-modulated precession, obliquity and eccentricity; three jumps
#' of the sedimentation rate by factors 1.7, 0.5 and 1.6 plus gradual changes; AR1 noise; scalogram from
#' \code{analyse_wavelet_sst(omega_nr = 12, n_iter = 4)}): the largest age error was about 2 percent of the duration with
#' noise half as strong as the signal, 3-4 percent with equally strong noise and 6-9 percent with noise twice as strong
#' (always about half of the error of picking the strongest precession power in a band). With noise three times as strong
#' as the signal the track only succeeded with a rough guide (about 7 percent); without one it locked onto the wrong
#' cycles. The credible bands then also become too narrow: they describe the uncertainty given that the track follows the
#' right cycles, which the method cannot guarantee at such noise levels. Check weak intervals (low \code{confidence}) by
#' eye and with independent methods (e.g. \code{\link{ecoco_classic}}, amplitude modulation tests).
#'
#' The astronomical periods must fit the age of the record: for deep time use precession and obliquity periods that
#' account for the recession of the Moon (e.g. from \code{\link{astro_terms}}, a k value from \code{\link{test_tuned}} or a
#' published solution).
#'
#' @param wavelet Scalogram object (class \code{analyze.wavelet}), see \code{\link{track_ridge_hmm}}.
#' @param astro_periods Periods of the astronomical cycles in time
#' \code{Default=c(405, 125, 95, 53.6, 41, 39.7, 23.7, 22.4, 19)} (kyr, Cenozoic values). Listing all the main terms
#' (including the weaker obliquity and precession terms) helps to tell the cycles apart.
#' @param weights Weights of the cycles in the evidence, best proportional to their expected amplitude. Default: for
#' the default periods \code{c(0.5, 0.3, 0.3, 0.3, 1, 0.4, 1, 0.8, 0.6)}, otherwise 1 for all.
#' @param sedrate Range \code{c(lower, upper)} of the sedimentation rate (depth units per time unit of
#' \code{astro_periods}) in which the track is searched.
#' @param whiten Divide the power of every period by its average along the record before tracking (removes the
#' red-noise slope) \code{Default=FALSE} (in tests this made the tracking worse, because noise at the long periods is amplified).
#' @param guide Optional rough sedimentation rate curve (data frame: depth, sedimentation rate), e.g. from a manual
#' tracking divided by the period of the tracked cycle; the search is restricted to \code{guide * (1 +- guide_tolerance)}.
#' In very noisy records this prevents the track from locking onto the wrong cycle.
#' @param guide_tolerance Relative width of the search band around the guide \code{Default=0.25}.
#' @param floor Value added to the normalised evidence (average one) before taking the logarithm. It limits the penalty
#' for a cycle that is absent at a position \code{Default=0.5}.
#' @param background \code{"none"} or \code{"ar1"}: evidence relative to an AR1 red-noise background \code{Default="none"}.
#' @param sd_change Standard deviation of the change of log2(sedimentation rate) over a distance of one cycle of the
#' shortest astronomical period \code{Default=0.15}.
#' @param jump_prob Probability of an abrupt change of the sedimentation rate between neighbouring positions
#' \code{Default=1e-4}; small, because a jump of the track by the ratio of two astronomical periods (e.g. from
#' precession to obliquity) is the typical failure of tracking.
#' @param smooth_period,smooth_depth,fallback,fallback_weight,sharpness,n_positions,n_samples,level,seed See
#' \code{\link{track_ridge_hmm}}; \code{smooth_depth} defaults to twice the depth period of the shortest cycle at the
#' (geometric) mean sedimentation rate. For this tracker a scalogram in which weak cycles are kept (the Morlet scalogram with omega_nr
#' around 12, or \code{\link{analyse_wavelet_sst}} with \code{n_iter = 4}) worked better in tests than the consensus of
#' \code{\link{analyse_superlet_sst}}, which removes weak cycles.
#'
#' @param jump_depths,jump_prob_at,constraints See \code{\link{track_ridge_hmm}}.
#' @param long_window,long_weight Long-window evidence (similar to eCOCO): the wavelet power averaged over windows of
#' \code{long_window} depth units (default eight cycles of the shortest period at the mean sedimentation rate), with the
#' red-noise slope removed, is correlated with the astronomical template and added with weight \code{long_weight}
#' \code{Default=0} (off). In stress tests this term kept the track on the right cycles when the noise was three times the
#' signal (largest age error about 400 instead of 700 kyr in 3000 kyr) but made the age model less accurate at lower
#' noise, because the averaged spectra are biased where the sedimentation rate changes within a window. Use it (e.g.
#' \code{long_weight = 1}) for very noisy records, or through \code{\link{track_astro_auto}} as an independent check.
#' @param local_weight Weight of the local (scalogram) evidence \code{Default=1}; 0 with \code{long_weight > 0} gives a
#' pure long-window estimate.
#' @param guide_penalty Penalty (log-likelihood per octave outside the band, plus 10 percent) for sedimentation rates
#' outside \code{guide * (1 +- guide_tolerance)} \code{Default=50}, i.e. practically a hard limit; smaller values make
#' the guide soft.
#' @param normalise How the scalogram is turned into evidence: \code{"median"} (default) divides the power by its running
#' median along the period axis (over +- \code{local_octaves}), which removes the red-noise slope without being pulled
#' up by the cycles themselves, and scores every cycle relative to the average evidence at that position, so that cycles
#' outside the period range of the scalogram neither help nor hurt; \code{"column"} (the behaviour before) divides by the
#' mean power at each position, \code{"local"} by a Gaussian average and \code{"loglog"} by a quadratic fit of log power
#' against log period. In synthetic tests \code{"column"} and \code{"median"} were best; in six real records
#' \code{"column"} was biased towards sedimentation rates that put the long cycles on the red-noise slope or the short
#' cycles outside the resolved range, and \code{"median"} gave the smallest errors when an approximate mean sedimentation
#' rate was given.
#' @param local_octaves Half width (octaves) of the background window of \code{normalise} \code{Default=1}.
#' @param consistent_cycles \code{TRUE}: use only the cycles that lie inside the period range of the scalogram at
#' every sedimentation rate of \code{sedrate}; \code{"partial"} (default): cycles that leave the period range at one
#' end of the rate range (typically precession below the sampling resolution at the lowest rates) are kept, but only
#' the changes of their evidence along the record are used (row means removed), so they show where the rate changes
#' without favouring the rates at which they are resolved; \code{FALSE}: all cycles, all evidence. In CIP Case 3
#' (precession at 3-5 samples per cycle) \code{TRUE} silently dropped the precession cycles and the 1.2-Myr
#' oscillation of the rate was lost; \code{"partial"} keeps them. If the set of usable cycles changed with the sedimentation
#' rate, rates at which more cycles fall inside the resolved range would be favoured: in the Klonk core, which resolves
#' only eccentricity, the ~110-kyr cycle was read as obliquity and the 405-kyr cycle as ~125 kyr. Give a realistic
#' \code{sedrate} range (or a guide) so that the cycles you expect stay in.
#' @param smooth_path Optional L1 smoothing of the tracked path (Viterbi path) with weights from the posterior
#' (majorisation-minimisation; Wang et al. 2023): \code{"tv"} penalises the first differences of the log sedimentation
#' rate (total variation: a piecewise constant rate with sharp jumps is favoured), \code{"trend"} the second
#' differences (piecewise linear). \code{Default="none"}.
#' @param smooth_lambda Strength of the smoothing: the penalty per octave of change of the log2 sedimentation rate
#' (\code{"tv"}) or per octave of change of its slope (\code{"trend"}), against the weighted squared distance from
#' the tracked path \code{Default=1}.
#' @param age_ma,k_yr Optional age (Ma) or precession constant: the astronomical periods and weights are then taken from
#' \code{\link{astro_periods_tracking}} (deep-time obliquity and precession periods), replacing \code{astro_periods}.
#' @return A data frame of class \code{wr_track} with \code{depth}, \code{sedrate} (Viterbi path), \code{sedrate_mean},
#' \code{sedrate_lower}, \code{sedrate_upper}, \code{confidence}, \code{time}, \code{time_lower}, \code{time_upper},
#' and a column \code{period_<P>} with the depth period of every astronomical cycle. The attribute
#' \code{astro_periods} holds the periods.
#'
#' @references
#' Rabiner, L.R., 1989. A tutorial on hidden Markov models and selected applications in speech recognition. Proceedings
#' of the IEEE 77, 257-286.
#'
#' Meyers, S.R., 2015. The evaluation of eccentricity-related amplitude modulation and bundling in paleoclimate data: an
#' inverse approach for astrochronologic testing and time scale optimization. Paleoceanography 30, 1625-1640.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1500, by = 0.5)
#' ecc <- 0.4 * cos(2 * pi * t / 405) + 0.25 * cos(2 * pi * t / 125)
#' x_t <- (1 + ecc) * (sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19)) +
#'   0.7 * sin(2 * pi * t / 41)
#' sr <- 2 + 0.8 * sin(2 * pi * t / 600)                 # sedimentation rate (cm/kyr)
#' zt <- c(0, cumsum(sr[-1] * diff(t)))
#' z <- seq(0, max(zt), by = 1)
#' x <- approx(zt, x_t, z)$y + rnorm(length(z), sd = 1)
#' s <- analyse_superlet_sst(data.frame(z, x), lowerPeriod = 15, upperPeriod = 400, dj = 1 / 50)
#' tr <- track_astro_hmm(s, astro_periods = c(405, 125, 41, 23.7, 19), sedrate = c(0.8, 5))
#' plot_track(tr, s)
#' }
#' @export
track_astro_hmm <- function(wavelet,
                            astro_periods = c(405, 125, 95, 53.6, 41, 39.7, 23.7, 22.4, 19),
                            weights = NULL,
                            sedrate = NULL,
                            guide = NULL,
                            guide_tolerance = 0.25,
                            floor = 0.5,
                            whiten = FALSE,
                            background = c("none", "ar1"),
                            sd_change = 0.15,
                            jump_prob = 1e-4,
                            smooth_period = 0.02,
                            smooth_depth = NULL,
                            fallback = "cwt",
                            fallback_weight = 0.2,
                            sharpness = 0.5,
                            n_positions = 1500,
                            n_samples = 200,
                            level = 0.95,
                            jump_depths = NULL,
                            jump_prob_at = 0.05,
                            long_window = NULL,
                            long_weight = 0,
                            local_weight = 1,
                            guide_penalty = 50,
                            normalise = c("median", "column", "local", "loglog"),
                            local_octaves = 1,
                            consistent_cycles = "partial",
                            constraints = NULL,
                            smooth_path = c("none", "tv", "trend"),
                            smooth_lambda = 1,
                            age_ma = NULL,
                            k_yr = NULL,
                            seed = 1) {
  background <- match.arg(background)
  smooth_path <- match.arg(smooth_path)
  normalise <- match.arg(normalise)
  if (!is.null(age_ma) || !is.null(k_yr)) {
    astro_periods <- astro_periods_tracking(if (is.null(age_ma)) 0 else age_ma, k_yr)
  }
  if (is.data.frame(astro_periods)) {
    if (is.null(weights) && !is.null(astro_periods$weight)) weights <- astro_periods$weight
    astro_periods <- astro_periods$period
  }
  per <- wavelet$Period
  z <- wavelet$axis.1
  n <- length(z)
  dj <- wavelet$dj
  if (is.null(weights)) {
    default_p <- c(405, 125, 95, 53.6, 41, 39.7, 23.7, 22.4, 19)
    default_w <- c(0.5, 0.3, 0.3, 0.3, 1, 0.4, 1, 0.8, 0.6)
    weights <- if (identical(astro_periods, default_p)) default_w else rep(1, length(astro_periods))
  }
  if (is.null(sedrate)) {
    # every state must place at least one cycle inside the scalogram
    sedrate <- c(min(per) / max(astro_periods), max(per) / min(astro_periods))
  }
  partial_mode <- identical(consistent_cycles, "partial")
  partial <- rep(FALSE, length(astro_periods))
  if (partial_mode) {
    # cycles inside the period range of the scalogram over the whole sedimentation-rate range contribute fully; cycles
    # that leave the range at one end (e.g. precession below the sampling resolution at the lowest rates) contribute
    # only the changes of their evidence along the record (row means removed), so they still show where the rate
    # changes but cannot favour the rates at which they happen to be resolved
    inside_all <- astro_periods * sedrate[1] >= min(per) & astro_periods * sedrate[2] <= max(per)
    inside_mid <- astro_periods * sqrt(prod(sedrate)) >= min(per) & astro_periods * sqrt(prod(sedrate)) <= max(per)
    keep <- inside_all | inside_mid
    partial <- keep & !inside_all
    if (any(partial)) message("track_astro_hmm: cycles used for changes of the rate only (not resolved over the whole ",
                              "sedimentation-rate range): ", paste(signif(astro_periods[partial], 4), collapse = ", "))
    if (any(!keep)) message("track_astro_hmm: cycles not used: ", paste(signif(astro_periods[!keep], 4), collapse = ", "))
    dropped <- astro_periods[!keep]
    astro_used <- astro_periods[keep]
    weights_used <- weights[keep]
    partial <- partial[keep]
  } else if (isTRUE(consistent_cycles)) {
    # use only cycles that lie inside the period range of the scalogram at every sedimentation rate of the search
    # range: otherwise states that move more (heavily weighted) cycles into the resolved range are favoured, e.g. a
    # sedimentation rate at which the 110-kyr cycle is read as obliquity in a record that cannot resolve obliquity
    keep <- astro_periods * sedrate[1] >= min(per) & astro_periods * sedrate[2] <= max(per)
    if (sum(keep) < 2) {
      warning("fewer than two astronomical cycles stay inside the period range of the scalogram over the whole ",
              "sedimentation-rate range; narrow the range (e.g. with a guide or duration) or set consistent_cycles = FALSE")
      keep <- astro_periods * sqrt(prod(sedrate)) >= min(per) & astro_periods * sqrt(prod(sedrate)) <= max(per)
    }
    if (any(!keep)) message("track_astro_hmm: cycles not used (not resolved over the whole sedimentation-rate range): ",
                            paste(signif(astro_periods[!keep], 4), collapse = ", "))
    dropped <- astro_periods[!keep]
    astro_used <- astro_periods[keep]
    weights_used <- weights[keep]
    partial <- partial[keep]
  } else {
    dropped <- numeric(0)
    astro_used <- astro_periods
    weights_used <- weights
  }
  ls <- seq(log2(sedrate[1]), log2(sedrate[2]), by = dj)
  sgrid <- 2^ls
  pref <- min(astro_used) * exp(mean(log(sedrate)))
  cols <- .wr_track_positions(n, min(n_positions, ceiling(diff(range(z)) / (pref / 4)) + 1))
  zc <- z[cols]
  if (is.null(smooth_depth)) smooth_depth <- 2 * pref
  smooth_depth <- max(smooth_depth, wavelet$dt)
  ev <- .wr_evidence(wavelet, cols, smooth_period, smooth_depth, fallback, fallback_weight, whiten, background, normalise, local_octaves)
  lper <- log2(per)
  logE <- matrix(0, length(ls), length(cols))
  # baseline: the average log evidence of a period at each position; a cycle contributes its evidence relative to
  # it, and a cycle outside the period range of the scalogram contributes nothing (no bias towards placing cycles
  # outside the range)
  base <- if (normalise != "column") colMeans(log(ev + floor)) else rep(log(1 + floor), ncol(ev))
  maps <- vector("list", length(astro_used))
  for (k in seq_along(astro_used)) {
    target <- log2(astro_used[k]) + ls
    inside <- target >= min(lper) & target <= max(lper)
    if (!any(inside)) next
    # interpolate the evidence at the depth period of cycle k for every sedimentation rate
    fi <- (target[inside] - lper[1]) / dj + 1
    i0 <- pmin(floor(fi), length(per) - 1)
    w1 <- fi - i0
    e <- (1 - w1) * ev[i0, , drop = FALSE] + w1 * ev[i0 + 1, , drop = FALSE]
    m <- matrix(0, length(ls), length(cols))
    m[inside, ] <- log(e + floor) - matrix(base, sum(inside), ncol(e), byrow = TRUE)
    if (partial[k]) m[inside, ] <- m[inside, ] - rowMeans(m[inside, , drop = FALSE])
    maps[[k]] <- m
  }
  for (k in seq_along(astro_used)) if (!is.null(maps[[k]])) logE <- logE + weights_used[k] * maps[[k]]
  step <- if (length(cols) > 1) stats::median(diff(zc)) else wavelet$dt
  sig <- sd_change * sqrt(step / pref) / dj
  logE <- local_weight * logE * sharpness * min(1, step / max(pref, smooth_depth))
  if (long_weight > 0) {
    if (is.null(long_window)) long_window <- 8 * pref
    C <- .wr_long_evidence(wavelet, zc, long_window, ls, astro_used, weights_used)
    logE <- logE + long_weight * 30 * C * min(1, step / long_window)
  }
  if (!is.null(guide)) {
    g <- stats::approx(guide[, 1], guide[, 2], zc, rule = 2)$y
    dist <- abs(outer(ls, log2(g), "-")) - log2(1 + guide_tolerance)
    outside <- dist > 0
    # penalty grows with the distance outside the band (guide_penalty per octave, plus a step)
    logE[outside] <- logE[outside] - guide_penalty * (0.1 + dist[outside])
  }
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  if (isTRUE(getOption("WaverideR.track_evidence"))) {
    return(list(logE = logE, ls = ls, zc = zc, astro_used = astro_used, weights_used = weights_used, partial = partial))
  }
  rho <- .wr_jump_rho(zc, jump_prob, jump_depths, jump_prob_at)
  hc <- .wr_hmm_constrained(logE, max(sig, 0.3), rho, n_samples, 1 / sgrid, zc, constraints)
  h <- hc$h
  vit <- ls[h$path]
  pm <- colSums(h$post * ls)
  if (smooth_path != "none") {
    w <- sapply(seq_along(cols), function(t) sum(h$post[abs(ls - vit[t]) <= log2(1.1), t]))
    vit <- .wr_l1_smooth(vit, w, smooth_lambda, order = if (smooth_path == "tv") 1 else 2)
  }
  conf <- sapply(seq_along(cols), function(t) sum(h$post[abs(ls - vit[t]) <= log2(1.02), t]))
  conf10 <- sapply(seq_along(cols), function(t) sum(h$post[abs(ls - vit[t]) <= log2(1.1), t]))
  out <- data.frame(depth = z,
                    sedrate = 2^stats::approx(zc, vit, z, rule = 2)$y,
                    sedrate_mean = 2^stats::approx(zc, pm, z, rule = 2)$y,
                    confidence = stats::approx(zc, conf, z, rule = 2)$y,
                    confidence10 = stats::approx(zc, conf10, z, rule = 2)$y)
  out <- .wr_track_age(out, out$sedrate, zc, sgrid[h$samples], level)
  out <- out[, c("depth", "sedrate", "sedrate_mean", "sedrate_lower", "sedrate_upper", "confidence", "confidence10",
                 "time", "time_median", "time_lower", "time_upper")]
  for (k in seq_along(astro_periods)) out[[paste0("period_", astro_periods[k])]] <- out$sedrate * astro_periods[k]
  attr(out, "astro_periods") <- astro_periods
  attr(out, "cycles_not_used") <- dropped
  attr(out, "cycle_weights") <- data.frame(period = astro_used, weight = weights_used, partial = partial)
  attr(out, "posterior") <- list(depth = zc, sedrate = sgrid, post = h$post, evidence = logE)
  attr(out, "constraints") <- hc$constraints
  attr(out, "jump_depths") <- jump_depths
  class(out) <- c("wr_track", "data.frame")
  out
}

#' @title Plot a continuous track on its scalogram
#'
#' @description Plots the scalogram with the track(s) of \code{\link{track_ridge_hmm}} or \code{\link{track_astro_hmm}}
#' and their credible bands, the confidence of the track and, when available, the sedimentation rate and the age model
#' with its uncertainty.
#'
#' @param track Output of \code{\link{track_ridge_hmm}} or \code{\link{track_astro_hmm}}.
#' @param wavelet The scalogram the track was made from (or another one to display).
#' @param show Periods of the astronomical cycles to draw (for \code{\link{track_astro_hmm}}) \code{Default=NULL}: all.
#' @param flags Optional data frame of intervals (columns \code{from}, \code{to}) to shade, e.g. the \code{flags} of
#' \code{\link{track_astro_auto}}.
#' @param palette Colour palette of \code{\link[grDevices]{hcl.colors}} \code{Default="viridis"}.
#' @param scale Colour scale as in all image plots of the package (see \code{\link{plot_wavelet}}) \code{Default="auto"}:
#' log for synchrosqueezed maps, quantiles of the power divided by its mean per period otherwise.
#' @param keep_editable Keep the plot settings after plotting \code{Default=FALSE}.
#' @return Invisibly the track.
#' @export
plot_track <- function(track, wavelet, show = NULL, flags = NULL, palette = "viridis", keep_editable = FALSE,
                       scale = "auto") {
  if (!keep_editable) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
  }
  has_age <- !is.null(track$time)
  graphics::layout(matrix(seq_len(if (has_age) 3 else 2), ncol = 1), heights = if (has_age) c(3, 1, 1.4) else c(3, 1))
  graphics::par(mar = c(1, 4.5, 1, 1), oma = c(3, 0, 0, 0), mgp = c(2.4, 0.7, 0))
  P <- wavelet$Power
  if (!isTRUE(wavelet$synchrosqueezed)) P <- P / rowMeans(P, na.rm = TRUE)       # flatten the red-noise slope
  .wr_image(wavelet$axis.1, log2(wavelet$Period), P, scale = scale, palette = palette,
            sparse = if (isTRUE(wavelet$synchrosqueezed)) TRUE else NULL, xlab = "", ylab = "period", yaxt = "n", xaxt = "n")
  at <- pretty(log2(wavelet$Period))
  graphics::axis(2, at, signif(2^at, 3), las = 1)
  addband <- function(lo, hi, mid, col) {
    graphics::polygon(c(track$depth, rev(track$depth)), log2(c(lo, rev(hi))), col = grDevices::adjustcolor(col, 0.25), border = NA)
    graphics::lines(track$depth, log2(mid), col = col, lwd = 1.5)
  }
  ap <- attr(track, "astro_periods")
  if (is.null(ap) && !is.null(show) && !is.null(track$sedrate)) ap <- show      # e.g. a chirp-matched track
  if (!is.null(ap) && !is.null(show)) {                # show: the periods to draw, matched within 10 %
    m <- vapply(show, function(p) { d <- abs(log(ap / p)); if (min(d) < log(1.1)) ap[which.min(d)] else p }, numeric(1))
    ap <- unique(m)
  }
  if (is.null(track$sedrate_lower) && !is.null(track$sedrate)) { track$sedrate_lower <- track$sedrate; track$sedrate_upper <- track$sedrate }
  if (!is.null(ap) && length(ap)) {
    cols <- rep(c("#FFFFFF", "#FF5A5F", "#FFB000", "#00D1FF", "#FF4DFF", "#7CFF4D"), length.out = length(ap))
    for (k in seq_along(ap)) addband(track$sedrate_lower * ap[k], track$sedrate_upper * ap[k], track$sedrate * ap[k], cols[k])
  } else if (!is.null(track$period)) {
    addband(track$period_lower, track$period_upper, track$period, "red")
  }
  if (isTRUE(flags)) flags <- attr(track, "flags")
  if (!is.data.frame(flags)) flags <- NULL
  shade <- function() if (!is.null(flags) && nrow(flags)) {
    u <- graphics::par("usr")
    graphics::rect(flags$from, u[3], flags$to, u[4], col = grDevices::adjustcolor("orange", 0.25), border = NA)
  }
  shade()
  graphics::plot(track$depth, track$confidence, type = "l", ylim = c(0, 1), xaxs = "i", xlab = "", ylab = "confidence",
                 xaxt = if (has_age) "n" else "s")
  shade()
  if (has_age) {
    graphics::plot(track$depth, track$time, type = "n", xaxs = "i", ylab = "time", xlab = "",
                   ylim = range(c(track$time_lower, track$time_upper)))
    shade()
    graphics::polygon(c(track$depth, rev(track$depth)), c(track$time_lower, rev(track$time_upper)),
                      col = "grey80", border = NA)
    graphics::lines(track$depth, track$time, lwd = 1.5)
  }
  graphics::mtext("depth", side = 1, line = 1.5, outer = TRUE)
  invisible(track)
}

#' @title Recommended workflow for continuous astronomical tracking of noisy records
#'
#' @description Runs the workflow that worked best in stress tests (changing sedimentation rates with jumps, red noise):
#' \enumerate{
#'   \item a scalogram that keeps weak cycles: \code{\link{analyse_wavelet_sst}} with \code{omega_nr = 12} and four
#'   reassignment iterations (\code{scalogram = "msst"}, default), or the sharper superlet synchrosqueezing
#'   (\code{\link{analyse_superlet_sst}}, \code{scalogram = "superlet_sst"}), which was more accurate in records with
#'   moderate noise but drops weak cycles in very noisy records;
#'   \item astronomically constrained tracking with all cycles (\code{\link{track_astro_hmm}}), with deep-time periods
#'   from \code{age_ma} or \code{k_yr}, a realistic sedimentation-rate range, and optionally a rough guide and the depths
#'   of known hiatuses or abrupt facies changes;
#'   \item an independent check: the sedimentation rate is also estimated from long windows only (an eCOCO-like
#'   correlation of red-noise-corrected spectra with the astronomical template), and intervals where the posterior of
#'   the track is spread are flagged. Check flagged intervals by eye, with \code{\link{ecoco_classic}} or with amplitude
#'   modulation tests.
#' }
#'
#' Tests on real records (ODP 1119, Orcadian Basin 12/13-1, KTG-01, Klonk, Bisciaro, Pod Barrandovem) showed that the
#' tracking is only reliable with an approximate mean sedimentation rate (e.g. from the approximate duration known from
#' biostratigraphy or radiometric ages), given as \code{guide} with a wide \code{guide_tolerance}; without it the track
#' often settled on a sedimentation rate that was about twice too high or too low (the wrong cycles). In records where only
#' the eccentricity cycles are resolved, the result followed such a prior closely, i.e. the record itself did not
#' constrain the mean rate. The \code{duration} argument sets such a prior directly.
#'
#' @param data Record: first column depth, second column proxy.
#' @param sedrate Range \code{c(lower, upper)} of the sedimentation rate (depth units per kyr).
#' @param age_ma,k_yr Age (Ma) or precession constant for the astronomical periods (see
#' \code{\link{astro_periods_tracking}}) \code{Default=0} (present-day periods).
#' @param astro_periods Optional data frame (term, family, period, weight) or vector of periods replacing the ones from
#' \code{age_ma}.
#' @param scalogram \code{"msst"} (default), \code{"superlet_sst"}, \code{"chirplet_sst"} (matched chirped wavelets,
#' \code{\link{analyse_chirplet_sst}}; sharper where the sedimentation rate changes quickly, about 4 times slower) or
#' \code{"cwt"} (a plain Morlet scalogram with \code{omega_nr}). In records in which only the eccentricity cycles are resolved, a CWT with a low \code{omega_nr}
#' (5-6) followed the ridges more continuously than the synchrosqueezed maps (Klonk, Pridoli).
#' @param omega_nr Morlet \code{omega_nr} for \code{scalogram = "cwt"} \code{Default=6}.
#' @param constraints Optional duration constraints (data frame \code{from}, \code{to}, \code{duration} in kyr and
#' optionally \code{tolerance}, default 0.025), e.g. from radiometric ages or well-dated stage boundaries, see
#' \code{\link{track_astro_hmm}}. Recommended whenever available: with the Pridoli duration (Klonk) the published age
#' model was reproduced within about 0.3 Myr; without it the track settled on other cycles.
#' @param guide,guide_tolerance,jump_depths Optional prior information, see \code{\link{track_astro_hmm}}.
#' @param duration Approximate duration of the record (kyr), e.g. from biostratigraphy or radiometric ages. It is turned
#' into a guide with the constant mean sedimentation rate and a tolerance of \code{duration_tolerance} (default 0.5,
#' i.e. the rate may differ by a factor 1.5 from the mean). Strongly recommended for real records.
#' @param duration_tolerance Relative tolerance of the \code{duration} guide \code{Default=0.5}.
#' @param long_weight Weight of the long-window evidence inside the main track \code{Default=0}; set to 1 for records in
#' which the cycles are barely visible.
#' @param flag_ratio Relative difference between the track and the long-window estimate above which an interval is
#' flagged when \code{flag_long = TRUE} \code{Default=0.25}.
#' @param min_confidence Intervals where less than this posterior probability lies within 10 percent of the tracked
#' sedimentation rate (\code{confidence10}) are flagged \code{Default=0.5}.
#' @param flag_long Also flag intervals where the long-window estimate disagrees \code{Default=FALSE}. In real records
#' the long-window estimate was itself often biased, so its disagreement is reported in \code{check} but not used for the
#' flags by default.
#' @param dj Spacing of the period axis \code{Default=1/50}.
#' @param method \code{"both"} (default): the scalogram is tracked with \code{\link{track_astro_hmm}}
#' (\code{track_hmm}) and, in addition, chirp-matched windowed spectra (\code{\link{analyse_chirp_eha}}, automatic
#' windows and gradient grid) are tracked with \code{\link{track_chirp_hmm}} with the scalogram evidence added for
#' locality; this chirp-matched track is returned as \code{track}, \code{check} holds both rates and their
#' disagreement, and \code{agreement} the fraction of the record where they are within 10 percent. \code{"chirp"}
#' and \code{"hmm"} run only one of them. In the tests the chirp-matched track was better or equal everywhere:
#' CIP Case 3 with the rate unknown over 1-20 cm/kyr, correlation with the true rate 0.97 against 0.08 for the HMM
#' track (92 against 40 percent within 10 percent); CIP Case 2, both 0.93-0.95; synthetic records with three abrupt
#' jumps, 0.86-0.94 against 0.00-0.92. Where the two agree both are probably right; where they disagree, look at
#' \code{\link{rate_gram}} and \code{\link{plot_astro_evidence}}. The chirp-matched tracking adds 1-2 minutes.
#' @param jump_prob Jump probability of the chirp-matched tracker, see \code{\link{track_chirp_hmm}}
#' \code{Default=0}. Tested on all records: with \code{1e-3} the age errors of two of the four synthetic records
#' with jumps halved (57 and 70 against 123 and 215 kyr) and those of the other two doubled, and Case 3 collapsed
#' (0.97 to 0.10), so it stays off unless hiatuses or abrupt facies changes are known.
#' @param ... Further arguments for \code{\link{track_astro_hmm}}.
#'
#' @return A list with \code{scalogram}, \code{track} (see \code{\link{track_astro_hmm}}), \code{check} (data frame:
#' depth, sedrate of the track, long-window sedrate, ratio, confidence, flagged) and \code{flags} (flagged intervals:
#' from, to). Plot with \code{plot_track(x$track, x$scalogram, flags = x$flags)}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1500, by = 0.5)
#' ecc <- 0.4 * cos(2 * pi * t / 405) + 0.25 * cos(2 * pi * t / 125)
#' x_t <- (1 + ecc) * (sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19)) +
#'   0.7 * sin(2 * pi * t / 41)
#' sr <- ifelse(t < 700, 1.5, 2.6)                       # a jump of the sedimentation rate
#' zt <- c(0, cumsum(sr[-1] * diff(t)))
#' z <- seq(0, max(zt), by = 1)
#' x <- approx(zt, x_t, z)$y + rnorm(length(z), sd = 1.5)
#' res <- track_astro_auto(data.frame(z, x), sedrate = c(0.8, 5))
#' plot_track(res$track, res$scalogram, show = c(405, 41, 23.7), flags = res$flags)
#' }
#' @export
track_astro_auto <- function(data,
                             sedrate,
                             age_ma = 0,
                             k_yr = NULL,
                             astro_periods = NULL,
                             scalogram = c("msst", "superlet_sst", "cwt", "chirplet_sst"),
                             omega_nr = 6,
                             guide = NULL,
                             guide_tolerance = 0.4,
                             duration = NULL,
                             duration_tolerance = 0.5,
                             constraints = NULL,
                             jump_depths = NULL,
                             long_weight = 0,
                             flag_ratio = 0.25,
                             min_confidence = 0.5,
                             flag_long = FALSE,
                             dj = 1 / 50,
                             method = c("both", "chirp", "hmm"),
                             jump_prob = 0,
                             ...) {
  scalogram <- match.arg(scalogram)
  method <- match.arg(method)
  if (!is.null(duration) && is.null(guide)) {
    # approximate duration -> constant mean sedimentation rate as a wide guide
    zr <- range(as.data.frame(data)[, 1], na.rm = TRUE)
    guide <- data.frame(zr, diff(zr) / duration)
    guide_tolerance <- duration_tolerance
  }
  ap <- if (is.null(astro_periods)) astro_periods_tracking(age_ma, k_yr) else astro_periods
  pv <- if (is.data.frame(ap)) ap$period else ap
  dat <- as.data.frame(data)[, 1:2]
  dat <- dat[stats::complete.cases(dat), ]
  dz <- stats::median(diff(sort(dat[, 1])))
  lo <- max(3 * dz, min(pv) * sedrate[1] * 0.95)
  hi <- min(diff(range(dat[, 1])) / 2, max(pv) * sedrate[2] * 1.05)
  sc <- if (scalogram == "msst") {
    analyse_wavelet_sst(dat, lowerPeriod = lo, upperPeriod = hi, dj = dj, omega_nr = 12, n_iter = 4)
  } else if (scalogram == "cwt") {
    analyse_wavelet(dat, lowerPeriod = lo, upperPeriod = hi, dj = dj, omega_nr = omega_nr, verbose = FALSE)
  } else if (scalogram == "chirplet_sst") {
    analyse_chirplet_sst(dat, lowerPeriod = lo, upperPeriod = hi, dj = dj)
  } else {
    analyse_superlet_sst(dat, lowerPeriod = lo, upperPeriod = hi, dj = dj)
  }
  tr <- track_astro_hmm(sc, astro_periods = ap, sedrate = sedrate, guide = guide, guide_tolerance = guide_tolerance,
                        jump_depths = jump_depths, long_weight = long_weight, constraints = constraints, ...)
  # independent long-window estimate
  lw <- track_astro_hmm(sc, astro_periods = ap, sedrate = sedrate, jump_depths = jump_depths, long_weight = 1,
                        local_weight = 0, n_samples = 1)
  ratio <- tr$sedrate / lw$sedrate
  # flag where the posterior is spread (less than min_confidence of the probability within 10 percent of the track)
  # and, optionally, where the long-window estimate disagrees
  flagged <- tr$confidence10 < min_confidence
  if (flag_long) flagged <- flagged | abs(log(ratio)) > log(1 + flag_ratio)
  # intervals (merge flagged runs separated by less than one shortest cycle)
  z <- tr$depth
  gap <- min(pv) * sedrate[1]
  flags <- data.frame(from = numeric(0), to = numeric(0))
  if (any(flagged)) {
    r <- rle(flagged)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    iv <- data.frame(from = z[starts[r$values]], to = z[ends[r$values]])
    for (i in seq_len(nrow(iv))) {
      if (nrow(flags) && iv$from[i] - flags$to[nrow(flags)] < gap) flags$to[nrow(flags)] <- iv$to[i]
      else flags <- rbind(flags, iv[i, ])
    }
  }
  check <- data.frame(depth = z, sedrate = tr$sedrate, sedrate_long = lw$sedrate, ratio = ratio,
                      confidence = tr$confidence, flagged = flagged)
  out <- list(scalogram = sc, track = tr, check = check, flags = flags, astro_periods = ap)
  if (method != "hmm") {
    # chirp-matched windowed spectra + gradient-steered HMM, with the scalogram evidence added for locality
    ce <- analyse_chirp_eha(dat, window = NULL, chirp = "auto", lowerPeriod = lo, upperPeriod = hi, dj = dj,
                            criterion = "template", astro_periods = ap, sedrate = sedrate)
    tc <- track_chirp_hmm(ce, jump_prob = jump_prob, scalogram = sc, scalogram_weight = 1,
                          scalogram_args = list(astro_periods = ap, jump_depths = jump_depths, guide = guide,
                                                guide_tolerance = guide_tolerance))
    out$chirp_eha <- ce
    out$track_chirp <- tc
    # agreement of the two trackers on the common depth range
    rc <- stats::approx(tc$depth, tc$sedrate, tr$depth, rule = 2)$y
    out$check$sedrate_chirp <- rc
    out$check$disagreement <- abs(log2(tr$sedrate / rc))
    out$agreement <- mean(out$check$disagreement < log2(1.1))
    out$track_hmm <- tr
    if (is.null(attr(tc, "astro_periods"))) attr(tc, "astro_periods") <- if (!is.null(out$track_hmm)) attr(out$track_hmm, "astro_periods") else ap
    out$track <- tc                                 # the chirp-matched track is the primary one (see 'method')
  }
  class(out) <- "wr_track_auto"
  out
}

# Weighted L1-penalised smoothing of a path (majorisation-minimisation): minimise sum w (u - v)^2 + lambda sum |D u|
# with D the first (order 1, total variation) or second (order 2, L1 trend filter) difference operator
.wr_l1_smooth <- function(v, w, lambda, order = 1, n_iter = 30, eps = 1e-4) {
  n <- length(v)
  if (n < 3 || lambda <= 0) return(v)
  w <- pmax(w, 1e-3)
  D <- diff(diag(n), differences = order)
  u <- v
  for (it in seq_len(n_iter)) {
    d <- abs(as.vector(D %*% u))
    Lam <- 1 / pmax(d, eps)
    A <- diag(w) + (lambda / 2) * t(D) %*% (Lam * D)
    u_new <- as.vector(solve(A, w * v))
    if (max(abs(u_new - u)) < 1e-5) { u <- u_new; break }
    u <- u_new
  }
  u
}

#' Map of the astronomical evidence against sedimentation rate and depth
#'
#' Shows what \code{\link{track_astro_hmm}} tracks: for every sedimentation rate of a (wide) range and every position,
#' the weighted sum of the log evidence of all astronomical cycles at the depth periods they would have at that rate.
#' Where the record carries the cycles, the map shows a ridge at the true rate, and changes of the rate (a ramp, a
#' jump, an oscillation) are visible as the shape of that ridge, with the aliases (rates at which other cycles fall on
#' the ridges) as weaker parallel bands. The map is the tool to look at before choosing \code{sedrate} for the
#' tracking: on CIP Case 3 (sedimentation rate oscillating between 3 and 6 cm/kyr with a 1.2-Myr period, very red
#' noise) the posterior map of a Morlet scalogram with \code{omega_nr = 16} over a rate range of 1-20 cm/kyr shows a
#' ridge varying between 3 and 6 cm/kyr, without assuming anything about the form of the variation, whereas the
#' synchrosqueezed maps show speckle and long horizontal streaks (the slow wander of the red noise reassigned to a few
#' periods) with a bias towards the edges of the range.
#'
#' @param wavelet A scalogram (\code{\link{analyse_wavelet}} with a high \code{omega_nr} works best for noisy
#' records; \code{\link{analyse_superlet_sst}}, \code{\link{analyse_chirplet_sst}} etc. for clean ones).
#' @param sedrate Range of the sedimentation rate to show (depth units per kyr); wide by default.
#' @param track Optional result of \code{\link{track_astro_hmm}} or \code{\link{track_astro_auto}} to draw on the map.
#' @param reference Optional data frame (depth, sedimentation rate) to draw, e.g. a known or published rate.
#' @param rate_unit Factor and label of the rate axis, e.g. \code{c(100, "cm/kyr")} for a record in metres
#' \code{Default=c(1, "depth units / kyr")}.
#' @param palette Colour palette name for \code{hcl.colors} \code{Default="viridis"}.
#' @param main Title.
#' @param type \code{"evidence"} (default) or \code{"posterior"}: the posterior of the sedimentation rate of
#' \code{track} (the attribute \code{posterior} of a \code{\link{track_astro_hmm}} result; its square root is shown).
#' The posterior combines the evidence with the continuity of the rate and is the cleaner picture of a weak ridge.
#' @param ... Arguments of \code{\link{track_astro_hmm}} that shape the evidence (\code{astro_periods},
#' \code{weights}, \code{age_ma}, \code{k_yr}, \code{smooth_depth}, \code{normalise}, \code{consistent_cycles},
#' ...).
#' @return Invisibly, a list with \code{depth}, \code{sedrate} and the matrix shown (rates in rows).
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 2000, by = 0.5)
#' x_t <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19) + 0.8 * sin(2 * pi * t / 41)
#' sr <- 2 + 0.7 * sin(2 * pi * t / 600)                   # oscillating sedimentation rate
#' zt <- c(0, cumsum(sr[-1] * diff(t)))
#' z <- seq(0, max(zt), by = 1)
#' x <- approx(zt, x_t, z)$y + rnorm(length(z), sd = 1.5)
#' w <- analyse_wavelet(data.frame(z, x), lowerPeriod = 10, upperPeriod = 500, dj = 1 / 50,
#'     omega_nr = 12)
#' plot_astro_evidence(w, sedrate = c(0.5, 8), reference = data.frame(z, approx(zt, sr, z)$y))
#' }
#' @export
plot_astro_evidence <- function(wavelet = NULL, sedrate = NULL, track = NULL, reference = NULL,
                                rate_unit = c(1, "depth units / kyr"), palette = "viridis", main = NULL,
                                type = c("evidence", "posterior"), ...) {
  if (inherits(track, "wr_track_auto")) track <- track$track
  type <- match.arg(type)
  dots <- list(...)
  f <- as.numeric(rate_unit[1])
  if (type == "posterior" || is.null(wavelet)) {
    P <- attr(track, "posterior")
    if (is.null(P)) stop("type = 'posterior' needs a track from track_astro_hmm() (attribute 'posterior')")
    E <- list(zc = P$depth, ls = log2(P$sedrate), logE = if (type == "posterior") sqrt(P$post) else P$evidence)
    if (is.null(main)) main <- if (type == "posterior") "Posterior of the sedimentation rate (square root)" else "Astronomical evidence against sedimentation rate"
  } else {
  if (is.null(sedrate)) {
    ap <- dots$astro_periods
    if (is.null(ap)) ap <- if (!is.null(dots$age_ma) || !is.null(dots$k_yr)) astro_periods_tracking(if (is.null(dots$age_ma)) 0 else dots$age_ma, dots$k_yr) else c(405, 125, 95, 53.6, 41, 39.7, 23.7, 22.4, 19)
    pv <- if (is.data.frame(ap)) ap$period else ap
    sedrate <- c(min(wavelet$Period) / min(pv), max(wavelet$Period) / max(pv))
  }
  old <- getOption("WaverideR.track_evidence")
  options(WaverideR.track_evidence = TRUE)
  on.exit(options(WaverideR.track_evidence = old), add = TRUE)
  E <- do.call(track_astro_hmm, c(list(wavelet, sedrate = sedrate, n_samples = 0), dots))
  if (is.null(main)) main <- "Astronomical evidence against sedimentation rate"
  }
  graphics::image(E$zc, E$ls, t(E$logE), col = grDevices::hcl.colors(64, palette), xlab = "Depth", yaxt = "n",
                  ylab = paste0("Sedimentation rate (", rate_unit[2], ")"), main = main, xaxs = "i")
  rng <- 2^range(E$ls) * f
  ticks <- c(1, 2, 3, 5) * rep(10^((floor(log10(rng[1])) - 1):(ceiling(log10(rng[2])) + 1)), each = 4)
  ticks <- ticks[ticks >= rng[1] & ticks <= rng[2]]
  graphics::axis(2, at = log2(ticks / f), labels = ticks)
  if (!is.null(reference)) graphics::lines(reference[, 1], log2(reference[, 2]), col = "red", lwd = 1.3, lty = 2)
  if (!is.null(track)) graphics::lines(track$depth, log2(track$sedrate), col = "white", lwd = 1.5)
  invisible(list(depth = E$zc, sedrate = 2^E$ls, map = E$logE, cycles = E$astro_used, weights = E$weights_used))
}
