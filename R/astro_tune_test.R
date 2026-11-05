#' @title Test a minimal tuning: multi-band rate tracking with surrogates run through the whole chain
#'
#' @description Estimates the sedimentation-rate history of one or more stratigraphic records
#' from all astronomical cycles at once, turns it into a (minimally) tuned age model with
#' uncertainty, and tests the tuned records for astronomical amplitude and frequency
#' modulation against surrogates that keep the power spectrum of the data.
#' The rate tracker here is a fast version of \code{\link{track_astro_hmm}} (Morlet power,
#' no synchrosqueezing), so that every surrogate can be tracked, tuned and tested in exactly the
#' same way as the data; this is what calibrates the gain that any tuning gives to any series.
#' For the best age model of a record use \code{\link{track_astro_auto}},
#' \code{\link{track_astro_hmm}} or \code{\link{track_astro_iterative}}; this function tests
#' whether such a (minimal) tuning is supported by astronomical structure beyond the power
#' spectrum, and which rate range (interpretation) is supported. Tests of an already tuned
#' record are in \code{\link{test_tuned}}.
#'
#' @details
#' \strong{1. Rate likelihood from all cycles.} The Morlet wavelet power of each record is
#' computed (in depth) at the frequencies that every trial sedimentation rate predicts for the
#' 405-kyr eccentricity, short eccentricity, obliquity and precession terms. Power is divided by
#' a robust background (median power over depth, smoothed with a running median over
#' frequency), and the weighted sum of the log power ratios over all terms gives the
#' log-likelihood of every rate at every depth. Combining all terms (rather than following one
#' cycle) is the multi-order approach of tacholess order tracking in machine diagnostics
#' (Leclere et al., 2016) and of multiband wavelet age modelling (Duesing et al., 2021).
#'
#' \strong{2. Smooth rate path (hidden Markov model).} The log rate follows a random walk with
#' standard deviation \code{sigma10} per 10 m (continuity prior, as the sedimentation-rate
#' transition model of HMM-Match; Lin et al., 2014). The Viterbi algorithm gives the most
#' probable rate path, the forward-backward algorithm the rate posterior at every depth, and
#' forward-filtering backward-sampling gives \code{n_paths} sampled rate paths, which give the
#' age-model uncertainty. \code{lambda} tempers the likelihood (neighbouring depths are not
#' independent). The marginal log-likelihood of the model is reported for comparing rate
#' ranges (e.g. two competing interpretations) through \code{sedmin} and \code{sedmax}.
#'
#' \strong{3. Tuning.} Time is the integral of 1 / rate along the most probable path. Because
#' the path is smooth (a few degrees of freedom per 405-kyr cycle), this is a minimal tuning:
#' it does not impose precession-scale frequency modulation (Zeeden et al., 2015).
#'
#' \strong{4. Modulation tests on the tuned records.} (a) precession amplitude versus the
#' record's own eccentricity band: wide precession filter (65\% bandwidth), Hilbert amplitude,
#' low-pass (< 1/75 kyr), Spearman correlation (Zeeden et al., 2015); (b) short-eccentricity
#' amplitude versus the record's 405-kyr band (wide filter, low-pass < 1/250 kyr); (c)
#' short-eccentricity frequency modulation (instantaneous frequency) versus the 405-kyr band
#' (Laurin et al., 2016). For an eccentricity-forced record (a) and (b) have the sign of the
#' proxy response and (c) the opposite sign; \code{S_am = |a + b| / 2} and
#' \code{S_all = |a + b - c| / 3} use this sign consistency. Wide filters and low-pass
#' envelopes avoid the amplitude modulation that narrow filters create from frequency
#' modulation (Huybers and Aharonson, 2010).
#'
#' \strong{5. Significance.} Phase-randomised surrogates of the depth records (Ebisuzaki,
#' 1997; with the same random phases for all records, so the correlation between records is
#' kept; Prichard and Theiler, 1994) are run through the whole chain (rate likelihood, path,
#' tuning and tests). The surrogates have exactly the power spectrum of the data, including the
#' power of any cycles, so the tests ask whether the data hold astronomical structure beyond
#' their power spectrum: a smooth rate path that aligns power at all astronomical frequency
#' ratios at the same depths (\code{p_loglik}), and the modulation patterns above.
#'
#' @param data Evenly spaced records (matrix or data frame): depth in metres in column 1, and one or more proxies
#' measured on the same depths in the further columns (no missing values). Every proxy is detrended (see
#' \code{detrend}) and standardised; the column names are used as record names. Required, no default
#' (\code{Default=NULL} gives an error).
#' @param sedmin,sedmax Lowest and highest sedimentation rate of the rate states of the hidden Markov model, in
#' cm/kyr (numeric, > 0, \code{sedmin} < \code{sedmax}). Comparing the marginal log-likelihood of different
#' ranges compares rate interpretations. \code{Default=2} and \code{Default=25}.
#' @param numsed Number of log-spaced rate states between \code{sedmin} and \code{sedmax} (integer, >= 2).
#' \code{Default=50}.
#' @param age_ma Age of the record in Ma (numeric) used to compute the target periods from AstroGeo22 (see
#' \code{\link{astro_sedrate_test}}); \code{NULL} uses present-day periods. \code{Default=NULL}.
#' @param targetE,targetP,targetO Your own eccentricity, precession and obliquity target periods in kyr (numeric
#' vectors); the longest value of \code{targetE} is taken as the 405-kyr term. They override the age-dependent
#' values; \code{NULL} keeps those. \code{Default=NULL}.
#' @param weights Weights of the four term groups in the rate likelihood (numeric vector of length 4, >= 0, in
#' the order 405-kyr eccentricity, short eccentricity, obliquity, precession); each group weight is shared equally
#' among its terms. \code{Default=c(0.5, 1.5, 1, 2)}.
#' @param sigma10 Standard deviation of the change in log sedimentation rate per 10 m of depth (numeric, > 0):
#' the continuity prior of the rate path; smaller values give smoother rate paths. \code{Default=0.05}.
#' @param lambda Tempering factor applied to the rate log-likelihood, to account for the dependence of
#' neighbouring depths (numeric, between 0 and 1; 1 = no tempering). \code{Default=0.3}.
#' @param step Depth step in metres between the states of the hidden Markov model (numeric, > 0; at least one
#' sample interval). \code{Default=0.5}.
#' @param omega0 Central angular frequency of the Morlet wavelet used for the wavelet power (numeric, > 0);
#' larger values give a better frequency and a poorer depth resolution. \code{Default=6}.
#' @param n_sur Number of phase-randomised surrogates run through the whole chain to compute the p-values
#' (integer, > 0). \code{Default=200}.
#' @param n_paths Number of rate paths sampled from the posterior (forward-filtering backward-sampling) for the
#' uncertainty of the age model and the rate (integer, > 0). \code{Default=200}.
#' @param detrend Remove a linear trend from every record before the analysis (logical, TRUE or FALSE).
#' \code{Default=TRUE}.
#' @param run_multicore Process the surrogates in parallel on several cores (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param n_cores Number of cores used for the surrogates (integer); \code{NULL} uses the option
#' \code{waverider.n_cores} or, with \code{run_multicore = TRUE}, all but one detected core (at most
#' \code{getOption("waverider.max_cores", 4)}). \code{Default=NULL}.
#' @param verbose Print progress messages (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @return A list of class \code{astro_tune_test} with
#' \describe{
#'   \item{age_model}{depth, time (kyr, median of the sampled paths), its 2.5 and 97.5\% quantiles,
#'   and the time of the most probable (Viterbi) path, which is used for the tests}
#'   \item{rate}{depth, most probable rate, posterior mean rate and its 2.5 and 97.5\% quantiles}
#'   \item{tests}{per record: r_prec, r_eam, r_efm, S_am, S_all and their probabilities}
#'   \item{summary}{mean rate and duration (median of the sampled paths) with its 95\% interval, the
#'   duration of the Viterbi path, marginal log-likelihood and its
#'   probability, mean S_am and S_all and their probabilities}
#'   \item{surrogates}{the surrogate statistics}
#'   \item{settings}{the settings and target periods}
#' }
#'
#' @references
#' Duesing, W. et al. (2021). Multiband wavelet age modeling for a ~293 m (~600 kyr) sediment
#' core from Chew Bahir Basin, southern Ethiopian Rift. Frontiers in Earth Science 9, 594047
#' <doi:10.3389/feart.2021.594047>
#'
#' Ebisuzaki, W. (1997). A method to estimate the statistical significance of a correlation
#' when the data are serially correlated. Journal of Climate 10, 2147-2153
#' <doi:10.1175/1520-0442(1997)010<2147:AMTETS>2.0.CO;2>
#'
#' Huybers, P. and Aharonson, O. (2010). Orbital tuning, eccentricity, and the frequency
#' modulation of climatic precession. Paleoceanography 25, PA4228 <doi:10.1029/2010PA001952>
#'
#' Laurin, J., Meyers, S.R., Galeotti, S. and Lanci, L. (2016). Frequency modulation reveals
#' the phasing of orbital eccentricity during Cretaceous Oceanic Anoxic Event II and the Eocene
#' hyperthermals. Earth and Planetary Science Letters 442, 143-156
#' <doi:10.1016/j.epsl.2016.02.047>
#'
#' Leclere, Q., Andre, H. and Antoni, J. (2016). A multi-order probabilistic approach for
#' instantaneous angular speed tracking: debriefing of the CMMNO'14 diagnosis contest.
#' Mechanical Systems and Signal Processing 81, 375-386 <doi:10.1016/j.ymssp.2016.02.053>
#'
#' Lin, L., Khider, D., Lisiecki, L.E. and Lawrence, C.E. (2014). Probabilistic sequence
#' alignment of stratigraphic records. Paleoceanography 29, 976-989
#' <doi:10.1002/2014PA002713>
#'
#' Prichard, D. and Theiler, J. (1994). Generating surrogate data for time series with several
#' simultaneously measured variables. Physical Review Letters 73, 951-954
#' <doi:10.1103/PhysRevLett.73.951>
#'
#' Zeeden, C., Meyers, S.R., Lourens, L.J. and Hilgen, F.J. (2015). Testing astronomically tuned
#' age models. Paleoceanography 30, 369-383 <doi:10.1002/2014PA002762>
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' res <- astro_tune_test(data = mag[1:600, ], sedmin = 0.3, sedmax = 2, numsed = 30,
#'     n_sur = 10, n_paths = 20)
#' res$summary
#' }
#' @export
#' @importFrom stats median approx fft runif cor quantile lm residuals sd
astro_tune_test <- function(data = NULL,
    sedmin = 2,
    sedmax = 25,
    numsed = 50,
    age_ma = NULL,
    targetE = NULL,
    targetP = NULL,
    targetO = NULL,
    weights = c(0.5, 1.5, 1, 2),
    sigma10 = 0.05,
    lambda = 0.3,
    step = 0.5,
    omega0 = 6,
    n_sur = 200,
    n_paths = 200,
    detrend = TRUE,
    run_multicore = FALSE,
    n_cores = NULL,
    verbose = FALSE) {
  dat <- unname(as.matrix(data))
  if (ncol(dat) < 2) {
    stop("astro_tune_test: data needs a depth column and at least one proxy", call. = FALSE)
  }
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  if (any(!is.finite(dat))) {
    stop("astro_tune_test: data must not contain missing values", call. = FALSE)
  }
  dd <- diff(dat[, 1])
  dx <- stats::median(dd)
  if (!all(abs(dd - dx) <= 1e-3 * abs(dx))) {
    stop("astro_tune_test: the depth values must be evenly spaced", call. = FALSE)
  }
  n <- nrow(dat)
  Y <- dat[, -1, drop = FALSE]
  rec_names <- colnames(data)[-1]
  if (is.null(rec_names) || length(rec_names) != ncol(Y)) {
    rec_names <- paste0("record", seq_len(ncol(Y)))
  }
  colnames(Y) <- rec_names
  for (j in seq_len(ncol(Y))) {
    y <- Y[, j]
    if (detrend) {
      y <- stats::residuals(stats::lm(y ~ dat[, 1]))
    }
    Y[, j] <- (y - mean(y)) / stats::sd(y)
  }
  tg <- wr_ast_targets(age_ma, targetE, targetP, targetO)
  terms <- wr_att_terms(tg, weights)
  rates <- exp(seq(log(sedmin), log(sedmax), length.out = numsed))
  set <- list(dx = dx, rates = rates, terms = terms, step = step, sigma10 = sigma10,
      lambda = lambda, omega0 = omega0, tg = tg)

  obs <- wr_att_chain(Y, set, keep = TRUE)
  # sampled rate paths -> age-model and rate uncertainty
  z <- (0:(n - 1)) * dx
  hm <- obs$hmm
  P <- wr_att_ffbs(hm, n_paths)
  TT <- matrix(NA_real_, n, n_paths)
  for (k in seq_len(n_paths)) {
    TT[, k] <- wr_att_time(z, obs$zs, rates[P[, k]])
  }
  tq <- matrix(NA_real_, n, 3)
  for (i in seq_len(n)) {
    tq[i, ] <- stats::quantile(TT[i, ], c(0.025, 0.5, 0.975))
  }
  rq <- matrix(NA_real_, nrow(P), 2)
  for (i in seq_len(nrow(P))) {
    rq[i, ] <- stats::quantile(rates[P[i, ]], c(0.025, 0.975))
  }
  durs <- TT[n, ]

  if (verbose) {
    message("astro_tune_test: ", n_sur, " surrogates")
  }
  sims <- wr_lapply(seq_len(n_sur), function(k) {
    wr_att_chain(wr_att_surrogate(Y), set, keep = FALSE)
  }, n_cores = n_cores, run_multicore = run_multicore, combine = "list", progress = verbose)
  L <- numeric(n_sur)
  Sa <- numeric(n_sur)
  Sl <- numeric(n_sur)
  R <- array(NA_real_, c(n_sur, ncol(Y), 3))
  for (k in seq_len(n_sur)) {
    L[k] <- sims[[k]]$loglik
    Sa[k] <- mean(sims[[k]]$mods$S_am)
    Sl[k] <- mean(sims[[k]]$mods$S_all)
    R[k, , 1] <- sims[[k]]$mods$r_prec
    R[k, , 2] <- sims[[k]]$mods$r_eam
    R[k, , 3] <- sims[[k]]$mods$r_efm
  }
  pv <- function(v0, v) {
    (1 + sum(v >= v0, na.rm = TRUE)) / (sum(is.finite(v)) + 1)
  }
  tests <- obs$mods
  tests$p_prec <- NA_real_
  tests$p_eam <- NA_real_
  tests$p_efm <- NA_real_
  tests$p_S_am <- NA_real_
  tests$p_S_all <- NA_real_
  for (j in seq_len(ncol(Y))) {
    tests$p_prec[j] <- pv(abs(tests$r_prec[j]), abs(R[, j, 1]))
    tests$p_eam[j] <- pv(abs(tests$r_eam[j]), abs(R[, j, 2]))
    tests$p_efm[j] <- pv(abs(tests$r_efm[j]), abs(R[, j, 3]))
    sa_j <- numeric(n_sur)
    sl_j <- numeric(n_sur)
    for (k in seq_len(n_sur)) {
      sa_j[k] <- sims[[k]]$mods$S_am[j]
      sl_j[k] <- sims[[k]]$mods$S_all[j]
    }
    tests$p_S_am[j] <- pv(tests$S_am[j], sa_j)
    tests$p_S_all[j] <- pv(tests$S_all[j], sl_j)
  }
  dmed <- stats::median(durs)
  summ <- data.frame(mean_rate = (z[n] - z[1]) * 100 / dmed, duration_kyr = dmed,
      duration_viterbi = obs$t[n], duration_lo = stats::quantile(durs, 0.025),
      duration_hi = stats::quantile(durs, 0.975),
      loglik = obs$loglik, p_loglik = pv(obs$loglik, L), S_am = mean(tests$S_am),
      p_S_am = pv(mean(tests$S_am), Sa), S_all = mean(tests$S_all), p_S_all = pv(mean(tests$S_all), Sl),
      row.names = NULL)
  age_model <- data.frame(depth = dat[, 1], time = tq[, 2], time_lo = tq[, 1], time_hi = tq[, 3],
      time_viterbi = obs$t)
  rate <- data.frame(depth = dat[1, 1] + obs$zs, rate = obs$path, rate_mean = hm$post_mean,
      rate_lo = rq[, 1], rate_hi = rq[, 2])
  out <- list(age_model = age_model, rate = rate, tests = tests, summary = summ,
      surrogates = list(loglik = L, S_am = Sa, S_all = Sl),
      settings = list(sedmin = sedmin, sedmax = sedmax, numsed = numsed, weights = weights,
      sigma10 = sigma10, lambda = lambda, step = step, omega0 = omega0, n_sur = n_sur,
      n_paths = n_paths, dx = dx, targets = tg, terms = terms))
  class(out) <- "astro_tune_test"
  out
}

# Internal helpers

wr_att_terms <- function(tg, weights) {
  ne <- length(tg$E) - 1
  per <- c(tg$E[1], tg$E[-1], tg$O, tg$P)
  wt <- c(weights[1], rep(weights[2] / ne, ne), rep(weights[3] / max(1, length(tg$O)),
      length(tg$O)), rep(weights[4] / length(tg$P), length(tg$P)))
  data.frame(period = per, weight = wt)
}

wr_att_morlet <- function(y, dx, f, omega0) {
  n <- length(y)
  yp <- c(rev(y), y, rev(y))
  N <- 2^ceiling(log2(length(yp)))
  yp <- c(yp, rep(0, N - length(yp)))
  Yf <- stats::fft(yp)
  k <- 0:(N - 1)
  w <- 2 * pi * k / (N * dx)
  ff <- 4 * pi / (omega0 + sqrt(2 + omega0^2))
  P <- matrix(NA_real_, n, length(f))
  for (j in seq_along(f)) {
    s <- 1 / (f[j] * ff)
    psi <- exp(-(s * w - omega0)^2 / 2)
    psi[k > N / 2 | k == 0] <- 0
    W <- stats::fft(Yf * psi, inverse = TRUE) / N
    P[, j] <- Mod(W[(n + 1):(2 * n)])^2
  }
  P
}

wr_att_background <- function(P, f, half_octave = 0.6) {
  m <- numeric(ncol(P))
  for (j in seq_len(ncol(P))) {
    m[j] <- stats::median(P[, j])
  }
  lf <- log2(f)
  b <- numeric(length(f))
  for (j in seq_along(f)) {
    b[j] <- stats::median(m[abs(lf - lf[j]) <= half_octave])
  }
  b
}

# emission log-likelihood (rows = HMM depth steps, columns = rate states)
wr_att_emission <- function(y, set, cap = 3) {
  n <- length(y)
  dx <- set$dx
  rates <- set$rates
  terms <- set$terms
  ff <- numeric(length(rates) * nrow(terms))
  idx <- matrix(NA_integer_, length(rates), nrow(terms))
  q <- 0
  for (k in seq_along(rates)) {
    for (i in seq_len(nrow(terms))) {
      q <- q + 1
      ff[q] <- 100 / (rates[k] * terms$period[i])
      idx[k, i] <- q
    }
  }
  nyq <- 0.5 / dx
  fmin <- 2 / (n * dx)
  ok_f <- ff < 0.8 * nyq & ff > fmin
  ffc <- pmin(pmax(ff, fmin), 0.8 * nyq)
  P <- wr_att_morlet(y, dx, ffc, set$omega0)
  B <- wr_att_background(P, ffc)
  rows <- unique(round(seq(1, n, by = max(1, round(set$step / dx)))))
  E <- matrix(0, length(rows), length(rates))
  for (k in seq_along(rates)) {
    for (i in seq_len(nrow(terms))) {
      j <- idx[k, i]
      if (ok_f[j]) {
        E[, k] <- E[, k] + terms$weight[i] * pmin(log(P[rows, j] / B[j]), cap)
      }
    }
  }
  list(E = E, zs = (rows - 1) * dx)
}

wr_att_lse <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

wr_att_hmm <- function(E, rates, step, sigma10, lambda) {
  ns <- nrow(E)
  K <- length(rates)
  lr <- log(rates)
  sig <- sigma10 * sqrt(step / 10)
  logT <- -outer(lr, lr, "-")^2 / (2 * sig^2)
  for (k in 1:K) {
    logT[k, ] <- logT[k, ] - wr_att_lse(logT[k, ])
  }
  lp0 <- rep(-log(K), K)
  em <- lambda * E
  V <- matrix(-Inf, ns, K)
  bp <- matrix(0L, ns, K)
  A <- matrix(-Inf, ns, K)
  V[1, ] <- lp0 + em[1, ]
  A[1, ] <- V[1, ]
  for (t in 2:ns) {
    for (k in 1:K) {
      cand <- V[t - 1, ] + logT[, k]
      b <- which.max(cand)
      bp[t, k] <- b
      V[t, k] <- cand[b] + em[t, k]
      A[t, k] <- wr_att_lse(A[t - 1, ] + logT[, k]) + em[t, k]
    }
  }
  path <- integer(ns)
  path[ns] <- which.max(V[ns, ])
  for (t in (ns - 1):1) {
    path[t] <- bp[t + 1, path[t + 1]]
  }
  Bm <- matrix(0, ns, K)
  for (t in (ns - 1):1) {
    for (k in 1:K) {
      Bm[t, k] <- wr_att_lse(logT[k, ] + em[t + 1, ] + Bm[t + 1, ])
    }
  }
  pm <- numeric(ns)
  for (t in 1:ns) {
    po <- A[t, ] + Bm[t, ]
    po <- exp(po - wr_att_lse(po))
    pm[t] <- exp(sum(po * lr))
  }
  list(path_idx = path, path = rates[path], post_mean = pm, loglik = wr_att_lse(A[ns, ]),
      A = A, logT = logT)
}

# forward-filtering backward-sampling of rate paths (state indices, rows = steps)
wr_att_ffbs <- function(hm, n_paths) {
  A <- hm$A
  ns <- nrow(A)
  K <- ncol(A)
  P <- matrix(0L, ns, n_paths)
  for (k in seq_len(n_paths)) {
    w <- exp(A[ns, ] - wr_att_lse(A[ns, ]))
    P[ns, k] <- sample.int(K, 1, prob = w)
    for (t in (ns - 1):1) {
      lw <- A[t, ] + hm$logT[, P[t + 1, k]]
      w <- exp(lw - wr_att_lse(lw))
      P[t, k] <- sample.int(K, 1, prob = w)
    }
  }
  P
}

wr_att_time <- function(z, zs, r) {
  rr <- stats::approx(zs, r, xout = z, rule = 2)$y
  c(0, cumsum(diff(z) * 100 / ((rr[-1] + rr[-length(rr)]) / 2)))
}

# Gaussian-tapered FFT bandpass (flow = 0: low-pass)
wr_att_bandpass <- function(x, dx, flow, fhigh, taper = 0.15) {
  n <- length(x)
  xm <- x - mean(x)
  xp <- c(rev(xm), xm, rev(xm))
  N <- length(xp)
  X <- stats::fft(xp)
  f <- (0:(N - 1)) / (N * dx)
  f[f > 0.5 / dx] <- f[f > 0.5 / dx] - 1 / dx
  af <- abs(f)
  g <- rep(1, N)
  lo <- af < flow
  if (flow > 0) {
    g[lo] <- exp(-((af[lo] - flow) / (taper * flow))^2 / 2)
  }
  hi <- af > fhigh
  g[hi] <- exp(-((af[hi] - fhigh) / (taper * fhigh))^2 / 2)
  out <- Re(stats::fft(X * g, inverse = TRUE)) / N
  out[(n + 1):(2 * n)]
}

wr_att_analytic <- function(x) {
  n <- length(x)
  h <- numeric(n)
  h[1] <- 1
  if (n %% 2 == 0) {
    h[n / 2 + 1] <- 1
    h[2:(n / 2)] <- 2
  } else {
    h[2:((n + 1) / 2)] <- 2
  }
  stats::fft(stats::fft(x) * h, inverse = TRUE) / n
}

wr_att_unwrap <- function(ph) {
  d <- diff(ph)
  d <- d - 2 * pi * round(d / (2 * pi))
  c(ph[1], ph[1] + cumsum(d))
}

wr_att_modulation <- function(yt, dt, tg) {
  E <- tg$E
  P <- tg$P
  sp <- function(a, b) {
    ok <- is.finite(a) & is.finite(b)
    suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))
  }
  n <- length(yt)
  edge <- seq_len(n) > round(0.05 * n) & seq_len(n) <= round(0.95 * n)
  e_band <- wr_att_bandpass(yt, dt, 1 / (E[1] * 1.6), 1 / (min(E[-1]) * 0.8))
  e405 <- wr_att_bandpass(yt, dt, 1 / (E[1] * 1.3), 1 / (E[1] * 0.75))
  fc <- mean(1 / P)
  pa <- Mod(wr_att_analytic(wr_att_bandpass(yt, dt, fc * 0.675, fc * 1.325)))
  pa <- wr_att_bandpass(pa, dt, 0, 1 / 75)
  r_prec <- sp(pa[edge], e_band[edge])
  fe <- mean(1 / E[-1])
  an <- wr_att_analytic(wr_att_bandpass(yt, dt, fe * 0.675, fe * 1.325))
  ea <- wr_att_bandpass(Mod(an), dt, 0, 1 / 250)
  r_eam <- sp(ea[edge], e405[edge])
  ifq <- c(diff(wr_att_unwrap(Arg(an))), NA) / (2 * pi * dt)
  ifq[!is.finite(ifq)] <- fe
  ifq <- wr_att_bandpass(ifq, dt, 0, 1 / 250)
  r_efm <- sp(ifq[edge], e405[edge])
  data.frame(r_prec = r_prec, r_eam = r_eam, r_efm = r_efm, S_am = abs(r_prec + r_eam) / 2,
      S_all = abs(r_prec + r_eam - r_efm) / 3, polarity = sign(r_prec + r_eam - r_efm))
}

# multivariate phase-randomised surrogate (same phases for all records)
wr_att_surrogate <- function(Y) {
  n <- nrow(Y)
  m <- floor((n - 1) / 2)
  ph <- stats::runif(m, 0, 2 * pi)
  sgn <- sample(c(-1, 1), 1)
  out <- Y
  for (j in seq_len(ncol(Y))) {
    Fy <- stats::fft(Y[, j] - mean(Y[, j]))
    Z <- Fy
    Z[2:(m + 1)] <- Fy[2:(m + 1)] * exp(1i * ph)
    Z[n - (1:m) + 1] <- Conj(Z[2:(m + 1)])
    if (n %% 2 == 0) {
      Z[n / 2 + 1] <- Fy[n / 2 + 1] * sgn
    }
    out[, j] <- Re(stats::fft(Z, inverse = TRUE)) / n
  }
  out
}

wr_att_chain <- function(Y, set, keep = FALSE) {
  n <- nrow(Y)
  z <- (0:(n - 1)) * set$dx
  E <- NULL
  for (j in seq_len(ncol(Y))) {
    em <- wr_att_emission(Y[, j], set)
    if (is.null(E)) {
      E <- em$E
    } else {
      E <- E + em$E
    }
  }
  E <- E / ncol(Y)
  hm <- wr_att_hmm(E, set$rates, set$step, set$sigma10, set$lambda)
  t <- wr_att_time(z, em$zs, hm$path)
  dt <- stats::median(diff(t))
  tg_ <- seq(0, t[n], by = dt)
  mods <- NULL
  for (j in seq_len(ncol(Y))) {
    yt <- stats::approx(t, Y[, j], xout = tg_)$y
    mods <- rbind(mods, data.frame(record = colnames(Y)[j], wr_att_modulation(yt, dt, set$tg)))
  }
  out <- list(loglik = hm$loglik, mods = mods)
  if (keep) {
    out$t <- t
    out$zs <- em$zs
    out$path <- hm$path
    out$hmm <- hm
  }
  out
}

#' @export
print.astro_tune_test <- function(x, ...) {
  s <- x$summary
  cat("astro_tune_test: rates", x$settings$sedmin, "-", x$settings$sedmax, "cm/kyr,",
      nrow(x$tests), "record(s),", x$settings$n_sur, "surrogates\n")
  cat(sprintf("mean rate %.2f cm/kyr; duration %.0f kyr (95%%: %.0f-%.0f; Viterbi path %.0f)\n",
      s$mean_rate,
      s$duration_kyr, s$duration_lo, s$duration_hi, s$duration_viterbi))
  cat(sprintf("rate path log-likelihood %.1f, p = %.3f (surrogates with the same spectrum)\n",
      s$loglik, s$p_loglik))
  cat(sprintf("amplitude modulation S_am %.3f, p = %.3f; with FM S_all %.3f, p = %.3f\n", s$S_am,
      s$p_S_am, s$S_all, s$p_S_all))
  print(x$tests, digits = 3, row.names = FALSE)
  invisible(x)
}

#' @title Plot an astro_tune_test result
#' @description Rate path with its 95\% band, the age model with its 95\% band, and the observed
#' statistics against the surrogates.
#' @param x Result to plot: an object of class \code{astro_tune_test}, as returned by
#' \code{\link{astro_tune_test}}. Required, no default.
#' @param ... Not used; present for compatibility with the generic \code{plot()}.
#' @return No return value, called for the plot.
#' @export
plot.astro_tune_test <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::layout(matrix(1:4, 2, 2, byrow = TRUE))
  graphics::par(mar = c(4, 4.2, 2.4, 1))
  r <- x$rate
  graphics::plot(r$depth, r$rate, type = "n", ylim = range(c(r$rate_lo, r$rate_hi)),
      xlab = "depth (m)", ylab = "sedimentation rate (cm/kyr)", main = "rate path", las = 1)
  graphics::polygon(c(r$depth, rev(r$depth)), c(r$rate_lo, rev(r$rate_hi)), col = "grey85",
      border = NA)
  graphics::lines(r$depth, r$rate, lwd = 2)
  a <- x$age_model
  graphics::plot(a$depth, a$time, type = "n", ylim = range(c(a$time_lo, a$time_hi)),
      xlab = "depth (m)", ylab = "time (kyr)", main = "age model", las = 1)
  graphics::polygon(c(a$depth, rev(a$depth)), c(a$time_lo, rev(a$time_hi)), col = "grey85",
      border = NA)
  graphics::lines(a$depth, a$time, lwd = 2)
  s <- x$summary
  graphics::hist(x$surrogates$loglik, breaks = 20, col = "grey80", border = "white",
      xlim = range(c(x$surrogates$loglik, s$loglik)), main = sprintf("rate path, p = %.3f",
      s$p_loglik), xlab = "log-likelihood (surrogates)")
  graphics::abline(v = s$loglik, lwd = 2, col = "#eb6834")
  graphics::hist(x$surrogates$S_am, breaks = 20, col = "grey80", border = "white",
      xlim = range(c(x$surrogates$S_am, s$S_am)), main = sprintf("amplitude modulation, p = %.3f",
      s$p_S_am), xlab = "S_am (surrogates)")
  graphics::abline(v = s$S_am, lwd = 2, col = "#eb6834")
  invisible(NULL)
}
