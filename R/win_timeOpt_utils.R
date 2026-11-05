# Internal helpers of win_timeOpt(): timeOpt settings per fit, timeOpt evaluation, noise-model fit,
# noise surrogates, the testable mask and the noise-level (null) tables. All functions are internal.

# Settings

# timeOpt settings of one interpretation: fit 1 (precession amplitude modulation), fit 2 (short
# eccentricity amplitude modulation) or fit 3 (obliquity amplitude modulation, run as timeOpt fit 1
# with obliquity carrier periods and the obliquity modulation periods).
wr_to_cfg <- function(fit, fitModPwr, flow, fhigh, roll, targetE, targetP, targetO, targetOmod,
    detrend, limit) {
  if (!(fit %in% c(1, 2, 3))) {
    stop("fit must be 1 (precession), 2 (short eccentricity) or 3 (obliquity)", call. = FALSE)
  }
  if (fit == 3) {
    cfg <- list(to_fit = 1, fitModPwr = FALSE,
        flow = if (is.null(flow)) { 1 / 46 } else { flow },
        fhigh = if (is.null(fhigh)) { 1 / 36 } else { fhigh },
        roll = roll, targetE = targetOmod, targetP = targetO,
        carrier_min = min(targetO), longest = max(targetOmod))
  } else if (fit == 1) {
    cfg <- list(to_fit = 1, fitModPwr = fitModPwr, flow = flow, fhigh = fhigh, roll = roll,
        targetE = targetE, targetP = targetP, carrier_min = min(targetP), longest = max(targetE))
  } else {
    cfg <- list(to_fit = 2, fitModPwr = fitModPwr, flow = flow, fhigh = fhigh, roll = roll,
        targetE = targetE, targetP = targetP, carrier_min = min(targetE[-1]), longest = targetE[1])
  }
  cfg$fit <- fit
  cfg$detrend <- detrend
  cfg$limit <- limit
  cfg
}

# text key of everything the noise tables depend on (used to check a supplied null_cache)
wr_to_cfg_key <- function(cfg, sedrate, dx, window_mode, window_duration) {
  paste(c(cfg$fit, cfg$to_fit, cfg$fitModPwr, format(cfg$flow), format(cfg$fhigh), format(cfg$roll),
      format(cfg$targetE), format(cfg$targetP), cfg$detrend, signif(sedrate, 8), signif(dx, 8),
      window_mode, if (window_mode == "time") { signif(window_duration, 8) } else { "" }),
      collapse = "|")
}

wr_sedrate_grid <- function(sedmin, sedmax, numsed, linLog) {
  if (numsed == 1) {
    return(sedmin)
  }
  if (linLog == 0) {
    sedinc <- (sedmax - sedmin) / (numsed - 1)
    return(sedmin + 0:(numsed - 1) * sedinc)
  }
  sedinc <- (log10(sedmax) - log10(sedmin)) / (numsed - 1)
  10^(log10(sedmin) + 0:(numsed - 1) * sedinc)
}

# TimeOpt evaluation

# envelope and spectral power r2 of one data window (matrix depth, value).
# sed: all rates of the grid (one timeOpt call, depth windows) or a single rate.
wr_to_run <- function(dat, sed, cfg, linLog = 1) {
  ns <- length(sed)
  out <- matrix(NA_real_, ns, 2)
  res <- tryCatch(astrochron::timeOpt(
      dat = dat,
      sedmin = min(sed),
      sedmax = max(sed),
      numsed = ns,
      linLog = linLog,
      limit = if (ns == 1) { FALSE } else { cfg$limit },
      fit = cfg$to_fit,
      fitModPwr = cfg$fitModPwr,
      flow = cfg$flow,
      fhigh = cfg$fhigh,
      roll = cfg$roll,
      targetE = cfg$targetE,
      targetP = cfg$targetP,
      detrend = cfg$detrend,
      output = 1,
      title = NULL,
      genplot = FALSE,
      check = FALSE,
      verbose = FALSE), error = function(e) { NULL })
  if (!is.null(res) && nrow(res) == ns) {
    out[, 1] <- res[, 2]
    out[, 2] <- res[, 3]
  }
  out
}

# Noise model

# Whittle fit of an AR(1) spectrum (flat at low frequency, bending to f^-2) and of a power-law
# spectrum (f^-beta) to the periodogram of a (detrended, standardised) series.
wr_noise_fit <- function(y, x = NULL, detrend = TRUE) {
  n <- length(y)
  if (detrend && !is.null(x)) {
    y <- stats::residuals(stats::lm(y ~ x))
  } else {
    y <- y - mean(y)
  }
  s <- stats::sd(y)
  if (!is.finite(s) || s <= 0) {
    return(list(rho = 0, beta = 0, ll_ar1 = 0, ll_pl = -Inf))
  }
  y <- y / s
  m <- floor((n - 1) / 2)
  if (m < 4) {
    return(list(rho = 0, beta = 0, ll_ar1 = 0, ll_pl = -Inf))
  }
  pg <- Mod(stats::fft(y))[2:(m + 1)]^2 / n
  fr <- (1:m) / n
  ll_ar1 <- function(r) {
    g <- 1 / (1 + r^2 - 2 * r * cos(2 * pi * fr))
    cc <- mean(pg / g)
    -sum(log(cc * g)) - m
  }
  ll_pl <- function(b) {
    g <- fr^(-b)
    cc <- mean(pg / g)
    -sum(log(cc * g)) - m
  }
  oa <- stats::optimize(ll_ar1, c(0, 0.99), maximum = TRUE)
  op <- stats::optimize(ll_pl, c(0, 3), maximum = TRUE)
  list(rho = oa$maximum, beta = op$maximum, ll_ar1 = oa$objective, ll_pl = op$objective)
}

# family (1 = AR(1), 2 = power law) and parameter of the chosen noise model
wr_noise_choose <- function(nf, noise_model) {
  if (noise_model == "ar1") {
    return(c(1, nf$rho))
  }
  if (noise_model == "powerlaw") {
    return(c(2, nf$beta))
  }
  if (nf$ll_pl > nf$ll_ar1) {
    c(2, nf$beta)
  } else {
    c(1, nf$rho)
  }
}

# stationary AR(1) series with unit variance
wr_ar1 <- function(n, rho) {
  e <- stats::rnorm(n)
  e[-1] <- e[-1] * sqrt(1 - rho^2)
  as.numeric(stats::filter(e, rho, method = "recursive"))
}

# power-law (f^-beta) noise, generated on twice the length to avoid wrap-around, unit variance
wr_powerlaw <- function(n, beta) {
  m <- 2 * n
  k <- 1:(m / 2)
  amp <- k^(-beta / 2)
  z <- complex(m)
  z[k + 1] <- complex(real = stats::rnorm(m / 2), imaginary = stats::rnorm(m / 2)) * amp
  kk <- k[-length(k)]
  z[m - kk + 1] <- Conj(z[kk + 1])
  z[m / 2 + 1] <- complex(real = Re(z[m / 2 + 1]))
  y <- Re(stats::fft(z, inverse = TRUE))[1:n]
  (y - mean(y)) / stats::sd(y)
}

wr_noise_series <- function(n, family, param) {
  if (family == 1) {
    wr_ar1(n, param)
  } else {
    wr_powerlaw(n, param)
  }
}

# Testable mask

# 0 = testable, 1 = fewer than min_samples per shortest carrier cycle, 2 = window spans fewer
# than min_cycles of the longest fitted period, 3 = window does not fit in the record
wr_to_testable <- function(sedrate, dx, npts, cfg, min_samples, min_cycles) {
  numsed <- length(sedrate)
  code <- matrix(0L, numsed, ncol(npts))
  for (j in seq_len(numsed)) {
    few_samples <- cfg$carrier_min * sedrate[j] / 100 < min_samples * dx
    dur <- (npts[j, ] - 1) * dx * 100 / sedrate[j]
    v <- integer(ncol(npts))
    v[dur < min_cycles * cfg$longest] <- 2L
    if (few_samples) {
      v[] <- 1L
    }
    v[is.na(npts[j, ])] <- 3L
    code[j, ] <- v
  }
  code
}

# Noise-level tables

wr_logit <- function(p, eps = 1e-6) {
  p[!is.finite(p)] <- NA
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

# moving average over neighbouring rates (half-width h), ignoring NA
wr_smooth_rates <- function(v, h) {
  if (h < 1) {
    return(v)
  }
  n <- length(v)
  out <- rep(NA_real_, n)
  for (j in seq_len(n)) {
    if (!is.finite(v[j])) {
      next
    }
    idx <- max(1, j - h):min(n, j + h)
    out[j] <- mean(v[idx], na.rm = TRUE)
  }
  out
}

wr_node_key <- function(family, param, npts) {
  sprintf("%d_%.4f_%d", as.integer(family), param, as.integer(npts))
}

# Build (or extend) the noise tables. nodes: data.frame(family, param, npts); npts = 0 for time
# windows, where the number of points per rate is npts_rate. Each surrogate is analysed with the
# same timeOpt settings as the data. Returns the cache (list with settings key and raw simulations).
wr_to_null_build <- function(nodes, cache, cfg_key, cfg, sedrate, dx, window_mode, npts_rate,
    linLog, n_null, n_cores, run_multicore, verbose) {
  if (is.null(cache) || !identical(cache$cfg_key, cfg_key)) {
    if (!is.null(cache) && verbose) {
      message("win_timeOpt: null_cache was made with other settings and is not used")
    }
    cache <- list(cfg_key = cfg_key, sims = list())
  }
  keys <- wr_node_key(nodes$family, nodes$param, nodes$npts)
  todo <- integer()
  for (k in seq_along(keys)) {
    have <- cache$sims[[keys[k]]]
    if (is.null(have) || dim(have)[3] < n_null) {
      todo <- c(todo, k)
    }
  }
  if (length(todo) == 0L) {
    return(cache)
  }
  if (verbose) {
    message("win_timeOpt: noise level for ", length(todo), " noise model(s) x ", n_null,
        " surrogates")
  }
  task_node <- rep(todo, each = n_null)
  numsed <- length(sedrate)
  sims <- wr_lapply(task_node, function(k) {
    fam <- nodes$family[k]
    par <- nodes$param[k]
    if (window_mode == "depth") {
      n <- nodes$npts[k]
      return(wr_to_run(cbind((0:(n - 1)) * dx, wr_noise_series(n, fam, par)), sedrate, cfg, linLog))
    }
    out <- matrix(NA_real_, numsed, 2)
    for (j in seq_len(numsed)) {
      n <- npts_rate[j]
      if (is.finite(n) && n >= 16) {
        out[j, ] <- wr_to_run(cbind((0:(n - 1)) * dx, wr_noise_series(n, fam, par)), sedrate[j],
            cfg, linLog)
      }
    }
    out
  }, n_cores = n_cores, run_multicore = run_multicore, combine = "list", progress = verbose)
  for (k in todo) {
    idx <- which(task_node == k)
    A <- array(NA_real_, c(numsed, 2, length(idx)))
    for (q in seq_along(idx)) {
      A[, , q] <- sims[[idx[q]]]
    }
    cache$sims[[keys[k]]] <- A
  }
  cache
}

# mean / sd of logit r2 (smoothed over rates) and the null maxima of z_opt over testable rates
wr_to_null_stats <- function(A, testable_rate, smooth_h) {
  numsed <- dim(A)[1]
  L <- wr_logit(A)
  mu <- matrix(NA_real_, numsed, 2)
  sdv <- mu
  for (q in 1:2) {
    m <- numeric(numsed)
    s <- numeric(numsed)
    for (j in seq_len(numsed)) {
      v <- L[j, q, ]
      if (sum(is.finite(v)) >= 3) {
        m[j] <- mean(v, na.rm = TRUE)
        s[j] <- stats::sd(v, na.rm = TRUE)
      } else {
        m[j] <- NA
        s[j] <- NA
      }
    }
    mu[, q] <- wr_smooth_rates(m, smooth_h)
    sdv[, q] <- pmax(wr_smooth_rates(s, smooth_h), 1e-3)
  }
  # null maxima with leave-one-out mean / sd, so that a noise series is not standardised with
  # statistics it helped to estimate (that would make the threshold too low)
  n_sim <- dim(A)[3]
  S <- matrix(0, numsed, 2)
  SS <- S
  cnt <- S
  for (q in 1:2) {
    for (j in seq_len(numsed)) {
      v <- L[j, q, ]
      v <- v[is.finite(v)]
      S[j, q] <- sum(v)
      SS[j, q] <- sum(v^2)
      cnt[j, q] <- length(v)
    }
  }
  zmax <- rep(NA_real_, n_sim)
  for (k in seq_len(n_sim)) {
    zq <- matrix(NA_real_, numsed, 2)
    for (q in 1:2) {
      x <- L[, q, k]
      ok <- is.finite(x) & cnt[, q] >= 4
      m_k <- rep(NA_real_, numsed)
      s_k <- m_k
      m_k[ok] <- (S[ok, q] - x[ok]) / (cnt[ok, q] - 1)
      v_k <- (SS[ok, q] - x[ok]^2 - (cnt[ok, q] - 1) * m_k[ok]^2) / (cnt[ok, q] - 2)
      s_k[ok] <- sqrt(pmax(v_k, 0))
      m_k <- wr_smooth_rates(m_k, smooth_h)
      s_k <- pmax(wr_smooth_rates(s_k, smooth_h), 1e-3)
      zq[, q] <- (x - m_k) / s_k
    }
    ze <- zq[, 1]
    zp <- zq[, 2]
    zo <- (ze + zp) / sqrt(2)
    zo <- zo[testable_rate & is.finite(zo)]
    if (length(zo)) {
      zmax[k] <- max(zo)
    }
  }
  list(mu = mu, sd = sdv, zmax = zmax[is.finite(zmax)])
}

# interpolation weights of value v on sorted grid g (indices and weights of the two neighbours)
wr_interp_w <- function(g, v) {
  if (length(g) == 1L) {
    return(list(i = c(1L, 1L), w = c(1, 0)))
  }
  v <- min(max(v, g[1]), g[length(g)])
  i <- max(1L, min(findInterval(v, g), length(g) - 1L))
  w2 <- (v - g[i]) / (g[i + 1] - g[i])
  list(i = c(i, i + 1L), w = c(1 - w2, w2))
}
