#' @title Significance of wavelet and superlet power: analytic and surrogate (Monte Carlo) levels
#'
#' @description
#' Calculates the significance level of the spectral power in a wavelet or superlet scalogram against a noise
#' background, and which cells of the scalogram are significant.
#'
#' Two ways are available:
#' \describe{
#' \item{\code{method = "analytic"}}{(wavelets only) the near instantaneous method of Torrence and Compo (1998): the
#' theoretical Fourier spectrum of a lag-1 autoregressive (red noise) process times a chi-square distribution with 2
#' degrees of freedom for the local wavelet power (their equations 16 and 18), and for the global (time averaged)
#' spectrum the chi-square distribution with the degrees of freedom of their equation 23. No simulations are needed.}
#' \item{\code{method = "surrogate"}}{(wavelets and superlets) noise surrogates of the record are transformed in exactly the
#' same way as the data. For every period the \code{siglvl} quantile of the power outside the cone of influence is
#' taken, and averaged over the surrogates (scale wise critical level, after Weesp et al.). The \code{siglvl} quantile of the
#' time averaged power of the surrogates gives the critical level of the global spectrum. The background can be white
#' noise, AR(1) red noise, AR(1) noise estimated after removing the periods of the candidate cycles (so that the cycles
#' that are tested do not inflate the noise background), or a power law (\eqn{1/f^\beta}) noise.}
#' }
#'
#' The result can be added to the plots with the argument \code{add_signif} of \code{\link{plot_wavelet}} and
#' \code{\link{plot_superlet}}.
#'
#' @param object Scalogram to test: a wavelet object (class \code{analyze.wavelet}, from \code{\link{analyse_wavelet}})
#' or a superlet object (class \code{analyze.superlet}, from \code{\link{analyse_superlet}}). Synchrosqueezed
#' transforms are not accepted. Required, no default.
#' @param method How the critical levels are computed (character, \code{"analytic"} or \code{"surrogate"}):
#' \code{"analytic"} the Torrence and Compo (1998) levels (wavelets only, instantaneous), \code{"surrogate"} Monte
#' Carlo levels from \code{n_sim} noise surrogates transformed like the data. \code{NULL} chooses \code{"analytic"} for
#' wavelets and \code{"surrogate"} for superlets. \code{Default=NULL}.
#' @param background Noise background of the null hypothesis (character, one of \code{"ar1"}, \code{"white"},
#' \code{"ar1_corrected"}, \code{"powerlaw"}): \code{"ar1"} red noise with lag-1 autocorrelation \code{lag1};
#' \code{"white"} white noise; \code{"ar1_corrected"} AR(1) noise whose coefficient is estimated after removing the
#' bands around \code{notch_periods}; \code{"powerlaw"} \eqn{1/f^\beta} noise with \eqn{\beta} fitted to the
#' periodogram (outside the \code{notch_periods} bands). The analytic method supports only \code{"ar1"} and
#' \code{"white"}. \code{Default="ar1"}.
#' @param siglvl Significance level (numeric, between 0 and 1, e.g. 0.95 for the 95\% level). \code{Default=0.95}.
#' @param n_sim Number of noise surrogates (integer, >= 1), only used with \code{method = "surrogate"}; more
#' surrogates give more stable levels and take longer. Use 1000 or more for publication. \code{Default=200}.
#' @param lag1 Lag-1 autocorrelation of the red noise (numeric, between 0 and 1), used for \code{"ar1"} and
#' \code{"ar1_corrected"}. \code{NULL} estimates it from the record: maximum likelihood (AR(1) fit) for the
#' surrogates, the lag-1 and lag-2 autocorrelation \eqn{(\alpha_1+\sqrt{\alpha_2})/2} of Torrence and Compo for the
#' analytic method. \code{Default=NULL}.
#' @param notch_periods Periods of the candidate cycles (numeric vector, in the depth/time unit of the record) that are
#' removed before the AR(1) coefficient (\code{"ar1_corrected"}) or the power law exponent (\code{"powerlaw"}) is
#' estimated, so that the cycles that are tested do not inflate the noise background; \code{NULL} removes nothing.
#' \code{Default=NULL}.
#' @param notch_tol Relative half width of the removed bands around \code{notch_periods} (numeric, between 0 and 1;
#' 0.2 removes periods within 20 percent of every notch period). \code{Default=0.2}.
#' @param candidates Optional periods of candidate cycles (numeric vector, optionally named, in the depth/time unit of
#' the record), only used with \code{method = "surrogate"}. For each, the p-value of the largest time averaged power
#' within \code{candidate_tol} of the period is calculated against the surrogates (the largest value in the band is
#' compared with the largest value of every surrogate, which corrects for searching a band). \code{Default=NULL}.
#' @param candidate_tol Relative half width of the candidate bands (numeric, between 0 and 1). \code{Default=0.2}.
#' @param gamma Decorrelation factor of the wavelet used in the degrees of freedom of the global spectrum of the
#' analytic method (numeric, > 0; Torrence and Compo 1998, eq. 23). \code{NULL} takes the value of Torrence and Compo
#' (1998, table 2) for the wavelet of the object: 2.32 (Morlet), 1.17 (Paul), 1.43 (DOG of order 2); 2.32 is used for
#' the other wavelets. \code{Default=NULL}.
#' @param run_multicore Run the surrogates in parallel on several cores (logical, TRUE or FALSE); with \code{TRUE} and
#' \code{n_cores = NULL} a safe number of cores is chosen automatically. \code{Default=FALSE}.
#' @param n_cores Number of cores for the surrogates (integer, >= 1); \code{NULL} uses the option
#' \code{waverider.n_cores} or the choice of \code{run_multicore}. The number is capped by the available cores and
#' memory (see \code{\link{analyse_wavelet}}). \code{Default=NULL}.
#' @param seed Seed of the random number generator for the surrogates (integer), so that the surrogate levels can be
#' reproduced; \code{NULL} derives it from the current random state (so \code{set.seed()} before the call also gives
#' reproducible results). Results do not depend on the number of cores, and the random state of the session is not
#' changed. \code{Default=NULL}.
#' @param verbose Show a progress bar of the surrogate simulations (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @return A list of class \code{"wr_significance"} with
#' \code{period}, \code{crit} (critical power level per period, in the units of \code{object$Power}),
#' \code{gws}, \code{gws_crit} (time averaged power outside the cone of influence and its critical level, per period),
#' \code{significant} (logical matrix with the same layout as \code{object$Power}: power above the critical level and
#' outside the cone of influence), \code{valid} (outside the cone of influence), \code{candidates}
#' (data frame of p-values), and the settings and estimated noise parameters (\code{lag1}, \code{beta}).
#'
#' @references
#' Torrence, C., and G. P. Compo, 1998. A practical guide to wavelet analysis.
#' Bulletin of the American Meteorological Society 79, 61-78.
#'
#' Weesp analysis: surrogate significance with white, AR(1), signal corrected AR(1) and power law backgrounds.
#'
#' Liu, Y., San Liang, X., Weisberg, R. H., 2007. Rectification of the bias in the wavelet power spectrum.
#' Journal of Atmospheric and Oceanic Technology 24, 2093-2102.
#'
#' @examples
#' set.seed(1)
#' z <- seq(0, 60, by = 0.1)
#' d <- data.frame(z, sin(2 * pi * z / 4) + as.numeric(arima.sim(list(ar = 0.5), length(z))))
#' wt <- analyse_wavelet(d, dj = 1 / 10, lowerPeriod = 0.4, upperPeriod = 30,
#'                        verbose = FALSE, omega_nr = 6)
#' sig <- wavelet_significance(wt)                    # analytic, instantaneous
#' sig$period[which(sig$gws > sig$gws_crit)]          # periods with significant global power
#' \donttest{
#' sig2 <- wavelet_significance(wt, method = "surrogate", background = "ar1", n_sim = 20, seed = 1)
#' }
#' @export
#' @importFrom stats qchisq quantile arima acf rnorm runif spec.pgram lm coef
wavelet_significance <- function(object,
    method = NULL,
    background = c("ar1", "white", "ar1_corrected", "powerlaw"),
    siglvl = 0.95,
    n_sim = 200,
    lag1 = NULL,
    notch_periods = NULL,
    notch_tol = 0.2,
    candidates = NULL,
    candidate_tol = 0.2,
    gamma = NULL,
    run_multicore = FALSE,
    n_cores = NULL,
    seed = NULL,
    verbose = FALSE) {
  is_superlet <- inherits(object, "analyze.superlet")
  if (!is_superlet && !inherits(object, "analyze.wavelet")) {
    stop("object must be the result of analyse_wavelet() or analyse_superlet()")
  }
  if (is.null(method)) {
    method <- if (is_superlet) { "surrogate" } else { "analytic" }
  }
  method <- match.arg(method, c("analytic", "surrogate"))
  background <- match.arg(background)
  if (method == "analytic" && is_superlet) {
    stop("the analytic (Torrence and Compo) levels are only valid for wavelets; use method = \"surrogate\" for superlets")
  }
  if (method == "analytic" && !background %in% c("ar1", "white")) {
    stop("method = \"analytic\" supports background = \"ar1\" or \"white\"")
  }

  if (isTRUE(object$synchrosqueezed)){
    stop("significance levels are not defined for a synchrosqueezed transform; use the wavelet transform (analyse_wavelet)")
  }
  wv <- if (is.null(object$wavelet)) { "morlet" } else { object$wavelet }
  if (is.null(gamma)){
    gamma <- switch(wv, morlet = 2.32, paul = 1.17, dog = 1.43, 2.32)
  }
  period <- as.numeric(object$Period)
  nper <- length(period)
  x <- as.numeric(object$x)
  y <- as.numeric(object$y)
  n <- length(x)
  dt <- if (is.null(object$dt)) { stats::median(diff(x)) } else { object$dt }
  Pobs <- .sig_power_pd(object)                              # period x depth
  valid <- .sig_valid(object)

  gws_obs <- numeric(nper)
  for (i in seq_len(nper)) {
    if (any(valid[i, ])) {
      gws_obs[i] <- mean(Pobs[i, valid[i, ]])
    } else {
      gws_obs[i] <- NA_real_
    }
  }
  out <- list(method = method, background = background, siglvl = siglvl, period = period, gws = gws_obs,
      n_sim = if (method == "surrogate") { n_sim } else { NA_integer_ }, lag1 = NA_real_, beta = NA_real_,
      candidates = NULL)

  if (method == "analytic") {
    if (is.null(lag1)) {
      if (background == "white") {
        lag1 <- 0
      } else {
        ac <- stats::acf(y, lag.max = 2, plot = FALSE)$acf[, 1, 1]
        lag1 <- (ac[2] + sqrt(max(ac[3], 0))) / 2
      }
    }
    omega0 <- object$omega_nr
    fac <- if (wv == "morlet") { 4 * pi / (omega0 + sqrt(2 + omega0^2)) } else { object$fourier_factor }   # Fourier period per unit scale (Torrence and Compo, table 1)
    fr <- dt / (fac * object$Scale)
    fft_theor <- (1 - lag1^2) / (1 + lag1^2 - 2 * lag1 * cos(2 * pi * fr))
    crit <- fft_theor * stats::qchisq(siglvl, 2) / 2 / object$Scale     # Power is |W|^2 / scale
    na <- rowSums(valid)
    dof <- 2 * sqrt(1 + (na * dt / (gamma * object$Scale))^2)
    gws_crit <- fft_theor * stats::qchisq(siglvl, dof) / dof / object$Scale
    # periods without a point outside the cone of influence get the instantaneous level (dof = 2), so that the
    # curve is continuous; the observed global spectrum 'gws' stays NA there
    out$lag1 <- lag1
    out$crit <- crit
    out$gws_crit <- gws_crit
  } else {
    pars <- .sig_background_pars(y, period, dt, background, lag1, notch_periods, notch_tol)
    out$lag1 <- pars$phi
    out$beta <- pars$beta
    transform_pd <- .sig_transform(object)
    sim_fun <- .sig_simulator(background)
    one <- function(k) {
      Ps <- transform_pd(sim_fun(n, pars))
      cr <- numeric(nper)
      for (i in seq_len(nper)) {
        v <- Ps[i, if (any(valid[i, ])) { valid[i, ] } else { TRUE }]
        cr[i] <- stats::quantile(v, siglvl, names = FALSE)
      }
      # periods without a point outside the cone of influence: use all points, so the curve stays continuous
      gw <- numeric(nper)
      for (i in seq_len(nper)) {
        gw[i] <- mean(Ps[i, if (any(valid[i, ])) { valid[i, ] } else { TRUE }])
      }
      c(cr, gw)
    }
    sims <- wr_lapply(seq_len(n_sim), one, n_cores = n_cores, run_multicore = run_multicore, combine = "cbind",
        progress = verbose, seed = seed)
    sims <- matrix(sims, nrow = 2 * nper)
    crit_sim <- sims[seq_len(nper), , drop = FALSE]
    gws_sim <- sims[nper + seq_len(nper), , drop = FALSE]
    out$crit <- rowMeans(crit_sim, na.rm = TRUE)
    gws_crit <- numeric(nrow(gws_sim))
    for (i in seq_len(nrow(gws_sim))) {
      gws_crit[i] <- stats::quantile(gws_sim[i, ], probs = siglvl, na.rm = TRUE, names = FALSE)
    }
    out$gws_crit <- gws_crit
    if (!is.null(candidates)) {
      cand_rows <- vector("list", length(candidates))
      names(cand_rows) <- names(candidates)
      for (ci in seq_along(candidates)) {
        p <- candidates[[ci]]
        rr <- which(period >= p * (1 - candidate_tol) & period <= p * (1 + candidate_tol))
        if (!length(rr) || all(is.na(gws_obs[rr]))) {
          next
        }      # the entry stays NULL
        obs_max <- max(gws_obs[rr], na.rm = TRUE)
        sim_max <- numeric(ncol(gws_sim))
        for (j in seq_len(ncol(gws_sim))) {
          sim_max[j] <- max(gws_sim[rr, j], na.rm = TRUE)
        }
        cand_rows[[ci]] <- data.frame(candidate = p, peak_period = period[rr][which.max(gws_obs[rr])], peak_power = obs_max,
            p_max_in_band = (1 + sum(sim_max >= obs_max)) / (1 + length(sim_max)),
            background = background, stringsAsFactors = FALSE)
      }
      out$candidates <- do.call(rbind, cand_rows)
    }
  }
  sig_pd <- (Pobs > out$crit) & valid
  sig_pd[is.na(sig_pd)] <- FALSE
  out$valid <- if (is_superlet) { t(valid) } else { valid }
  out$significant <- if (is_superlet) { t(sig_pd) } else { sig_pd }
  class(out) <- "wr_significance"
  out
}

#' @export
print.wr_significance <- function(x, ...) {
  cat("Significance of wavelet/superlet power (", x$method, ", ", x$background, " background, level ", x$siglvl, ")\n", sep = "")
  if (is.finite(x$lag1)) {
    cat("  lag-1 autocorrelation:", round(x$lag1, 3), "\n")
  }
  if (is.finite(x$beta)) {
    cat("  power law exponent:", round(x$beta, 3), "\n")
  }
  if (is.finite(x$n_sim)) {
    cat("  surrogates:", x$n_sim, "\n")
  }
  cat("  significant cells:", sum(x$significant), "of", sum(x$valid), "outside the cone of influence\n")
  if (!is.null(x$candidates)) {
    cat("  candidate cycles:\n")
    print(x$candidates, row.names = FALSE)
  }
  invisible(x)
}

# Helpers

# power as period x depth for wavelets and superlets
.sig_power_pd <- function(object) {
  if (inherits(object, "analyze.superlet")) {
    t(object$Power)
  } else {
    object$Power
  }
}

# TRUE for cells outside the cone of influence, period x depth
.sig_valid <- function(object) {
  n <- length(object$x)
  edge <- pmin(seq_len(n) - 1, n - seq_len(n))
  if (inherits(object, "analyze.superlet")) {
    dt <- if (is.null(object$dt)) { stats::median(diff(object$x)) } else { object$dt }
    return(outer(object$coi, edge * dt, "<="))
  }
  edge[edge == 0] <- 1e-5
  if (is.null(object$fourier_factor)){
    coi <- (2 * pi / object$omega_nr) * sqrt(2) * object$dt * edge       # as in analyse_wavelet()
  } else {
    wd <- .wr_wavelet(object$wavelet, object$wavelet_order, object$omega_nr)
    coi <- object$fourier_factor * .wr_coi_factor(wd) * object$dt * edge
  }
  outer(as.numeric(object$Period), coi, "<=")
}

# function that returns the period x depth power of a new record, using the settings of the object
.sig_transform <- function(object) {
  x <- as.numeric(object$x)
  if (inherits(object, "analyze.superlet")) {
    function(y) {
      s <- analyse_superlet(cbind(x, y), upperPeriod = max(object$Period), lowerPeriod = min(object$Period),
          Nf = object$nr, c1 = object$c1, o = object$o,
          mult = if (is.null(object$mult)) { TRUE } else { object$mult },
          order_scaling = object$order_scaling, order_alpha = object$order_alpha,
          detrend = if (is.null(object$detrend)) { "linear" } else { object$detrend },
          standardize = if (is.null(object$standardize)) { TRUE } else { object$standardize })
      t(s$Power)
    }
  } else {
    function(y) {
      w <- analyse_wavelet(cbind(x, y), dj = object$dj, lowerPeriod = min(object$Period),
          upperPeriod = max(object$Period), verbose = FALSE, omega_nr = object$omega_nr,
          wavelet = if (is.null(object$wavelet)) { "morlet" } else { object$wavelet },
          wavelet_order = object$wavelet_order)
      w$Power
    }
  }
}

# AR(1) coefficient (maximum likelihood, autocorrelation if that fails)
.sig_phi <- function(y) {
  phi <- tryCatch(as.numeric(stats::arima(y, order = c(1, 0, 0), include.mean = TRUE, method = "ML")$coef["ar1"]),
      error = function(e) NA_real_)
  if (!is.finite(phi)) {
    phi <- stats::acf(y, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
  }
  max(min(phi, 0.99), -0.99)
}

# remove narrow period bands from a record in the Fourier domain
.sig_remove_bands <- function(x, periods, dt, tol) {
  n <- length(x)
  mu <- mean(x)
  f <- stats::fft(x - mu)
  freq <- (0:(n - 1)) / n                      # cycles per sample
  drop <- rep(FALSE, n)
  for (p in periods / dt) {
    lo <- 1 / (p * (1 + tol))
    hi <- 1 / (p * (1 - tol))
    drop <- drop | (freq >= lo & freq <= hi) | (freq >= 1 - hi & freq <= 1 - lo)
  }
  f[drop] <- 0
  Re(stats::fft(f, inverse = TRUE) / n) + mu
}

.sig_background_pars <- function(y, period, dt, background, lag1, notch_periods, notch_tol) {
  pars <- list(mean = mean(y), sd = stats::sd(y), phi = NA_real_, beta = NA_real_)
  if (background == "ar1") {
    pars$phi <- if (is.null(lag1)) { .sig_phi(y) } else { lag1 }
  }
  if (background == "white") {
    pars$phi <- 0
  }
  if (background == "ar1_corrected") {
    yc <- if (length(notch_periods)) { .sig_remove_bands(y, notch_periods, dt, notch_tol) } else { y }
    pars$phi <- if (is.null(lag1)) { .sig_phi(yc) } else { lag1 }
  }
  if (background == "powerlaw") {
    sp <- stats::spec.pgram(y, spans = c(3, 3), taper = 0.1, detrend = TRUE, plot = FALSE)
    per <- dt / sp$freq
    keep <- is.finite(sp$spec) & sp$spec > 0 & per >= min(period) & per <= max(period)
    for (p in notch_periods) {
      keep <- keep & !(per >= p * (1 - notch_tol) & per <= p * (1 + notch_tol))
    }
    if (sum(keep) < 5) {
      stop("too few frequencies to fit the power law exponent")
    }
    pars$beta <- -as.numeric(stats::coef(stats::lm(log(sp$spec[keep]) ~ log(sp$freq[keep])))[2])
  }
  pars
}

.sig_simulator <- function(background) {
  if (background == "powerlaw") {
    return(function(n, par) {
      nn <- 2^ceiling(log2(n))
      fr <- (1:(nn / 2)) / nn
      sp <- complex(modulus = fr^(-par$beta / 2), argument = stats::runif(nn / 2, 0, 2 * pi))
      yy <- Re(stats::fft(c(0, sp, Conj(rev(sp))[-1]), inverse = TRUE))[1:n]
      par$mean + par$sd * (yy - mean(yy)) / stats::sd(yy)
    })
  }
  function(n, par) {
    if (par$phi == 0) {
      stats::rnorm(n, par$mean, par$sd)
    } else {
      wr_coloured_noise(n, par$mean, par$sd, par$phi)
    }
  }
}

# Plotting helpers

# contour of the significant cells on a scalogram that has just been drawn
.wr_add_signif <- function(signif, object, horizontal = TRUE, col = "black", lwd = 1) {
  if (is.null(signif)) {
    return(invisible(NULL))
  }
  if (!inherits(signif, "wr_significance")) {
    stop("add_signif must be the result of wavelet_significance()")
  }
  m <- signif$significant
  if (inherits(object, "analyze.superlet")) {
    m <- t(m)
  }                       # period x depth
  if (nrow(m) != length(object$axis.2) || ncol(m) != length(object$x)) {
    warning("add_signif does not match the scalogram (different size): not drawn")
    return(invisible(NULL))
  }
  m <- matrix(as.numeric(m), nrow(m), ncol(m))
  if (!any(m > 0) || all(m > 0)) {
    return(invisible(NULL))
  }
  if (horizontal) {
    graphics::contour(x = object$x, y = object$axis.2, z = t(m), levels = 0.5, drawlabels = FALSE, add = TRUE,
        col = col, lwd = lwd)
  } else {
    graphics::contour(x = object$axis.2, y = object$x, z = m, levels = 0.5, drawlabels = FALSE, add = TRUE,
        col = col, lwd = lwd)
  }
  invisible(NULL)
}

# axis limits of the average power panel so that the critical level is visible
.wr_gws_lim <- function(signif, power_avg) {
  if (is.null(signif)) {
    return(NULL)
  }
  r <- range(c(power_avg, signif$gws_crit), finite = TRUE)
  if (!all(is.finite(r))) {
    NULL
  } else {
    r
  }
}

# critical level of the global spectrum on the average power panel
.wr_add_signif_gws <- function(signif, horizontal_panel) {
  if (is.null(signif)) {
    return(invisible(NULL))
  }
  if (horizontal_panel) {         # power on the x-axis, period on the y-axis
    graphics::lines(x = signif$gws_crit, y = signif$period, col = "red", lty = 2, lwd = 1.5)
  } else {
    graphics::lines(x = signif$period, y = signif$gws_crit, col = "red", lty = 2, lwd = 1.5)
  }
  invisible(NULL)
}
