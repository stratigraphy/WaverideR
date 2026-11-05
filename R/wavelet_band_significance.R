#' Significance of scale-averaged wavelet power (Torrence and Compo, section 4d)
#'
#' Averages the wavelet power over a band of periods (a "scale-averaged wavelet power" time series, eq. 24 of
#' Torrence and Compo, 1998) and tests it against a white or AR(1) red noise background. This is the test to
#' use when one asks "is there significant power in this Milankovitch band at this depth?" without looking at
#' every scale separately: the band average is much less noisy than the power at one scale.
#'
#' Two ways of obtaining the null distribution are offered. Both are analytic and run instantaneously.
#'
#' \describe{
#' \item{\code{method = "exact"}}{(default) The scale-average of the squared modulus of a linear transform of a
#' Gaussian process is a weighted sum of independent exponential variables. The weights are the eigenvalues of a
#' small (number of scales x number of scales) matrix that is computed from the wavelet and the background spectrum
#' with the same Gaussian integrals that define the Morlet wavelet. This is exact for any \code{omega_nr}, any
#' \code{dj}, any band and any background spectrum, and it needs no tabulated constants. The distribution of the
#' weighted sum is approximated by a three-moment (shifted gamma) fit, which was checked against Monte Carlo
#' simulations with the package's own \code{\link{analyse_wavelet}} (see the package tests).}
#' \item{\code{method = "torrence"}}{The published closed form: the band average is a scaled chi-square variable with
#' \eqn{\nu = 2 n_a (S_{avg}/S_{mid}) \sqrt{1 + (n_a \delta j / \delta j_0)^2}} degrees of freedom (eq. 28), where
#' \eqn{n_a} is the number of scales, \eqn{S_{avg}} the harmonic-type mean scale and \eqn{S_{mid}} the
#' geometric mid-scale of the band, and \eqn{\delta j_0} an empirical decorrelation length. Torrence and Compo
#' only tabulated \eqn{\delta j_0 = 0.60} for a Morlet wavelet with \code{omega_nr = 6}. Because the Morlet
#' wavelet becomes narrower in log-scale when \code{omega_nr} increases, that constant under-estimates the degrees
#' of freedom for larger values, so the test becomes too conservative (at \code{omega_nr = 10} only about 2\% instead of 5\% of white noise points exceed the 95\% level). Here the default is rescaled as
#' \eqn{\delta j_0 = 0.60 \cdot 6 / \omega_0}, with \eqn{\omega_0} the value of \code{omega_nr},, which was fitted against the exact method for
#' \code{omega_nr} between 4 and 12 and is within about 10\% of the exact degrees of freedom for bands wider than
#' one octave. Use \code{method = "exact"} when the exact answer matters.}
#' }
#'
#' Only points whose whole period band lies outside the cone of influence are tested.
#'
#' @param object Wavelet to test (list of class \code{analyze.wavelet}), the output of \code{\link{analyse_wavelet}}
#'   (Morlet wavelet; its \code{omega_nr}, \code{dj}, \code{dt}, \code{Scale}, \code{Power} and \code{y} are used).
#'   Required, no default.
#' @param period_range Band of periods to average (numeric vector of length 2, \code{c(lowest, highest)}, in the units
#'   of the depth/time axis; at least one scale of the wavelet must fall inside), for example \code{c(0.3, 0.5)} for an
#'   obliquity band. Required, no default.
#' @param background Noise background of the null hypothesis (character, one of \code{"ar1"}, \code{"white"}):
#'   \code{"ar1"} red noise with the lag-1 autocorrelation \code{lag1}, \code{"white"} white noise (\code{lag1} is then
#'   set to 0). \code{Default="ar1"}.
#' @param siglvl Significance level of the test (numeric, between 0 and 1, e.g. 0.95 for the 95\% level).
#'   \code{Default=0.95}.
#' @param method How the null distribution of the band average is computed (character, one of \code{"exact"},
#'   \code{"torrence"}): \code{"exact"} the weighted sum of exponential variables with a three-moment gamma fit,
#'   \code{"torrence"} the scaled chi-square distribution of Torrence and Compo (1998, eq. 28); see Details.
#'   \code{Default="exact"}.
#' @param lag1 Lag-1 autocorrelation of the AR(1) background (numeric, between 0 and 1). \code{NULL} estimates it from
#'   the record as in \code{\link{wavelet_significance}} (mean of the lag-1 autocorrelation and the square root of the
#'   lag-2 autocorrelation). Ignored for \code{background = "white"}. \code{Default=NULL}.
#' @param delta_j0 Decorrelation length of the scales, only used by \code{method = "torrence"} (numeric, > 0);
#'   \code{NULL} uses \code{0.6 * 6 / omega_nr}, the value 0.6 of Torrence and Compo for \code{omega_nr = 6} rescaled to
#'   the \code{omega_nr} of \code{object}. \code{Default=NULL}.
#' @param genplot Plot the band-averaged power with the significance level (logical, TRUE or FALSE), using the
#'   \code{plot} method. \code{Default=FALSE}.
#' @return An object of class \code{wr_band_significance}, a list with: \code{x} (depth/time axis), \code{power}
#'   (the scale-averaged power, in units of variance when the record is standardised), \code{crit} (significance
#'   level of the band average), \code{pvalue}, \code{significant} (logical, only points outside the cone of
#'   influence), \code{valid}, \code{dof} (equivalent chi-square degrees of freedom), \code{background_mean}
#'   (expected band average under the null), \code{period_range}, \code{n_scales}, \code{lag1},
#'   \code{Cdelta}, \code{delta_j0}, \code{method}, \code{background}, \code{siglvl}.
#'   There are \code{print} and \code{plot} methods.
#' @details The band average is \eqn{\delta j \delta t / C_\delta \sum_j |W_j|^2 / s_j} (eq. 24). The factor
#' \eqn{C_\delta} is evaluated numerically for the \code{omega_nr} and \code{dj} of the object, normalised so that
#' the average over all scales of unit-variance white noise is 1. It only sets the units of \code{power}; the
#' significance decision does not depend on it.
#'
#' The test is for a single point in depth. The band average is strongly autocorrelated along depth, so the
#' number of significant points expected by chance is 5\% of the points, but they occur in runs.
#'
#' @references
#' Torrence, C., Compo, G.P., 1998. A practical guide to wavelet analysis. Bulletin of the American
#' Meteorological Society 79, 61-78.
#' @seealso \code{\link{wavelet_significance}}
#' @examples
#' set.seed(1)
#' z <- seq(0, 60, by = 0.1)
#' d <- data.frame(z, sin(2 * pi * z / 4) + as.numeric(arima.sim(list(ar = 0.5), length(z))))
#' wt <- analyse_wavelet(d, dj = 1 / 10, lowerPeriod = 0.4, upperPeriod = 30,
#'                        verbose = FALSE, omega_nr = 8)
#' bs <- wavelet_band_significance(wt, period_range = c(3, 5))
#' bs
#' @export
wavelet_band_significance <- function(object,
    period_range,
    background = c("ar1", "white"),
    siglvl = 0.95,
    method = c("exact", "torrence"),
    lag1 = NULL,
    delta_j0 = NULL,
    genplot = FALSE) {
  if (!inherits(object, "analyze.wavelet")) {
    stop("object must be the result of analyse_wavelet()")
  }
  background <- match.arg(background)
  method <- match.arg(method)
  if (length(period_range) != 2 || !all(is.finite(period_range)) || period_range[1] >= period_range[2]) {
    stop("period_range must be c(lowest period, highest period)")
  }
  period <- as.numeric(object$Period)
  sel <- which(period >= min(period_range) & period <= max(period_range))
  if (length(sel) < 1) {
    stop("no scales of the wavelet fall inside period_range")
  }
  omega0 <- object$omega_nr
  dt <- if (is.null(object$dt)) { stats::median(diff(object$x)) } else { object$dt }
  dj <- object$dj
  scales <- as.numeric(object$Scale)[sel]
  y <- as.numeric(object$y)

  if (is.null(lag1)) {
    lag1 <- .sig_lag1(y, background)
  }
  if (background == "white") {
    lag1 <- 0
  }
  spec <- function(w) (1 - lag1^2) / (1 + lag1^2 - 2 * lag1 * cos(w * dt))        # unit-variance AR(1) spectrum

  Cd <- .wr_cdelta(omega0, dj)
  wt_avg <- dj * dt / (Cd * scales)                                                 # weights of |W|^2
  Pband <- object$Power[sel, , drop = FALSE]
  power <- as.numeric(dj * dt / Cd * colSums(Pband))                                 # Power is |W|^2 / scale

  valid <- as.logical(.sig_valid(object)[sel[which.max(period[sel])], ])
  dj0 <- NA_real_

  if (method == "exact") {
    ev <- .wr_band_eigen(scales, omega0, dt, wt_avg, spec)
    fit <- .wr_sumexp_fit(ev)
    crit <- .wr_sumexp_q(fit, siglvl)
    pval <- .wr_sumexp_p(fit, power)
    dof <- fit$dof
    bg_mean <- fit$mean
  } else {
    if (is.null(delta_j0)) {
      delta_j0 <- 0.6 * 6 / omega0
    }
    dj0 <- delta_j0
    na <- length(scales)
    fac <- 4 * pi / (omega0 + sqrt(2 + omega0^2))
    fr <- dt / (fac * scales)
    ft <- (1 - lag1^2) / (1 + lag1^2 - 2 * lag1 * cos(2 * pi * fr))
    bg_mean <- sum(wt_avg * ft)
    Savg <- 1 / sum(1 / scales)
    Smid <- exp(mean(range(log(scales))))
    dof <- 2 * na * Savg / Smid * sqrt(1 + (na * dj / dj0)^2)
    crit <- bg_mean * stats::qchisq(siglvl, dof) / dof
    pval <- stats::pchisq(power * dof / bg_mean, dof, lower.tail = FALSE)
  }
  significant <- valid & is.finite(power) & power > crit
  out <- list(x = as.numeric(object$axis.1), power = power, crit = crit, pvalue = pval, significant = significant,
      valid = valid, dof = dof, background_mean = bg_mean, period_range = range(period[sel]),
      n_scales = length(sel), lag1 = lag1, Cdelta = Cd, delta_j0 = dj0, method = method,
      background = background, siglvl = siglvl)
  class(out) <- "wr_band_significance"
  if (genplot) {
    plot(out)
  }
  out
}

#' @export
print.wr_band_significance <- function(x, ...) {
  cat("Significance of scale-averaged wavelet power (", x$method, ", ", x$background, " background, level ",
      x$siglvl, ")\n", sep = "")
  cat("  band: periods ", signif(x$period_range[1], 4), " to ", signif(x$period_range[2], 4), " (",
      x$n_scales, " scales)\n", sep = "")
  cat("  equivalent degrees of freedom:", round(x$dof, 2), "\n")
  cat("  expected band average under the null:", signif(x$background_mean, 4),
      "  critical value:", signif(x$crit, 4), "\n")
  cat("  significant points:", sum(x$significant), "of", sum(x$valid), "outside the cone of influence\n")
  invisible(x)
}

#' @export
#' @param x Result of \code{wavelet_band_significance} (object of class \code{wr_band_significance}) to plot.
#'   Required, no default.
#' @param ... Further graphical arguments passed to \code{\link[graphics]{plot}} (e.g. \code{xlim}, \code{las};
#'   not \code{xlab}, \code{ylab}, \code{ylim} or \code{main}, which are set by the method).
#' @rdname wavelet_band_significance
plot.wr_band_significance <- function(x, ...) {
  ok <- x$valid
  yl <- range(c(x$power[ok], x$crit), na.rm = TRUE)
  graphics::plot(x$x, x$power, type = "n", ylim = yl, xlab = "", ylab = "scale-averaged power",
      main = paste0("Band ", signif(x$period_range[1], 3), " - ", signif(x$period_range[2], 3)), ...)
  if (any(!ok)) {
    r <- rle(!ok)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    for (i in which(r$values)) {
      graphics::rect(x$x[starts[i]], yl[1] - diff(yl), x$x[ends[i]], yl[2] + diff(yl), col = "grey90", border = NA)
    }
  }
  graphics::lines(x$x, x$power, lwd = 1.5)
  graphics::abline(h = x$crit, lty = 2, col = "red")
  graphics::abline(h = x$background_mean, lty = 3, col = "blue")
  sg <- which(x$significant)
  if (length(sg)) {
    graphics::points(x$x[sg], x$power[sg], pch = 16, col = "red", cex = 0.5)
  }
  graphics::legend("topright", c("band average", paste0(100 * x$siglvl, "% level"), "null mean"),
      lty = c(1, 2, 3), col = c("black", "red", "blue"), bty = "n", cex = 0.8)
  invisible(x)
}

# Helpers

# lag-1 autocorrelation of the background, as in wavelet_significance()
.sig_lag1 <- function(y, background) {
  if (background == "white") {
    return(0)
  }
  ac <- stats::acf(y, lag.max = 2, plot = FALSE)$acf[, 1, 1]
  (ac[2] + sqrt(max(ac[3], 0))) / 2
}

# log-spaced angular frequency grid with trapezoid weights, from the lowest relevant frequency up to Nyquist
.wr_omega_grid <- function(smax, dt, ratio = 1.03) {
  lo <- 0.05 / smax
  hi <- pi / dt
  w <- exp(seq(log(lo), log(hi), length.out = max(50, ceiling(log(hi / lo) / log(ratio)))))
  h <- diff(log(w))
  lw <- c(h[1] / 2, (h[-1] + h[-length(h)]) / 2, h[length(h)] / 2)               # trapezoid in log(omega)
  list(w = w, dw = lw * w)
}

# eigenvalues of the quadratic form  sum_j wt_j |W_j|^2  for a Gaussian process with spectrum spec(omega)
# W_j = sum_k xhat_k * psi_hat(s_j omega_k), xhat_k with variance S_k / N  (Torrence and Compo eq. 4-6)
.wr_band_eigen <- function(scales, omega0, dt, wt, spec) {
  g <- .wr_omega_grid(max(scales), dt)
  psi <- matrix(NA_real_, length(g$w), length(scales))
  for (j in seq_along(scales)) {
    psi[, j] <- sqrt(2 * pi * scales[j] / dt) * pi^(-1 / 4) * exp(-(scales[j] * g$w - omega0)^2 / 2)
  }
  psi <- matrix(psi, ncol = length(scales))
  a <- psi * sqrt(spec(g$w) * g$dw * dt / (2 * pi))
  G <- crossprod(a) * outer(sqrt(wt), sqrt(wt))
  ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > max(ev) * 1e-12]
  ev
}

# three moment (shifted gamma) fit of  sum_m lambda_m E_m,  E_m ~ Exp(1)
.wr_sumexp_fit <- function(ev) {
  k1 <- sum(ev)
  k2 <- sum(ev^2)
  k3 <- sum(ev^3)
  theta <- k3 / k2
  alpha <- k2 / theta^2
  list(mean = k1, shape = alpha, scale = theta, shift = k1 - alpha * theta, dof = 2 * k1^2 / k2)
}
.wr_sumexp_q <- function(fit, p) fit$shift + stats::qgamma(p, shape = fit$shape, scale = fit$scale)
.wr_sumexp_p <- function(fit, v) stats::pgamma(pmax(v - fit$shift, 0), shape = fit$shape, scale = fit$scale, lower.tail = FALSE)

# normalising constant of the scale average (eq. 24): the band average over all scales of unit-variance white
# noise equals one.  Evaluated with the same Nyquist-truncated Gaussian integrals as the wavelet transform.
.wr_cdelta <- function(omega0, dj) {
  J <- ceiling(14 / dj)
  s <- 2 * 2^((0:J) * dj)                                                         # dt = 1
  g <- .wr_omega_grid(max(s), 1)
  ew <- numeric(length(s))                                                        # E|W|^2, white noise
  for (i in seq_along(s)) {
    ew[i] <- sum(s[i] / sqrt(pi) * exp(-(s[i] * g$w - omega0)^2) * g$dw)
  }
  dj * sum(ew / s)
}
