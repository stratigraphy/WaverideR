#' @title Significance of synchrosqueezed (and other) scalograms against autocorrelated noise
#'
#' @description Tests which parts of a synchrosqueezed transform (\code{\link{analyse_wavelet_sst}},
#' \code{\link{analyse_superlet_sst}}, \code{\link{analyse_conceft}}) are unlikely to be produced by noise. The noise of
#' the record is modelled by an autoregressive process (order chosen by AIC, \code{ar_order_max}); optionally its variance
#' may change along the record (a moving standard deviation of the innovations, as in the time-varying AR bootstrap of
#' Chen and Wu, 2025). \code{n_sim} noise records are simulated and transformed with exactly the same settings as the
#' original. For every period the (1 - alpha) quantile of the noise power, pooled over all positions, gives the
#' threshold; with \code{local_variance = TRUE} it is scaled by the local noise variance. Coefficients below the
#' threshold are set to zero in the returned object, so only significant ridges remain in plots and tracking.
#'
#' Under the null hypothesis (no cycles) the record itself is used as the noise model; this is conservative when strong
#' cycles are present, because they increase the variance attributed to noise.
#'
#' @param wavelet Object made by \code{\link{analyse_wavelet_sst}}, \code{\link{analyse_superlet_sst}} or
#' \code{\link{analyse_conceft}}.
#' @param n_sim Number of simulated noise records \code{Default=50}.
#' @param alpha Significance level \code{Default=0.05}.
#' @param ar_order_max Maximum order of the AR noise model \code{Default=5}; 1 gives the classical red noise (AR1).
#' @param local_variance Scale the threshold with the local variance of the noise \code{Default=FALSE}.
#' @param window Length (in depth/time units) of the moving window of the local variance \code{Default=NULL}: twice
#' the longest period of the transform.
#' @param seed Seed of the simulations \code{Default=1}.
#' @param verbose Print messages.
#'
#' @return The input object with \code{Power}, \code{Wave}, \code{Ampl} and \code{Power.avg} restricted to the
#' significant coefficients, and the extra elements \code{significant} (logical matrix), \code{Power_all} (the
#' unmasked power), \code{threshold_period} (noise power threshold per period) and \code{alpha}.
#'
#' @references
#' Chen, Z., Wu, H.-T., 2025. Uncertainty quantification of synchrosqueezing transform under complicated nonstationary
#' noise. arXiv:2506.00779.
#'
#' Torrence, C., Compo, G.P., 1998. A practical guide to wavelet analysis. Bulletin of the American
#' Meteorological Society 79, 61-78.
#'
#' @examples
#' \donttest{
#' z <- 0:999
#' x <- sin(2 * pi * z / 40) + as.numeric(arima.sim(list(ar = 0.7), 1000))
#' s <- analyse_wavelet_sst(data.frame(z, x), lowerPeriod = 10, upperPeriod = 200, dj = 1 / 40)
#' sig <- sst_significance(s, n_sim = 20)
#' plot_wavelet(sig)
#' }
#' @export
sst_significance <- function(wavelet,
                             n_sim = 50,
                             alpha = 0.05,
                             ar_order_max = 5,
                             local_variance = FALSE,
                             window = NULL,
                             seed = 1,
                             verbose = FALSE) {
  if (is.null(wavelet$transform) || is.null(wavelet$settings)) {
    stop("wavelet must be made by analyse_wavelet_sst(), analyse_superlet_sst() or analyse_conceft()")
  }
  fun <- get(wavelet$transform, mode = "function")
  x <- wavelet$y
  z <- wavelet$x
  n <- length(x)
  xs <- (x - mean(x)) / stats::sd(x)
  fit <- stats::ar(xs, aic = TRUE, order.max = ar_order_max, method = "yule-walker")
  innov <- fit$resid
  innov[is.na(innov)] <- 0
  dt <- wavelet$dt
  if (is.null(window)) window <- 2 * max(wavelet$Period)
  L <- max(5, 2 * floor(window / dt / 2) + 1)
  L <- min(L, 2 * floor((n - 1) / 2) + 1)
  local_sd <- sqrt(as.numeric(stats::filter(innov^2, rep(1 / L, L), sides = 2)))
  ok <- which(is.finite(local_sd))
  local_sd <- stats::approx(ok, local_sd[ok], seq_len(n), rule = 2)$y
  local_sd <- local_sd / sqrt(mean(local_sd^2))
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  ns <- length(wavelet$Period)
  # pool the noise power per period: a fixed random subset of positions of every simulation keeps memory small
  take <- unique(round(seq(1, n, length.out = min(n, 400))))
  pool <- matrix(NA_real_, ns, n_sim * length(take))
  sd_innov <- sqrt(fit$var.pred)
  for (i in seq_len(n_sim)) {
    e <- stats::rnorm(n + 200) * sd_innov
    sim <- if (fit$order > 0) as.numeric(stats::filter(e, fit$ar, method = "recursive"))[-(1:200)] else e[-(1:200)]
    if (local_variance) sim <- sim * local_sd
    args <- c(list(data = data.frame(z, sim)), wavelet$settings)
    w <- do.call(fun, args)
    P <- w$Power
    if (local_variance) P <- P / matrix(local_sd^2, ns, n, byrow = TRUE)
    pool[, (i - 1) * length(take) + seq_along(take)] <- P[, take]
    if (verbose) cat("simulation", i, "of", n_sim, "\n")
  }
  thr <- apply(pool, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
  # synchrosqueezed noise power is sparse: when more than a fraction (1 - alpha) of the noise pixels of a period
  # are zero the quantile is zero and every non-zero pixel is significant at that level (use a smaller alpha for a
  # stricter test)
  Tm <- matrix(thr, ns, n)
  if (local_variance) Tm <- Tm * matrix(local_sd^2, ns, n, byrow = TRUE)
  sig <- wavelet$Power > Tm
  out <- wavelet
  out$Power_all <- wavelet$Power
  out$Power <- wavelet$Power * sig
  out$Wave <- wavelet$Wave * sig
  if (!is.null(out$Ampl)) out$Ampl <- out$Ampl * sig
  out$Power.avg <- rowMeans(out$Power)
  out$significant <- sig
  out$threshold_period <- thr
  out$alpha <- alpha
  out$noise_model <- list(ar = fit$ar, order = fit$order, local_variance = local_variance)
  invisible(out)
}
