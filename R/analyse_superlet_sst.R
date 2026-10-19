#' @title Superlet synchrosqueezing: consensus of synchrosqueezed transforms with short and long wavelets
#'
#' @description Combines the idea of superlets (Moca et al. 2021) with synchrosqueezing. A short wavelet
#' (few cycles) follows fast changes of the period, for instance where the sedimentation rate changes, but
#' smears the energy over a broad range of periods; a long wavelet (many cycles) separates close periods but
#' blurs changes along the record. Superlets combine a set of wavelets with an increasing number of cycles by the
#' geometric mean of their amplitudes, which keeps the energy only where all of them agree. Here the wavelet
#' transforms are first synchrosqueezed (first or second order, see \code{\link{analyse_wavelet_sst}}), slightly
#' smoothed along the period axis, and then combined by the geometric mean. A cycle that is present in the
#' record is found by every order at (nearly) the same period and survives; noise is squeezed to different,
#' random periods by the different orders and is strongly suppressed; ridges of the long wavelets that lag behind a
#' change of the period are removed by the short wavelets and vice versa. This gives sharp ridges that follow
#' changes of the sedimentation rate and separate close cycles in noisy records.
#'
#' The defaults (orders with omega_0 = 8, 12 and 16, second-order synchrosqueezing with 4 iterations of the
#' reassignment, smoothing of 0.01 octave) gave the most concentrated and best separated ridges in synthetic tests with
#' eccentricity-modulated precession, obliquity, changing sedimentation rates and red noise (see NEWS).
#'
#' @param data Input is a time series with the first column being depth or time and the second column being a proxy.
#' @param lowerPeriod,upperPeriod Shortest and longest period.
#' @param dj Spacing of the period axis in octaves \code{Default=1/100}.
#' @param c1,o,mult Superlet style definition of the orders, used when \code{omega_nr = NULL}: \code{o} orders with
#' \code{i * c1} cycles (\code{mult = TRUE}) or \code{c1 + i - 1} cycles (\code{mult = FALSE}); cycles are converted to
#' omega_0 = 2 pi cycles / 5, as in \code{\link{analyse_superlet}}.
#' @param omega_nr omega_0 (Morlet equivalent central angular frequencies) of the orders \code{Default=c(8, 12, 16)}.
#' Set to \code{NULL} to use \code{c1}, \code{o} and \code{mult}.
#' @param wavelet,wavelet_order Mother wavelet, see \code{\link{analyse_wavelet}} \code{Default="morlet0"} (the exactly
#' zero-mean Morlet wavelet, which stays admissible for few cycles). The orders of other wavelets get the same relative
#' bandwidth as the Morlet wavelet with the same number of cycles.
#' @param sst_order Order of the synchrosqueezing (1 or 2) \code{Default=2}.
#' @param n_iter Iterations of the frequency reassignment (multisynchrosqueezing) of every order, see
#' \code{\link{analyse_wavelet_sst}} \code{Default=4}.
#' @param threshold Noise threshold of the synchrosqueezing, see \code{\link{analyse_wavelet_sst}}.
#' @param smooth Standard deviation (in octaves) of the Gaussian smoothing along the period axis before the orders are
#' combined. It absorbs the small differences in the period estimated by the different orders \code{Default=0.01}.
#' @param combine \code{"geometric"} (default, as superlets), \code{"min"} (strictest consensus) or \code{"mean"}.
#' @param verbose Print messages.
#'
#' @return A list of class \code{analyze.wavelet} (see \code{\link{analyse_wavelet_sst}}) that can be plotted with
#' \code{\link{plot_wavelet}} and used by the tracking functions, with the extra elements \code{superlet = TRUE},
#' \code{omega_nr_all} (omega_0 of every order) and \code{Power_orders} (the squeezed power of every order, a list).
#'
#' @references
#' Moca, V.V., Barzan, H., Nagy-Dabacan, A., Muresan, R.C., 2021. Time-frequency super-resolution with superlets.
#' Nature Communications 12, 337.
#'
#' Daubechies, I., Lu, J., Wu, H.-T., 2011. Synchrosqueezed wavelet transforms: an empirical mode
#' decomposition-like tool. Applied and Computational Harmonic Analysis 30, 243-261.
#'
#' Oberlin, T., Meignen, S., Perrier, V., 2015. Second-order synchrosqueezing transform or invertible reassignment?
#' Towards ideal time-frequency representations. IEEE Transactions on Signal Processing 63, 1335-1344.
#'
#' @examples
#' \donttest{
#' z <- 0:1499
#' x <- ifelse(z < 750, sin(2 * pi * z / 40), sin(2 * pi * z / 44)) + rnorm(length(z), sd = 1)
#' s <- analyse_superlet_sst(data.frame(z, x), lowerPeriod = 20, upperPeriod = 100, dj = 1 / 50)
#' plot_wavelet(s)
#' }
#' @export
analyse_superlet_sst <- function(data = NULL,
                                 lowerPeriod = 2,
                                 upperPeriod = 1024,
                                 dj = 1 / 100,
                                 c1 = 4,
                                 o = 5,
                                 mult = TRUE,
                                 omega_nr = c(8, 12, 16),
                                 wavelet = "morlet0",
                                 wavelet_order = NULL,
                                 sst_order = 2,
                                 n_iter = 4,
                                 threshold = "riskshrink",
                                 smooth = 0.01,
                                 combine = c("geometric", "min", "mean"),
                                 verbose = FALSE) {
  settings <- mget(setdiff(names(formals(sys.function())), c("data", "verbose")))
  combine <- match.arg(combine)
  dat <- as.matrix(data)
  dat <- dat[stats::complete.cases(dat[, 1:2]), 1:2, drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  dt <- stats::median(diff(dat[, 1]))
  xout <- seq(dat[1, 1], dat[nrow(dat), 1], by = dt)
  xi <- stats::approx(dat[, 1], dat[, 2], xout)$y
  n <- length(xi)
  xs <- (xi - mean(xi)) / stats::sd(xi)
  if (is.null(omega_nr)) {
    cycles <- if (mult) c1 * seq_len(o) else c1 + seq_len(o) - 1
    omega_nr <- 2 * pi * cycles / 5
  }
  J <- as.integer(log2(upperPeriod / lowerPeriod) / dj)
  periods <- lowerPeriod * 2^((0:J) * dj)
  ns <- length(periods)
  # Gaussian smoothing kernel along the period axis
  sb <- smooth / dj
  kh <- if (sb > 0) max(1, ceiling(3 * sb)) else 0
  kern <- if (kh > 0) exp(-((-kh):kh)^2 / (2 * sb^2)) else 1
  kern <- kern / sum(kern)
  smooth_rows <- function(M) {
    if (kh == 0) return(M)
    Mp <- rbind(matrix(0, kh, ncol(M)), M, matrix(0, kh, ncol(M)))
    out <- matrix(0, nrow(M), ncol(M))
    for (i in seq_along(kern)) out <- out + kern[i] * Mp[i:(i + nrow(M) - 1), , drop = FALSE]
    out
  }
  amps <- vector("list", length(omega_nr))
  phasor <- matrix(0 + 0i, ns, n)
  wdef1 <- NULL
  for (i in seq_along(omega_nr)) {
    wdef <- .wr_wavelet_equiv(wavelet, omega_nr[i], wavelet_order, ngrid = 8001)
    if (!wdef$analytic) stop("superlet synchrosqueezing needs an analytic wavelet")
    if (is.null(wdef1)) wdef1 <- wdef
    scales <- periods / wdef$fourier_factor
    r <- .wr_sst_core(xs, dt, wdef, scales, dj, sst_order, "sst", threshold, n_iter)
    T_i <- r$tx / sqrt(scales)                      # reconstruction contributions, comparable between orders
    amps[[i]] <- smooth_rows(Mod(T_i))
    phasor <- phasor + smooth_rows(Re(T_i)) + 1i * smooth_rows(Im(T_i))
    if (verbose) cat("order", i, "omega_0 =", round(omega_nr[i], 2), "done\n")
  }
  eps <- 1e-6 * max(sapply(amps, max))
  G <- switch(combine,
    geometric = exp(Reduce(`+`, lapply(amps, function(a) log(a + eps))) / length(amps)) - eps,
    min = Reduce(pmin, amps),
    mean = Reduce(`+`, amps) / length(amps))
  G[G < 0] <- 0
  unit <- ifelse(Mod(phasor) > 0, phasor / Mod(phasor), 0)
  scales <- periods / wdef1$fourier_factor
  smat <- matrix(rep(sqrt(scales), n), nrow = ns)
  Wave <- G * unit * smat
  Power <- Mod(Wave)^2 / smat^2
  lp <- log2(periods)
  coi <- wdef1$fourier_factor * .wr_coi_factor(wdef1) * dt * c(1e-5, 1:((n + 1) / 2 - 1), rev(1:(n / 2 - 1)), 1e-5)
  coi.x <- c(xout[c(1, 1)] - dt * 0.5, xout, xout[c(n, n)] + dt * 0.5)
  logyint <- lp[2] - lp[1]
  yl <- c(lp[ns] + 0.5 * logyint, lp[1] - 0.5 * logyint)
  coi.y <- c(yl, log2(coi), rev(yl))
  out <- list(Wave = Wave, Phase = Arg(Wave), Ampl = G, Power = Power, dt = dt, dj = dj,
              Power.avg = rowMeans(Power), Period = periods, Scale = scales,
              coi.1 = coi.x, coi.2 = coi.y, nc = n, nr = ns, axis.1 = xout, axis.2 = lp,
              omega_nr = omega_nr[1], wavelet = wdef1$name, wavelet_order = wavelet_order,
              fourier_factor = wdef1$fourier_factor, psi0 = wdef1$psi0, x = xout, y = xi,
              Power.avg.pval = NULL, Power.pval = NULL, synchrosqueezed = TRUE, superlet = TRUE,
              sst_order = sst_order, combine = combine, omega_nr_all = omega_nr,
              Power_orders = lapply(amps, function(a) a^2))
  out$settings <- settings
  out$transform <- "analyse_superlet_sst"
  class(out) <- "analyze.wavelet"
  invisible(out)
}
