#' @title Multitaper synchrosqueezing (ConceFT) of a record
#'
#' @description Computes the "Concentration of Frequency and Time" (ConceFT) representation of Daubechies, Wang and Wu
#' (2016), a multitaper version of the wavelet synchrosqueezing transform that is robust against noise. The synchrosqueezed
#' transform of a noisy record contains spurious short ridges because noise also has an instantaneous frequency.
#' ConceFT uses a set of K orthonormal generalised Morse wavelets (Olhede and Walden 2002) and computes the
#' synchrosqueezed transform for many random linear combinations of them (random unit vectors). The noise in these
#' transforms is (nearly) independent while the ridges of the cycles are the same, so the average keeps the cycles and
#' averages the noise away, with much less broadening than averaging the scalograms of the separate wavelets
#' (\code{method = "multitaper"}, which is also available).
#'
#' @param data Input record: a matrix or data frame with depth/time in column 1 and the proxy value in column 2
#' (further columns are ignored). Incomplete rows are removed, the rows are sorted by depth/time, the record is
#' linearly interpolated to the median sampling interval and standardised. Required (\code{Default=NULL} gives an
#' error).
#' @param lowerPeriod,upperPeriod Shortest and longest period of the transform, in the same unit as the depth/time
#' column of \code{data} (numeric, 0 < \code{lowerPeriod} < \code{upperPeriod}; \code{lowerPeriod} at least twice
#' the sampling interval). Defaults: \code{lowerPeriod}: \code{Default=2}; \code{upperPeriod}: \code{Default=1024}.
#' @param dj Spacing of the period axis in octaves (numeric, > 0; \code{1/x} gives x periods per doubling of the
#' period). \code{Default=1/100}.
#' @param gamma,beta Parameters gamma and beta of the generalised Morse wavelets (numeric, > 0); the defaults are
#' Lilly's default family (Daubechies et al. 2016 used gamma = 9, beta = 30). A larger beta gives longer wavelets
#' (better period resolution, worse depth/time resolution). Defaults: \code{gamma}: \code{Default=3}; \code{beta}:
#' \code{Default=20}.
#' @param K Number of orthonormal Morse wavelets (tapers) of orders 0 to K-1 (integer, >= 1). More tapers suppress
#' noise better but broaden the ridges. \code{Default=2}.
#' @param n_comb Number of random combinations of the tapers whose synchrosqueezed transforms are averaged
#' (integer, >= 1); only used when \code{method="conceft"}. The first combination is always the zero-order Morse
#' wavelet itself. More combinations average out more noise at a higher computation time. \code{Default=20}.
#' @param method Averaging method (character), one of \code{"conceft"} (average over \code{n_comb} random unit-vector
#' combinations of the K tapers) or \code{"multitaper"} (average of the synchrosqueezed transforms of the K
#' orthonormal wavelets themselves). \code{Default="conceft"} (the first element of
#' \code{c("conceft", "multitaper")}).
#' @param sst_order Order of the synchrosqueezing of every combination (numeric, \code{1} or \code{2}), see
#' \code{\link{analyse_wavelet_sst}}. \code{Default=1}.
#' @param n_iter Number of iterations of the frequency reassignment of every combination (integer, >= 1), see
#' \code{\link{analyse_wavelet_sst}}. \code{Default=1}.
#' @param threshold Noise threshold of the synchrosqueezing: \code{"riskshrink"}, \code{"none"} or a number (numeric,
#' in standard deviations of the standardised record), see \code{\link{analyse_wavelet_sst}}.
#' \code{Default="riskshrink"}.
#' @param seed Seed of the random number generator used to draw the random combinations (integer), so that results
#' are reproducible; the previous random state is restored afterwards. Only used when \code{method="conceft"}.
#' \code{Default=1}.
#' @param verbose Print a message after each combination (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @return A list of class \code{analyze.wavelet} (see \code{\link{analyse_wavelet_sst}}) with the averaged power in
#' \code{Power}; \code{Wave} has this power with the phase of the synchrosqueezed transform of the first (zero order)
#' Morse wavelet, so the extraction functions can use it.
#'
#' @references
#' Daubechies, I., Wang, Y., Wu, H.-T., 2016. ConceFT: concentration of frequency and time via a multitapered
#' synchrosqueezed transform. Philosophical Transactions of the Royal Society A 374, 20150193.
#'
#' Olhede, S.C., Walden, A.T., 2002. Generalized Morse wavelets. IEEE Transactions on Signal Processing 50, 2661-2670.
#'
#' Lilly, J.M., Olhede, S.C., 2012. Generalized Morse wavelets as a superfamily of analytic wavelets. IEEE Transactions
#' on Signal Processing 60, 6036-6041.
#'
#' @examples
#' \donttest{
#' z <- 0:999
#' x <- sin(2 * pi * z / 40) + sin(2 * pi * z / 50) + rnorm(length(z), sd = 1.5)
#' cf <- analyse_conceft(data.frame(z, x), lowerPeriod = 20, upperPeriod = 100, dj = 1 / 50,
#'     n_comb = 10)
#' plot_wavelet(cf)
#' }
#' @export
analyse_conceft <- function(data = NULL,
    lowerPeriod = 2,
    upperPeriod = 1024,
    dj = 1 / 100,
    gamma = 3,
    beta = 20,
    K = 2,
    n_comb = 20,
    method = c("conceft", "multitaper"),
    sst_order = 1,
    n_iter = 1,
    threshold = "riskshrink",
    seed = 1,
    verbose = FALSE) {
  settings <- mget(setdiff(names(formals(sys.function())), c("data", "verbose")))
  method <- match.arg(method)
  dat <- as.matrix(data)
  dat <- dat[stats::complete.cases(dat[, 1:2]), 1:2, drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  dt <- stats::median(diff(dat[, 1]))
  xout <- seq(dat[1, 1], dat[nrow(dat), 1], by = dt)
  xi <- stats::approx(dat[, 1], dat[, 2], xout)$y
  n <- length(xi)
  xs <- (xi - mean(xi)) / stats::sd(xi)
  J <- as.integer(log2(upperPeriod / lowerPeriod) / dj)
  periods <- lowerPeriod * 2^((0:J) * dj)
  ns <- length(periods)
  tapers <- vector("list", K)
  for (k in seq_len(K)) {
    tapers[[k]] <- .wr_morse_taper(k - 1, gamma, beta)
  }
  if (method == "multitaper") {
    combos <- diag(K)
  } else {
    old <- if (exists(".Random.seed", envir = globalenv())) { get(".Random.seed", envir = globalenv()) } else { NULL }
    on.exit(if (!is.null(old)) { assign(".Random.seed", old, envir = globalenv()) }, add = TRUE)
    set.seed(seed)
    combos <- matrix(NA_real_, n_comb, K)
    for (i in seq_len(n_comb)) {
      r <- stats::rnorm(K)
      r <- r / sqrt(sum(r^2))
      if (r[1] < 0) {
        r <- -r
      }        # half sphere: the combination is positively correlated with the first taper
      combos[i, ] <- r
    }
    if (K == 1) {
      combos <- matrix(1, n_comb, 1)
    }
    combos[1, ] <- c(1, rep(0, K - 1))   # the first combination is the zero-order wavelet itself (phase reference)
  }
  Pacc <- matrix(0, ns, n)
  ref <- NULL
  for (i in seq_len(nrow(combos))) {
    rr <- combos[i, ]
    f <- function(w) {
      out <- 0
      for (k in seq_len(K)) {
        if (rr[k] != 0) {
          out <- out + rr[k] * tapers[[k]](w)
        }
      }
      out
    }
    wdef <- .wr_wavelet(f, ngrid = 8001)
    scales <- periods / wdef$fourier_factor
    r <- .wr_sst_core(xs, dt, wdef, scales, dj, sst_order, "sst", threshold, n_iter)
    Pacc <- Pacc + Mod(r$tx)^2 / scales
    if (is.null(ref)) {
      ref <- list(r = r, wdef = wdef, scales = scales)
    }
    if (verbose) {
      cat("combination", i, "of", nrow(combos), "\n")
    }
  }
  Power <- Pacc / nrow(combos)
  scales <- ref$scales
  wdef <- ref$wdef
  smat <- matrix(rep(sqrt(scales), n), nrow = ns)
  ph <- ifelse(Mod(ref$r$tx) > 0, ref$r$tx / Mod(ref$r$tx), 0)
  Wave <- sqrt(Power) * smat * ph
  lp <- log2(periods)
  coi <- wdef$fourier_factor * .wr_coi_factor(wdef) * dt * c(1e-5, 1:((n + 1) / 2 - 1), rev(1:(n / 2 - 1)), 1e-5)
  coi.x <- c(xout[c(1, 1)] - dt * 0.5, xout, xout[c(n, n)] + dt * 0.5)
  logyint <- lp[2] - lp[1]
  yl <- c(lp[ns] + 0.5 * logyint, lp[1] - 0.5 * logyint)
  coi.y <- c(yl, log2(coi), rev(yl))
  out <- list(Wave = Wave, Phase = Arg(Wave), Ampl = sqrt(Power), Power = Power, dt = dt, dj = dj,
      Power.avg = rowMeans(Power), Period = periods, Scale = scales,
      coi.1 = coi.x, coi.2 = coi.y, nc = n, nr = ns, axis.1 = xout, axis.2 = lp,
      omega_nr = NA, wavelet = "morse", wavelet_order = c(gamma, beta),
      fourier_factor = wdef$fourier_factor, psi0 = wdef$psi0, x = xout, y = xi,
      Power.avg.pval = NULL, Power.pval = NULL, synchrosqueezed = TRUE, conceft = TRUE,
      method = method, K = K, n_comb = nrow(combos), sst_order = sst_order, n_iter = n_iter)
  out$settings <- settings
  out$transform <- "analyse_conceft"
  class(out) <- "analyze.wavelet"
  invisible(out)
}

# Orthonormal generalised Morse wavelet of order k (Olhede and Walden 2002): psi_hat_k(w) is proportional to
# w^beta exp(-w^gamma) L_k^c(2 w^gamma), c = (2 beta + 1) / gamma - 1, with the generalised Laguerre polynomial L
.wr_morse_taper <- function(k, gamma, beta) {
  cc <- (2 * beta + 1) / gamma - 1
  lognorm <- 0.5 * (log(gamma) + cc * log(2) + lgamma(k + 1) - lgamma(k + cc + 1)) + 0.5 * log(pi) * 0
  function(w) {
    wp <- pmax(w, 0)
    x <- 2 * wp^gamma
    L0 <- rep(1, length(w))
    if (k == 0) {
      L <- L0
    } else {
      L1 <- 1 + cc - x
      Lm <- L0
      L <- L1
      if (k > 1) {
        for (m in 1:(k - 1)) {
          Ln <- ((2 * m + 1 + cc - x) * L - (m + cc) * Lm) / (m + 1)
          Lm <- L
          L <- Ln
        }
      }
    }
    v <- exp(lognorm + beta * log(pmax(wp, 1e-300)) - wp^gamma) * L
    ifelse(w > 0, v, 0)
  }
}
