#' Synchrosqueezed wavelet-chirplet transform
#'
#' Scalogram with chirped Morlet wavelets: at every scale and position the wavelet whose linear chirp (the change of
#' its frequency along the record) matches the signal best is selected, and the coefficient is squeezed to the
#' instantaneous period estimated with that matched wavelet. Where the period of a cycle changes quickly along the
#' record (a changing sedimentation rate), the ordinary wavelet sees a frequency that sweeps across its bandwidth and
#' the ridge blurs; the matched chirped wavelet follows the sweep, so the ridge stays sharp and its instantaneous
#' period estimate is unbiased for a locally linear chirp (synchrosqueezed chirplet transforms: Chen et al. 2023; Li
#' et al. 2026). The selected chirp gives, as a by-product, a map of the local rate of change of the period, i.e. of
#' the sedimentation rate.
#'
#' Several \code{omega_nr} values are combined as in \code{\link{analyse_superlet_sst}} (geometric mean of the
#' squeezed amplitudes). The chirp parameter is dimensionless: \code{chirp = c} means that the frequency of the
#' wavelet changes by \code{c / omega_nr} (relative) over one scale length (the period times omega_nr / 2 pi).
#'
#' @param data Input record: a matrix or data frame with depth/time in column 1 and the proxy value in column 2
#' (further columns are ignored). Incomplete rows are removed, the rows are sorted by depth/time, the record is
#' linearly interpolated to the median sampling interval and standardised. Required (\code{Default=NULL} gives an
#' error).
#' @param lowerPeriod,upperPeriod,dj Period axis: shortest and longest period, in the same unit as the depth/time
#' column of \code{data} (numeric, 0 < \code{lowerPeriod} < \code{upperPeriod}; \code{lowerPeriod} at least twice
#' the sampling interval), and the spacing of the periods in octaves (numeric, > 0; \code{1/x} gives x periods per
#' doubling of the period). Defaults: \code{lowerPeriod}: \code{Default=2}; \code{upperPeriod}:
#' \code{Default=1024}; \code{dj}: \code{Default=1/100}.
#' @param omega_nr Central angular frequencies omega_0 (number of cycles) of the chirped Morlet wavelets of the
#' orders that are combined (numeric vector, values > 0); a single value gives a single-order transform.
#' \code{Default=c(8, 12, 16)}.
#' @param chirp Dimensionless chirp values that are tried at every scale and position (numeric vector; 0 is the
#' ordinary Morlet wavelet, negative and positive values a decreasing or increasing frequency, see the
#' description). The best-matching value (largest modulus) is used; a wider and denser set follows faster changes
#' of the period at a higher computation time. \code{Default=seq(-2, 2, by = 0.25)}.
#' @param n_iter Number of iterations of the frequency reassignment (multisynchrosqueezing, integer, >= 1), see
#' \code{\link{analyse_wavelet_sst}}. \code{Default=4}.
#' @param threshold Noise threshold of the synchrosqueezing: \code{"riskshrink"}, \code{"none"} or a number (numeric,
#' in standard deviations of the standardised record), as in \code{\link{analyse_wavelet_sst}}.
#' \code{Default="riskshrink"}.
#' @param smooth Standard deviation, in octaves, of the Gaussian smoothing along the period axis applied to every
#' order before the orders are combined (numeric, >= 0; 0 switches it off). \code{Default=0.01}.
#' @param combine How the squeezed amplitudes of the orders are combined (character), one of \code{"geometric"}
#' (geometric mean), \code{"min"} (minimum: strictest consensus) or \code{"mean"} (arithmetic mean).
#' \code{Default="geometric"} (the first element of \code{c("geometric", "min", "mean")}).
#' @param verbose Print a message when each order is done (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @return A list with the structure of \code{\link{analyse_wavelet}} (so that \code{\link{plot_wavelet}},
#' \code{\link{track_astro_hmm}} etc. work on it) with in addition \code{Chirp}: the relative rate of change of the
#' instantaneous period along the record (per depth unit; negative where the period shortens, i.e. where the
#' sedimentation rate falls) of the matched wavelet at every period and position, averaged over the orders.
#'
#' @references Chen, Z., Jiang, Q., Lin, J. and Shu, C. (2023). Synchrosqueezed chirplet transforms. Also:
#' Li, S. et al. (2026). High-order synchrosqueezed wavelet-chirplet transform for instantaneous frequency and
#' chirprate estimation. J. Comput. Appl. Math.
#'
#' @examples
#' \donttest{
#' z <- seq(0, 600, by = 0.5)
#' sr <- 1 + z / 600                      # sedimentation rate doubles along the record
#' t <- cumsum(c(0, diff(z) / sr[-1]))
#' x <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 41) + rnorm(length(z), sd = 0.5)
#' cs <- analyse_chirplet_sst(data.frame(z, x), lowerPeriod = 10, upperPeriod = 200)
#' plot_wavelet(cs)
#' }
#' @export
analyse_chirplet_sst <- function(data = NULL,
    lowerPeriod = 2,
    upperPeriod = 1024,
    dj = 1 / 100,
    omega_nr = c(8, 12, 16),
    chirp = seq(-2, 2, by = 0.25),
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
  J <- as.integer(log2(upperPeriod / lowerPeriod) / dj)
  periods <- lowerPeriod * 2^((0:J) * dj)
  ns <- length(periods)
  sb <- smooth / dj
  kh <- if (sb > 0) { max(1, ceiling(3 * sb)) } else { 0 }
  kern <- if (kh > 0) { exp(-((-kh):kh)^2 / (2 * sb^2)) } else { 1 }
  kern <- kern / sum(kern)
  smooth_rows <- function(M) {
    if (kh == 0) {
      return(M)
    }
    Mp <- rbind(matrix(0, kh, ncol(M)), M, matrix(0, kh, ncol(M)))
    out <- matrix(0, nrow(M), ncol(M))
    for (i in seq_along(kern)) {
      out <- out + kern[i] * Mp[i:(i + nrow(M) - 1), , drop = FALSE]
    }
    out
  }
  amps <- vector("list", length(omega_nr))
  chirps <- vector("list", length(omega_nr))
  phasor <- matrix(0 + 0i, ns, n)
  wdef1 <- NULL
  for (i in seq_along(omega_nr)) {
    wdef <- .wr_wavelet_equiv("morlet", omega_nr[i], NULL, ngrid = 8001)
    if (is.null(wdef1)) {
      wdef1 <- wdef
    }
    scales <- periods / wdef$fourier_factor
    r <- .wr_chirplet_core(xs, dt, omega_nr[i], scales, dj, chirp, threshold, n_iter)
    T_i <- r$tx / sqrt(scales)
    amps[[i]] <- smooth_rows(Mod(T_i))
    chirps[[i]] <- r$chirp
    phasor <- phasor + smooth_rows(Re(T_i)) + 1i * smooth_rows(Im(T_i))
    if (verbose) {
      cat("order", i, "omega_0 =", round(omega_nr[i], 2), "done\n")
    }
  }
  amp_max <- numeric(length(amps))
  for (i in seq_along(amps)) {
    amp_max[i] <- max(amps[[i]])
  }
  eps <- 1e-6 * max(amp_max)
  G <- switch(combine,
      geometric = {
        acc <- log(amps[[1]] + eps)
        for (i in seq_along(amps)[-1]) {
          acc <- acc + log(amps[[i]] + eps)
        }
        exp(acc / length(amps)) - eps
      },
      min = {
        acc <- amps[[1]]
        for (i in seq_along(amps)[-1]) {
          acc <- pmin(acc, amps[[i]])
        }
        acc
      },
      mean = {
        acc <- amps[[1]]
        for (i in seq_along(amps)[-1]) {
          acc <- acc + amps[[i]]
        }
        acc / length(amps)
      })
  chirp_sum <- chirps[[1]]
  for (i in seq_along(chirps)[-1]) {
    chirp_sum <- chirp_sum + chirps[[i]]
  }
  power_orders <- vector("list", length(amps))
  for (i in seq_along(amps)) {
    power_orders[[i]] <- amps[[i]]^2
  }
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
      omega_nr = omega_nr[1], wavelet = "morlet", wavelet_order = NULL,
      fourier_factor = wdef1$fourier_factor, psi0 = wdef1$psi0, x = xout, y = xi,
      Power.avg.pval = NULL, Power.pval = NULL, synchrosqueezed = TRUE, superlet = length(omega_nr) > 1,
      chirplet = TRUE, combine = combine, omega_nr_all = omega_nr,
      Chirp = chirp_sum / length(chirps),
      Power_orders = power_orders)
  out$settings <- settings
  out$transform <- "analyse_chirplet_sst"
  class(out) <- "analyze.wavelet"
  invisible(out)
}

# Chirped Morlet in the frequency domain: psi(t) = pi^(-1/4) exp(i w0 t) exp(-t^2 (1 - i c) / 2), so that
# psi_hat(w) = pi^(-1/4) (1 - i c)^(-1/2) exp(-(w - w0)^2 / (2 (1 - i c))); for c = 0 the Morlet wavelet.
.wr_chirplet_core <- function(xs, dt, w0, scales, dj, chirp, threshold, n_iter) {
  n <- length(xs)
  ns <- length(scales)
  pad <- 2^(trunc(log2(n) + 0.5) + 1) - n
  N <- n + pad
  omega.k <- 1:floor(N / 2) * (2 * pi) / (N * dt)
  omega.k <- c(0, omega.k, -omega.k[floor((N - 1) / 2):1])
  fx <- fft(c(xs, rep(0, pad)))
  inv <- function(v) (fft(fx * v, inverse = TRUE) / N)[1:n]
  periods <- scales * 2 * pi / w0
  lp <- log2(periods)
  daughter <- function(s, cc) {
    a <- 1 - 1i * cc
    sqrt(2 * pi * s / dt) * pi^(-1 / 4) * a^(-1 / 2) * exp(-(s * omega.k - w0)^2 / (2 * a))
  }
  nc <- length(chirp)
  wave <- matrix(0 + 0i, ns, n)          # matched coefficient
  wder <- matrix(0 + 0i, ns, n)          # its derivative along the record
  cmap <- matrix(0, ns, n)               # matched chirp (dimensionless)
  for (j in seq_len(ns)) {
    best <- rep(-1, n)
    wb <- rep(0 + 0i, n)
    wd <- rep(0 + 0i, n)
    cb <- rep(0, n)
    for (k in seq_len(nc)) {
      dau <- daughter(scales[j], chirp[k])
      w <- inv(dau)
      m <- Mod(w)
      better <- m > best
      if (any(better)) {
        w1 <- inv(dau * 1i * omega.k)
        best[better] <- m[better]
        wb[better] <- w[better]
        wd[better] <- w1[better]
        cb[better] <- chirp[k]
      }
    }
    wave[j, ] <- wb
    wder[j, ] <- wd
    cmap[j, ] <- cb
  }
  amod <- Mod(wave)
  if (is.character(threshold)) {
    threshold <- match.arg(threshold, c("riskshrink", "none"))
    gam <- if (threshold == "riskshrink") {
      nfine <- min(ns, max(1, round(1 / dj)))
      sqrt(2 * log(n)) * 1.4826 * stats::mad(as.numeric(amod[seq_len(nfine), ]), constant = 1)
    } else {
      sqrt(.Machine$double.eps)
    }
  } else {
    gam <- threshold
  }
  if (is.character(threshold)) {
    gam <- min(gam, stats::quantile(amod, 0.9, names = FALSE))
  }   # see .wr_sst_core
  tx <- matrix(0 + 0i, ns, n)
  bmat <- matrix(NA_integer_, ns, n)
  rel_chirp <- matrix(NA_real_, ns, n)
  for (j in seq_len(ns)) {
    valid <- amod[j, ] > gam
    wv <- ifelse(valid, wave[j, ], NA)
    om <- Im(wder[j, ] / wv)
    pj <- ifelse(is.finite(om) & om > 0, 2 * pi / om, NA_real_)
    b <- round((log2(pj) - lp[1]) / dj) + 1
    b[!is.finite(b) | b < 1 | b > ns] <- NA
    bmat[j, ] <- b
    # chirp c / s^2 is d omega / d b; d ln(period) / d b = -c / (s^2 omega)
    rel_chirp[j, ] <- ifelse(is.finite(om) & om > 0, -cmap[j, ] / (scales[j]^2 * om), NA_real_)
  }
  if (n_iter > 1) {
    cols <- rep(seq_len(n), each = ns)
    for (it in 2:n_iter) {
      cur <- as.vector(bmat)
      nxt <- cur
      okk <- !is.na(cur)
      nxt[okk] <- bmat[cbind(cur[okk], cols[okk])]
      nxt[okk & is.na(nxt)] <- cur[okk & is.na(nxt)]
      bmat <- matrix(nxt, ns, n)
    }
  }
  chirp_out <- matrix(0, ns, n)
  wsum <- matrix(0, ns, n)
  for (j in seq_len(ns)) {
    b <- bmat[j, ]
    ok <- which(!is.na(b))
    if (!length(ok)) {
      next
    }
    idx <- cbind(b[ok], ok)
    tx[idx] <- tx[idx] + wave[j, ok] / sqrt(scales[j])
    wgt <- amod[j, ok]
    rc <- rel_chirp[j, ok]
    good <- is.finite(rc)
    chirp_out[idx[good, , drop = FALSE]] <- chirp_out[idx[good, , drop = FALSE]] + (wgt * rc)[good]
    wsum[idx[good, , drop = FALSE]] <- wsum[idx[good, , drop = FALSE]] + wgt[good]
  }
  chirp_out <- ifelse(wsum > 0, chirp_out / wsum, 0)
  list(tx = tx * sqrt(scales), wave = wave, chirp = chirp_out, gam = gam)
}
