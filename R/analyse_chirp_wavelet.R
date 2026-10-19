#' Chirp-matched wavelet scalogram with the astronomical template over rate and rate gradient
#'
#' The wavelet counterpart of \code{\link{analyse_chirp_eha}}. At every scale the record is transformed with a
#' family of chirped Morlet wavelets whose frequency sweep corresponds to a relative gradient of the sedimentation
#' rate g = d ln(rate) / d depth (a cycle of fixed duration has d ln(period) / d depth = g, i.e. the chirp of the
#' matched wavelet at scale s is c = -g s omega_0). The astronomical template is then scored, for every position, rate
#' and gradient, on the scalogram of the wavelets matched to that gradient: the score is the weighted mean log ratio
#' of the power at the depth periods of the cycles to the local continuum (running median along the period axis).
#' \code{\link{track_chirp_hmm}} tracks the rate and the gradient jointly through this score. With \code{steer}, the
#' scalogram of the gradient of a track is returned (every cell with the wavelet matched to the local gradient), as
#' power and, with \code{squeeze = TRUE}, synchrosqueezed with the matched wavelets' instantaneous period.
#'
#' Compared with \code{\link{analyse_chirplet_sst}}, which picks the chirp of the strongest coefficient cell by cell
#' (and therefore follows noise in a noisy record), the chirp here is a property of the sedimentation rate shared by
#' all cycles at a position, so it is decided by all cycles together and by the continuity of the track.
#'
#' @param data Record: first column depth, second column proxy.
#' @param lowerPeriod,upperPeriod,dj Period grid.
#' @param omega_nr Morlet omega_0 \code{Default=10} (long wavelets: the de-chirping pays off most for them).
#' @param chirp Relative rate gradients g (per depth unit) \code{Default=seq(-0.05, 0.05, by = 0.005)}.
#' @param astro_periods Astronomical periods (vector, or data frame with \code{period} and \code{weight}).
#' @param sedrate Rate range searched \code{c(lower, upper)} (depth units per kyr).
#' @param step Spacing of the positions at which the template is scored \code{Default=5 dz}.
#' @param smooth_depth Smoothing of the power along the record before the template is scored (depth units)
#' \code{Default=2 x shortest period x geometric mean rate}.
#' @param local_octaves Half width (octaves) of the running median that defines the local continuum
#' \code{Default=0.75}.
#' @param chirp_max Largest dimensionless chirp of a wavelet (\code{|g| x scale x omega_0}) \code{Default=1.5}: a
#' wavelet that would have to sweep by more than this within its own length is not matched (the plain wavelet is
#' used), because the linear-chirp model does not hold over a wavelet longer than the scale of the rate changes
#' (the 405-kyr wavelet spans a whole 1.2-Myr oscillation in CIP Case 3).
#' @param continuum Continuum estimate of a position's spectrum, see \code{\link{analyse_chirp_eha}}.
#' @param steer Optional track (\code{\link{track_chirp_hmm}}, \code{\link{track_astro_hmm}}): the scalogram of the
#' wavelets matched to its gradient is returned in \code{Power}.
#' @param squeeze With \code{steer}: also synchrosqueeze the steered scalogram \code{Default=FALSE}.
#' @param verbose Print the progress.
#' @return A list of class \code{analyze.wavelet} with \code{Power} (the steered scalogram if \code{steer} is given,
#' otherwise the g = 0 Morlet scalogram), \code{Power_sst} (if squeezed), \code{Period}, \code{axis.1},
#' \code{axis.2}, and for \code{\link{track_chirp_hmm}}: \code{template_score_all} (rates x gradients x positions),
#' \code{sedrate_grid}, \code{chirp_grid}, \code{positions}, \code{dt} (the position step).
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 3000, by = 1)
#' x_t <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19) + 0.8 * sin(2 * pi * t / 41)
#' sr <- 2 + 0.7 * sin(2 * pi * t / 600)
#' zt <- c(0, cumsum(sr[-1] * diff(t)))
#' z <- seq(0, max(zt), by = 1)
#' x <- approx(zt, x_t, z)$y + rnorm(length(z), sd = 1.5)
#' cw <- analyse_chirp_wavelet(data.frame(z, x), lowerPeriod = 10, upperPeriod = 300,
#'     chirp = seq(-0.004, 0.004,
#'                             by = 0.0005), astro_periods = c(41, 23.7, 19), sedrate = c(0.8, 5))
#' tr <- track_chirp_hmm(cw)
#' plot(tr$depth, tr$sedrate, type = "l"); lines(z, approx(zt, sr, z)$y, col = 2)
#' }
#' @export
analyse_chirp_wavelet <- function(data, lowerPeriod, upperPeriod, dj = 1 / 50, omega_nr = 10,
                                  chirp = seq(-0.05, 0.05, by = 0.005), astro_periods, sedrate, step = NULL,
                                  smooth_depth = NULL, local_octaves = 0.75, chirp_max = 1.5, steer = NULL,
                                  squeeze = FALSE, continuum = c("lowess", "median"), verbose = FALSE) {
  continuum <- match.arg(continuum)
  if (inherits(steer, "wr_track_auto")) steer <- steer$track
  tw <- if (is.data.frame(astro_periods)) astro_periods$weight else rep(1, length(astro_periods))
  tp <- if (is.data.frame(astro_periods)) astro_periods$period else astro_periods
  if (is.null(tw)) tw <- rep(1, length(tp))
  dat <- as.matrix(data)
  dat <- dat[stats::complete.cases(dat[, 1:2]), 1:2, drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  dz <- stats::median(diff(dat[, 1]))
  z <- seq(dat[1, 1], dat[nrow(dat), 1], by = dz)
  x <- stats::approx(dat[, 1], dat[, 2], z)$y
  x <- (x - mean(x)) / stats::sd(x)
  n <- length(x)
  J <- as.integer(log2(upperPeriod / lowerPeriod) / dj)
  periods <- lowerPeriod * 2^((0:J) * dj); lp <- log2(periods); np <- length(periods)
  w0 <- omega_nr
  ff <- 4 * pi / (w0 + sqrt(2 + w0^2))
  scales <- periods / ff
  pad <- 2^(trunc(log2(n) + 0.5) + 1) - n
  N <- n + pad
  omega.k <- 1:floor(N / 2) * (2 * pi) / (N * dz)
  omega.k <- c(0, omega.k, -omega.k[floor((N - 1) / 2):1])
  fx <- stats::fft(c(x, rep(0, pad)))
  inv <- function(v) (stats::fft(fx * v, inverse = TRUE) / N)[1:n]
  daughter <- function(s, cc) {
    a <- 1 - 1i * cc
    sqrt(2 * pi * s / dz) * pi^(-1 / 4) * a^(-1 / 2) * exp(-(s * omega.k - w0)^2 / (2 * a))
  }
  chirp <- sort(unique(c(0, chirp))); ng <- length(chirp); g0 <- which(chirp == 0)
  # rates and cycles
  lr <- seq(log2(sedrate[1]), log2(sedrate[2]), by = dj); rgrid <- 2^lr; nr <- length(lr)
  if (is.null(step)) step <- 5 * dz
  cols <- unique(round(seq(1, n, by = max(1, step / dz))))
  nc <- length(cols); positions <- z[cols]
  if (is.null(smooth_depth)) smooth_depth <- 2 * min(tp) * sqrt(prod(sedrate))
  ks <- max(1, round(smooth_depth / dz)); ks <- 2 * floor(ks / 2) + 1
  smooth_cols <- function(v) { if (ks < 3) return(v); f <- as.numeric(stats::filter(v, rep(1 / ks, ks), sides = 2)); f[is.na(f)] <- v[is.na(f)]; f }
  # steering: gradient per position of the record
  g_steer <- NULL
  if (!is.null(steer)) {
    lrt <- log(stats::approx(steer$depth, steer$sedrate, z, rule = 2)$y)
    kw <- max(3, 2 * floor(ks / 2) + 1)
    lrs <- as.numeric(stats::filter(lrt, rep(1 / kw, kw), sides = 2)); lrs[is.na(lrs)] <- lrt[is.na(lrs)]
    g_steer <- c(diff(lrs)[1], diff(lrs)) / dz
    gi_steer <- vapply(g_steer, function(g) which.min(abs(chirp - g)), integer(1))
  }
  # power of the chirped family at the scored positions: periods x gradients x positions (smoothed along depth)
  PWg <- array(0, c(np, ng, nc))
  P0 <- matrix(0, np, n); Pst <- if (!is.null(steer)) matrix(0, np, n) else NULL
  Wst <- if (!is.null(steer) && squeeze) matrix(0 + 0i, np, n) else NULL
  Dst <- if (!is.null(steer) && squeeze) matrix(0 + 0i, np, n) else NULL
  for (j in seq_len(np)) {
    s <- scales[j]
    for (gi in seq_len(ng)) {
      cc <- -chirp[gi] * s * w0                      # chirp of the wavelet matched to the gradient at this scale
      # a wavelet longer than the scale over which the rate changes cannot be matched by a linear chirp (the
      # frequency would sweep by more than a factor e within the wavelet); such cells keep the plain wavelet
      if (abs(cc) > chirp_max) cc <- 0
      dau <- daughter(s, cc)
      w <- inv(dau)
      pw <- Mod(w)^2 / s
      PWg[j, gi, ] <- smooth_cols(pw)[cols]
      if (gi == g0) P0[j, ] <- pw
      if (!is.null(steer)) {
        sel <- gi_steer == gi
        if (any(sel)) {
          Pst[j, sel] <- pw[sel]
          if (squeeze) { Wst[j, sel] <- w[sel]; Dst[j, sel] <- inv(dau * 1i * omega.k)[sel] }
        }
      }
    }
    if (verbose && j %% 50 == 0) cat("scale", j, "of", np, "\n")
  }
  # template score: log ratio of the power to the local continuum (running median along the period axis)
  kb <- max(3, 2 * floor(local_octaves / dj) + 1)
  tgt <- outer(lr, log2(tp), "+")
  fi <- (tgt - lp[1]) / dj + 1
  t_ok <- fi >= 1 & fi <= np                        # a cycle outside the grid contributes the continuum level
  i0 <- pmin(pmax(floor(fi), 1), np - 1); w1 <- fi - i0
  Wm <- matrix(tw, nr, length(tp), byrow = TRUE) * t_ok; wsum <- sum(tw)
  TS <- array(NA_real_, c(nr, ng, nc))
  for (ci in seq_len(nc)) for (gi in seq_len(ng)) {
    pg <- PWg[, gi, ci]
    pn <- .wr_line_ratio(pg, lp, continuum)
    v <- (1 - w1) * pn[i0] + w1 * pn[i0 + 1]
    TS[, gi, ci] <- rowSums(Wm * v) / wsum
  }
  Power <- if (!is.null(steer)) Pst else P0
  out <- list(Wave = sqrt(Power) + 0i, Phase = matrix(0, np, n), Ampl = sqrt(Power), Power = Power, Power0 = P0,
              dt = step, dj = dj, Power.avg = rowMeans(Power), Period = periods, Scale = scales,
              coi.1 = NULL, coi.2 = NULL, nc = n, nr = np, axis.1 = z, axis.2 = lp, omega_nr = w0,
              wavelet = "morlet", fourier_factor = ff, psi0 = pi^(-1 / 4), x = z, y = x,
              Power.avg.pval = NULL, Power.pval = NULL, synchrosqueezed = FALSE,
              template_score_all = TS, sedrate_grid = rgrid, chirp_grid = chirp, positions = positions,
              window = smooth_depth, cycles = tp, weights = tw, chirp = if (!is.null(steer)) g_steer else NULL)
  if (!is.null(steer) && squeeze) {
    # synchrosqueezing of the steered scalogram with the instantaneous period of the matched wavelets
    tx <- matrix(0 + 0i, np, n)
    for (j in seq_len(np)) {
      om <- Im(Dst[j, ] / Wst[j, ])
      pj <- ifelse(is.finite(om) & om > 0, 2 * pi / om, NA_real_)
      b <- round((log2(pj) - lp[1]) / dj) + 1
      ok <- which(is.finite(b) & b >= 1 & b <= np)
      idx <- cbind(b[ok], ok)
      tx[idx] <- tx[idx] + Wst[j, ok] / sqrt(scales[j])
    }
    out$Power_sst <- Mod(tx)^2
  }
  out$transform <- "analyse_chirp_wavelet"
  class(out) <- "analyze.wavelet"
  invisible(out)
}
