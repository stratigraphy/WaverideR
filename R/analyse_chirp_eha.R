#' Chirp-matched evolutive harmonic analysis (de-chirped windowed multitaper spectra)
#'
#' Evolutive harmonic analysis (sliding-window multitaper spectra with the harmonic F-test, Thomson 1982; Meyers
#' 2012) in which, in every window, the depth axis is first warped by a linear change of the sedimentation rate so
#' that cycles whose period changes along the window become lines again. Within a window of length \code{window} the
#' rate is modelled as r(z) = r0 (1 + g (z - z0)); the warped coordinate u = ln(1 + g (z - z0)) / g makes a cycle of
#' fixed duration stationary in u, with the period it has at the window centre z0. For every window a grid of
#' relative rate gradients \code{chirp} (g, per depth unit) is tried, and the gradient that concentrates the spectrum
#' most (lowest Renyi entropy of the whitened power in the analysed period band, or the largest total F-test) is
#' kept. The result is a period-depth map in which each cycle appears as a sharp line that moves with the
#' sedimentation rate, plus the gradient itself (\code{chirp}), an independent estimate of where and how fast the
#' rate changes.
#'
#' The ordinary EHA smears a cycle whose period changes by more than the spectral resolution within the window
#' (with a 20 m window and a rate oscillating by a factor two over 54 m, the precession period changes by up to 60
#' percent inside a window), the wavelet transform smears it over the wavelet length, and the synchrosqueezed
#' transforms only sharpen what the wavelet already resolves. The de-chirping here is the windowed-spectrum
#' counterpart of the matched chirped wavelets of \code{\link{analyse_chirplet_sst}}; with the multitaper estimate it
#' is far more robust to red noise than the squeezed maps.
#'
#' @param data Record: first column depth, second column proxy (uniformly sampled; it is interpolated otherwise).
#' @param window Window length (depth units), or several lengths, whose template scores are then averaged (short
#' windows resolve abrupt changes, long ones resolve the lines). \code{NULL} with \code{criterion = "template"}:
#' windows of about 12 and 25 cycles of the shortest cycle at the geometric mean of \code{sedrate}. A window should
#' hold 15-25 cycles of the cycles that carry the evidence and be shorter than the scale over which the rate changes.
#' @param step Step between window centres \code{Default=window / 20}.
#' @param chirp Relative rate gradients tried (per depth unit): g = (dr/dz) / r at the window centre
#' \code{Default=seq(-0.05, 0.05, by = 0.005)}; \code{"auto"} (with the template): 25 values up to a factor 3 change
#' of the rate within the longest window.
#' @param k Number of Slepian tapers; by default \code{2 * nw - 1}.
#' @param nw Multitaper time-bandwidth product \code{Default=3}; \code{k} tapers \code{Default=2 nw - 1}.
#' @param lowerPeriod,upperPeriod,dj Log2 period grid of the output.
#' @param band Period range (depth units) within which the concentration criterion is evaluated
#' \code{Default=c(lowerPeriod, min(upperPeriod, window / 3))}.
#' @param criterion How the gradient of a window is chosen: \code{"ftest"} (default; sum of the harmonic F-test
#' values above 2 in \code{band}, i.e. the strength of the lines), \code{"template"} (the match with the
#' astronomical template, see \code{astro_periods}; pools all cycles and is the most robust in noisy records),
#' \code{"peak"} (the highest whitened peak) or
#' \code{"entropy"} (Renyi entropy of order 2 of the power in \code{band}; biased towards strong warps in noisy
#' records).
#' @param steer Optional tracked sedimentation rate (a \code{\link{track_astro_hmm}} or \code{\link{track_astro_auto}}
#' result): every window is then de-chirped with the gradient of the track (nearest value of \code{chirp}) instead
#' of the gradient chosen by \code{criterion}. Where the track is right, all cycles become sharp lines that follow
#' the changes of the rate; where it is wrong, they do not, which makes the map a check of the track that is
#' independent of the power the track was fitted to.
#' @param astro_periods,sedrate For \code{criterion = "template"}: the astronomical periods (vector, or data frame
#' with \code{period} and \code{weight}) and the sedimentation-rate range searched; the de-chirped spectrum of every
#' window is then matched to the astronomical template (summed log line-to-continuum ratio at the periods of all
#' cycles, over a grid of rates, as in eCOCO), the gradient with the best match is chosen, and the best rate per
#' window is returned as \code{sedrate} with the score matrix \code{template_score} (windows x rates).
#' @param continuum How the continuum against which a line is measured is estimated from the spectrum of a window:
#' \code{"lowess"} (default; a smooth fit of log power against log period over the whole window spectrum, which is
#' not pulled up by clusters of lines such as the precession triplet) or \code{"median"} (running median over
#' 0.5 octave; it was pulled up inside line clusters, so that isolated aliases scored better than the true rate).
#' @param smooth_chirp The criterion is averaged over this length (depth units) along the record before the gradient
#' is chosen, because the gradient changes slowly and the choice from one noisy window is unreliable
#' \code{Default=window / 2}.
#' @param pad Zero-padding factor of the window spectra \code{Default=8}.
#' @param detrend Remove a linear trend from every window \code{Default=TRUE}.
#' @param verbose Print the progress \code{Default=FALSE}.
#'
#' @return A list of class \code{analyze.wavelet} (usable in \code{\link{plot_wavelet}},
#' \code{\link{track_astro_hmm}}, \code{\link{plot_astro_evidence}}) with \code{Power} (periods x windows, de-chirped
#' multitaper power), \code{Power_norm} (the power divided by the running median of its own window along the period
#' axis, i.e. lines against their local continuum: the map to look at), \code{Ftest}, \code{Period}, \code{axis.1} (window centres), \code{axis.2} (log2 periods),
#' \code{chirp} (the chosen gradient per window), \code{chirp_grid}, \code{criterion} (windows x gradients), the
#' power without de-chirping \code{Power0} and the record.
#'
#' @references Thomson, D.J. (1982). Spectrum estimation and harmonic analysis. Proc. IEEE 70, 1055-1096.
#'
#' Meyers, S.R. (2012). Seeing red in cyclic stratigraphy: spectral noise estimation for astrochronology.
#' Paleoceanography 27, PA3228.
#'
#' Sinnesael, M. et al. (2019). The Cyclostratigraphy Intercomparison Project (CIP): consistency, merits and pitfalls.
#' Earth-Science Reviews 199, 102965.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 3000, by = 1)
#' x_t <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19) + 0.8 * sin(2 * pi * t / 41)
#' sr <- 2 + 0.7 * sin(2 * pi * t / 600)
#' zt <- c(0, cumsum(sr[-1] * diff(t)))
#' z <- seq(0, max(zt), by = 1)
#' x <- approx(zt, x_t, z)$y + rnorm(length(z), sd = 1.5)
#' ce <- analyse_chirp_eha(data.frame(z, x), window = 300, lowerPeriod = 10, upperPeriod = 300)
#' plot_wavelet(ce, whiten = TRUE)
#' }
#' @export
analyse_chirp_eha <- function(data, window = NULL, step = NULL, chirp = seq(-0.05, 0.05, by = 0.005), nw = 3, k = NULL,
                              lowerPeriod = NULL, upperPeriod = NULL, dj = 1 / 50, band = NULL,
                              criterion = c("ftest", "template", "peak", "entropy"), astro_periods = NULL,
                              sedrate = NULL, smooth_chirp = NULL, steer = NULL, pad = 8, detrend = TRUE,
                              continuum = c("lowess", "median"), verbose = FALSE) {
  continuum <- match.arg(continuum)
  criterion <- match.arg(criterion)
  args <- as.list(environment())
  if (criterion == "template" && !is.null(astro_periods) && !is.null(sedrate)) {
    tp0 <- if (is.data.frame(astro_periods)) astro_periods$period else astro_periods
    r0 <- sqrt(prod(sedrate))
    if (is.null(window)) {
      # windows holding about 12 and 25 cycles of the shortest cycle at the geometric mean of the rate range: the
      # short one resolves abrupt changes, the long one resolves the lines; their scores are combined
      window <- c(12, 25) * min(tp0) * r0
      message("analyse_chirp_eha: windows ", paste(signif(window, 3), collapse = " and "))
    }
    if (identical(chirp, "auto")) {
      # the gradient grid: up to a factor 3 change within the longest window
      gmax <- log(3) / max(window)
      chirp <- seq(-gmax, gmax, length.out = 25)
    }
  }
  if (is.null(window)) stop("give a window length")
  if (length(window) > 1) {
    # several windows: the template scores are averaged on the positions of the shortest window
    window <- sort(window)
    if (is.null(step)) step <- window[1] / 20
    res <- lapply(window, function(w) {
      a <- args; a$window <- w; a$step <- step; a$chirp <- chirp
      a <- a[names(a) %in% names(formals(analyse_chirp_eha))]
      do.call(analyse_chirp_eha, a)
    })
    out <- res[[length(res)]]                     # the longest window gives the sharpest lines for display
    pos <- res[[1]]$axis.1
    if (!is.null(out$template_score_all)) {
      TS <- 0
      for (r in res) {
        ts <- r$template_score_all
        TSr <- array(NA_real_, c(dim(ts)[1:2], length(pos)))
        for (gi in seq_len(dim(ts)[2])) TSr[, gi, ] <- t(apply(ts[, gi, , drop = FALSE], 1, function(v) stats::approx(r$axis.1, as.vector(v), pos, rule = 2)$y))
        TS <- TS + TSr
      }
      out$template_score_all <- TS / length(res)
      out$positions <- pos
      out$dt <- step
      out$windows <- window
    }
    return(invisible(out))
  }
  if (criterion == "template") {
    if (is.null(astro_periods)) stop("criterion = 'template' needs astro_periods")
    tw <- if (is.data.frame(astro_periods)) astro_periods$weight else rep(1, length(astro_periods))
    tp <- if (is.data.frame(astro_periods)) astro_periods$period else astro_periods
    if (is.null(tw)) tw <- rep(1, length(tp))
  }
  if (inherits(steer, "wr_track_auto")) steer <- steer$track
  dat <- as.matrix(data)
  dat <- dat[stats::complete.cases(dat[, 1:2]), 1:2, drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  dz <- stats::median(diff(dat[, 1]))
  z <- seq(dat[1, 1], dat[nrow(dat), 1], by = dz)
  x <- stats::approx(dat[, 1], dat[, 2], z)$y
  x <- (x - mean(x)) / stats::sd(x)
  n <- length(x)
  if (is.null(step)) step <- window / 20
  if (is.null(k)) k <- 2 * nw - 1
  nwin <- max(3, round(window / dz))
  if (nwin %% 2 == 0) nwin <- nwin + 1
  half <- (nwin - 1) / 2
  if (is.null(lowerPeriod)) lowerPeriod <- 3 * dz
  if (is.null(upperPeriod)) upperPeriod <- window / 2
  if (is.null(band)) band <- c(lowerPeriod, min(upperPeriod, window / 3))
  if (is.null(smooth_chirp)) smooth_chirp <- window / 2
  J <- as.integer(log2(upperPeriod / lowerPeriod) / dj)
  periods <- lowerPeriod * 2^((0:J) * dj)
  lp <- log2(periods)
  centres <- seq(z[1] + half * dz, z[n] - half * dz, by = step)
  nc <- length(centres)
  # tapers
  E <- .wr_dpss(nwin, nw, k)
  nfft <- 2^ceiling(log2(nwin * pad))
  freq <- (0:(nfft / 2)) / (nfft * dz)
  fsel <- freq > 0
  fper <- 1 / freq[fsel]
  in_band <- fper >= band[1] & fper <= band[2]
  ng <- length(chirp)
  Pw <- matrix(0, length(periods), nc); Ft <- matrix(0, length(periods), nc); P0 <- matrix(0, length(periods), nc)
  crit <- matrix(NA_real_, nc, ng); gbest <- numeric(nc)
  PW <- array(0, c(sum(fsel), ng, nc)); FT <- array(0, c(sum(fsel), ng, nc))
  if (criterion == "template") {
    if (is.null(sedrate)) sedrate <- c(band[1] / min(tp), band[2] / max(tp))
    lr <- seq(log2(sedrate[1]), log2(sedrate[2]), by = dj); rgrid <- 2^lr
    # a cycle outside the period grid at a rate contributes the continuum level (zero) at that rate, and the score
    # is always divided by the full weight: rates at which cycles are not resolved are penalised rather than
    # favoured (averaging over the resolved cycles only favoured the edges of the range)
    RB <- matrix(NA_real_, nc, ng)                           # best rate per window and gradient
    TS <- array(NA_real_, c(length(lr), ng, nc))             # template score per rate
    # interpolation tables of the template: for every rate the fractional index of every cycle on the period grid
    tgt <- outer(lr, log2(tp), "+")                          # rates x cycles (log2 depth period)
    fi <- (tgt - lp[1]) / dj + 1
    t_ok <- fi >= 1 & fi <= length(lp)
    t_i0 <- pmin(pmax(floor(fi), 1), length(lp) - 1); t_w1 <- fi - t_i0
    t_w <- matrix(tw, nrow(tgt), length(tp), byrow = TRUE) * t_ok
    t_wsum <- rep(sum(tw), nrow(tgt))
  }
  kb <- max(3, 2 * floor(0.5 / dj) + 1)
  taper_spec <- function(seg) {
    if (detrend) seg <- stats::residuals(stats::lm(seg ~ seq_along(seg))) else seg <- seg - mean(seg)
    yk <- matrix(0 + 0i, sum(fsel), k)
    for (j in seq_len(k)) {
      f <- stats::fft(c(seg * E$E[, j], rep(0, nfft - nwin)))
      yk[, j] <- f[1:(nfft / 2 + 1)][fsel]
    }
    pw <- rowMeans(Mod(yk)^2)
    ft <- .wr_ftest(E$E, yk)
    list(pw = pw, ft = ft)
  }
  to_grid <- function(v, scale = 1) {
    # interpolate a spectrum on the frequency axis to the log2 period grid; scale = spacing of the (warped)
    # samples relative to dz, so that the periods are in the units of the window centre
    stats::approx(log2(fper * scale), v, lp, rule = 2)$y
  }
  SC <- matrix(1, ng, nc)
  rel <- (-half:half) * dz
  g_steer <- NULL
  if (!is.null(steer)) {
    # the gradient of every window from the tracked sedimentation rate (log rate smoothed over the window length)
    lr <- log(stats::approx(steer$depth, steer$sedrate, z, rule = 2)$y)
    kw <- max(3, 2 * floor(nwin / 2) + 1)
    lrs <- as.numeric(stats::filter(lr, rep(1 / kw, kw), sides = 2)); lrs[is.na(lrs)] <- lr[is.na(lrs)]
    grad <- c(diff(lrs)[1], diff(lrs)) / dz
    g_steer <- stats::approx(z, grad, centres, rule = 2)$y
    gl <- 0.95 / (half * dz)                      # a gradient at which the rate would vanish inside the window
    g_steer <- pmin(pmax(g_steer, -gl), gl)
    chirp <- chirp[abs(chirp) < gl]
    chirp <- sort(unique(c(0, chirp))); ng <- length(chirp)
    crit <- matrix(NA_real_, nc, ng); PW <- array(0, c(sum(fsel), ng, nc)); FT <- array(0, c(sum(fsel), ng, nc))
  }
  for (ci in seq_len(nc)) {
    i0 <- which.min(abs(z - centres[ci]))
    idx <- (i0 - half):(i0 + half)
    seg0 <- x[idx]
    for (gi in seq_len(ng)) {
      g <- chirp[gi]
      if (!is.null(g_steer) && abs(g) > 1e-9 && abs(g - g_steer[ci]) > min(diff(chirp)) / 2 + 1e-12 &&
          which.min(abs(chirp - g_steer[ci])) != gi) next          # steered: only the gradient of the track (and 0)
      if (abs(g) < 1e-9) {
        seg <- seg0; sc <- 1
      } else {
        if (abs(g) * half * dz >= 0.999) { crit[ci, gi] <- NA; next }   # the rate would reach zero inside the window
        u <- log(1 + g * rel) / g                      # warped coordinate, du/dz = 1 at the centre
        ug <- seq(min(u), max(u), length.out = nwin)
        sc <- (ug[2] - ug[1]) / dz                     # spacing of the warped samples relative to dz
        seg <- stats::approx(u, seg0, ug)$y
      }
      s <- taper_spec(seg)
      PW[, gi, ci] <- s$pw; FT[, gi, ci] <- s$ft; SC[gi, ci] <- sc
      if (criterion == "entropy") { p <- s$pw[in_band]; p <- p / sum(p); crit[ci, gi] <- -log(sum(p^2)) }  # Renyi entropy
      if (abs(g) < 1e-9) P0[, ci] <- to_grid(s$pw)
      if (criterion %in% c("ftest", "peak")) {
        # evaluate the criterion in the band of the window centre
        pb <- fper * sc >= band[1] & fper * sc <= band[2]
        crit[ci, gi] <- if (criterion == "ftest") -sum(s$ft[pb][s$ft[pb] > 2]) else -max(s$pw[pb] / mean(s$pw[pb]))
      } else if (criterion == "template") {
        # astronomical template (eCOCO-like) on the de-chirped spectrum: the summed log of the line-to-continuum
        # ratio at the periods of all cycles for every sedimentation rate; the best rate gives the criterion
        pg <- to_grid(s$pw, sc)
        pn <- .wr_line_ratio(pg, lp, continuum)
        v <- (1 - t_w1) * pn[t_i0] + t_w1 * pn[t_i0 + 1]
        sc_r <- rowSums(t_w * v) / t_wsum
        TS[, gi, ci] <- sc_r
        ri <- which.max(sc_r); RB[ci, gi] <- rgrid[ri]; crit[ci, gi] <- -sc_r[ri]
      }
    }
    if (verbose && ci %% 25 == 0) cat("window", ci, "of", nc, "\n")
  }
  # the gradient changes slowly along the record: the criterion is averaged over neighbouring windows before the
  # gradient is chosen (the chirp of a single noisy window is unreliable)
  crit_s <- crit
  if (smooth_chirp > 0 && nc > 2) {
    kw <- max(1, round(smooth_chirp / step)); kw <- 2 * floor(kw / 2) + 1
    if (kw > 1) for (gi in seq_len(ng)) {
      v <- crit[, gi]; v[!is.finite(v)] <- max(crit, na.rm = TRUE)
      f <- stats::filter(v, rep(1 / kw, kw), sides = 2)
      f[is.na(f)] <- v[is.na(f)]
      crit_s[, gi] <- as.numeric(f)
    }
  }
  for (ci in seq_len(nc)) {
    gi <- if (is.null(g_steer)) which.min(crit_s[ci, ]) else which.min(abs(chirp - g_steer[ci]))
    gbest[ci] <- chirp[gi]
    Pw[, ci] <- to_grid(PW[, gi, ci], SC[gi, ci]); Ft[, ci] <- to_grid(FT[, gi, ci], SC[gi, ci])
  }
  # background-normalised power: every window spectrum divided by its own running median along the period axis
  # (half width 0.5 octave), so that a line is measured against the local continuum of its own window (the
  # continuum of a warped window differs from that of the unwarped one, which otherwise shows as a shift of all rows)
  kb <- max(3, 2 * floor(0.5 / dj) + 1)
  Pn <- apply(Pw, 2, function(v) v / pmax(stats::runmed(v, kb, endrule = "median"), 1e-12 * max(v)))
  smat <- matrix(1, length(periods), nc)
  out <- list(Wave = sqrt(Pw) + 0i, Phase = matrix(0, length(periods), nc), Ampl = sqrt(Pw), Power = Pw,
              Power0 = P0, Ftest = Ft, Power_norm = Pn, dt = step, dj = dj, Power.avg = rowMeans(Pw), Period = periods,
              Scale = periods / (4 * pi / (6 + sqrt(2 + 6^2))), coi.1 = NULL, coi.2 = NULL, nc = nc,
              nr = length(periods), axis.1 = centres, axis.2 = lp, omega_nr = NA, wavelet = "chirp_eha",
              fourier_factor = 4 * pi / (6 + sqrt(2 + 6^2)), psi0 = NA, x = z, y = x,
              Power.avg.pval = NULL, Power.pval = NULL, synchrosqueezed = FALSE,
              chirp = gbest, chirp_grid = chirp, criterion = crit, window = window, nw = nw, k = k)
  if (criterion == "template") {
    gi <- vapply(seq_len(nc), function(ci) which(chirp == gbest[ci]), integer(1))
    out$sedrate <- RB[cbind(seq_len(nc), gi)]
    out$sedrate_grid <- rgrid
    out$template_score <- t(sapply(seq_len(nc), function(ci) TS[, gi[ci], ci]))
    out$template_score_all <- TS                             # rates x gradients x windows, for track_chirp_hmm()
  }
  out$args <- args[setdiff(names(args), "data")]
  out$transform <- "analyse_chirp_eha"
  class(out) <- "analyze.wavelet"
  invisible(out)
}

#' Significance of the lines of a (de-chirped) EHA map against warped noise
#'
#' The de-chirped map of \code{\link{analyse_chirp_eha}} warps the noise inside every window as well as the cycles,
#' so wavy structure in it is not by itself evidence. This function runs the same analysis, with the same windows,
#' gradients (or \code{steer}) and settings, on AR1 surrogates of the record and compares the observed
#' line-to-continuum ratios with those of the surrogates cell by cell. Lines that are only artefacts of the warping
#' appear in the surrogates as well and are not significant.
#'
#' @param ce Result of \code{\link{analyse_chirp_eha}}.
#' @param data The record it was computed from.
#' @param n_sim Number of surrogates \code{Default=20}.
#' @param rho Lag-1 autocorrelation of the surrogates \code{Default=NULL}: estimated from the record.
#' @param seed Seed.
#' @return \code{ce} with \code{Power_sig} (observed \code{Power_norm} divided by the 95th percentile of the
#' surrogates: above 1 is significant at that level), \code{Power_pval} (fraction of surrogates at or above the
#' observed value) and \code{surrogate_q95}.
#' @export
chirp_eha_significance <- function(ce, data, n_sim = 20, rho = NULL, seed = 1) {
  if (is.null(ce$args)) stop("ce must come from analyse_chirp_eha()")
  dat <- as.matrix(data)
  dat <- dat[stats::complete.cases(dat[, 1:2]), 1:2, drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  x <- dat[, 2] - mean(dat[, 2]); n <- length(x)
  if (is.null(rho)) rho <- stats::cor(x[-1], x[-n])
  set.seed(seed)
  a <- ce$args; a$verbose <- FALSE
  if (!is.null(a$steer)) a$chirp <- 0                      # steered: only the gradient of the track is computed
  a <- a[names(a) %in% names(formals(analyse_chirp_eha))]
  cnt <- matrix(0, nrow(ce$Power_norm), ncol(ce$Power_norm))
  sims <- array(NA_real_, c(dim(ce$Power_norm), n_sim))
  for (i in seq_len(n_sim)) {
    e <- as.numeric(stats::arima.sim(list(ar = rho), n = n, sd = sqrt(1 - rho^2)))
    e <- e / stats::sd(e) * stats::sd(x)
    r <- do.call(analyse_chirp_eha, c(list(data = data.frame(dat[, 1], e)), a))
    pn <- r$Power_norm
    if (!all(dim(pn) == dim(cnt))) {
      # the surrogate grid can differ by one column: put it on the observed positions
      pn <- t(apply(pn, 1, function(v) stats::approx(r$axis.1, v, ce$axis.1, rule = 2)$y))
    }
    pn[!is.finite(pn)] <- 0
    sims[, , i] <- pn
    cnt <- cnt + (pn >= ce$Power_norm)
  }
  q95 <- apply(sims, c(1, 2), stats::quantile, probs = 0.95, na.rm = TRUE)
  cnt[!is.finite(ce$Power_norm)] <- n_sim
  ce$Power_pval <- cnt / n_sim
  ce$Power_sig <- ce$Power_norm / pmax(q95, 1e-12)
  ce$surrogate_q95 <- q95
  ce
}

# Discrete prolate spheroidal sequences (Slepian tapers) by the tridiagonal eigenproblem, n small
.wr_dpss <- function(n, nw, k) {
  w <- nw / n
  i <- 0:(n - 1)
  dg <- ((n - 1 - 2 * i) / 2)^2 * cos(2 * pi * w)
  od <- (1:(n - 1)) * (n - (1:(n - 1))) / 2
  Tm <- diag(dg)
  Tm[cbind(1:(n - 1), 2:n)] <- od; Tm[cbind(2:n, 1:(n - 1))] <- od
  ev <- eigen(Tm, symmetric = TRUE)
  E <- ev$vectors[, 1:k, drop = FALSE]
  for (j in 1:k) {
    E[, j] <- E[, j] / sqrt(sum(E[, j]^2))
    if (j %% 2 == 1 && sum(E[, j]) < 0) E[, j] <- -E[, j]
  }
  list(E = E)
}

# Harmonic F-test of Thomson (1982) from the eigencoefficients yk (frequencies x tapers)
.wr_ftest <- function(E, yk) {
  k <- ncol(E)
  vn0 <- colSums(E)
  vn0sq <- sum(vn0^2)
  mu <- as.vector(yk %*% vn0) / vn0sq
  fup <- (k - 1) * Mod(mu)^2 * vn0sq
  resid <- yk - outer(mu, vn0)
  fdw <- rowSums(Mod(resid)^2)
  fup / pmax(fdw, 1e-300)
}

# log ratio of a spectrum (on a log2 period grid) to its continuum: a lowess fit of log power against log period
# (robust to narrow lines, follows the red-noise slope) or a running median; + 0.5 keeps the log finite
.wr_line_ratio <- function(pg, lp, continuum = "lowess") {
  lpg <- log(pmax(pg, 1e-12 * max(pg)))
  bg <- if (continuum == "lowess") stats::lowess(lp, lpg, f = 0.4, iter = 2)$y
        else log(pmax(stats::runmed(pg, max(3, 2 * floor(0.5 / (lp[2] - lp[1])) + 1), endrule = "median"), 1e-12 * max(pg)))
  log(exp(lpg - bg) + 0.5)
}
