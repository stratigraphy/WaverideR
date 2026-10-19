#' Build a track (depth, sedimentation rate, time) from a tracked period curve, a rate curve or an age model
#'
#' The tracking functions return \code{wr_track} objects; this builds one from the curves older workflows produce,
#' so that they can be given to \code{\link{analyse_dechirped}} and \code{\link{check_agemodel}}: a tracked period
#' curve (\code{\link{track_period_wavelet}}: depth and the period of one astronomical cycle in depth units), a
#' sedimentation-rate curve (e.g. from an EHA/sedrate workflow) or an age-depth model.
#'
#' @param x Data frame or matrix: depth in the first column and, in the second, the tracked period (depth units),
#' the sedimentation rate (depth units per kyr, times \code{rate_factor}) or the age (kyr, times \code{rate_factor}).
#' @param type \code{"period"} (default), \code{"sedrate"} or \code{"age"}.
#' @param cycle Period (kyr) of the tracked cycle for \code{type = "period"} \code{Default=405}.
#' @param rate_factor Multiplier that brings the second column to depth units per kyr (rate) or to kyr (age):
#' \code{0.01} for a rate in cm/kyr with depth in m, \code{1000} for an age in Ma \code{Default=1}.
#' @param step Depth step of the output; \code{NULL} keeps the input positions \code{Default=NULL}.
#' @return A \code{wr_track} data frame with \code{depth}, \code{sedrate} and \code{time} (kyr from the top of
#' the curve, increasing with depth).
#' @examples
#' curve <- data.frame(depth = 0:100, period = 12 + 3 * sin(0:100 / 10))
#' tr <- as_wr_track(curve, "period", cycle = 405)
#' head(tr)
#' @export
as_wr_track <- function(x, type = c("period", "sedrate", "age"), cycle = 405, rate_factor = 1, step = NULL) {
  type <- match.arg(type)
  x <- as.data.frame(x)[, 1:2]; names(x) <- c("depth", "v")
  x <- x[is.finite(x$depth) & is.finite(x$v), ]; x <- x[order(x$depth), ]
  x <- stats::aggregate(v ~ depth, data = x, FUN = mean)
  if (nrow(x) < 2) stop("need at least two finite points")
  if (type == "age") {
    t <- x$v * rate_factor; t <- abs(t - t[1])
    if (any(diff(t) <= 0)) stop("the age model must be strictly monotonic")
    r <- diff(x$depth) / diff(t); r <- c(r[1], r)
    r <- (c(r[-1], r[length(r)]) + r) / 2
    out <- data.frame(depth = x$depth, sedrate = r, time = t)
  } else {
    r <- if (type == "period") x$v / cycle else x$v * rate_factor
    if (any(r <= 0)) stop("rates must be positive")
    out <- data.frame(depth = x$depth, sedrate = r)
  }
  if (!is.null(step)) {
    z <- seq(min(out$depth), max(out$depth), by = step)
    out <- data.frame(depth = z, sedrate = exp(stats::approx(out$depth, log(out$sedrate), z)$y))
  }
  if (is.null(out$time) || !is.null(step)) {
    dz <- diff(out$depth); rm <- 2 / (1 / out$sedrate[-1] + 1 / out$sedrate[-nrow(out)])
    out$time <- c(0, cumsum(dz / rm))
  }
  class(out) <- c("wr_track", "data.frame")
  out
}

#' Check an age model against the de-chirped record: where do the astronomical cycles bend away?
#'
#' An existing track (a tracked period curve, a sedimentation-rate curve or an age model, see
#' \code{\link{as_wr_track}}) is used to de-chirp the record with \code{\link{analyse_dechirped}}: in the time of the
#' track every astronomical cycle should be a horizontal line at its period. Wherever the track is wrong the lines bend
#' away from the expected periods, and the local offset of the ridge gives the correction of the sedimentation rate
#' that the record itself asks for: an observed period \code{P_obs} where \code{P} was expected means the true rate is
#' \code{P_obs / P} times the tracked rate (the track ran through that interval too fast or too slowly by that factor). The offset is measured per cycle and per transform (the smooth Morlet
#' wavelet, the multitaper EHA, the superlet-SST, the chirplet-SST and the combined map each see the ridges
#' differently). The offset is found by a comb search: a wrong rate shifts every cycle by the same factor, so the
#' log-power at all expected periods is summed for every common shift of the comb (as the harmonic product spectrum
#' of pitch detection does) and the best shift is the correction; the cycles are then checked one by one at that
#' shift. The corrections are combined by a weighted median and an interval is flagged when the combined
#' correction exceeds \code{tolerance} with at least \code{min_cycles} cycles and \code{min_transforms} transforms
#' agreeing on its direction. Optionally a blind chirp-matched re-track (\code{\link{track_astro_auto}}) is run for
#' comparison.
#'
#' The measured correction is only as good as the ridges: a cycle that is not expressed in the record gives no
#' offset (its strength falls below \code{min_strength} and it is dropped), and an interval flagged by one weak
#' cycle alone is not flagged. The flagged intervals are a list of places to look at in the maps of
#' \code{\link{plot_agemodel_check}}, not a new age model; use \code{\link{track_astro_auto}} or
#' \code{\link{track_chirp_hmm}} to re-track.
#'
#' @param data Data frame or matrix with depth and the proxy.
#' @param track The existing track (\code{wr_track} or \code{\link{track_astro_auto}} result).
#' @param age_ma Age (Ma) used to pick the astronomical periods with \code{\link{astro_periods_tracking}} when
#' \code{astro_periods} is not given \code{Default=0}.
#' @param astro_periods Periods (kyr) of the cycles to measure; by default the 405, 124 and 95 kyr eccentricity
#' cycles, the main obliquity term and the two strongest precession terms for \code{age_ma}.
#' @param transforms Which transforms of \code{\link{analyse_dechirped}} to use
#' \code{Default=c("wavelet", "eha", "superlet_sst", "chirplet_sst")}.
#' @param band Half-width (octaves) of the band around each (shifted) expected period in which a cycle's own ridge is
#' looked for; it is narrowed where two cycles are closer than that \code{Default=0.3}.
#' @param max_shift Largest common shift (octaves) of the comb that is searched, i.e. a rate error of up to a factor
#' \code{2^max_shift}; a best shift at the edge of the grid is discarded \code{Default=1.2}.
#' @param min_contrast Minimum contrast of the best shift (its score above the median over all shifts, in MADs) for a
#' column to count \code{Default=1.5}.
#' @param tolerance Relative rate correction above which an interval is flagged \code{Default=0.1}.
#' @param min_cycles,min_transforms Minimum number of cycles and transforms supporting a flag \code{Default=2}.
#' @param min_strength Minimum ridge prominence (ratio of the band maximum to the median of the column) for an
#' offset to count; squeezed transforms use a tenth of it on their sparse maps \code{Default=1.5}.
#' @param eha_window Window (kyr) of the multitaper EHA in the time of the track; by default 2.5 times the longest
#' cycle, so that the 405-kyr cycle is resolved \code{Default=NULL}.
#' @param harmonics Add the 203 and 135 kyr harmonics of the 405-kyr cycle to the comb (weight 0.4), so that a
#' non-sinusoidal long eccentricity cycle is not read as shifted 124 and 95 kyr cycles \code{Default=TRUE}.
#' @param clip Restrict the record to the depth range of the track (beyond it the track is extrapolated with its
#' end rates, which can dominate the time axis) \code{Default=TRUE}.
#' @param min_length Shortest flagged interval kept and shortest gap merged, in kyr of the track; by default half the
#' longest cycle \code{Default=NULL}.
#' @param n_max Records longer than this are resampled in depth before the transforms \code{Default=6000}.
#' @param retrack Run a blind chirp-matched re-track for comparison \code{Default=TRUE}.
#' @param sedrate Rate range for the blind re-track; by default half to twice the range of the existing track.
#' @param dj Period resolution (octaves) \code{Default=1/50}.
#' @param verbose Print progress \code{Default=TRUE}.
#' @param ... Passed to \code{\link{analyse_dechirped}}.
#' @return A list of class \code{wr_agemodel_check}: \code{profile} (data frame with \code{depth}, \code{time},
#' \code{sedrate} of the track, one column of log2 rate correction per transform, \code{correction} (the combined
#' factor), \code{n_cycles}, \code{n_transforms}, \code{flag}), \code{intervals} (data frame of the flagged
#' intervals: \code{from}, \code{to}, \code{correction}, \code{cycles}, \code{transforms}), \code{offsets} (per
#' transform: \code{shift}, the comb shift per time column; \code{offset}, the log2 offset of every cycle at that
#' shift; \code{offset_single}, the offset of every cycle searched on its own around its expected period), \code{dechirped} (the \code{wr_dechirped} objects),
#' \code{retrack} (the blind track, if run), \code{astro_periods}, \code{track}, \code{data}.
#' @seealso \code{\link{plot_agemodel_check}}, \code{\link{analyse_dechirped}}, \code{\link{as_wr_track}}
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 3000, 1)
#' x <- cos(2 * pi * t / 405) + 0.8 * cos(2 * pi * t / 100) + 0.6 * cos(2 * pi * t / 41) +
#'   rnorm(length(t), 0, 0.8)
#' rate <- ifelse(t < 1500, 0.05, 0.08)      # a step in the rate the age model does not know
#' z <- cumsum(rate); data <- data.frame(z, x)
#' track <- as_wr_track(data.frame(z, 0.065), "sedrate")   # a constant-rate age model
#' chk <- check_agemodel(data, track, age_ma = 0, transforms = c("wavelet", "superlet_sst"),
#'     retrack = FALSE)
#' chk$intervals
#' plot_agemodel_check(chk)
#' }
#' @export
check_agemodel <- function(data, track, age_ma = 0, astro_periods = NULL,
                           transforms = c("wavelet", "eha", "superlet_sst", "chirplet_sst"),
                           band = 0.3, max_shift = 1.2, tolerance = 0.1, min_cycles = 2, min_transforms = 2, min_strength = 1.5,
                           min_contrast = 1.5, eha_window = NULL, n_max = 6000, clip = TRUE, min_length = NULL, harmonics = TRUE, retrack = TRUE, sedrate = NULL, dj = 1 / 50, verbose = TRUE, ...) {
  if (inherits(track, "wr_track_auto")) track <- track$track
  dat <- as.data.frame(data)[, 1:2]; names(dat) <- c("depth", "x")
  dat <- dat[stats::complete.cases(dat), ]; dat <- dat[order(dat$depth), ]
  if (clip) dat <- dat[dat$depth >= min(track$depth) & dat$depth <= max(track$depth), ]
  if (nrow(dat) > n_max) {
    z <- seq(min(dat$depth), max(dat$depth), length.out = n_max)
    k <- ceiling(nrow(dat) / n_max)
    if (k > 1) dat$x <- stats::filter(dat$x, rep(1 / k, k), sides = 2)
    ok <- is.finite(dat$x)
    dat <- data.frame(depth = z, x = stats::approx(dat$depth[ok], dat$x[ok], z, rule = 2)$y)
  }
  if (is.null(astro_periods)) {
    at <- astro_periods_tracking(age_ma)
    fam <- if (!is.null(at$family)) at$family else rep(c("eccentricity", "obliquity", "precession"), c(5, 3, 4))
    pr <- at$period[fam == "precession"]; pw <- at$weight[fam == "precession"]
    pr <- pr[order(-pw)][1:2]
    ob <- at$period[fam == "obliquity"][which.max(at$weight[fam == "obliquity"])]
    astro_periods <- c(405.7, 123.9, 94.9, ob, sort(pr, decreasing = TRUE))
  }
  P <- sort(astro_periods, decreasing = TRUE)
  dtt <- diff(range(stats::approx(track$depth, track$time, dat$depth, rule = 2)$y)) / (nrow(dat) - 1)
  if (any(P < 2.5 * dtt)) {
    if (verbose) message("check_agemodel: cycles below the resolution of the record in the time of the track dropped: ", paste(round(P[P < 2.5 * dtt], 1), collapse = ", "), " kyr (sample spacing ", round(dtt, 1), " kyr)")
    P <- P[P >= 2.5 * dtt]
  }
  if (length(P) < 2) stop("fewer than two astronomical cycles are resolvable in the time of the track")
  lP <- log2(P)
  # comb weights: the 405-kyr cycle, obliquity and the main precession term count fully, the shorter eccentricity
  # terms and the second precession term less; optional harmonics of the 405-kyr cycle (a non-sinusoidal long
  # eccentricity cycle puts ridges at 203 and 135 kyr, which must not be read as shifted 124 and 95 kyr cycles)
  cw <- rep(1, length(P))
  cw[abs(P - 123.9) < 5 | abs(P - 94.9) < 5] <- 0.6
  if (sum(P < 25) >= 2) cw[which(P < 25)[-1]] <- 0.8
  Ph <- P; wh <- cw
  if (harmonics && any(abs(P - 405.7) < 5)) { Ph <- c(Ph, 405.7 / 2, 405.7 / 3); wh <- c(wh, 0.4, 0.4) }
  lPh <- log2(Ph)
  # band of every cycle, narrowed to half the distance to its neighbours
  bw <- rep(band, length(P))
  if (length(P) > 1) {
    d <- abs(diff(lP)) / 2
    bw <- pmin(bw, c(d, Inf), c(Inf, d))
  }
  if (is.null(eha_window)) eha_window <- min(2.5 * max(P), diff(range(track$time)) / 3)
  lower <- min(P) * 2^(-max_shift) * 0.9; upper <- max(P) * 2^max_shift * 1.1
  sparse_tf <- c("superlet_sst", "wavelet_sst", "chirplet_sst")
  dcs <- list(); offs <- list(); corr <- list()
  for (tf in transforms) {
    if (verbose) message("check_agemodel: ", tf)
    dc <- analyse_dechirped(dat, track, tf, lowerPeriod = lower, upperPeriod = upper, dj = dj, window = eha_window, ...)
    dcs[[tf]] <- dc
    Pt <- dc$tuned$Power; lp <- dc$tuned$axis.2; ax <- dc$tuned$axis.1
    sparse <- tf %in% sparse_tf
    Pt[!is.finite(Pt)] <- 0
    if (!sparse) {                                             # lines against the continuum of every column
      Pt <- apply(Pt, 2, function(v) { ok <- v > 0; if (sum(ok) < 10) return(v * 0)
        bg <- rep(NA_real_, length(v)); bg[ok] <- exp(stats::lowess(lp[ok], log(v[ok]), f = 0.5, iter = 2)$y)
        r <- v / bg; r[!is.finite(r)] <- 0; r })
    }
    if (sparse) {
      # smooth the squeezed map over a tenth of an octave so that a ridge a bin away still scores
      kk <- max(3, 2 * floor(0.05 / dj) + 1)
      Pt <- apply(Pt, 2, function(v) stats::filter(v, rep(1 / kk, kk), sides = 2, circular = FALSE))
      Pt[!is.finite(Pt)] <- 0
      Pt <- Pt / max(Pt, 1e-12)
    }
    # the comb search: a wrong rate shifts every cycle by the same factor, so the log-power at all expected periods
    # is summed for every common shift of the comb (the harmonic product spectrum of pitch detection) and the best
    # shift is the correction; cycles are then checked individually at that shift
    shifts <- seq(-max_shift, max_shift, by = dj)
    eps <- if (sparse) 1e-3 else 0.05
    logP <- log(Pt + eps)
    ws <- wh; S <- matrix(NA_real_, length(shifts), ncol(Pt))
    for (si in seq_along(shifts)) {
      v <- 0; wsum <- 0
      for (i in seq_along(Ph)) {
        r <- (lPh[i] + shifts[si] - lp[1]) / dj + 1
        if (r < 1 || r > nrow(logP)) next
        r0 <- floor(r); fr <- r - r0; r1 <- min(r0 + 1, nrow(logP))
        v <- v + ws[i] * ((1 - fr) * logP[r0, ] + fr * logP[r1, ]); wsum <- wsum + ws[i]
      }
      S[si, ] <- if (wsum >= 0.5 * sum(ws)) v / wsum else NA
    }
    S[!is.finite(S)] <- -Inf
    best <- apply(S, 2, which.max); sh <- shifts[best]
    for (j in seq_len(ncol(S))) {                               # parabolic refinement
      b <- best[j]
      if (b > 1 && b < length(shifts)) { a <- S[b - 1, j]; m <- S[b, j]; cc <- S[b + 1, j]; den <- a - 2 * m + cc
        if (is.finite(den) && den < 0) sh[j] <- sh[j] + 0.5 * (a - cc) / den * dj }
    }
    # contrast: the best shift against the median over all shifts, per column
    contrast <- apply(S, 2, function(v) { v <- v[is.finite(v)]; if (length(v) < 5) return(NA); (max(v) - stats::median(v)) / max(1e-9, stats::mad(v)) })
    sh[!is.finite(contrast) | contrast < min_contrast | best <= 1 | best >= length(shifts)] <- NA
    # cycles individually: is there a ridge within band/2 of the shifted period?
    O <- matrix(NA_real_, length(P), ncol(Pt)); Sg <- O
    colmed <- apply(Pt, 2, function(v) { v <- v[v > 0]; if (length(v)) stats::median(v) else NA })
    for (i in seq_along(P)) for (j in seq_len(ncol(Pt))) {
      if (!is.finite(sh[j])) next
      rows <- which(abs(lp - (lP[i] + sh[j])) <= bw[i])
      if (length(rows) < 3) next
      v <- Pt[rows, j]; r <- which.max(v)
      if (r == 1 || r == length(rows)) next
      strength <- if (sparse) v[r] / max(1e-12, stats::quantile(Pt[, j][Pt[, j] > 0], 0.5, na.rm = TRUE)) else v[r] / colmed[j]
      if (!is.finite(strength) || strength < (if (sparse) min_strength / 10 else min_strength)) next
      O[i, j] <- lp[rows[r]] - lP[i]; Sg[i, j] <- strength
    }
    # every cycle on its own (no comb): the ridge nearest its expected period within the band, for the per-cycle view
    OI <- matrix(NA_real_, length(P), ncol(Pt))
    for (i in seq_along(P)) {
      rows <- which(abs(lp - lP[i]) <= max(bw[i], 0.2))
      if (length(rows) < 3) next
      sub <- Pt[rows, , drop = FALSE]; jm <- apply(sub, 2, which.max); mx <- sub[cbind(jm, seq_len(ncol(sub)))]
      # only a real ridge counts here: a tenth of the column maximum on a squeezed map, min_strength over the continuum otherwise
      strength <- if (sparse) mx / pmax(1e-12, apply(Pt, 2, max)) else mx / colmed
      ok <- is.finite(strength) & strength >= (if (sparse) 0.1 else min_strength) & jm > 1 & jm < length(rows)
      oi <- lp[rows[jm]] - lP[i]; oi[!ok] <- NA
      OI[i, ] <- .wr_runmed_na(oi, max(3, 2 * floor(0.75 * P[i] / stats::median(diff(ax))) + 1))
    }
    # smooth the shift along time over half the longest cycle (a rate error cannot be resolved finer than that)
    w <- max(3, 2 * floor(0.5 * max(P) / stats::median(diff(ax))) + 1)
    shs <- .wr_runmed_na(sh, w)
    nc <- colSums(is.finite(O)); nc <- round(.wr_runmed_na(nc, w)); nc[!is.finite(nc)] <- 0
    offs[[tf]] <- list(time = ax, shift = sh, offset = O, strength = Sg, offset_single = OI, score = S, shifts = shifts, contrast = contrast)
    corr[[tf]] <- list(time = ax, log2 = shs, n_cycles = nc)
  }
  # common depth grid: the depth axis of the first de-chirped map
  z <- dcs[[1]]$depth$axis.1; tz <- dcs[[1]]$time; rz <- dcs[[1]]$rate
  prof <- data.frame(depth = z, time = tz, sedrate = rz)
  NC <- matrix(0, length(z), length(transforms)); LC <- matrix(NA_real_, length(z), length(transforms))
  for (k in seq_along(transforms)) {
    cr <- corr[[transforms[k]]]
    LC[, k] <- stats::approx(cr$time, cr$log2, tz, rule = 2)$y
    NC[, k] <- round(stats::approx(cr$time, cr$n_cycles, tz, rule = 2)$y)
    prof[[transforms[k]]] <- LC[, k]
  }
  comb <- apply(LC, 1, function(v) if (any(is.finite(v))) stats::median(v, na.rm = TRUE) else NA_real_)
  nt <- rowSums(is.finite(LC) & abs(LC) > log2(1 + tolerance) & sign(LC) == sign(matrix(comb, length(z), length(transforms))))
  prof$correction <- 2^comb
  prof$n_cycles <- apply(NC, 1, max)
  prof$n_transforms <- nt
  prof$flag <- is.finite(comb) & abs(comb) > log2(1 + tolerance) & prof$n_cycles >= min_cycles & nt >= min_transforms
  # intervals: gaps shorter than min_length (kyr) between flags of the same sign are closed, short flags dropped
  if (is.null(min_length)) min_length <- 0.5 * max(P)
  sg <- sign(comb); sg[!prof$flag] <- 0
  r <- rle(sg); e <- cumsum(r$lengths); s0 <- e - r$lengths + 1
  for (i in seq_along(r$values)) if (r$values[i] == 0 && i > 1 && i < length(r$values) && r$values[i - 1] == r$values[i + 1] &&
                                     tz[e[i]] - tz[s0[i]] < min_length) sg[s0[i]:e[i]] <- r$values[i - 1]
  r <- rle(sg); e <- cumsum(r$lengths); s0 <- e - r$lengths + 1
  for (i in seq_along(r$values)) if (r$values[i] != 0 && tz[e[i]] - tz[s0[i]] < min_length) sg[s0[i]:e[i]] <- 0
  prof$flag <- sg != 0
  r <- rle(prof$flag); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  iv <- NULL
  for (i in which(r$values)) {
    s <- starts[i]:ends[i]
    cyc <- max(prof$n_cycles[s]); trf <- max(prof$n_transforms[s])
    iv <- rbind(iv, data.frame(from = z[s[1]], to = z[s[length(s)]], correction = round(exp(mean(log(prof$correction[s]))), 3),
                               cycles = cyc, transforms = trf))
  }
  if (is.null(iv)) iv <- data.frame(from = numeric(0), to = numeric(0), correction = numeric(0), cycles = integer(0), transforms = integer(0))
  rt <- NULL
  if (retrack) {
    if (is.null(sedrate)) sedrate <- range(track$sedrate) * c(0.5, 2)
    if (verbose) message("check_agemodel: blind chirp-matched re-track")
    rt <- tryCatch(track_astro_auto(dat, sedrate = sedrate, age_ma = age_ma, astro_periods = astro_periods, method = "chirp", dj = dj),
                   error = function(e) { if (verbose) message("re-track failed: ", conditionMessage(e)); NULL })
    if (!is.null(rt)) {
      rtr <- if (inherits(rt, "wr_track_auto")) rt$track else rt
      prof$sedrate_retrack <- exp(stats::approx(rtr$depth, log(rtr$sedrate), z, rule = 2)$y)
    }
  }
  out <- list(profile = prof, intervals = iv, offsets = offs, dechirped = dcs, retrack = rt, astro_periods = P,
              band = bw, tolerance = tolerance, track = track, data = dat)
  class(out) <- "wr_agemodel_check"
  out
}

.wr_runmed_na <- function(v, k) {
  k <- max(3, 2 * floor(k / 2) + 1); n <- length(v); out <- rep(NA_real_, n); h <- (k - 1) / 2
  for (j in seq_len(n)) {
    s <- v[max(1, j - h):min(n, j + h)]; s <- s[is.finite(s)]
    if (length(s) >= max(2, k / 3)) out[j] <- stats::median(s)
  }
  out
}
.wr_wmedian <- function(v, w) {
  o <- order(v); v <- v[o]; w <- w[o]; cw <- cumsum(w) / sum(w)
  v[which(cw >= 0.5)[1]]
}

#' @export
print.wr_agemodel_check <- function(x, ...) {
  cat("Age-model check:", nrow(x$profile), "positions,", length(x$dechirped), "transforms,",
      length(x$astro_periods), "cycles (", paste(round(x$astro_periods, 1), collapse = ", "), "kyr )\n")
  fl <- mean(x$profile$flag, na.rm = TRUE)
  cat(sprintf("Flagged: %.0f %% of the record in %d interval(s) (tolerance %.0f %%)\n", 100 * fl, nrow(x$intervals), 100 * x$tolerance))
  if (nrow(x$intervals)) print(x$intervals, row.names = FALSE)
  if (!is.null(x$profile$sedrate_retrack)) {
    rr <- x$profile$sedrate_retrack / x$profile$sedrate
    cat(sprintf("Blind re-track / existing track: median %.2f, range %.2f-%.2f\n", stats::median(rr), min(rr), max(rr)))
  }
  invisible(x)
}

#' Plot an age-model check: the de-chirped maps, the rate corrections the cycles ask for and the flagged intervals
#'
#' @param x Result of \code{\link{check_agemodel}}.
#' @param transforms Which de-chirped maps to show (default all that were computed).
#' @param view \code{"tuned"} (maps in the time of the track, cycles should be horizontal) or \code{"depth"}
#' (mapped back to depth) \code{Default="tuned"}.
#' @param rate_unit Factor and label of the rate axis, e.g. \code{c(100, "cm/kyr")}.
#' @param main Title.
#' @param ... Passed to \code{\link{plot_dechirped}}.
#' @return Invisibly \code{x}.
#' @export
plot_agemodel_check <- function(x, transforms = names(x$dechirped), view = c("tuned", "depth"),
                                rate_unit = c(1, "depth units / kyr"), main = NULL, ...) {
  view <- match.arg(view)
  pr <- x$profile; f <- as.numeric(rate_unit[1])
  n <- length(transforms)
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::layout(matrix(seq_len(n + 3), ncol = 1), heights = c(1, rep(1.6, n), 1.6, 1.3))
  graphics::par(mar = c(2.5, 4.2, 2, 1), mgp = c(2.2, 0.6, 0))
  xl <- if (view == "tuned") "Time of the track (kyr)" else "Depth"
  xx <- if (view == "tuned") pr$time else pr$depth
  graphics::plot(xx, stats::approx(x$data$depth, x$data$x, pr$depth)$y, type = "l", xlab = "", ylab = "Proxy", xaxs = "i",
                 main = if (is.null(main)) "Age-model check: cycles should be horizontal in the time of the track" else main)
  shade <- function() { if (!any(pr$flag)) return(invisible()); fl <- rle(pr$flag); e <- cumsum(fl$lengths); s <- e - fl$lengths + 1
    u <- graphics::par("usr"); for (i in which(fl$values)) graphics::rect(xx[s[i]], u[3], xx[e[i]], u[4], col = grDevices::adjustcolor("red", 0.12), border = NA) }
  shade()
  for (tf in transforms) {
    dc <- x$dechirped[[tf]]
    if (view == "tuned") plot_dechirped(dc, "tuned", show = x$astro_periods, main = paste(tf, "in the time of the track (periods in kyr)"), ...)
    else plot_dechirped(dc, "plain", show = x$astro_periods, main = paste(tf, "mapped back to depth"), xlab = "Depth", ...)
    if (any(pr$flag)) {
      fl <- rle(pr$flag); e <- cumsum(fl$lengths); s <- e - fl$lengths + 1
      for (i in which(fl$values)) graphics::segments(xx[s[i]], graphics::par("usr")[4], xx[e[i]], graphics::par("usr")[4], col = "red", lwd = 4, xpd = NA)
    }
  }
  # corrections
  L <- as.matrix(pr[, transforms, drop = FALSE])
  yl <- range(c(0.5, 2, 2^L, pr$correction), na.rm = TRUE)
  graphics::plot(xx, pr$correction, type = "n", log = "y", ylim = yl, xlab = xl, xaxs = "i",
                 ylab = "Rate correction factor", main = "Sedimentation-rate correction the cycles ask for (1 = track is right); flagged intervals shaded")
  shade()
  graphics::abline(h = 1, lwd = 1); graphics::abline(h = c(1 + x$tolerance, 1 / (1 + x$tolerance)), lty = 3)
  cols <- grDevices::hcl.colors(max(3, n), "Dark 3")
  for (k in seq_len(n)) graphics::lines(xx, 2^L[, k], col = cols[k], lwd = 0.8)
  graphics::lines(xx, pr$correction, lwd = 2.2)
  leg <- transforms; lc <- cols[seq_len(n)]; lw <- rep(0.8, n)
  if (!is.null(pr$sedrate_retrack)) {
    graphics::lines(xx, pr$sedrate_retrack / pr$sedrate, col = "blue", lwd = 1.5, lty = 2)
    leg <- c(leg, "combined", "blind re-track / track"); lc <- c(lc, "black", "blue"); lw <- c(lw, 2.2, 1.5)
  } else { leg <- c(leg, "combined"); lc <- c(lc, "black"); lw <- c(lw, 2.2) }
  graphics::legend("topleft", legend = leg, col = lc, lwd = lw, lty = c(rep(1, n + 1), if (!is.null(pr$sedrate_retrack)) 2), bty = "n", cex = 0.8, ncol = 2)
  # every cycle on its own, from the sharpest transform available
  tfc <- intersect(c("superlet_sst", "chirplet_sst", "wavelet_sst", "eha", "wavelet", "combined"), transforms)[1]
  o <- x$offsets[[tfc]]; OI <- o$offset_single; xo <- if (view == "tuned") o$time else stats::approx(pr$time, pr$depth, o$time, rule = 2)$y
  graphics::plot(range(xx), c(0.5, 2), type = "n", log = "y", xlab = xl, xaxs = "i", ylab = "Observed / expected period",
                 main = paste0("Every cycle on its own (", tfc, "): where its ridge sits relative to its expected period"))
  graphics::abline(h = 1); graphics::abline(h = c(1 + x$tolerance, 1 / (1 + x$tolerance)), lty = 3)
  cc <- grDevices::hcl.colors(nrow(OI) + 1, "Zissou 1")
  for (i in seq_len(nrow(OI))) graphics::lines(xo, 2^OI[i, ], col = cc[i], lwd = 1.2)
  graphics::legend("topleft", legend = paste(round(x$astro_periods, 1), "kyr"), col = cc[seq_len(nrow(OI))], lwd = 1.2, bty = "n", cex = 0.8, ncol = min(6, nrow(OI)))
  invisible(x)
}
