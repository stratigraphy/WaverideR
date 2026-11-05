#' Iterative (demodulated) astronomical tracking of the sedimentation rate
#'
#' Tracks the sedimentation rate, converts the record to approximate time with the result, transforms it again in
#' that (nearly stationary) time domain with sharper wavelets, tracks the residual stretch of the time axis, composes
#' the age models and repeats. Demodulating a signal by its own estimated instantaneous frequency before the
#' synchrosqueezing makes the ridges much sharper than any higher-order synchrosqueezing of the raw signal (Meignen et
#' al. 2017; Jiang et al. 2022), and in cyclostratigraphy the demodulation is the age model itself: in the time
#' domain the astronomical cycles have (nearly) fixed periods, so wavelets with many more cycles can be used, the
#' eccentricity quartet and the precession bands separate, and the tracked residual converges to the true age model.
#'
#' Every iteration is checked with the posterior of the tracking (\code{confidence10}, the probability within 10
#' percent of the track): when the confidence falls, the previous iteration is kept and the loop stops. The
#' demodulation sharpens the ridges but cannot repair a first track that follows the wrong cycles over a long
#' interval (the warped record then shows the same ambiguity at a different stretch), so the prior information of
#' \code{\link{track_astro_auto}} (\code{duration}, \code{constraints}, \code{jump_depths}) matters more than
#' the iterations. The credible band of an iteration is the band of the residual widened in quadrature by the band of
#' the previous iteration.
#'
#' @param data Record to analyse (required, no default): matrix or data frame with depth in column 1 and the proxy
#' value in column 2; further columns are ignored, rows with missing values are removed and the rows are sorted by depth.
#' @param sedrate Range of the sedimentation rate in which the first (depth) track is searched (required, no default):
#' numeric vector \code{c(lower, upper)}, > 0, in depth units per kyr (e.g. m/kyr for a record in metres).
#' @param n_iter Maximum number of warped (time-domain) iterations (integer, >= 0; 0 returns the depth track only).
#' In synthetic tests (three jumps of the sedimentation rate, AR1 noise) the first warped iteration lowered the largest
#' age error by about 10-15 percent and a second one gave nothing more (it sometimes regressed), so one iteration is the
#' default. \code{Default=1}.
#' @param residual Largest residual stretch of the time axis searched in the warped iterations (numeric, > 1): the
#' residual rate is searched in \code{c(1/residual, residual)} in the first warped iteration and in
#' \code{c(1/sqrt(residual), sqrt(residual))} in later ones. \code{Default=1.6}.
#' @param scalogram Scalogram of the first (depth) iteration (character): \code{"chirplet_sst"} (in the synthetic tests
#' its tracks had about 35 percent lower age errors than those of \code{"superlet_sst"}, at four times the cost),
#' \code{"superlet_sst"}, \code{"msst"} or \code{"cwt"} (best for records that resolve only the eccentricity, e.g. Klonk);
#' see \code{\link{track_astro_auto}}. \code{Default=c("chirplet_sst", "superlet_sst", "msst", "cwt")}, i.e.
#' \code{"chirplet_sst"}.
#' @param omega_nr Number of oscillations of the Morlet wavelet for \code{scalogram = "cwt"} (numeric, > 0).
#' \code{Default=6}.
#' @param omega_nr_warped Morlet omega_0 values of the superlet-SST orders (\code{omega_nr} of
#' \code{\link{analyse_superlet_sst}}) used in the warped (time) iterations (numeric vector, > 0); higher values give
#' sharper period resolution, which the nearly stationary time domain allows. \code{Default=c(12, 18, 24)}.
#' @param age_ma Age of the record in Ma (numeric, >= 0), used to compute the astronomical periods with
#' \code{\link{astro_periods_tracking}}; ignored when \code{astro_periods} is given. \code{Default=0}.
#' @param k_yr Precession constant (numeric, the period of the axial precession in years) passed to
#' \code{\link{astro_periods_tracking}} instead of the AstroGeo22 value for \code{age_ma}; \code{NULL} uses AstroGeo22.
#' \code{Default=NULL}.
#' @param astro_periods Astronomical periods in kyr replacing the ones from \code{age_ma}: a numeric vector of periods or
#' a data frame in the format of \code{\link{astro_periods_tracking}}; \code{NULL} computes them. \code{Default=NULL}.
#' @param guide Rough sedimentation rate curve for the first track (two-column data frame or matrix: depth, rate in
#' depth units per kyr, or \code{NULL}); see \code{\link{track_astro_auto}}. \code{Default=NULL}.
#' @param guide_tolerance Relative half width of the band around \code{guide} (numeric, > 0). \code{Default=0.4}.
#' @param duration Approximate duration of the record in kyr (numeric, or \code{NULL}); turned into a guide when
#' \code{guide = NULL}, see \code{\link{track_astro_auto}}. Strongly recommended for real records. \code{Default=NULL}.
#' @param duration_tolerance Relative tolerance of the guide built from \code{duration} (numeric, > 0).
#' \code{Default=0.5}.
#' @param constraints Duration constraints (data frame with the columns \code{from}, \code{to} (depths), \code{duration}
#' (kyr) and optionally \code{tolerance}, or \code{NULL}); see \code{\link{track_astro_auto}} and
#' \code{\link{track_astro_hmm}}. Given in depth and converted to model time for the warped iterations.
#' \code{Default=NULL}.
#' @param jump_depths Depths where an abrupt change of the sedimentation rate is plausible (numeric vector in depth
#' units, or \code{NULL}); see \code{\link{track_astro_hmm}}. Converted to model time for the warped iterations.
#' \code{Default=NULL}.
#' @param dj Spacing of the period axis of the scalograms, in octaves (numeric, > 0; smaller values give a finer grid
#' and a slower run). \code{Default=1/50}.
#' @param n_samples Number of paths drawn from the posterior for the credible bands (integer, > 0). \code{Default=100}.
#' @param verbose Print the confidence and duration of every iteration (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param ... Further arguments passed to \code{\link{track_astro_auto}} (first iteration) and to
#' \code{\link{track_astro_hmm}} (warped iterations), e.g. \code{sd_change}, \code{seed}.
#'
#' @return A \code{wr_track} data frame in depth as from \code{\link{track_astro_hmm}} (columns depth, sedrate,
#' confidence, confidence10, time, time_lower, time_upper, ...), with attributes \code{iterations} (the per-iteration
#' summaries: mean confidence10, duration), \code{tracks} (the track of every iteration in depth), \code{warped}
#' (the last warped scalogram with its period axis in kyr, i.e. the tuned scalogram, and the record in time) and
#' \code{first} (the result of the depth iteration, a \code{\link{track_astro_auto}} object).
#'
#' @references Meignen, S., Pham, D.-H. and McLaughlin, S. (2017). On demodulation, ridge detection, and
#' synchrosqueezing for multicomponent signals. IEEE Trans. Signal Process. 65, 2093-2103.
#'
#' Jiang, Q., Prater-Bennette, A., Suter, B.W. and Zeyani, A. (2022). Instantaneous frequency-embedded
#' synchrosqueezing transform for signal separation. Front. Appl. Math. Stat. 8, 830530.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1500, by = 0.5)
#' ecc <- 0.4 * cos(2 * pi * t / 405) + 0.25 * cos(2 * pi * t / 125)
#' x_t <- (1 + ecc) * (sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19)) +
#'   0.7 * sin(2 * pi * t / 41)
#' sr <- ifelse(t < 700, 1.5, 2.6)
#' zt <- c(0, cumsum(sr[-1] * diff(t)))
#' z <- seq(0, max(zt), by = 1)
#' x <- approx(zt, x_t, z)$y + rnorm(length(z), sd = 1.5)
#' res <- track_astro_iterative(data.frame(z, x), sedrate = c(0.8, 5))
#' plot(res$depth, res$sedrate, type = "l")
#' lines(z, approx(zt, sr, z)$y, col = 2)
#' }
#' @export
track_astro_iterative <- function(data,
    sedrate,
    n_iter = 1,
    residual = 1.6,
    scalogram = c("chirplet_sst", "superlet_sst", "msst", "cwt"),
    omega_nr = 6,
    omega_nr_warped = c(12, 18, 24),
    age_ma = 0,
    k_yr = NULL,
    astro_periods = NULL,
    guide = NULL,
    guide_tolerance = 0.4,
    duration = NULL,
    duration_tolerance = 0.5,
    constraints = NULL,
    jump_depths = NULL,
    dj = 1 / 50,
    n_samples = 100,
    verbose = TRUE,
    ...) {
  scalogram <- match.arg(scalogram)
  dat <- as.data.frame(data)[, 1:2]
  dat <- dat[stats::complete.cases(dat), ]
  dat <- dat[order(dat[, 1]), ]
  ap <- if (is.null(astro_periods)) { astro_periods_tracking(age_ma, k_yr) } else { astro_periods }
  pv <- if (is.data.frame(ap)) { ap$period } else { ap }
  dots <- list(...)
  first <- track_astro_auto(dat, sedrate = sedrate, astro_periods = ap, scalogram = scalogram, omega_nr = omega_nr, guide = guide,
      guide_tolerance = guide_tolerance, duration = duration,
      duration_tolerance = duration_tolerance, constraints = constraints,
      jump_depths = jump_depths, dj = dj, n_samples = n_samples, method = "hmm", ...)
  tr <- first$track
  tracks <- list(tr)
  its <- data.frame(iteration = 0, confidence10 = mean(tr$confidence10), duration = diff(range(tr$time)))
  if (verbose) {
    cat("iteration 0 (depth): confidence10", round(its$confidence10, 3), " duration", round(its$duration), "\n")
  }
  z <- tr$depth
  warped <- NULL
  for (it in seq_len(n_iter)) {
    f <- if (it == 1) { residual } else { sqrt(residual) }
    # the record on a regular grid of model time
    tz <- tr$time
    if (any(diff(tz) <= 0)) {
      tz <- cummax(tz + seq_along(tz) * 1e-9)
    }
    xz <- stats::approx(dat[, 1], dat[, 2], z, rule = 2)$y
    dtau <- diff(range(tz)) / (length(z) - 1)
    tau <- seq(min(tz), max(tz), by = dtau)
    xt <- stats::approx(tz, xz, tau)$y
    lo <- max(3 * dtau, min(pv) / f * 0.95)
    hi <- min(diff(range(tau)) / 2, max(pv) * f * 1.05)
    sc <- analyse_superlet_sst(data.frame(tau, xt), lowerPeriod = lo, upperPeriod = hi, dj = dj, omega_nr = omega_nr_warped)
    # prior information carried into time
    to_tau <- function(d) stats::approx(z, tz, d, rule = 2)$y
    cons_t <- NULL
    if (!is.null(constraints)) {
      cons_t <- as.data.frame(constraints)
      cons_t[, 1] <- to_tau(cons_t[, 1])
      cons_t[, 2] <- to_tau(cons_t[, 2])
    }
    jd_t <- if (!is.null(jump_depths)) { to_tau(jump_depths) } else { NULL }
    args <- c(list(sc, astro_periods = ap, sedrate = c(1 / f, f), n_samples = n_samples, constraints = cons_t,
        jump_depths = jd_t), dots)
    rt <- do.call(track_astro_hmm, args)
    # compose: the residual stretch s(tau) = d(model time) / d(true time); true time = rt$time
    tnew <- stats::approx(rt$depth, rt$time, tz, rule = 2)$y
    tlo <- stats::approx(rt$depth, rt$time_lower, tz, rule = 2)$y
    thi <- stats::approx(rt$depth, rt$time_upper, tz, rule = 2)$y
    tmed <- stats::approx(rt$depth, rt$time_median, tz, rule = 2)$y
    sres <- stats::approx(rt$depth, rt$sedrate, tz, rule = 2)$y
    new <- tr
    new$sedrate <- tr$sedrate * sres
    new$sedrate_mean <- tr$sedrate_mean * stats::approx(rt$depth, rt$sedrate_mean, tz, rule = 2)$y
    new$sedrate_lower <- tr$sedrate_lower * stats::approx(rt$depth, rt$sedrate_lower, tz, rule = 2)$y
    new$sedrate_upper <- tr$sedrate_upper * stats::approx(rt$depth, rt$sedrate_upper, tz, rule = 2)$y
    new$confidence <- stats::approx(rt$depth, rt$confidence, tz, rule = 2)$y
    new$confidence10 <- stats::approx(rt$depth, rt$confidence10, tz, rule = 2)$y
    # credible band: the band of the residual (from the warped posterior) widened in quadrature by the band of the
    # previous iteration, since the warped posterior is conditional on the previous age model
    hw_prev <- (tr$time_upper - tr$time_lower) / 2
    hw_lo <- sqrt((tnew - tlo)^2 + hw_prev^2)
    hw_hi <- sqrt((thi - tnew)^2 + hw_prev^2)
    new$time <- tnew
    new$time_lower <- tnew - hw_lo
    new$time_upper <- tnew + hw_hi
    new$time_median <- tmed
    for (k in seq_along(pv)) {
      new[[paste0("period_", pv[k])]] <- new$sedrate * pv[k]
    }
    attr(new, "constraints") <- attr(rt, "constraints")
    c10 <- mean(new$confidence10)
    if (verbose) {
      cat("iteration", it, "(warped, residual", round(1 / f, 2), "-", round(f, 2), "): confidence10",
          round(c10, 3), " duration", round(diff(range(new$time))), "\n")
    }
    its <- rbind(its, data.frame(iteration = it, confidence10 = c10, duration = diff(range(new$time))))
    tracks[[it + 1]] <- new
    if (c10 < its$confidence10[it] - 0.02 && it > 1) {
      if (verbose) {
        cat("confidence fell: keeping iteration", it - 1, "\n")
      }
      break
    }
    tr <- new
    warped <- list(scalogram = sc, track = rt, time = tz, record = data.frame(time = tau, x = xt))
  }
  attr(tr, "iterations") <- its
  attr(tr, "tracks") <- tracks
  attr(tr, "warped") <- warped
  attr(tr, "first") <- first
  attr(tr, "astro_periods") <- pv
  class(tr) <- c("wr_track", "data.frame")
  tr
}
