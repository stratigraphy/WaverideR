#' @title Windowed timeOpt sedimentation rate estimation
#'
#' @description The \code{\link{win_timeOpt}} function conducts a windowed
#' timeOpt sedimentation rate estimation.
#' This function is based on \code{\link[astrochron:timeOpt]{timeOpt}} and
#' \code{\link[astrochron:eTimeOpt]{eTimeOpt}} of the 'astrochron' package, but allows for
#' multithreaded analysis, windows of fixed duration, a correction of the r-squared values for
#' their bias towards high sedimentation rates, a mask of rates that cannot be tested and an
#' obliquity version of the fit.
#'
#' @details
#' \strong{Bias towards high sedimentation rates.} In a window of fixed thickness, a higher
#' trial sedimentation rate means a shorter time span, so fewer cycles of the target periods fit
#' in the window. The envelope fit and the spectral power fit then have fewer independent
#' data points per fitted sinusoid, and for red noise the target frequencies move to the strong
#' low-frequency end of the spectrum. Both raise the r-squared values that pure noise reaches,
#' so the raw r-squared surface has a rising baseline, and a genuine peak at a low sedimentation
#' rate can be smaller than the noise level at high sedimentation rates.
#'
#' \strong{Windows of fixed duration} (\code{window_mode = "time"}). For every trial
#' sedimentation rate the window is cut to \code{window_duration} kyr (its thickness is
#' \code{window_duration * sedrate / 100} m), so every rate is tested on the same number of
#' cycles. This removes most of the bias at its source. Rates whose window does not fit in the
#' record around a given depth are marked as untestable there. timeOpt is run once per rate.
#' With \code{window_mode = "depth"} (default) the window has a fixed thickness
#' \code{window_size} (the behaviour of earlier versions).
#'
#' \strong{Noise-level correction} (\code{bias_correction = "null"}). For each window a noise
#' model is fitted to the periodogram (Whittle likelihood; for time windows on the window of the
#' middle rate that fits in the record): an AR(1) spectrum (flat at low
#' frequencies, bending to f^-2) and a power law (f^-beta, e.g. pink noise for beta = 1).
#' With \code{noise_model = "auto"} the better-fitting one is used. Noise series of the same model,
#' length and sample spacing are analysed with exactly the same timeOpt settings. This gives,
#' for every trial sedimentation rate, the mean and standard deviation of the noise r-squared
#' values (on a logit scale, because r-squared is bounded by 0 and 1). These are smoothed over
#' neighbouring rates (\code{null_smooth}) and the observed values are expressed as z-scores:
#' \deqn{z_{env} = (logit(r^2_{env}) - \mu_{env}) / \sigma_{env}}
#' and likewise for the spectral power fit; both are combined (Stouffer's method) as
#' \eqn{z_{opt} = (z_{env} + z_{pwr}) / \sqrt{2}}. A z-score of 0 is the noise level at
#' that sedimentation rate. To keep the computation affordable, the noise level is computed on a
#' grid of model parameters (steps \code{null_rho_step} and \code{null_beta_step}; for depth windows
#' also three window lengths) and interpolated. The noise tables can be reused with
#' \code{null_cache}.
#'
#' \strong{Probabilities.} For each grid node the maximum z_opt over all testable rates is
#' recorded for every noise series. This gives a family-wise false-alarm threshold per window
#' (\code{z_threshold}, the 1 - \code{alpha} quantile) and a family-wise probability per cell
#' (\code{p_fw}): the probability that noise reaches this z-score anywhere on the rate grid.
#'
#' \strong{Untestable rates} (\code{untestable_code}): 1 = fewer than \code{min_samples} samples
#' per shortest carrier cycle (precession for fit 1, short eccentricity for fit 2, obliquity for
#' fit 3); 2 = the window spans fewer than \code{min_cycles} of the longest fitted period;
#' 3 = the window does not fit in the record. Untestable cells get no z-score.
#'
#' \strong{Obliquity version} (\code{fit = 3}): timeOpt fit 1 with obliquity carrier periods
#' (\code{targetO}) bandpassed between \code{flow} and \code{fhigh} (default 1/46-1/36 cycles/kyr)
#' and the obliquity modulation periods (\code{targetOmod}) in the envelope fit. A precession
#' record read as obliquity gives a peak at half the true rate, an obliquity record read as
#' precession (fit 1) a peak at twice the true rate; compare the fits with
#' \code{\link{win_timeOpt_multi}}.
#'
#' @param data Record to analyse (matrix or data frame with two columns: depth in metres in column 1, proxy value in
#' column 2), evenly sampled (the depth steps must be equal within 0.1 percent; use e.g. \code{linterp} of
#' 'astrochron' to resample). Must be given (the \code{NULL} default fails). \code{Default=NULL}.
#' @param window_size Thickness of the moving window in metres (numeric, > 0), only used with
#' \code{window_mode = "depth"}. \code{Default=10}.
#' @param sedmin Lowest sedimentation rate of the grid (numeric, > 0, in cm/kyr). \code{Default=0.5}.
#' @param sedmax Highest sedimentation rate of the grid (numeric, > \code{sedmin}, in cm/kyr). \code{Default=2}.
#' @param numsed Number of sedimentation rates of the grid between \code{sedmin} and \code{sedmax} (integer, >= 1;
#' 1 tests \code{sedmin} only). \code{Default=100}.
#' @param limit Limit the evaluated sedimentation rates to the range in which the full target signal can be
#' recovered (logical, TRUE or FALSE; passed to \code{\link[astrochron]{timeOpt}}). Only used with
#' \code{window_mode = "depth"}; the untestable mask (\code{min_samples}, \code{min_cycles}) is usually the better
#' choice. \code{Default=FALSE}.
#' @param fit Which amplitude modulation is tested (numeric, one of 1, 2, 3): 1 precession modulated by eccentricity,
#' 2 short eccentricity modulated by the 405-kyr cycle, 3 obliquity modulated by the obliquity modulation periods
#' (\code{targetO}, \code{targetOmod}; see Details). \code{Default=2}.
#' @param fitModPwr Include the modulation periods in the spectral power fit (logical, TRUE or FALSE; passed to
#' \code{\link[astrochron]{timeOpt}}; not used for \code{fit = 3}). \code{Default=TRUE}.
#' @param flow Low-frequency cut-off of the Taner bandpass filter that isolates the carrier (numeric, half power point
#' in cycles/kyr). \code{NULL} uses the default of \code{\link[astrochron]{timeOpt}} for \code{fit = 1} and
#' \code{fit = 2}, and 1/46 for \code{fit = 3}. \code{Default=NULL}.
#' @param fhigh High-frequency cut-off of the Taner bandpass filter (numeric, half power point in cycles/kyr).
#' \code{NULL} uses the default of \code{\link[astrochron]{timeOpt}} for \code{fit = 1} and \code{fit = 2}, and 1/36
#' for \code{fit = 3}. \code{Default=NULL}.
#' @param roll Roll-off rate of the Taner filter (numeric, > 0, in dB/octave); larger values give steeper filter
#' edges. \code{Default=10^6}.
#' @param targetE Eccentricity periods to evaluate (numeric vector, in kyr), in order of decreasing period with 405
#' kyr (or 405.7) first. \code{Default=c(405.7, 130.7, 123.8, 98.9, 94.9)}.
#' @param targetP Precession periods to evaluate (numeric vector, in kyr), in order of decreasing period; used for
#' \code{fit = 1} and in the spectral fit. \code{Default=c(20.9, 19.9, 17.1, 17.2)}.
#' @param targetO Obliquity carrier periods (numeric vector, in kyr), only used for \code{fit = 3}.
#' \code{Default=c(53.6, 41.0, 39.6)}.
#' @param targetOmod Obliquity modulation periods (numeric vector, in kyr) used in the envelope fit of \code{fit = 3}.
#' \code{Default=173}.
#' @param detrend Remove a linear trend from every window before the analysis (logical, TRUE or FALSE); also used for
#' the noise model fit. \code{Default=TRUE}.
#' @param normalize Divide the raw r-squared curve of every window by its maximum, so that it peaks at 1 (logical,
#' TRUE or FALSE). The bias-corrected z-scores are not normalised. \code{Default=TRUE}.
#' @param linLog Spacing of the sedimentation rate grid (numeric, 0 or 1): 0 linear, 1 logarithmic.
#' \code{Default=1}.
#' @param window_mode How the window is defined (character, one of \code{"depth"}, \code{"time"}): \code{"depth"} a
#' window of fixed thickness \code{window_size} (the behaviour of earlier versions); \code{"time"} a window of fixed
#' duration \code{window_duration} for every trial sedimentation rate, which removes most of the bias towards high
#' rates (see Details). \code{Default="depth"}.
#' @param window_duration Duration of the window in kyr (numeric, > 0), only used with \code{window_mode = "time"}; its
#' thickness at a trial rate is \code{window_duration * sedrate / 100} m. \code{Default=1500}.
#' @param step Spacing of the window centres in samples: every \code{step}-th data point is a window centre (integer,
#' >= 1); larger values are faster. Ignored when \code{centres} is given. \code{Default=1}.
#' @param centres Depths of the window centres to evaluate (numeric vector, in metres; each is moved to the nearest
#' sample); overrides \code{step}. \code{NULL} uses \code{step}. \code{Default=NULL}.
#' @param bias_correction Correction of the r-squared values for their bias towards high sedimentation rates
#' (character, one of \code{"null"}, \code{"none"}): \code{"null"} standardises them against noise series of a fitted
#' noise model and gives z-scores and family-wise probabilities (see Details); \code{"none"} only returns the raw
#' r-squared values (much faster). \code{Default="null"}.
#' @param noise_model Noise model of the correction (character, one of \code{"auto"}, \code{"ar1"},
#' \code{"powerlaw"}): \code{"ar1"} AR(1) red noise, \code{"powerlaw"} \eqn{f^{-\beta}} noise, \code{"auto"} the one of
#' the two that fits the periodogram of the window best (Whittle likelihood). Only used with
#' \code{bias_correction = "null"}. \code{Default="auto"}.
#' @param n_null Number of noise series per grid node of the noise model (integer, >= 1); more series give more stable
#' z-scores and probabilities and take longer (the family-wise probabilities need at least 5). \code{Default=50}.
#' @param null_rho_step Step of the grid of AR(1) coefficients on which the noise level is computed and interpolated
#' (numeric, between 0 and 1); smaller steps are more precise and slower. \code{Default=0.05}.
#' @param null_beta_step Step of the grid of power-law exponents on which the noise level is computed and interpolated
#' (numeric, > 0). \code{Default=0.1}.
#' @param null_smooth Half width (in grid steps of the sedimentation rate, integer >= 0) of the moving average of the
#' noise mean and standard deviation over neighbouring rates; 0 means no smoothing. \code{Default=2}.
#' @param alpha False-alarm level of the family-wise threshold \code{z_threshold} and of the column
#' \code{significant} (numeric, between 0 and 1). \code{Default=0.05}.
#' @param min_samples Minimum number of samples per shortest carrier cycle (precession for fit 1, short eccentricity
#' for fit 2, obliquity for fit 3) for a rate to be testable (numeric, > 0; in tests timeOpt still worked with 3.4
#' samples per shortest precession cycle). \code{Default=3}.
#' @param min_cycles Minimum number of cycles of the longest fitted period that the window must span for a rate to be
#' testable (numeric, > 0). \code{Default=2}.
#' @param null_cache Noise tables to reuse: the \code{null_cache} element of an earlier result of
#' \code{win_timeOpt} with the same settings (list); missing grid nodes are added. \code{NULL} computes all tables.
#' \code{Default=NULL}.
#' @param run_multicore Run the windows and the noise simulations in parallel on several cores (logical, TRUE or
#'   FALSE). \code{Default=FALSE}.
#' @param n_cores Number of cores used when \code{run_multicore = TRUE} (integer; \code{NULL}: an automatic, safe
#'   number: physical cores minus one, at most 4, and limited by the free memory). Set
#'   \code{options(waverider.n_cores = ...)} to change the default. A value > 1 also switches on parallel computing.
#'   Results do not depend on the number of cores. \code{Default=NULL}.
#' @param verbose Show progress bars of the window and noise computations (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param normalise British spelling of \code{normalize}; if given (not \code{NULL}) it overrides \code{normalize}.
#' \code{Default=NULL}.
#'
#' @author
#' Based on the \code{\link[astrochron:timeOpt]{timeOpt}} and
#' \code{\link[astrochron:eTimeOpt]{eTimeOpt}}
#' functions of the 'astrochron' R package.
#'
#' @references
#' Meyers, S.R. (2019). Cyclostratigraphy and the problem of astrochronologic testing.
#' Earth-Science Reviews 190, 190-223 <doi:10.1016/j.earscirev.2018.11.015>
#'
#' @examples
#' \donttest{
#' #Conduct a windowed timeOpt on the magnetic susceptibility record
#' #of the Sullivan core of Pas et al., (2018).
#' set.seed(1)
#' mag_win_timeOpt <- win_timeOpt(
#' data = mag[1:250, ],
#' window_size = 6,
#' sedmin = 0.1,
#' sedmax = 1,
#' numsed = 5,
#' limit = FALSE,
#' fit = 2,
#' fitModPwr = TRUE,
#' flow = NULL,
#' fhigh = NULL,
#' roll = 10 ^ 6,
#' targetE = c(405.7, 130.7, 123.8, 98.9, 94.9),
#' targetP = c(20.9, 19.9, 17.1, 17.2),
#' detrend = TRUE,
#' normalize = TRUE,
#' linLog = 1,
#' bias_correction = "null",
#' n_null = 6,
#' noise_model = "ar1",
#' step = 10,
#' run_multicore = FALSE,
#' verbose = FALSE)
#' }
#'
#' @return
#' Returns a list of class \code{win_timeOpt.result} with the elements
#' \describe{
#'   \item{r_2_envelope, r_2_power, r_2_opt}{matrices (sedimentation rate x window) of the
#'   raw r-squared values (normalised when \code{normalize = TRUE})}
#'   \item{r_2_envelope_avg, r_2_power_avg, r_2_opt_avg}{their averages over all windows}
#'   \item{sedrate}{the sedimentation rates (cm/ka)}
#'   \item{x, y}{depth and proxy value of the window centres}
#'   \item{linLog}{the linLog value}
#'   \item{z_envelope, z_power, z_opt}{z-scores against noise (\code{NULL} when
#'   \code{bias_correction = "none"}; \code{NA} for untestable cells)}
#'   \item{z_envelope_avg, z_power_avg, z_opt_avg}{their averages over all windows}
#'   \item{p_fw}{family-wise probability that noise reaches z_opt anywhere on the rate grid}
#'   \item{z_threshold}{family-wise 1 - alpha threshold of z_opt per window}
#'   \item{best_sedrate, best_z, p_best, significant}{per window: the rate with the highest
#'   z_opt among the testable rates (raw r2_opt when there is no correction), its z-score and
#'   family-wise probability, and whether \code{p_best < alpha}}
#'   \item{untestable_code, testable}{see Details}
#'   \item{noise_family, noise_param}{fitted noise model per cell (1 = AR(1) coefficient, 2 =
#'   power-law exponent)}
#'   \item{npts}{number of points of each window}
#'   \item{data}{the input data}
#'   \item{window_mode, window_size, window_duration, fit, bias_correction, noise_model, alpha}{
#'   the settings}
#'   \item{null_cache}{the noise tables, for reuse}
#' }
#'
#' @export
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom astrochron timeOpt
#' @importFrom stats rnorm lm residuals sd quantile median optimize fft filter

win_timeOpt <- function(data = NULL,
    window_size = 10,
    sedmin = 0.5,
    sedmax = 2,
    numsed = 100,
    limit = FALSE,
    fit = 2,
    fitModPwr = TRUE,
    flow = NULL,
    fhigh = NULL,
    roll = 10 ^ 6,
    targetE = c(405.7, 130.7, 123.8, 98.9, 94.9),
    targetP = c(20.9, 19.9, 17.1, 17.2),
    targetO = c(53.6, 41.0, 39.6),
    targetOmod = 173,
    detrend = TRUE,
    normalize = TRUE,
    linLog = 1,
    window_mode = c("depth", "time"),
    window_duration = 1500,
    step = 1,
    centres = NULL,
    bias_correction = c("null", "none"),
    noise_model = c("auto", "ar1", "powerlaw"),
    n_null = 50,
    null_rho_step = 0.05,
    null_beta_step = 0.1,
    null_smooth = 2,
    alpha = 0.05,
    min_samples = 3,
    min_cycles = 2,
    null_cache = NULL,
    run_multicore = FALSE,
    n_cores = NULL,
    verbose = FALSE,
    normalise = NULL) {
  if (!is.null(normalise)) {
    normalize <- normalise
  }  # British spelling synonym
  bias_correction <- match.arg(bias_correction)
  window_mode <- match.arg(window_mode)
  noise_model <- match.arg(noise_model)

  proxy <- unname(as.matrix(data))
  n <- nrow(proxy)
  dd <- diff(proxy[, 1])
  dx <- stats::median(dd)
  if (!all(abs(dd - dx) <= 1e-3 * abs(dx))) {
    stop("win_timeOpt: the depth values must be evenly spaced", call. = FALSE)
  }
  dx <- abs(dx)
  cfg <- wr_to_cfg(fit, fitModPwr, flow, fhigh, roll, targetE, targetP, targetO, targetOmod,
      detrend, limit)
  sedrate <- wr_sedrate_grid(sedmin, sedmax, numsed, linLog)
  numsed <- length(sedrate)
  if (is.null(centres)) {
    centres <- seq(1L, n, by = max(1L, as.integer(step)))
  } else {
    idx <- integer(length(centres))
    for (k in seq_along(centres)) {
      idx[k] <- wr_closest(proxy[, 1], centres[k], which = TRUE)[1]
    }
    centres <- idx
  }
  n_cent <- length(centres)

  # Windows: first row and number of points for every rate and centre
  start <- matrix(NA_integer_, numsed, n_cent)
  npts <- matrix(NA_integer_, numsed, n_cent)
  npts_rate <- rep(NA_real_, numsed)
  if (window_mode == "depth") {
    for (k in seq_len(n_cent)) {
      op <- centres[k]
      r1 <- wr_closest(proxy[, 1], proxy[op, 1] - (window_size / 2), which = TRUE)[1]
      r2 <- wr_closest(proxy[, 1], proxy[op, 1] + (window_size / 2), which = TRUE)[1]
      start[, k] <- min(r1, r2)
      npts[, k] <- abs(r2 - r1) + 1L
    }
  } else {
    for (j in seq_len(numsed)) {
      m <- round(window_duration * sedrate[j] / 100 / dx) + 1
      npts_rate[j] <- m
      for (k in seq_len(n_cent)) {
        s0 <- centres[k] - floor((m - 1) / 2)
        if (s0 >= 1 && s0 + m - 1 <= n) {
          start[j, k] <- s0
          npts[j, k] <- m
        }
      }
    }
  }
  code <- wr_to_testable(sedrate, dx, npts, cfg, min_samples, min_cycles)

  # Raw r2 and the noise model of every window
  win_res <- wr_lapply(seq_len(n_cent), function(k) {
    out <- matrix(NA_real_, numsed, 4)
    if (window_mode == "depth") {
      rows <- start[1, k]:(start[1, k] + npts[1, k] - 1L)
      d_sub <- proxy[rows, , drop = FALSE]
      out[, 1:2] <- wr_to_run(d_sub, sedrate, cfg, linLog)
      if (bias_correction == "null") {
        nf <- wr_noise_fit(d_sub[, 2], d_sub[, 1], detrend)
        out[, 3:4] <- rep(wr_noise_choose(nf, noise_model), each = numsed)
      }
    } else {
      for (j in seq_len(numsed)) {
        if (code[j, k] == 3L) {
          next
        }
        rows <- start[j, k]:(start[j, k] + npts[j, k] - 1L)
        out[j, 1:2] <- wr_to_run(proxy[rows, , drop = FALSE], sedrate[j], cfg, linLog)
      }
      # one noise model per centre, fitted on the window of the middle complete rate
      ok <- which(code[, k] != 3L)
      if (bias_correction == "null" && length(ok)) {
        j_ref <- ok[ceiling(length(ok) / 2)]
        rows <- start[j_ref, k]:(start[j_ref, k] + npts[j_ref, k] - 1L)
        nf <- wr_noise_fit(proxy[rows, 2], proxy[rows, 1], detrend)
        fam_par <- wr_noise_choose(nf, noise_model)
        out[ok, 3] <- fam_par[1]
        out[ok, 4] <- fam_par[2]
      }
    }
    out
  }, n_cores = n_cores, run_multicore = run_multicore, combine = "list", progress = verbose)

  r_2_envelope <- matrix(NA_real_, numsed, n_cent)
  r_2_power <- r_2_envelope
  noise_family <- r_2_envelope
  noise_param <- r_2_envelope
  for (k in seq_len(n_cent)) {
    r_2_envelope[, k] <- win_res[[k]][, 1]
    r_2_power[, k] <- win_res[[k]][, 2]
    noise_family[, k] <- win_res[[k]][, 3]
    noise_param[, k] <- win_res[[k]][, 4]
  }
  r_2_opt <- r_2_envelope * r_2_power

  # Noise-level correction
  z_envelope <- NULL
  z_power <- NULL
  z_opt <- NULL
  p_fw <- NULL
  z_threshold <- rep(NA_real_, n_cent)
  if (bias_correction == "null") {
    nodes <- wr_to_nodes(noise_family, noise_param, npts, code, window_mode, null_rho_step,
        null_beta_step)
    cfg_key <- wr_to_cfg_key(cfg, sedrate, dx, window_mode, window_duration)
    null_cache <- wr_to_null_build(nodes, null_cache, cfg_key, cfg, sedrate, dx, window_mode,
        npts_rate, linLog, n_null, n_cores, run_multicore, verbose)
    # statistics per node
    keys <- wr_node_key(nodes$family, nodes$param, nodes$npts)
    stats_list <- list()
    for (q in seq_along(keys)) {
      A <- null_cache$sims[[keys[q]]]
      A <- A[, , seq_len(min(n_null, dim(A)[3])), drop = FALSE]
      if (window_mode == "depth") {
        np_node <- matrix(nodes$npts[q], numsed, 1)
      } else {
        np_node <- matrix(npts_rate, numsed, 1)
      }
      t_rate <- wr_to_testable(sedrate, dx, np_node, cfg, min_samples, min_cycles)[, 1] == 0L
      stats_list[[keys[q]]] <- wr_to_null_stats(A, t_rate, null_smooth)
    }
    z_envelope <- matrix(NA_real_, numsed, n_cent)
    z_power <- z_envelope
    p_fw <- z_envelope
    zthr_cell <- z_envelope
    lg_env <- wr_logit(r_2_envelope)
    lg_pwr <- wr_logit(r_2_power)
    for (k in seq_len(n_cent)) {
      for (j in seq_len(numsed)) {
        if (code[j, k] != 0L || !is.finite(lg_env[j, k]) || !is.finite(lg_pwr[j, k])) {
          next
        }
        nb <- wr_to_node_weights(nodes, noise_family[j, k], noise_param[j, k], npts[j, k],
            window_mode)
        m_e <- 0
        s_e <- 0
        m_p <- 0
        s_p <- 0
        wsum <- 0
        for (q in seq_along(nb$key)) {
          st <- stats_list[[nb$key[q]]]
          if (is.null(st) || !is.finite(st$mu[j, 1]) || !is.finite(st$mu[j, 2])) {
            next
          }
          m_e <- m_e + nb$w[q] * st$mu[j, 1]
          s_e <- s_e + nb$w[q] * st$sd[j, 1]
          m_p <- m_p + nb$w[q] * st$mu[j, 2]
          s_p <- s_p + nb$w[q] * st$sd[j, 2]
          wsum <- wsum + nb$w[q]
        }
        if (wsum <= 0) {
          next
        }
        z_envelope[j, k] <- (lg_env[j, k] - m_e / wsum) / (s_e / wsum)
        z_power[j, k] <- (lg_pwr[j, k] - m_p / wsum) / (s_p / wsum)
        zo <- (z_envelope[j, k] + z_power[j, k]) / sqrt(2)
        pf <- 0
        th <- 0
        wsum2 <- 0
        for (q in seq_along(nb$key)) {
          st <- stats_list[[nb$key[q]]]
          if (is.null(st) || length(st$zmax) < 5) {
            next
          }
          pf <- pf + nb$w[q] * (1 + sum(st$zmax >= zo)) / (length(st$zmax) + 1)
          th <- th + nb$w[q] * stats::quantile(st$zmax, 1 - alpha, names = FALSE)
          wsum2 <- wsum2 + nb$w[q]
        }
        if (wsum2 > 0) {
          p_fw[j, k] <- min(1, pf / wsum2)
          zthr_cell[j, k] <- th / wsum2
        }
      }
    }
    z_opt <- (z_envelope + z_power) / sqrt(2)
  }

  # Best rate per window
  best_sedrate <- rep(NA_real_, n_cent)
  best_z <- rep(NA_real_, n_cent)
  p_best <- rep(NA_real_, n_cent)
  for (k in seq_len(n_cent)) {
    if (is.null(z_opt)) {
      v <- r_2_opt[, k]
    } else {
      v <- z_opt[, k]
    }
    v[code[, k] != 0L] <- NA
    if (any(is.finite(v))) {
      jb <- which.max(v)
      best_sedrate[k] <- sedrate[jb]
      best_z[k] <- v[jb]
      if (!is.null(p_fw)) {
        p_best[k] <- p_fw[jb, k]
        z_threshold[k] <- zthr_cell[jb, k]
      }
    }
  }

  if (normalize == TRUE) {
    r_2_envelope <- wr_col_normalise(r_2_envelope)
    r_2_power <- wr_col_normalise(r_2_power)
    r_2_opt <- wr_col_normalise(r_2_opt)
  }

  output <- list(
      r_2_envelope = r_2_envelope,
      r_2_power = r_2_power,
      r_2_opt = r_2_opt,
      r_2_envelope_avg = rowMeans(r_2_envelope, na.rm = TRUE),
      r_2_power_avg = rowMeans(r_2_power, na.rm = TRUE),
      r_2_opt_avg = rowMeans(r_2_opt, na.rm = TRUE),
      sedrate = sedrate,
      x = proxy[centres, 1],
      y = proxy[centres, 2],
      linLog = linLog,
      z_envelope = z_envelope,
      z_power = z_power,
      z_opt = z_opt,
      z_envelope_avg = if (is.null(z_envelope)) { NULL } else { rowMeans(z_envelope, na.rm = TRUE) },
      z_power_avg = if (is.null(z_power)) { NULL } else { rowMeans(z_power, na.rm = TRUE) },
      z_opt_avg = if (is.null(z_opt)) { NULL } else { rowMeans(z_opt, na.rm = TRUE) },
      p_fw = p_fw,
      z_threshold = z_threshold,
      best_sedrate = best_sedrate,
      best_z = best_z,
      p_best = p_best,
      significant = if (is.null(p_fw)) { NULL } else { is.finite(p_best) & p_best < alpha },
      untestable_code = code,
      testable = code == 0L,
      noise_family = if (bias_correction == "null") { noise_family } else { NULL },
      noise_param = if (bias_correction == "null") { noise_param } else { NULL },
      npts = npts,
      data = proxy,
      window_mode = window_mode,
      window_size = window_size,
      window_duration = window_duration,
      fit = fit,
      bias_correction = bias_correction,
      noise_model = noise_model,
      alpha = alpha,
      null_cache = if (bias_correction == "null") { null_cache } else { NULL })

  class(output) <- "win_timeOpt.result"
  return(invisible(output))
}

# Internal helpers

# divide every column by its maximum (non-finite values stay NA)
wr_col_normalise <- function(M) {
  for (j in seq_len(ncol(M))) {
    m <- suppressWarnings(max(M[, j], na.rm = TRUE))
    if (is.finite(m) && m > 0) {
      M[, j] <- M[, j] / m
    }
  }
  M
}

# grid nodes (family, parameter, window length) needed for the noise models of all testable cells
wr_to_nodes <- function(noise_family, noise_param, npts, code, window_mode, rho_step, beta_step) {
  ok <- code == 0L & is.finite(noise_param)
  out <- NULL
  if (window_mode == "depth") {
    np <- npts[1, ]
    np <- np[is.finite(np)]
    if (min(np) == max(np)) {
      n_nodes <- min(np)
    } else {
      n_nodes <- unique(round(exp(seq(log(min(np)), log(max(np)), length.out = 3))))
    }
  } else {
    n_nodes <- 0
  }
  for (fam in c(1, 2)) {
    sel <- ok & noise_family == fam
    if (!any(sel)) {
      next
    }
    st <- if (fam == 1) { rho_step } else { beta_step }
    hi_lim <- if (fam == 1) { 0.99 } else { 3 }
    lo <- max(0, floor(min(noise_param[sel]) / st - 1e-9) * st)
    hi <- min(hi_lim, ceiling(max(noise_param[sel]) / st + 1e-9) * st)
    pn <- seq(lo, hi, by = st)
    if (hi > max(pn) + 1e-9) {
      pn <- c(pn, hi)
    }
    pn <- round(pn, 4)
    for (nn in n_nodes) {
      out <- rbind(out, data.frame(family = fam, param = pn, npts = nn))
    }
  }
  if (is.null(out)) {
    out <- data.frame(family = 1, param = 0, npts = if (window_mode == "depth") { n_nodes[1] } else { 0 })
  }
  out
}

# keys and weights of the grid nodes around one cell (linear in the parameter; for depth windows
# bilinear in parameter and log window length)
wr_to_node_weights <- function(nodes, family, param, npts_cell, window_mode) {
  sel <- nodes$family == family
  pg <- sort(unique(nodes$param[sel]))
  wp <- wr_interp_w(pg, param)
  if (window_mode == "depth") {
    ng <- sort(unique(nodes$npts[sel]))
    wn <- wr_interp_w(log(ng), log(npts_cell))
  } else {
    ng <- 0
    wn <- list(i = c(1L, 1L), w = c(1, 0))
  }
  key <- character(4)
  w <- numeric(4)
  q <- 0
  for (a in 1:2) {
    for (b in 1:2) {
      q <- q + 1
      key[q] <- wr_node_key(family, pg[wp$i[a]], ng[wn$i[b]])
      w[q] <- wp$w[a] * wn$w[b]
    }
  }
  list(key = key, w = w)
}
