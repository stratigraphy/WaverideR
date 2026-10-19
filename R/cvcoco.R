#' @title Cross-validated correlation coefficient (cvCOCO) analysis
#'
#' @description The \code{\link{cvcoco}} function performs the bidirectional
#' held-out ("cross-validated") correlation coefficient (cvCOCO) test of
#' Acycle for astronomical forcing in a depth series. The cleaned record is
#' split into two folds: two contiguous halves at the depth midpoint
#' (\code{method = "blocked"}, Blocked cvCOCO, the Acycle public default) or
#' the odd and even observations (\code{method = "interleaved"}, Interleaved
#' cvCOCO). Each fold is linearly interpolated to its own median spacing.
#' In each direction one fold is used for training: at every tested
#' sedimentation rate the energies in four one-Rayleigh union bands (long
#' eccentricity, short eccentricity, obliquity, precession) are
#' leakage-corrected with an exact nonnegative least-squares solve of the
#' finite-record 4x4 leakage matrix, and the rate with the highest
#' spectral correlation gives four trained group amplitudes. These
#' amplitudes are frozen and the held-out fold searches its own
#' sedimentation rate (validation). The same-rate minimum of both
#' directional validation curves is the consensus curve. A single outer
#' stationary AR(1) Monte Carlo reproduces the complete bidirectional
#' train/validate pipeline and gives global (rate-search corrected) and
#' local p-values, the robust directional p-value max(pA, pB) and pCOCO.
#' Two target designs are available: \code{"four-group-coherent-nine"} (the
#' four group amplitudes expanded to nine zero-phase sinusoids summed
#' coherently; Acycle cvcoco9B.m / interleavedcvcoco.m, the GUI default) and
#' \code{"four-group"} (phase-averaged noncoherent four-group target;
#' Acycle cvcoco2.m, for which Acycle issues the confirmatory conclusion
#' report of cocoConclusionReport.m).
#' The COCO method was introduced by Li et al. (2018). Translated from
#' Acycle v3.1 cvcoco.m, cvcoco2.m, cvcoco9B.m, interleavedcvcoco.m,
#' cocoConclusionReport.m and helpers (Mingsong Li).
#'
#' This function calls the following helpers of the cocoR port (which must be
#' loaded): \code{coco_orbit_groups}, \code{detrend_linear},
#' \code{ml_periodogram} and \code{theored_ar1_ml};
#' and, for \code{red = 3}, \code{swa_background}.
#'
#' @param data Matrix with two columns: depth (m) and proxy value.
#' Unsorted, duplicated and irregularly spaced depths are allowed.
#' @param orbit9 Nine astronomical periods (kyr): one long eccentricity, four
#' short eccentricity, one obliquity and three precession periods (groups
#' are identified by descending period rank).
#' \code{Default=c(405.6912, 130.6979, 123.8532, 98.8517, 94.8856, 40.9897,
#' 23.6820, 22.3758, 18.9519)}
#' @param pad Periodogram length (zero padding); each fold uses
#' \code{max(pad, n_fold)}. \code{Default=0} (no padding, as in the Acycle
#' GUI when padding is unchecked)
#' @param sr1 Minimum sedimentation rate (cm/kyr). \code{Default=1}
#' @param sr2 Maximum sedimentation rate (cm/kyr). \code{Default=20}
#' @param srstep Sedimentation-rate step (cm/kyr). \code{Default=0.1}
#' @param red Red-noise background removal: 0 none, 1 classical AR(1),
#' 2 robust AR(1), 3 smoothed-window average (SWA). \code{Default=0}
#' @param n_sim Number of Monte Carlo simulations (0 gives only the
#' deterministic curves). \code{Default=1000}
#' @param cor_method Correlation method, "Pearson" or "Spearman".
#' \code{Default="Pearson"}
#' @param method Split design: "blocked" (depth midpoint halves) or
#' "interleaved" (odd/even observations). \code{Default="blocked"}
#' @param target_model "four-group-coherent-nine" or "four-group" (only for
#' \code{method = "blocked"}). \code{Default="four-group-coherent-nine"}
#' @param max_frequency Maximum temporal frequency (cycle/kyr) used in all
#' spectral correlations. \code{NULL} uses 1.2 times the highest orbital
#' frequency. \code{Default=NULL}
#' @param batch_size Monte Carlo realisations processed simultaneously.
#' \code{Default=100}
#' @param seed Seed for the random-number generator; \code{NULL} leaves the
#' generator untouched. \code{Default=NULL}
#' @param interleaved_phase 0: the first cleaned observation belongs to the
#' Odd fold; 1: it belongs to the Even fold. \code{Default=0}
#' @param genplot Plot the results with \code{\link{plot_cvcoco}}.
#' \code{Default=FALSE}
#' @param verbose Print preprocessing information and a progress bar.
#' \code{Default=FALSE}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li.
#'
#' @references
#' Li, M., Kump, L.R., Hinnov, L.A., Mann, M.E., 2018. Tracking variable
#' sedimentation rates and astronomical forcing in Phanerozoic
#' paleoclimate proxy series with evolutionary correlation coefficients and
#' hypothesis testing. Earth and Planetary Science Letters 501, 165-179.
#' \doi{10.1016/j.epsl.2018.08.041}
#'
#' Li, M., Hinnov, L., Kump, L., 2019. Acycle: Time-series analysis software
#' for paleoclimate research and education. Computers & Geosciences 127,
#' 12-22. \doi{10.1016/j.cageo.2019.02.011}
#'
#' @examples
#' \donttest{
#' # synthetic ETP-like series in depth (5 cm/kyr) with AR(1) noise
#' set.seed(2)
#' t_kyr <- seq(0, 1600, by = 2)
#' etp <- sin(2 * pi * t_kyr / 405) + 0.6 * sin(2 * pi * t_kyr / 125) +
#'   0.5 * sin(2 * pi * t_kyr / 95) + 0.8 * sin(2 * pi * t_kyr / 41) +
#'   0.5 * sin(2 * pi * t_kyr / 23.7) + 0.4 * sin(2 * pi * t_kyr / 19)
#' noise <- as.numeric(stats::filter(stats::rnorm(length(t_kyr)), 0.5,
#'   method = "recursive"))
#' dat <- cbind(t_kyr * 5 / 100, etp + 0.5 * noise)
#' res <- cvcoco(data = dat, sr1 = 2, sr2 = 10, srstep = 0.5,
#'   n_sim = 50, seed = 1)
#' res$consensus$best_rate
#' res$p_robust
#' plot_cvcoco(res)
#' }
#'
#' @return A list with elements
#' \code{name} ("Blocked cvCOCO" or "Interleaved cvCOCO"),
#' \code{target_model}, \code{sr_grid} (tested rates),
#' \code{split_depth}, \code{fold_labels},
#' \code{data_clean}, \code{data_A}, \code{data_B} (regularized folds),
#' \code{raw_index_A}, \code{raw_index_B}, \code{interpolation_A},
#' \code{interpolation_B}, \code{sampling_interval_A},
#' \code{sampling_interval_B},
#' \code{train_A}, \code{train_B} (training curve, best rate, amplitudes,
#' group amplitudes, frozen weights, leakage matrix at best rate),
#' \code{validate_AtoB}, \code{validate_BtoA} (validation curve, best rate,
#' score, directional p, Wilson interval, global and local p curves,
#' participating periods), \code{consensus} (same-rate minimum curve, best
#' rate, global p, p curves, pCOCO curve and maximum),
#' \code{score_symmetric}, \code{score_consensus}, \code{score_mean},
#' \code{p_sym}, \code{p_consensus}, \code{p_robust}, \code{p_A}
#' (B-to-A, fold A held out), \code{p_B} (A-to-B, fold B held out),
#' their confidence intervals, \code{p_curve_AtoB}, \code{p_curve_BtoA},
#' \code{p_local_curve_AtoB}, \code{p_local_curve_BtoA},
#' \code{p_curve_consensus}, \code{p_local_curve_consensus}, \code{pcoco},
#' null distributions (\code{null_symmetric}, \code{null_consensus},
#' \code{null_AtoB}, \code{null_BtoA}, \code{null_best_rate_AtoB},
#' \code{null_best_rate_BtoA}), AR(1) coefficients (\code{rho_A},
#' \code{rho_B}, \code{rho_full}), geometry diagnostics
#' (\code{orbit_count_A/B}, \code{resolved_group_count_A/B},
#' \code{resolved_orbit_mask_A/B}, \code{valid_rate_mask_A/B},
#' \code{training_rate_mask_A/B}, \code{group_leakage_rcond_A/B},
#' \code{all_nine_rate_range_A/B/shared}, \code{active_orbit_count_AtoB/BtoA/
#' consensus}, \code{active_group_count_AtoB/BtoA}), \code{degraded_mode},
#' \code{training_completeness}, \code{spectra} (data and frozen-target
#' temporal spectra at the best rates), \code{n_sim_completed},
#' \code{config}, and \code{conclusion} (confirmatory report; only for
#' \code{target_model = "four-group"}, otherwise \code{NULL}).
#'
#' @export
#' @importFrom stats rnorm
#' @importFrom stats approx
#' @importFrom stats median
#' @importFrom stats cor
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
cvcoco <- function(data = NULL,
                   orbit9 = c(405.6912, 130.6979, 123.8532, 98.8517, 94.8856,
                              40.9897, 23.6820, 22.3758, 18.9519),
                   pad = 0,
                   sr1 = 1,
                   sr2 = 20,
                   srstep = 0.1,
                   red = 0,
                   n_sim = 1000,
                   cor_method = "Pearson",
                   method = "blocked",
                   target_model = "four-group-coherent-nine",
                   max_frequency = NULL,
                   batch_size = 100,
                   seed = NULL,
                   interleaved_phase = 0,
                   genplot = FALSE,
                   verbose = FALSE){
  # ---- internal helper functions (local to this function) ----

  # =========================================================================
  # internal helpers
  # =========================================================================

  # @keywords internal
  # @noRd
  cvcoco_eps <- function(x){
    # MATLAB eps(x): spacing of doubles at abs(x)
    x <- abs(x)
    if (!is.finite(x)){
      return(NaN)
    }
    if (x < 2^-1022){
      return(2^-1074)
    }
    e <- floor(log2(x))
    if (2^e > x){
      e <- e - 1
    }
    if (2^(e + 1) <= x){
      e <- e + 1
    }
    2^(e - 52)
  }

  # @keywords internal
  # @noRd
  cvcoco_colon <- function(a, step, b){
    # sr1:srstep:sr2 as a + k*step (Octave range convention; identical to
    # MATLAB colon for exactly representable grids)
    n <- floor((b - a) / step + 1e-10)
    if (n < 0){
      return(numeric(0))
    }
    a + (0:n) * step
  }

  # @keywords internal
  # @noRd
  cvcoco_prepare_data <- function(data, method, interleaved_phase, analysis_name,
                                  verbose){
    keep <- is.finite(data[, 1]) & is.finite(data[, 2])
    clean <- data[keep, , drop = FALSE]
    if (nrow(clean) < 8){
      stop("At least eight finite observations are required")
    }
    clean <- clean[order(clean[, 1]), , drop = FALSE]
    depth_unique <- unique(clean[, 1])
    grp <- match(clean[, 1], depth_unique)
    sums <- rowsum(clean[, 2], grp, reorder = TRUE)
    counts <- tabulate(grp, nbins = length(depth_unique))
    clean <- cbind(depth_unique, as.numeric(sums) / counts)
    colnames(clean) <- NULL
    if (nrow(clean) < 8){
      stop("At least eight unique depth levels are required after de-duplication")
    }
    n <- nrow(clean)
    split_depth <- NA
    if (method == "interleaved"){
      if (interleaved_phase == 0){
        index_A <- seq(1, n, by = 2)
        index_B <- seq(2, n, by = 2)
      } else {
        index_A <- seq(2, n, by = 2)
        index_B <- seq(1, n, by = 2)
      }
      if (length(index_A) < 4 || length(index_B) < 4){
        stop("Each odd/even fold must contain at least four unique observations")
      }
      label_A <- "Odd fold"
      label_B <- "Even fold"
    } else {
      split_depth <- clean[1, 1] / 2 + clean[n, 1] / 2
      index_A <- which(clean[, 1] <= split_depth)
      index_B <- which(clean[, 1] > split_depth)
      if (length(index_A) < 4 || length(index_B) < 4){
        stop("Each midpoint half must contain at least four unique observations")
      }
      label_A <- "Segment A"
      label_B <- "Segment B"
    }
    raw_A <- clean[index_A, , drop = FALSE]
    raw_B <- clean[index_B, , drop = FALSE]
    A <- cvcoco_regularize(raw_A, label_A)
    B <- cvcoco_regularize(raw_B, label_B)
    dz_A <- stats::median(diff(A[, 1]))
    dz_B <- stats::median(diff(B[, 1]))
    info_A <- cvcoco_interp_info(raw_A, A, dz_A)
    info_B <- cvcoco_interp_info(raw_B, B, dz_B)
    if (verbose == TRUE){
      if (info_A$applied == TRUE){
        message(sprintf("%s: %s linearly interpolated from %d to %d points at %.6g m (largest gap %.3g x median)",
                        analysis_name, label_A, info_A$original_points,
                        info_A$interpolated_points, dz_A, info_A$max_gap_ratio))
      }
      if (info_B$applied == TRUE){
        message(sprintf("%s: %s linearly interpolated from %d to %d points at %.6g m (largest gap %.3g x median)",
                        analysis_name, label_B, info_B$original_points,
                        info_B$interpolated_points, dz_B, info_B$max_gap_ratio))
      }
    }
    list(clean = clean, A = A, B = B, split_depth = split_depth, dz_A = dz_A,
         dz_B = dz_B, info_A = info_A, info_B = info_B, index_A = index_A,
         index_B = index_B)
  }

  # @keywords internal
  # @noRd
  cvcoco_regularize <- function(raw, label){
    spacing <- diff(raw[, 1])
    spacing <- spacing[is.finite(spacing) & spacing > 0]
    if (length(spacing) == 0){
      stop(paste(label, "has no positive finite depth spacing"))
    }
    dz <- stats::median(spacing)
    if (!is.finite(dz) || dz <= 0){
      stop(paste(label, "has no positive finite depth spacing"))
    }
    tolerance <- 1e-5 * abs(dz)
    if (max(abs(spacing - dz)) <= tolerance){
      return(raw)
    }
    n_raw <- nrow(raw)
    count_exact <- (raw[n_raw, 1] - raw[1, 1]) / dz
    count_rounded <- cvcoco_round(count_exact)
    if (abs(count_exact - count_rounded) <= 1e-10 * max(1, abs(count_exact))){
      interval_count <- count_rounded
    } else {
      interval_count <- floor(count_exact)
    }
    if (!is.finite(interval_count + 1) || interval_count + 1 < 4 ||
        interval_count + 1 > 1e6){
      stop(paste(label, "median-spacing interpolation would create an invalid number of points"))
    }
    grid <- raw[1, 1] + (0:interval_count) * dz
    endpoint_tolerance <- 16 * cvcoco_eps(max(1, max(abs(raw[c(1, n_raw), 1])))) *
      max(1, length(grid))
    if (abs(grid[length(grid)] - raw[n_raw, 1]) <= endpoint_tolerance){
      grid[length(grid)] <- raw[n_raw, 1]
    }
    values <- stats::approx(raw[, 1], raw[, 2], xout = grid, method = "linear",
                            rule = 1)$y
    ok <- is.finite(grid) & is.finite(values)
    out <- cbind(grid[ok], values[ok])
    if (nrow(out) < 4){
      stop(paste(label, "contains fewer than four finite interpolated points"))
    }
    out
  }

  # @keywords internal
  # @noRd
  cvcoco_interp_info <- function(raw, out, dz){
    spacing <- diff(raw[, 1])
    tolerance <- 1e-5 * abs(dz)
    deviation <- max(abs(spacing - dz))
    list(applied = deviation > tolerance,
         original_points = nrow(raw),
         interpolated_points = nrow(out),
         original_depth_range = raw[c(1, nrow(raw)), 1],
         interpolated_depth_range = out[c(1, nrow(out)), 1],
         spacing_min = min(spacing), spacing_median = dz,
         spacing_max = max(spacing), max_gap_ratio = max(spacing) / dz,
         tolerance = tolerance, max_spacing_deviation = deviation)
  }

  # @keywords internal
  # @noRd
  cvcoco_resolved_variance <- function(original, detrended){
    n <- length(detrended)
    tolerance <- 256 * .Machine$double.eps
    reference_scale <- max(abs(original))
    residual_scale <- max(abs(detrended))
    out <- list(resolved = FALSE, std = 0)
    if (n < 2 || reference_scale <= 0 || residual_scale <= 0){
      return(out)
    }
    scaled <- detrended / residual_scale
    energy <- sum(scaled^2)
    if (!is.finite(energy) || energy <= 0){
      return(out)
    }
    relative_rms <- (residual_scale / reference_scale) * sqrt(energy / n)
    residual_std <- residual_scale * sqrt(energy / (n - 1))
    out$resolved <- is.finite(relative_rms) && relative_rms > tolerance &&
      is.finite(residual_std) && residual_std > 0
    out$std <- residual_std
    out
  }

  # @keywords internal
  # @noRd
  cvcoco_rho <- function(values){
    # conditional least-squares AR(1) slope on demeaned values (estimateRhoM)
    x <- as.numeric(values)
    x <- x - mean(x)
    scale <- max(abs(x))
    if (length(x) < 4 || any(!is.finite(x)) || !is.finite(scale) || scale <= 0){
      stop("Each fold needs at least four finite, variable detrended values")
    }
    scaled <- x / scale
    n <- length(scaled)
    previous <- scaled[1:(n - 1)]
    nxt <- scaled[2:n]
    denominator <- sum(previous * previous)
    if (!is.finite(denominator) || denominator <= 0){
      stop("The conditional AR(1) denominator is zero or nonfinite")
    }
    rho <- sum(previous * nxt) / denominator
    min(max(rho, -0.999), 0.999)
  }

  # @keywords internal
  # @noRd
  cvcoco_geometry <- function(depth, periods, sr_grid, pad, group_index,
                              target_model, max_frequency){
    n <- length(depth)
    dz <- stats::median(diff(depth))
    fs <- 1 / dz
    nfft <- max(pad, n)
    nf <- floor(nfft / 2) + 1
    frequency <- (0:(nf - 1)) * (fs / nfft)
    rayleigh <- fs / n
    nyquist <- fs / 2
    n_rate <- length(sr_grid)
    is_coherent <- target_model == "four-group-coherent-nine"
    df <- fs / nfft

    geom <- list(n = n, dz = dz, fs = fs, nfft = nfft, frequency = frequency,
                 rayleigh = rayleigh, nyquist = nyquist, df = df,
                 band_lo = matrix(0, n_rate, 9), band_hi = matrix(0, n_rate, 9),
                 basis = vector("list", n_rate),
                 group_template = vector("list", n_rate),
                 group_band_index = vector("list", n_rate),
                 group_unit_energy = matrix(NA_real_, n_rate, 4),
                 group_leakage_matrix = vector("list", n_rate),
                 group_leakage_rcond = rep(NA_real_, n_rate),
                 group_leakage_solve_matrix = vector("list", n_rate),
                 partial_group_leakage_rcond = rep(NA_real_, n_rate),
                 minimum_leakage_rcond = 1e-10,
                 corr_index = vector("list", n_rate),
                 valid_rate_mask = rep(FALSE, n_rate),
                 training_rate_mask = rep(FALSE, n_rate),
                 strict_training_rate_mask = rep(FALSE, n_rate),
                 partial_training_rate_mask = rep(FALSE, n_rate),
                 training_active_group_mask = matrix(FALSE, n_rate, 4),
                 orbit_count = rep(0, n_rate),
                 resolved_orbit_mask = matrix(FALSE, n_rate, 9),
                 resolved_group_count = matrix(0, n_rate, 4),
                 cross_group_band_overlap = rep(FALSE, n_rate),
                 group_index = group_index, target_model = target_model,
                 all_nine_rate_range = c(max(200 * dz / periods),
                                         min(100 * n * dz / periods)),
                 sr_grid = sr_grid, spectrum_scale = 1 / (fs * n))
    one_sided <- rep(2, nf)
    one_sided[1] <- 1
    if (nfft %% 2 == 0){
      one_sided[nf] <- 1
    }
    geom$one_sided_weight <- one_sided
    spectrum_scale <- 1 / (fs * n)

    for (m in seq_len(n_rate)){
      sr <- sr_grid[m]
      spatial_orbit <- 100 / (periods * sr)
      is_resolved <- spatial_orbit >= rayleigh & spatial_orbit < nyquist
      geom$orbit_count[m] <- sum(is_resolved)
      geom$resolved_orbit_mask[m, ] <- is_resolved
      for (g in 1:4){
        geom$resolved_group_count[m, g] <- sum(is_resolved & group_index == g)
      }
      temporal_frequency <- frequency * sr / 100
      corr_index <- which(temporal_frequency <= max_frequency +
                            16 * cvcoco_eps(max(1, max_frequency)))
      geom$corr_index[[m]] <- corr_index
      if (length(corr_index) < 3 || !any(is_resolved)){
        next
      }
      lo <- numeric(9)
      hi <- numeric(9)
      for (j in 1:9){
        in_band <- which(abs(frequency - spatial_orbit[j]) <= rayleigh)
        if (length(in_band) == 0){
          in_band <- which.min(abs(frequency - spatial_orbit[j]))
        }
        lo[j] <- in_band[1]
        hi[j] <- in_band[length(in_band)]
      }
      time <- (0:(n - 1)) * (100 * dz / sr)
      phase <- matrix(0, n, 9)
      for (j in 1:9){
        phase[, j] <- 2 * pi * time / periods[j]
      }
      H <- cvcoco_fft_half(detrend_linear(sin(phase)), nfft, nf)
      if (any(!is.finite(H))){
        next
      }
      geom$band_lo[m, ] <- lo
      geom$band_hi[m, ] <- hi
      Hc <- cvcoco_fft_half(detrend_linear(cos(phase)), nfft, nf)
      if (any(!is.finite(Hc))){
        next
      }
      unit_power <- 0.5 * (Mod(H)^2 + Mod(Hc)^2) * one_sided * spectrum_scale
      template <- matrix(0, nf, 4)
      upper <- min(nyquist, max_frequency * 100 / sr)
      for (j in 1:9){
        if (is_resolved[j] == TRUE){
          template[, group_index[j]] <- template[, group_index[j]] +
            unit_power[, j]
        }
      }
      band_list <- vector("list", 4)
      for (g in 1:4){
        members <- which(group_index == g & is_resolved)
        band_mask <- rep(FALSE, nf)
        for (member in members){
          band_mask[lo[member]:hi[member]] <- TRUE
        }
        band_mask <- band_mask & frequency <= upper + 64 * cvcoco_eps(max(1, upper))
        idx_group <- which(band_mask)
        band_list[[g]] <- idx_group
        if (length(idx_group) > 0){
          geom$group_unit_energy[m, g] <- sum(template[idx_group, g]) * df
        }
      }
      geom$group_band_index[[m]] <- band_list
      geom$cross_group_band_overlap[m] <- cvcoco_band_overlap(
        spatial_orbit, rayleigh, upper, is_resolved, group_index)
      leakage <- matrix(0, 4, 4)
      for (observed in 1:4){
        band <- band_list[[observed]]
        if (length(band) > 0){
          leakage[observed, ] <- colSums(template[band, , drop = FALSE]) * df
        }
      }
      geom$group_leakage_matrix[[m]] <- leakage
      if (all(is.finite(leakage))){
        geom$group_leakage_rcond[m] <- cvcoco_rcond(leakage)
      }
      geom$group_template[[m]] <- template
      if (is_coherent == TRUE){
        geom$basis[[m]] <- H
      }
      geom$valid_rate_mask[m] <- any(template > 0)
      active <- geom$resolved_group_count[m, ] > 0 &
        is.finite(geom$group_unit_energy[m, ]) & geom$group_unit_energy[m, ] > 0
      active[is.na(active)] <- FALSE
      geom$training_active_group_mask[m, ] <- active
      solve_matrix <- leakage
      for (g in which(!active)){
        solve_matrix[g, ] <- 0
        solve_matrix[, g] <- 0
        solve_matrix[g, g] <- 1
      }
      geom$group_leakage_solve_matrix[[m]] <- solve_matrix
      if (any(active)){
        active_matrix <- leakage[active, active, drop = FALSE]
        if (all(is.finite(active_matrix)) && all(diag(active_matrix) > 0)){
          geom$partial_group_leakage_rcond[m] <- cvcoco_rcond(active_matrix)
        }
      }
      prc <- geom$partial_group_leakage_rcond[m]
      geom$partial_training_rate_mask[m] <- geom$valid_rate_mask[m] &&
        any(active) && is.finite(prc) && prc >= geom$minimum_leakage_rcond &&
        !geom$cross_group_band_overlap[m]
      energy <- geom$group_unit_energy[m, ]
      frc <- geom$group_leakage_rcond[m]
      geom$strict_training_rate_mask[m] <- all(is_resolved) &&
        all(is.finite(energy) & energy > 0) && is.finite(frc) &&
        frc >= geom$minimum_leakage_rcond && !geom$cross_group_band_overlap[m]
      geom$training_rate_mask[m] <- geom$strict_training_rate_mask[m]
    }
    geom
  }

  # @keywords internal
  # @noRd
  cvcoco_fft_half <- function(x, nfft, nf){
    x <- as.matrix(x)
    if (nfft > nrow(x)){
      x <- rbind(x, matrix(0, nfft - nrow(x), ncol(x)))
    }
    stats::mvfft(x)[seq_len(nf), , drop = FALSE]
  }

  # @keywords internal
  # @noRd
  cvcoco_rcond <- function(M){
    # MATLAB rcond: LAPACK 1-norm reciprocal condition estimate (0 if singular)
    out <- tryCatch(rcond(M, norm = "O"), error = function(e){
      0
    })
    out
  }

  # @keywords internal
  # @noRd
  cvcoco_band_overlap <- function(center, rayleigh, upper, is_resolved,
                                  group_index){
    active <- which(is_resolved)
    if (length(active) < 2){
      return(FALSE)
    }
    lo <- pmax(0, center - rayleigh)
    hi <- pmin(upper, center + rayleigh)
    scale <- max(c(1, abs(lo), abs(hi), abs(upper)))
    tolerance <- 128 * cvcoco_eps(scale)
    for (aa in 1:(length(active) - 1)){
      first <- active[aa]
      for (bb in (aa + 1):length(active)){
        second <- active[bb]
        if (group_index[first] == group_index[second]){
          next
        }
        width <- min(hi[first], hi[second]) - max(lo[first], lo[second])
        if (width > tolerance){
          return(TRUE)
        }
      }
    }
    FALSE
  }

  # @keywords internal
  # @noRd
  cvcoco_spectra <- function(values, geom, red){
    values <- detrend_linear(as.matrix(values))
    if (any(!is.finite(values))){
      stop("Linear detrending overflowed; rescale the proxy values")
    }
    sp <- ml_periodogram(values, geom$nfft, geom$fs)
    power <- as.matrix(sp$p)
    f <- sp$f
    if (any(!is.finite(power))){
      stop("The periodogram returned nonfinite power; rescale the proxy values")
    }
    negative_tolerance <- 64 * cvcoco_eps(max(1, max(abs(power))))
    if (any(power < -negative_tolerance)){
      stop("The periodogram returned materially negative power")
    }
    power[power < 0] <- 0
    if (red == 3 && nrow(power) < 33){
      stop("SWA background removal requires at least 33 spectral bins; increase pad")
    }
    if (red == 0){
      return(list(power = power, f = f))
    }
    for (k in seq_len(ncol(power))){
      p <- power[, k]
      x <- values[, k]
      if (red == 1){
        background <- theored_ar1_ml(x, f, mean(p), geom$dz)
      } else if (red == 2){
        background <- cvcoco_redconf_any(2 * pi * f * geom$dz, p, geom$dz,
                                         0.25, 2)
      } else {
        positive <- p[is.finite(p) & p > 0]
        if (length(positive) == 0){
          stop("SWA background removal requires a positive periodogram ordinate")
        }
        log_floor <- max(.Machine$double.xmin,
                         max(positive) * .Machine$double.eps)
        p_log <- p
        p_log[p_log < log_floor] <- log_floor
        background <- swa_background(freq = f, log_power = log10(p_log),
                                     n_data = length(x), verbose = FALSE)$swa
      }
      background <- as.numeric(background)
      if (length(background) != length(p) || any(!is.finite(background)) ||
          any(background <= 0)){
        stop("The red-noise method returned an invalid spectral background")
      }
      p <- p - background
      p[p < 0] <- 0
      power[, k] <- p
    }
    list(power = power, f = f)
  }

  # @keywords internal
  # @noRd
  cvcoco_round <- function(x){
    # MATLAB round: halves away from zero (R's round() rounds half to even)
    sign(x) * floor(abs(x) + 0.5)
  }

  # @keywords internal
  # @noRd
  cvcoco_move_median <- function(x, w){
    # moveMedian.m: w-point moving median with shrinking end windows
    x <- as.numeric(x)
    n <- length(x)
    halfw <- floor(w / 2)
    if (w %% 2 == 1){
      back <- halfw
      ahead <- halfw
    } else {
      back <- max(0, halfw - 1)
      ahead <- halfw
    }
    y <- numeric(n)
    for (i in seq_len(n)){
      y[i] <- stats::median(x[max(1, i - back):min(n, i + ahead)])
    }
    y
  }

  # @keywords internal
  # @noRd
  cvcoco_rho_s0_distance <- function(rho_grid, s0_grid, cospara, power, linlog){
    n_rho <- length(rho_grid)
    n_s0 <- length(s0_grid)
    shape <- matrix(0, length(cospara), n_rho)
    for (i in seq_len(n_rho)){
      shape[, i] <- (1 - rho_grid[i]^2) /
        (1 - 2 * cospara * rho_grid[i] + rho_grid[i]^2)
    }
    if (any(!is.finite(shape) | shape <= 0) || any(!is.finite(s0_grid) |
                                                    s0_grid <= 0)){
      stop("The rho/S0 search produced a nonpositive theoretical spectrum")
    }
    if (linlog == 1){
      a <- colSums(shape^2)
      b <- as.numeric(crossprod(shape, power))
      distance <- matrix(a, n_rho, n_s0) *
        matrix(s0_grid^2, n_rho, n_s0, byrow = TRUE) -
        2 * matrix(b, n_rho, n_s0) * matrix(s0_grid, n_rho, n_s0, byrow = TRUE) +
        sum(power^2)
    } else {
      base <- log(shape) - log(power)
      rs <- colSums(base)
      rss <- colSums(base^2)
      ls <- log(s0_grid)
      distance <- matrix(rss, n_rho, n_s0) +
        2 * matrix(rs, n_rho, n_s0) * matrix(ls, n_rho, n_s0, byrow = TRUE) +
        length(power) * matrix(ls^2, n_rho, n_s0, byrow = TRUE)
    }
    tol <- 256 * cvcoco_eps(max(1, max(abs(distance))))
    distance[distance < 0 & distance >= -tol] <- 0
    if (any(!is.finite(distance)) || any(distance < 0)){
      stop("The rho/S0 search produced an invalid squared distance")
    }
    distance
  }

  # @keywords internal
  # @noRd
  cvcoco_minirhos0 <- function(s0, fn, ft, pxxsmooth, linlog = 2){
    # minirhos0.m: two-stage grid search for the AR(1) rho and S0
    cospara <- cos(pi * as.numeric(ft) / fn)
    power <- as.numeric(pxxsmooth)
    if (any(!is.finite(power)) || any(power <= 0)){
      stop("The median-smoothed spectrum must be finite and positive (red = 2)")
    }
    nn <- 50
    rhoi <- seq(0.001, 0.999, length.out = nn)
    s0i <- seq(0.2 * s0, 5 * s0, length.out = nn)
    d <- cvcoco_rho_s0_distance(rhoi, s0i, cospara, power, linlog)
    k <- which.min(d)
    x <- (k - 1) %% nrow(d) + 1
    y <- (k - 1) %/% nrow(d) + 1
    mm <- nn / 2
    for (it in 1:3){
      s0max <- s0i[y] + (s0i[2] - s0i[1])
      s0min <- s0i[y] - (s0i[2] - s0i[1])
      rhomax <- rhoi[x] + (rhoi[2] - rhoi[1])
      rhomin <- rhoi[x] - (rhoi[2] - rhoi[1])
      if (rhomax >= 1){
        rhomax <- 0.9999
      }
      if (rhomin <= 0){
        rhomin <- 0.0001
      }
      rhoi <- seq(rhomin, rhomax, length.out = mm)
      s0i <- seq(s0min, s0max, length.out = mm)
      d <- cvcoco_rho_s0_distance(rhoi, s0i, cospara, power, linlog)
      k <- which.min(d)
      x <- (k - 1) %% nrow(d) + 1
      y <- (k - 1) %/% nrow(d) + 1
    }
    list(rho = rhoi[x], s0 = s0i[y])
  }

  # @keywords internal
  # @noRd
  cvcoco_redconf_any <- function(f, pxx, dt, smoothwin = 0.25, linlog = 2){
    # redconf_any.m: robust AR(1) background of a median-smoothed spectrum;
    # f in rad/sample
    fn <- 1 / (2 * dt)
    ft <- f / pi * fn
    smoothn <- cvcoco_round(smoothwin * length(pxx))
    pxxsmooth <- cvcoco_move_median(pxx, smoothn)
    s0 <- mean(pxxsmooth)
    fit <- cvcoco_minirhos0(s0, fn, ft, pxxsmooth, linlog)
    fit$s0 * (1 - fit$rho^2) / (1 - (2 * fit$rho * cos(pi * ft / fn)) +
                                  fit$rho^2)
  }

  # @keywords internal
  # @noRd
  cvcoco_col_max <- function(x){
    x <- as.matrix(x)
    out <- x[1, ]
    if (nrow(x) > 1){
      for (i in 2:nrow(x)){
        out <- pmax(out, x[i, ])
      }
    }
    out
  }

  # @keywords internal
  # @noRd
  cvcoco_rank <- function(M){
    d <- svd(M, nu = 0, nv = 0)$d
    if (length(d) == 0){
      return(0)
    }
    tol <- max(dim(M)) * cvcoco_eps(max(d))
    sum(d > tol)
  }

  # @keywords internal
  # @noRd
  cvcoco_nnls4 <- function(M, E){
    # exact nonnegative least squares min ||M q - E||, q >= 0, by enumerating
    # all 15 active sets (cocoNonnegativeLeakageSolve.m)
    E <- as.matrix(E)
    n_series <- ncol(E)
    col_scale <- cvcoco_col_max(E)
    col_scale[!is.finite(col_scale) | col_scale <= 0] <- 1
    S <- E / matrix(col_scale, 4, n_series, byrow = TRUE)
    G <- matrix(0, 4, n_series)
    best <- colSums(S^2)
    for (code in 1:15){
      active <- bitwAnd(code, c(1L, 2L, 4L, 8L)) > 0
      n_active <- sum(active)
      D <- M[, active, drop = FALSE]
      if (cvcoco_rank(D) < n_active){
        next
      }
      if (n_active == 4){
        cand <- solve(D, S)
      } else {
        cand <- qr.coef(qr(D), S)
      }
      cand <- matrix(cand, nrow = n_active)
      tol <- 256 * .Machine$double.eps * pmax(1, cvcoco_col_max(abs(cand)))
      ok <- cand >= -matrix(tol, n_active, n_series, byrow = TRUE)
      ok[is.na(ok)] <- FALSE
      feasible <- colSums(!ok) == 0
      if (!any(feasible)){
        next
      }
      cand <- pmax(cand, 0)
      residual <- colSums((D %*% cand - S)^2)
      better <- feasible & is.finite(residual) & residual < best
      if (any(better)){
        G[, better] <- 0
        G[active, better] <- cand[, better, drop = FALSE]
        best[better] <- residual[better]
      }
    }
    if (any(!is.finite(best))){
      stop("No finite nonnegative least-squares residual was found")
    }
    G * matrix(col_scale, 4, n_series, byrow = TRUE)
  }

  # @keywords internal
  # @noRd
  cvcoco_target_power <- function(geom, m, group_amplitude){
    group_amplitude <- as.matrix(group_amplitude)
    if (geom$target_model == "four-group-coherent-nine"){
      expanded <- group_amplitude[geom$group_index, , drop = FALSE] *
        as.numeric(geom$resolved_orbit_mask[m, ])
      target_fft <- geom$basis[[m]] %*% expanded
      target <- Mod(target_fft)^2 * geom$one_sided_weight * geom$spectrum_scale
    } else {
      target <- geom$group_template[[m]] %*% (group_amplitude^2)
    }
    target
  }

  # @keywords internal
  # @noRd
  cvcoco_train <- function(power, geom, cor_method, keep_curve){
    power <- as.matrix(power)
    n_rate <- length(geom$sr_grid)
    n_series <- ncol(power)
    curve <- NULL
    if (keep_curve == TRUE){
      curve <- matrix(NA_real_, n_rate, n_series)
    }
    best_score <- rep(-Inf, n_series)
    best_rate <- rep(NA_real_, n_series)
    amplitude <- matrix(NA_real_, 9, n_series)
    tol2 <- (1e-8)^2
    for (m in seq_len(n_rate)){
      if (geom$training_rate_mask[m] == FALSE){
        next
      }
      active <- geom$training_active_group_mask[m, ]
      energy <- matrix(0, 4, n_series)
      for (g in 1:4){
        if (active[g] == FALSE){
          next
        }
        band <- geom$group_band_index[[m]][[g]]
        if (length(band) == 0){
          stop("An eligible training group has no frequency-integration band")
        }
        energy[g, ] <- colSums(power[band, , drop = FALSE]) * geom$df
      }
      if (any(!is.finite(energy)) || any(energy < 0)){
        stop("A four-group integration band produced nonfinite or negative energy")
      }
      energy[!active, ] <- 0
      group_power <- cvcoco_nnls4(geom$group_leakage_solve_matrix[[m]], energy)
      group_power[!active, ] <- 0
      max_power <- cvcoco_col_max(group_power)
      negligible <- group_power <= matrix(max_power * tol2, 4, n_series,
                                          byrow = TRUE)
      group_power[negligible] <- 0
      group_amplitude <- sqrt(group_power)
      target <- cvcoco_target_power(geom, m, group_amplitude)
      idx <- geom$corr_index[[m]]
      rho <- cvcoco_column_correlation(power[idx, , drop = FALSE],
                                       target[idx, , drop = FALSE], cor_method)
      if (keep_curve == TRUE){
        curve[m, ] <- rho
      }
      better <- is.finite(rho) & rho > best_score
      if (any(better)){
        best_score[better] <- rho[better]
        best_rate[better] <- geom$sr_grid[m]
        amplitude[, better] <- group_amplitude[geom$group_index, better,
                                               drop = FALSE]
      }
    }
    best_score[!is.finite(best_rate)] <- NA
    list(curve = curve, best_score = best_score, best_rate = best_rate,
         amplitude = amplitude)
  }

  # @keywords internal
  # @noRd
  cvcoco_group_amplitudes <- function(amplitude, group_index){
    amplitude <- as.matrix(amplitude)
    n_series <- ncol(amplitude)
    raw <- matrix(NA_real_, 4, n_series)
    for (g in 1:4){
      a <- amplitude[group_index == g, , drop = FALSE]
      raw[g, ] <- sqrt(colMeans(a^2))
    }
    scale <- cvcoco_col_max(raw)
    normalized <- raw / matrix(scale, 4, n_series, byrow = TRUE)
    finite_col <- colSums(!is.finite(raw)) == 0
    zero_target <- is.finite(scale) & scale == 0 & finite_col
    normalized[, zero_target] <- 0
    bad <- !is.finite(scale) | scale < 0
    bad[is.na(bad)] <- TRUE
    normalized[, bad] <- NA
    small <- !is.na(normalized) & normalized >= 0 & normalized <= 1e-8
    normalized[small] <- 0
    list(raw = raw, normalized = normalized)
  }

  # @keywords internal
  # @noRd
  cvcoco_validate <- function(power, geom, weights, cor_method){
    power <- as.matrix(power)
    weights <- as.matrix(weights)
    n_rate <- length(geom$sr_grid)
    n_series <- ncol(power)
    if (ncol(weights) == 1 && n_series > 1){
      weights <- matrix(weights, 4, n_series)
    }
    if (nrow(weights) != 4 || ncol(weights) != n_series){
      stop("Fixed validation weights have an invalid size")
    }
    curve <- matrix(NA_real_, n_rate, n_series)
    best_score <- rep(-Inf, n_series)
    best_rate <- rep(NA_real_, n_series)
    for (m in seq_len(n_rate)){
      if (geom$valid_rate_mask[m] == FALSE){
        next
      }
      target <- cvcoco_target_power(geom, m, weights)
      idx <- geom$corr_index[[m]]
      rho <- cvcoco_column_correlation(power[idx, , drop = FALSE],
                                       target[idx, , drop = FALSE], cor_method)
      curve[m, ] <- rho
      better <- is.finite(rho) & rho > best_score
      if (any(better)){
        best_score[better] <- rho[better]
        best_rate[better] <- geom$sr_grid[m]
      }
    }
    best_score[!is.finite(best_rate)] <- NA
    list(curve = curve, best_score = best_score, best_rate = best_rate)
  }

  # @keywords internal
  # @noRd
  cvcoco_column_correlation <- function(x, y, cor_method){
    x <- as.matrix(x)
    y <- as.matrix(y)
    n_series <- ncol(x)
    tiny <- .Machine$double.xmin
    rel <- 256 * .Machine$double.eps
    if (cor_method == "Spearman"){
      rho <- rep(NA_real_, n_series)
      for (k in seq_len(n_series)){
        ok <- is.finite(x[, k]) & is.finite(y[, k])
        if (sum(ok) >= 2){
          xk <- x[ok, k]
          yk <- y[ok, k]
          if (cvcoco_effectively_constant(xk) || cvcoco_effectively_constant(yk)){
            rho[k] <- 0
          } else {
            rho[k] <- stats::cor(xk, yk, method = "spearman")
          }
        }
      }
      return(rho)
    }
    ok <- is.finite(x) & is.finite(y)
    n <- colSums(ok)
    x0 <- x
    y0 <- y
    x0[!ok] <- 0
    y0[!ok] <- 0
    nr <- nrow(x)
    mx <- colSums(x0) / pmax(n, 1)
    my <- colSums(y0) / pmax(n, 1)
    xc <- (x0 - matrix(mx, nr, n_series, byrow = TRUE)) * ok
    yc <- (y0 - matrix(my, nr, n_series, byrow = TRUE)) * ok
    numerator <- colSums(xc * yc)
    ssx <- colSums(xc^2)
    ssy <- colSums(yc^2)
    scale_x <- cvcoco_col_max(abs(x0))
    scale_y <- cvcoco_col_max(abs(y0))
    tol_x <- n * (rel * pmax(scale_x, tiny))^2
    tol_y <- n * (rel * pmax(scale_y, tiny))^2
    rho <- numerator / sqrt(ssx * ssy)
    degenerate <- n >= 2 & (ssx <= tol_x | ssy <= tol_y)
    rho[degenerate] <- 0
    rho[n < 2 | !is.finite(rho)] <- NA
    rho
  }

  # @keywords internal
  # @noRd
  cvcoco_effectively_constant <- function(x){
    x <- x[is.finite(x)]
    if (length(x) < 2){
      return(TRUE)
    }
    centered <- x - mean(x)
    scale <- max(abs(x))
    tolerance <- length(x) * (256 * .Machine$double.eps *
                                max(scale, .Machine$double.xmin))^2
    !all(is.finite(centered)) || sum(centered^2) <= tolerance
  }

  # @keywords internal
  # @noRd
  cvcoco_ar1_filter <- function(innovations, rho){
    # MATLAB filter(1,[1,-rho],innovations,[],1): y(t) = e(t) + rho*y(t-1)
    y <- as.matrix(innovations)
    n <- nrow(y)
    if (n > 1){
      for (t in 2:n){
        y[t, ] <- y[t, ] + rho * y[t - 1, ]
      }
    }
    y
  }

  # @keywords internal
  # @noRd
  cvcoco_min_pair <- function(x, y){
    z <- x
    z[] <- NA
    ok <- is.finite(x) & is.finite(y)
    z[ok] <- pmin(x[ok], y[ok])
    z
  }

  # @keywords internal
  # @noRd
  cvcoco_max_pair <- function(x, y){
    z <- x
    z[] <- NA
    ok <- is.finite(x) & is.finite(y)
    z[ok] <- pmax(x[ok], y[ok])
    z
  }

  # @keywords internal
  # @noRd
  cvcoco_curve_max <- function(curve, sr_grid){
    curve <- as.numeric(curve)
    finite_index <- which(is.finite(curve))
    if (length(finite_index) == 0){
      return(list(score = NA_real_, rate = NA_real_, index = NA_integer_))
    }
    relative <- which.max(curve[finite_index])
    index <- finite_index[relative]
    list(score = curve[index], rate = sr_grid[index], index = index)
  }

  # @keywords internal
  # @noRd
  cvcoco_value_at <- function(x, index){
    if (length(index) == 0 || !is.finite(index)){
      return(NA_real_)
    }
    x[index]
  }

  # @keywords internal
  # @noRd
  cvcoco_pcoco <- function(rho, p){
    p_safe <- p
    p_safe[!is.finite(p_safe) | p_safe <= 0] <- NA
    p_safe[!is.na(p_safe) & p_safe > 1] <- 1
    rho * abs(log10(p_safe))
  }

  # @keywords internal
  # @noRd
  cvcoco_accumulate <- function(simulated, observed, exceed, valid){
    simulated <- as.matrix(simulated)
    n_series <- ncol(simulated)
    obs <- matrix(observed, length(observed), n_series)
    finite <- is.finite(simulated) & is.finite(obs)
    ge <- finite & simulated >= obs
    ge[is.na(ge)] <- FALSE
    list(exceed = exceed + rowSums(ge), valid = valid + rowSums(finite))
  }

  # @keywords internal
  # @noRd
  cvcoco_check_counts <- function(valid, observed, n_sim_completed, label){
    if (n_sim_completed > 0 &&
        any(valid[is.finite(observed)] != n_sim_completed)){
      stop(sprintf("At least one %s null realization produced a nonfinite same-rate statistic",
                   label))
    }
    invisible(TRUE)
  }

  # @keywords internal
  # @noRd
  cvcoco_pointwise_p <- function(exceed, valid, observed){
    p <- rep(NA_real_, length(observed))
    ok <- is.finite(observed) & valid > 0
    p[ok] <- (exceed[ok] + 1) / (valid[ok] + 1)
    p
  }

  # @keywords internal
  # @noRd
  cvcoco_wilson <- function(k, n){
    if (n <= 0){
      return(c(NA_real_, NA_real_))
    }
    z <- 1.95996398454005
    phat <- k / n
    denominator <- 1 + z^2 / n
    center <- (phat + z^2 / (2 * n)) / denominator
    half <- z / denominator * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2))
    c(max(0, center - half), min(1, center + half))
  }

  # @keywords internal
  # @noRd
  cvcoco_mc_p <- function(null_values, observed){
    null_values <- null_values[is.finite(null_values)]
    n_valid <- length(null_values)
    if (n_valid == 0 || !is.finite(observed)){
      return(list(p = NA_real_, n_valid = n_valid, n_exceed = 0,
                  ci = c(NA_real_, NA_real_)))
    }
    n_exceed <- sum(null_values >= observed)
    list(p = (n_exceed + 1) / (n_valid + 1), n_valid = n_valid,
         n_exceed = n_exceed, ci = cvcoco_wilson(n_exceed, n_valid))
  }

  # @keywords internal
  # @noRd
  cvcoco_max_stat_p_curve <- function(null_maximum, observed){
    p <- rep(NA_real_, length(observed))
    null_maximum <- sort(null_maximum[is.finite(null_maximum)])
    n_null <- length(null_maximum)
    if (n_null == 0){
      return(p)
    }
    for (i in which(is.finite(observed))){
      n_exceed <- sum(null_maximum >= observed[i])
      p[i] <- (1 + n_exceed) / (n_null + 1)
    }
    p
  }

  # @keywords internal
  # @noRd
  cvcoco_active_counts <- function(geom, weights){
    active_group <- as.numeric(weights) > 1e-8
    active_orbit <- active_group[geom$group_index]
    mask <- geom$resolved_orbit_mask &
      matrix(active_orbit, nrow(geom$resolved_orbit_mask), 9, byrow = TRUE)
    group_counts <- matrix(0, nrow(mask), 4)
    for (g in 1:4){
      group_counts[, g] <- rowSums(mask[, geom$group_index == g, drop = FALSE])
    }
    list(mask = mask, count = rowSums(mask), group_counts = group_counts,
         active_orbit = active_orbit)
  }

  # @keywords internal
  # @noRd
  cvcoco_train_result <- function(tr, grp, sr_grid, geom){
    best_rate <- tr$best_rate[1]
    best_index <- which(sr_grid == best_rate)[1]
    curve <- as.numeric(tr$curve)
    out <- list(curve = curve, best_rate = best_rate, best_index = best_index,
                best_correlation = cvcoco_value_at(curve, best_index),
                amplitudes9 = as.numeric(tr$amplitude),
                group_raw = as.numeric(grp$raw),
                group_normalized = as.numeric(grp$normalized),
                frozen_validation_weights = as.numeric(grp$normalized),
                training_rate_mask = geom$training_rate_mask,
                group_leakage_matrix = matrix(NA_real_, 4, 4),
                group_leakage_rcond = NA_real_,
                resolved_group_mask = rep(FALSE, 4),
                resolvable_period_count = NA_real_,
                training_completeness = "unresolved")
    if (is.finite(best_index)){
      out$group_leakage_matrix <- geom$group_leakage_matrix[[best_index]]
      out$group_leakage_rcond <- geom$partial_group_leakage_rcond[best_index]
      out$resolved_group_mask <- geom$training_active_group_mask[best_index, ]
      out$resolvable_period_count <- geom$orbit_count[best_index]
      if (out$resolvable_period_count < 9){
        out$training_completeness <- "partial-orbit"
      } else {
        out$training_completeness <- "complete-nine"
      }
    }
    out
  }

  # @keywords internal
  # @noRd
  cvcoco_validation_result <- function(curve, score, best_rate, sr_grid, geom,
                                       weights){
    best_index <- which(sr_grid == best_rate)[1]
    out <- list(curve = curve, best_rate = best_rate, best_index = best_index,
                score = score, valid_rate_mask = geom$valid_rate_mask,
                resolvable_period_count = NA_real_,
                resolvable_group_counts = rep(NA_real_, 4),
                participating_period_count = NA_real_,
                participating_group_counts = rep(NA_real_, 4),
                all_nine_periods_participate = FALSE)
    if (is.finite(best_index)){
      act <- cvcoco_active_counts(geom, weights)
      out$resolvable_period_count <- geom$orbit_count[best_index]
      out$resolvable_group_counts <- geom$resolved_group_count[best_index, ]
      out$participating_period_count <- act$count[best_index]
      out$participating_group_counts <- act$group_counts[best_index, ]
      out$all_nine_periods_participate <- act$count[best_index] == 9
    }
    out
  }

  # @keywords internal
  # @noRd
  cvcoco_spectrum_diagnostic <- function(power, geom, rate, weights, group_index,
                                         mode){
    out <- list(rate = rate, frequency = numeric(0), data_power = numeric(0),
                target_power = numeric(0), mode = mode)
    if (!is.finite(rate)){
      return(out)
    }
    m <- which(abs(geom$sr_grid - rate) <= 16 * cvcoco_eps(max(1, abs(rate))))[1]
    if (is.na(m) || geom$valid_rate_mask[m] == FALSE){
      return(out)
    }
    if (mode == "adaptive"){
      if (length(weights) != 9 || any(!is.finite(weights))){
        return(out)
      }
      group_weights <- numeric(4)
      for (g in 1:4){
        group_weights[g] <- sqrt(mean(weights[group_index == g]^2))
      }
    } else {
      if (length(weights) != 4 || any(!is.finite(weights))){
        return(out)
      }
      group_weights <- weights
    }
    target <- cvcoco_target_power(geom, m, matrix(group_weights, 4, 1))
    idx <- geom$corr_index[[m]]
    jacobian <- 100 / rate
    out$frequency <- geom$frequency[idx] * rate / 100
    out$data_power <- as.matrix(power)[idx, 1] * jacobian
    out$target_power <- as.numeric(target[idx, 1]) * jacobian
    out
  }

  # @keywords internal
  # @noRd
  cvcoco_conclusion <- function(res){
    # numerical core of cocoConclusionReport('confirmatory', cv)
    alpha <- 0.05
    p_AtoB <- suppressWarnings(min(res$p_curve_AtoB[is.finite(res$p_curve_AtoB)]))
    p_BtoA <- suppressWarnings(min(res$p_curve_BtoA[is.finite(res$p_curve_BtoA)]))
    if (!is.finite(p_AtoB) || !is.finite(p_BtoA)){
      return(list(available = FALSE,
                  conclusion = "Monte Carlo p-values are required for a confirmatory report"))
    }
    i_AtoB <- res$validate_AtoB$best_index
    i_BtoA <- res$validate_BtoA$best_index
    period_AtoB <- res$active_orbit_count_AtoB[i_AtoB]
    period_BtoA <- res$active_orbit_count_BtoA[i_BtoA]
    resolvable_AtoB <- res$orbit_count_B[i_AtoB]
    resolvable_BtoA <- res$orbit_count_A[i_BtoA]
    all_nine <- period_AtoB == 9 && period_BtoA == 9
    p_robust <- max(p_AtoB, p_BtoA)
    p_sym <- res$p_sym
    directional_pass <- p_robust < alpha
    symmetric_pass <- is.finite(p_sym) && p_sym < alpha
    ci_A <- res$p_A_confidence_interval
    ci_B <- res$p_B_confidence_interval
    ci_sym <- res$p_sym_confidence_interval
    resolved_A <- (res$p_A < alpha && ci_A[2] < alpha) ||
      (res$p_A >= alpha && ci_A[1] >= alpha)
    resolved_B <- (res$p_B < alpha && ci_B[2] < alpha) ||
      (res$p_B >= alpha && ci_B[1] >= alpha)
    mc_resolved <- resolved_A && resolved_B
    sym_resolved <- (p_sym < alpha && ci_sym[2] < alpha) ||
      (p_sym >= alpha && ci_sym[1] >= alpha)
    degraded <- isTRUE(res$degraded_mode)
    if (degraded == TRUE){
      if (directional_pass && symmetric_pass){
        classification <- "Partial-orbit exploratory; directional and symmetric criteria met"
      } else if (directional_pass){
        classification <- "Partial-orbit exploratory; directional criterion met"
      } else if (symmetric_pass){
        classification <- "Partial-orbit exploratory; joint-only / inconclusive"
      } else {
        classification <- "Partial-orbit exploratory; negative"
      }
    } else if (directional_pass && symmetric_pass){
      classification <- "Robust positive"
    } else if (directional_pass){
      classification <- "Directional robust; symmetric test not significant"
    } else if (symmetric_pass){
      classification <- "Joint-only / inconclusive"
    } else {
      classification <- "Negative"
    }
    if (!all_nine && directional_pass){
      classification <- paste0(classification, "; partial-period target")
    }
    if (!mc_resolved){
      classification <- paste0(classification, "; Monte Carlo threshold uncertain")
    }
    if (!isTRUE(sym_resolved)){
      classification <- paste0(classification, "; p_sym threshold uncertain")
    }
    list(available = TRUE,
         mode = ifelse(degraded, "partial-orbit-exploratory", "confirmatory"),
         confirmatory_eligible = !degraded,
         alpha = alpha,
         p_global_AtoB = p_AtoB, p_global_BtoA = p_BtoA,
         p_A = p_BtoA, p_B = p_AtoB, p_robust = p_robust, p_sym = p_sym,
         best_rate_AtoB = res$validate_AtoB$best_rate,
         best_rate_BtoA = res$validate_BtoA$best_rate,
         participating_periods_AtoB = period_AtoB,
         participating_periods_BtoA = period_BtoA,
         resolvable_periods_AtoB = resolvable_AtoB,
         resolvable_periods_BtoA = resolvable_BtoA,
         all_nine_at_both_best_rates = all_nine,
         monte_carlo_threshold_resolved = mc_resolved,
         p_sym_threshold_resolved = sym_resolved,
         directional_pass = directional_pass,
         symmetric_pass = symmetric_pass,
         pass = (!degraded) && directional_pass,
         classification = classification)
  }

  detrend_linear <- function(x){
    # Least-squares removal of a straight line from every column (MATLAB detrend(x,1)).
    xm <- as.matrix(x)
    n <- nrow(xm)
    if (n < 2){
      if (is.matrix(x)){
        return(xm * 0)
      } else {
        return(as.vector(xm) * 0)
      }
    }
    t <- seq_len(n)
    X <- cbind(1, t)
    q <- qr(X)
    res <- qr.resid(q, xm)
    if (is.matrix(x)){
      res
    } else {
      as.vector(res)
    }
  }

  # --- MATLAB periodogram (rectangular window, one-sided PSD) -------------------
  ml_periodogram <- function(x, nfft, fs){
    # x: vector or samples-by-series matrix. Every column is zero padded to nfft.
    xm <- as.matrix(x)
    n <- nrow(xm)
    if (nfft < n){
      stop("pad (nfft) must be >= number of samples")
    }
    if (nfft > n){
      xm <- rbind(xm, matrix(0, nfft - n, ncol(xm)))
    }
    X <- mvfft(xm)
    if (nfft %% 2 == 0){
      nf <- nfft / 2 + 1
    } else {
      nf <- (nfft + 1) / 2
    }
    P <- (Mod(X[seq_len(nf), , drop = FALSE])^2) / (fs * n)
    if (nf > 2){
      if (nfft %% 2 == 0){
        P[2:(nf - 1), ] <- 2 * P[2:(nf - 1), ]
      } else {
        P[2:nf, ] <- 2 * P[2:nf, ]
      }
    } else if (nf == 2 && nfft %% 2 == 1){
      P[2, ] <- 2 * P[2, ]
    }
    f <- (0:(nf - 1)) * fs / nfft
    if (is.matrix(x)){
      p_out <- P
    } else {
      p_out <- as.vector(P)
    }
    list(p = p_out, f = f)
  }

  # --- AR(1) red-noise backgrounds ---------------------------------------------
  rho_ar1_ml <- function(x){
    # Mudelsee (2002) lag-1 estimator on demeaned data (Husson 2012 / Li 2017).
    x <- as.vector(x) - mean(x)
    n <- length(x)
    sum(x[2:n] * x[1:(n - 1)]) / sum(x[1:(n - 1)]^2)
  }

  theored_ar1_ml <- function(datax, f, meanp, dt){
    rho <- rho_ar1_ml(datax)
    fnyq <- 1 / (2 * dt)
    theored <- (1 - rho^2) / (1 - 2 * rho * cos(pi * f / fnyq) + rho^2)
    tmp <- theored
    tmp[1] <- 0
    art <- mean(tmp)
    theored * meanp / art
  }

  # Physical groups of a 9-period target by descending rank (cocoOrbitGroups.m)
  coco_orbit_groups <- function(orbit9){
    if (length(orbit9) != 9 || length(unique(orbit9)) != 9){
      stop("orbit9 must contain nine distinct periods")
    }
    o <- order(orbit9, decreasing = TRUE)
    idx <- integer(9)
    idx[o[1]] <- 1L
    idx[o[2:5]] <- 2L
    idx[o[6]] <- 3L
    idx[o[7:9]] <- 4L
    idx
  }

  # ---- main code ----

  # ---------------------------------------------------------------- checks
  if (is.null(data)){
    stop("data must be a two-column matrix [depth, value]")
  }
  data <- as.matrix(data)
  if (ncol(data) < 2 || nrow(data) < 1){
    stop("data must be a two-column matrix [depth, value]")
  }
  data <- data[, 1:2, drop = FALSE]
  method <- match.arg(method, c("blocked", "interleaved"))
  target_model <- match.arg(target_model,
                            c("four-group-coherent-nine", "four-group"))
  cor_method <- match.arg(cor_method, c("Pearson", "Spearman"))
  is_interleaved <- method == "interleaved"
  if (is_interleaved == TRUE && target_model != "four-group-coherent-nine"){
    stop("Interleaved cvCOCO fixes its target design to 'four-group-coherent-nine'")
  }
  is_coherent <- target_model == "four-group-coherent-nine"
  orbit9 <- as.numeric(orbit9)
  if (length(orbit9) != 9 || any(!is.finite(orbit9)) || any(orbit9 <= 0)){
    stop("orbit9 must contain nine finite positive periods")
  }
  if (length(unique(orbit9)) != 9){
    stop("orbit9 must contain nine distinct orbital periods")
  }
  if (length(pad) != 1 || !is.finite(pad) || pad < 0 || pad != round(pad)){
    stop("pad must be a nonnegative integer")
  }
  if (!is.finite(sr1) || sr1 <= 0 || !is.finite(sr2) || sr2 < sr1 ||
      !is.finite(srstep) || srstep <= 0){
    stop("sr1 and srstep must be positive and sr2 >= sr1")
  }
  if (!(red %in% c(0, 1, 2, 3))){
    stop("red must be 0, 1, 2 or 3")
  }
  if (!is.finite(n_sim) || n_sim < 0 || n_sim != round(n_sim) || n_sim > 1e6){
    stop("n_sim must be an integer between 0 and 1e6")
  }
  if (!is.finite(batch_size) || batch_size < 1 || batch_size != round(batch_size)){
    stop("batch_size must be a positive integer")
  }
  if (!(interleaved_phase %in% c(0, 1))){
    stop("interleaved_phase must be 0 or 1")
  }
  if (red == 3 && !exists("swa_background", mode = "function")){
    stop("red = 3 requires swa_background() (R/swa_background.R) to be loaded")
  }

  highest_frequency <- max(1 / orbit9)
  max_frequency_default <- is.null(max_frequency)
  if (max_frequency_default == TRUE){
    max_frequency <- 1.2 * highest_frequency
  }
  if (max_frequency < highest_frequency -
      64 * cvcoco_eps(max(1, highest_frequency))){
    stop(sprintf(paste0("max_frequency (%.12g cycle/kyr) is below the highest ",
                        "orbital frequency (%.12g cycle/kyr)"),
                 max_frequency, highest_frequency))
  }
  analysis_name <- "Blocked cvCOCO"
  fold_labels <- c("Segment A", "Segment B")
  if (is_interleaved == TRUE){
    analysis_name <- "Interleaved cvCOCO"
    fold_labels <- c("Odd", "Even")
  }
  group_index <- coco_orbit_groups(orbit9)
  group_names <- c("Long eccentricity", "Short eccentricity", "Obliquity",
                   "Precession")

  n_rate_estimate <- floor((sr2 - sr1) / srstep) + 1
  if (!is.finite(n_rate_estimate) || n_rate_estimate > 10000){
    stop("The sedimentation-rate grid is too large; increase srstep")
  }
  sr_grid <- cvcoco_colon(sr1, srstep, sr2)
  n_rate <- length(sr_grid)

  # --------------------------------------------------------- preprocessing
  prep <- cvcoco_prepare_data(data = data, method = method,
                              interleaved_phase = interleaved_phase,
                              analysis_name = analysis_name,
                              verbose = verbose)
  data_A <- prep$A
  data_B <- prep$B
  detrended_A <- as.numeric(detrend_linear(data_A[, 2]))
  detrended_B <- as.numeric(detrend_linear(data_B[, 2]))
  var_A <- cvcoco_resolved_variance(data_A[, 2], detrended_A)
  var_B <- cvcoco_resolved_variance(data_B[, 2], detrended_B)
  if (var_A$resolved == FALSE || var_B$resolved == FALSE){
    stop("Each held-out fold must retain resolved variance after linear detrending")
  }
  rho_A <- cvcoco_rho(detrended_A)
  rho_B <- cvcoco_rho(detrended_B)
  std_full <- NA
  if (is_interleaved == TRUE){
    detrended_full <- as.numeric(detrend_linear(prep$clean[, 2]))
    var_full <- cvcoco_resolved_variance(prep$clean[, 2], detrended_full)
    if (var_full$resolved == FALSE){
      stop("The complete cleaned series must retain resolved variance after detrending")
    }
    std_full <- var_full$std
    rho_full <- cvcoco_rho(detrended_full)
  } else {
    rho_full <- mean(c(rho_A, rho_B))
  }

  nf_A <- floor(max(pad, nrow(data_A)) / 2) + 1
  nf_B <- floor(max(pad, nrow(data_B)) / 2) + 1
  per_bin <- 4
  if (is_coherent == TRUE){
    per_bin <- 22
  }
  geometry_bytes <- 8 * per_bin * n_rate * (nf_A + nf_B)
  if (geometry_bytes > 512 * 1024^2){
    stop("The requested pad/rate grid needs too much memory; reduce pad or the rate grid")
  }

  geom_A <- cvcoco_geometry(depth = data_A[, 1], periods = orbit9,
                            sr_grid = sr_grid, pad = pad,
                            group_index = group_index,
                            target_model = target_model,
                            max_frequency = max_frequency)
  geom_B <- cvcoco_geometry(depth = data_B[, 1], periods = orbit9,
                            sr_grid = sr_grid, pad = pad,
                            group_index = group_index,
                            target_model = target_model,
                            max_frequency = max_frequency)
  if (!any(geom_A$valid_rate_mask) || !any(geom_B$valid_rate_mask)){
    stop("No tested sedimentation rate has usable frequency geometry in both folds")
  }

  degraded_mode <- FALSE
  completeness_A <- "complete-nine"
  completeness_B <- "complete-nine"
  warning_message <- ""
  strict_A <- geom_A$strict_training_rate_mask
  strict_B <- geom_B$strict_training_rate_mask
  if (!any(strict_A) || !any(strict_B)){
    if (any(geom_A$partial_training_rate_mask) &&
        any(geom_B$partial_training_rate_mask)){
      degraded_mode <- TRUE
      if (!any(strict_A)){
        geom_A$training_rate_mask <- geom_A$partial_training_rate_mask
        completeness_A <- "partial-orbit"
      }
      if (!any(strict_B)){
        geom_B$training_rate_mask <- geom_B$partial_training_rate_mask
        completeness_B <- "partial-orbit"
      }
      warning_message <- sprintf(paste0(
        "%s has no tested rate that supports complete all-nine training in ",
        "at least one fold. Continuing with a partial-orbit four-group target ",
        "(unresolved group weights fixed to zero); the result is exploratory. ",
        "Theoretical all-nine ranges (cm/kyr): %s (%.6g, %.6g], %s (%.6g, %.6g]."),
        analysis_name, fold_labels[1], geom_A$all_nine_rate_range[1],
        geom_A$all_nine_rate_range[2], fold_labels[2],
        geom_B$all_nine_rate_range[1], geom_B$all_nine_rate_range[2])
      warning(warning_message)
    } else {
      stop(sprintf(paste0(
        "%s must train the complete four-group target at a tested rate where ",
        "all nine periods are resolvable; no such rate exists in at least one ",
        "fold. Theoretical all-nine ranges (cm/kyr): %s (%.6g, %.6g], ",
        "%s (%.6g, %.6g]."), analysis_name, fold_labels[1],
        geom_A$all_nine_rate_range[1], geom_A$all_nine_rate_range[2],
        fold_labels[2], geom_B$all_nine_rate_range[1],
        geom_B$all_nine_rate_range[2]))
    }
  }

  # ------------------------------------------------------ observed pipeline
  spec_A <- cvcoco_spectra(values = data_A[, 2], geom = geom_A, red = red)
  spec_B <- cvcoco_spectra(values = data_B[, 2], geom = geom_B, red = red)
  power_A <- spec_A$power
  power_B <- spec_B$power

  tr_A <- cvcoco_train(power = power_A, geom = geom_A,
                       cor_method = cor_method, keep_curve = TRUE)
  tr_B <- cvcoco_train(power = power_B, geom = geom_B,
                       cor_method = cor_method, keep_curve = TRUE)
  grp_A <- cvcoco_group_amplitudes(tr_A$amplitude, group_index)
  grp_B <- cvcoco_group_amplitudes(tr_B$amplitude, group_index)
  frozen_A <- grp_A$normalized
  frozen_B <- grp_B$normalized
  if (any(!is.finite(frozen_A)) || any(!is.finite(frozen_B))){
    stop("At least one training fold produced nonfinite frozen target weights")
  }
  va_AtoB <- cvcoco_validate(power = power_B, geom = geom_B,
                             weights = frozen_A, cor_method = cor_method)
  va_BtoA <- cvcoco_validate(power = power_A, geom = geom_A,
                             weights = frozen_B, cor_method = cor_method)
  curve_AtoB <- as.numeric(va_AtoB$curve)
  curve_BtoA <- as.numeric(va_BtoA$curve)
  score_AtoB <- va_AtoB$best_score[1]
  score_BtoA <- va_BtoA$best_score[1]

  curve_consensus <- cvcoco_min_pair(curve_AtoB, curve_BtoA)
  cons_max <- cvcoco_curve_max(curve_consensus, sr_grid)
  score_consensus <- cons_max$score
  score_symmetric <- cvcoco_min_pair(score_AtoB, score_BtoA)
  score_mean <- NA
  if (is.finite(score_AtoB) && is.finite(score_BtoA)){
    score_mean <- (score_AtoB + score_BtoA) / 2
  }

  # ------------------------------------------------------------ Monte Carlo
  null_symmetric <- rep(NA_real_, n_sim)
  null_consensus <- rep(NA_real_, n_sim)
  null_AtoB <- rep(NA_real_, n_sim)
  null_BtoA <- rep(NA_real_, n_sim)
  null_rate_AtoB <- rep(NA_real_, n_sim)
  null_rate_BtoA <- rep(NA_real_, n_sim)
  exceed_AtoB <- numeric(n_rate)
  exceed_BtoA <- numeric(n_rate)
  exceed_cons <- numeric(n_rate)
  valid_AtoB <- numeric(n_rate)
  valid_BtoA <- numeric(n_rate)
  valid_cons <- numeric(n_rate)
  n_sim_completed <- 0
  actual_batch <- 0
  if (n_sim > 0){
    if (!is.null(seed)){
      set.seed(seed)
    }
    n_A <- nrow(data_A)
    n_B <- nrow(data_B)
    if (is_interleaved == TRUE){
      n_total <- nrow(prep$clean)
      raw_depth_A <- prep$clean[prep$index_A, 1]
      raw_depth_B <- prep$clean[prep$index_B, 1]
      workspace <- n_total + n_A + n_B
      innovation_std_full <- sqrt(max(0, 1 - rho_full^2))
    } else {
      n_total <- n_A + n_B
      workspace <- n_total
      innovation_std_A <- sqrt(max(0, 1 - rho_A^2))
      innovation_std_B <- sqrt(max(0, 1 - rho_B^2))
    }
    bytes_per_sim <- 64 * max(1, workspace + nf_A + nf_B + 2 * n_rate)
    memory_batch <- max(1, floor(256 * 1024^2 / bytes_per_sim))
    actual_batch <- min(c(batch_size, n_sim, memory_batch))
    if (verbose == TRUE){
      pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
    }
    first <- 1
    while (first <= n_sim){
      last <- min(first + actual_batch - 1, n_sim)
      n_batch <- last - first + 1
      innovations <- matrix(stats::rnorm(n_total * n_batch), n_total, n_batch)
      if (is_interleaved == TRUE){
        if (n_total > 1){
          innovations[2:n_total, ] <- innovation_std_full *
            innovations[2:n_total, , drop = FALSE]
        }
        surrogate_full <- std_full * cvcoco_ar1_filter(innovations, rho_full)
        surrogate_A <- matrix(0, n_A, n_batch)
        surrogate_B <- matrix(0, n_B, n_batch)
        for (k in seq_len(n_batch)){
          surrogate_A[, k] <- stats::approx(raw_depth_A,
                                            surrogate_full[prep$index_A, k],
                                            xout = data_A[, 1],
                                            method = "linear", rule = 1)$y
          surrogate_B[, k] <- stats::approx(raw_depth_B,
                                            surrogate_full[prep$index_B, k],
                                            xout = data_B[, 1],
                                            method = "linear", rule = 1)$y
        }
        if (any(!is.finite(surrogate_A)) || any(!is.finite(surrogate_B))){
          stop("The odd/even interpolation map produced nonfinite Monte Carlo values")
        }
      } else {
        innov_A <- innovations[seq_len(n_A), , drop = FALSE]
        innov_B <- innovations[n_A + seq_len(n_B), , drop = FALSE]
        if (n_A > 1){
          innov_A[2:n_A, ] <- innovation_std_A * innov_A[2:n_A, , drop = FALSE]
        }
        if (n_B > 1){
          innov_B[2:n_B, ] <- innovation_std_B * innov_B[2:n_B, , drop = FALSE]
        }
        surrogate_A <- var_A$std * cvcoco_ar1_filter(innov_A, rho_A)
        surrogate_B <- var_B$std * cvcoco_ar1_filter(innov_B, rho_B)
      }
      p_A <- cvcoco_spectra(values = surrogate_A, geom = geom_A, red = red)$power
      p_B <- cvcoco_spectra(values = surrogate_B, geom = geom_B, red = red)$power

      sim_tr_A <- cvcoco_train(power = p_A, geom = geom_A,
                               cor_method = cor_method, keep_curve = FALSE)
      sim_frozen_A <- cvcoco_group_amplitudes(sim_tr_A$amplitude,
                                              group_index)$normalized
      sim_AtoB <- cvcoco_validate(power = p_B, geom = geom_B,
                                  weights = sim_frozen_A,
                                  cor_method = cor_method)
      acc <- cvcoco_accumulate(sim_AtoB$curve, curve_AtoB, exceed_AtoB,
                               valid_AtoB)
      exceed_AtoB <- acc$exceed
      valid_AtoB <- acc$valid

      sim_tr_B <- cvcoco_train(power = p_B, geom = geom_B,
                               cor_method = cor_method, keep_curve = FALSE)
      sim_frozen_B <- cvcoco_group_amplitudes(sim_tr_B$amplitude,
                                              group_index)$normalized
      sim_BtoA <- cvcoco_validate(power = p_A, geom = geom_A,
                                  weights = sim_frozen_B,
                                  cor_method = cor_method)
      acc <- cvcoco_accumulate(sim_BtoA$curve, curve_BtoA, exceed_BtoA,
                               valid_BtoA)
      exceed_BtoA <- acc$exceed
      valid_BtoA <- acc$valid

      sim_cons <- cvcoco_min_pair(sim_AtoB$curve, sim_BtoA$curve)
      acc <- cvcoco_accumulate(sim_cons, curve_consensus, exceed_cons,
                               valid_cons)
      exceed_cons <- acc$exceed
      valid_cons <- acc$valid
      for (k in seq_len(n_batch)){
        col_values <- sim_cons[, k]
        col_values <- col_values[is.finite(col_values)]
        if (length(col_values) > 0){
          null_consensus[first + k - 1] <- max(col_values)
        }
      }
      null_AtoB[first:last] <- sim_AtoB$best_score
      null_BtoA[first:last] <- sim_BtoA$best_score
      null_rate_AtoB[first:last] <- sim_AtoB$best_rate
      null_rate_BtoA[first:last] <- sim_BtoA$best_rate
      null_symmetric[first:last] <- cvcoco_min_pair(sim_AtoB$best_score,
                                                    sim_BtoA$best_score)
      n_sim_completed <- last
      if (verbose == TRUE){
        utils::setTxtProgressBar(pb, last)
      }
      first <- last + 1
    }
    if (verbose == TRUE){
      close(pb)
    }
  }

  # ------------------------------------------------------------- p-values
  mc_sym <- cvcoco_mc_p(null_symmetric, score_symmetric)
  mc_cons <- cvcoco_mc_p(null_consensus, score_consensus)
  mc_AtoB <- cvcoco_mc_p(null_AtoB, score_AtoB)
  mc_BtoA <- cvcoco_mc_p(null_BtoA, score_BtoA)
  if (n_sim_completed > 0 &&
      (mc_sym$n_valid != n_sim_completed || mc_cons$n_valid != n_sim_completed ||
       mc_AtoB$n_valid != n_sim_completed || mc_BtoA$n_valid != n_sim_completed)){
    stop("At least one Monte Carlo realization produced a nonfinite statistic")
  }
  cvcoco_check_counts(valid_AtoB, curve_AtoB, n_sim_completed, "A-to-B")
  cvcoco_check_counts(valid_BtoA, curve_BtoA, n_sim_completed, "B-to-A")
  cvcoco_check_counts(valid_cons, curve_consensus, n_sim_completed, "consensus")
  p_curve_AtoB <- cvcoco_max_stat_p_curve(null_AtoB, curve_AtoB)
  p_curve_BtoA <- cvcoco_max_stat_p_curve(null_BtoA, curve_BtoA)
  p_curve_cons <- cvcoco_max_stat_p_curve(null_consensus, curve_consensus)
  p_local_AtoB <- cvcoco_pointwise_p(exceed_AtoB, valid_AtoB, curve_AtoB)
  p_local_BtoA <- cvcoco_pointwise_p(exceed_BtoA, valid_BtoA, curve_BtoA)
  p_local_cons <- cvcoco_pointwise_p(exceed_cons, valid_cons, curve_consensus)
  pcoco <- cvcoco_pcoco(curve_consensus, p_curve_cons)
  pcoco_max <- cvcoco_curve_max(pcoco, sr_grid)

  # --------------------------------------------------------- result lists
  train_A <- cvcoco_train_result(tr_A, grp_A, sr_grid, geom_A)
  train_B <- cvcoco_train_result(tr_B, grp_B, sr_grid, geom_B)
  validate_AtoB <- cvcoco_validation_result(curve_AtoB, score_AtoB,
                                            va_AtoB$best_rate[1], sr_grid,
                                            geom_B, frozen_A)
  validate_BtoA <- cvcoco_validation_result(curve_BtoA, score_BtoA,
                                            va_BtoA$best_rate[1], sr_grid,
                                            geom_A, frozen_B)
  validate_AtoB$p_directional <- mc_AtoB$p
  validate_AtoB$p_confidence_interval <- mc_AtoB$ci
  validate_AtoB$p_global_curve <- p_curve_AtoB
  validate_AtoB$p_local_curve <- p_local_AtoB
  validate_AtoB$p_local_at_best <- cvcoco_value_at(p_local_AtoB,
                                                   validate_AtoB$best_index)
  validate_BtoA$p_directional <- mc_BtoA$p
  validate_BtoA$p_confidence_interval <- mc_BtoA$ci
  validate_BtoA$p_global_curve <- p_curve_BtoA
  validate_BtoA$p_local_curve <- p_local_BtoA
  validate_BtoA$p_local_at_best <- cvcoco_value_at(p_local_BtoA,
                                                   validate_BtoA$best_index)
  active_AtoB <- cvcoco_active_counts(geom_B, frozen_A)
  active_BtoA <- cvcoco_active_counts(geom_A, frozen_B)
  active_consensus <- rowSums(active_AtoB$mask & active_BtoA$mask)

  consensus <- list(curve = curve_consensus,
                    best_rate = cons_max$rate,
                    best_index = cons_max$index,
                    best_correlation = score_consensus,
                    valid_rate_mask = geom_A$valid_rate_mask &
                      geom_B$valid_rate_mask,
                    p_global = mc_cons$p,
                    p_confidence_interval = mc_cons$ci,
                    p_global_curve = p_curve_cons,
                    p_local_curve = p_local_cons,
                    p_global_at_best = cvcoco_value_at(p_curve_cons,
                                                       cons_max$index),
                    p_local_at_best = cvcoco_value_at(p_local_cons,
                                                      cons_max$index),
                    p_directional_max = cvcoco_max_pair(mc_AtoB$p, mc_BtoA$p),
                    p_directional_max_curve = cvcoco_max_pair(p_curve_AtoB,
                                                              p_curve_BtoA),
                    active_orbit_count_curve = active_consensus,
                    active_orbit_count_at_best =
                      cvcoco_value_at(active_consensus, cons_max$index),
                    null_maximum = null_consensus[seq_len(n_sim_completed)],
                    local_exceedance_count = exceed_cons,
                    local_valid_count = valid_cons,
                    pcoco = pcoco,
                    best_pcoco = pcoco_max$score,
                    best_pcoco_rate = pcoco_max$rate,
                    best_pcoco_index = pcoco_max$index)

  spectra <- list(
    train_A = cvcoco_spectrum_diagnostic(power_A, geom_A, tr_A$best_rate[1],
                                         tr_A$amplitude[, 1], group_index,
                                         "adaptive"),
    train_B = cvcoco_spectrum_diagnostic(power_B, geom_B, tr_B$best_rate[1],
                                         tr_B$amplitude[, 1], group_index,
                                         "adaptive"),
    validate_AtoB = cvcoco_spectrum_diagnostic(power_B, geom_B,
                                               va_AtoB$best_rate[1],
                                               frozen_A[, 1], group_index,
                                               "fixed"),
    validate_BtoA = cvcoco_spectrum_diagnostic(power_A, geom_A,
                                               va_BtoA$best_rate[1],
                                               frozen_B[, 1], group_index,
                                               "fixed"))

  if (is_interleaved == TRUE){
    split_depth <- NA
  } else {
    split_depth <- prep$split_depth
  }
  keep <- seq_len(n_sim_completed)
  result <- list(
    name = analysis_name,
    method = method,
    target_model = target_model,
    degraded_mode = degraded_mode,
    training_completeness = ifelse(degraded_mode, "partial-orbit",
                                   "complete-nine"),
    training_completeness_A = completeness_A,
    training_completeness_B = completeness_B,
    warning_message = warning_message,
    sr_grid = sr_grid,
    split_depth = split_depth,
    interleaved_phase = interleaved_phase,
    fold_labels = fold_labels,
    data_clean = prep$clean,
    raw_index_A = prep$index_A,
    raw_index_B = prep$index_B,
    data_A = data_A,
    data_B = data_B,
    sampling_interval_A = prep$dz_A,
    sampling_interval_B = prep$dz_B,
    interpolation_A = prep$info_A,
    interpolation_B = prep$info_B,
    train_A = train_A,
    train_B = train_B,
    validate_AtoB = validate_AtoB,
    validate_BtoA = validate_BtoA,
    consensus = consensus,
    spectra = spectra,
    score_symmetric = score_symmetric,
    score_consensus = score_consensus,
    score_mean = score_mean,
    p_sym = mc_sym$p,
    p_sym_confidence_interval = mc_sym$ci,
    p_consensus = mc_cons$p,
    p_consensus_confidence_interval = mc_cons$ci,
    p_robust = cvcoco_max_pair(mc_AtoB$p, mc_BtoA$p),
    p_AtoB = mc_AtoB$p,
    p_BtoA = mc_BtoA$p,
    p_A = mc_BtoA$p,
    p_B = mc_AtoB$p,
    p_A_confidence_interval = mc_BtoA$ci,
    p_B_confidence_interval = mc_AtoB$ci,
    p_curve_AtoB = p_curve_AtoB,
    p_curve_BtoA = p_curve_BtoA,
    p_local_curve_AtoB = p_local_AtoB,
    p_local_curve_BtoA = p_local_BtoA,
    p_curve_consensus = p_curve_cons,
    p_local_curve_consensus = p_local_cons,
    pcoco = pcoco,
    best_pcoco = pcoco_max$score,
    best_pcoco_rate = pcoco_max$rate,
    null_symmetric = null_symmetric[keep],
    null_consensus = null_consensus[keep],
    null_AtoB = null_AtoB[keep],
    null_BtoA = null_BtoA[keep],
    null_best_rate_AtoB = null_rate_AtoB[keep],
    null_best_rate_BtoA = null_rate_BtoA[keep],
    local_exceedance_count_AtoB = exceed_AtoB,
    local_exceedance_count_BtoA = exceed_BtoA,
    local_valid_count_AtoB = valid_AtoB,
    local_valid_count_BtoA = valid_BtoA,
    rho_A = rho_A,
    rho_B = rho_B,
    rho_full = rho_full,
    rho_used_in_null = ifelse(is_interleaved, "rho_full", "rho_A and rho_B"),
    n_sim_requested = n_sim,
    n_sim_completed = n_sim_completed,
    n_sim_valid = mc_sym$n_valid,
    n_sim_valid_consensus = mc_cons$n_valid,
    n_sim_valid_AtoB = mc_AtoB$n_valid,
    n_sim_valid_BtoA = mc_BtoA$n_valid,
    seed = seed,
    orbit_periods = orbit9,
    group_names = group_names,
    group_index = group_index,
    valid_rate_mask_A = geom_A$valid_rate_mask,
    valid_rate_mask_B = geom_B$valid_rate_mask,
    training_rate_mask_A = geom_A$training_rate_mask,
    training_rate_mask_B = geom_B$training_rate_mask,
    strict_training_rate_mask_A = strict_A,
    strict_training_rate_mask_B = strict_B,
    partial_training_rate_mask_A = geom_A$partial_training_rate_mask,
    partial_training_rate_mask_B = geom_B$partial_training_rate_mask,
    training_active_group_mask_A = geom_A$training_active_group_mask,
    training_active_group_mask_B = geom_B$training_active_group_mask,
    orbit_count_A = geom_A$orbit_count,
    orbit_count_B = geom_B$orbit_count,
    resolved_group_count_A = geom_A$resolved_group_count,
    resolved_group_count_B = geom_B$resolved_group_count,
    resolved_orbit_mask_A = geom_A$resolved_orbit_mask,
    resolved_orbit_mask_B = geom_B$resolved_orbit_mask,
    active_orbit_count_AtoB = active_AtoB$count,
    active_orbit_count_BtoA = active_BtoA$count,
    active_group_count_AtoB = active_AtoB$group_counts,
    active_group_count_BtoA = active_BtoA$group_counts,
    active_orbit_count_consensus = active_consensus,
    cross_group_band_overlap_A = geom_A$cross_group_band_overlap,
    cross_group_band_overlap_B = geom_B$cross_group_band_overlap,
    group_leakage_rcond_A = geom_A$partial_group_leakage_rcond,
    group_leakage_rcond_B = geom_B$partial_group_leakage_rcond,
    full_group_leakage_rcond_A = geom_A$group_leakage_rcond,
    full_group_leakage_rcond_B = geom_B$group_leakage_rcond,
    all_nine_rate_range_A = geom_A$all_nine_rate_range,
    all_nine_rate_range_B = geom_B$all_nine_rate_range,
    all_nine_rate_range_shared = c(max(geom_A$all_nine_rate_range[1],
                                       geom_B$all_nine_rate_range[1]),
                                   min(geom_A$all_nine_rate_range[2],
                                       geom_B$all_nine_rate_range[2])),
    config = list(pad = pad, nfft_A = geom_A$nfft, nfft_B = geom_B$nfft,
                  sr1 = sr1, sr2 = sr2, srstep = srstep, red = red,
                  n_sim = n_sim, cor_method = cor_method,
                  batch_size = actual_batch, seed = seed, method = method,
                  target_model = target_model,
                  max_frequency = max_frequency,
                  max_frequency_default = max_frequency_default,
                  minimum_leakage_rcond = geom_A$minimum_leakage_rcond,
                  active_group_weight_tolerance = 1e-8),
    conclusion = NULL)
  if (target_model == "four-group"){
    result$conclusion <- cvcoco_conclusion(result)
  }
  if (genplot == TRUE){
    plot_cvcoco(result)
  }
  result
}
