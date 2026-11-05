#' @title Blocked eCOCO (cross-fitted evolutionary correlation coefficient)
#'
#' @description Blocked eCOCO evaluates complete sliding validation windows
#' with astronomical targets trained on strictly non-overlapping anchor
#' windows. Four leakage-corrected orbital-group amplitudes (long
#' eccentricity, short eccentricity, obliquity, precession; groups identified
#' by descending period rank) are fitted at regularly spaced anchor windows
#' (spacing \code{anchor_fraction * window}) using the complete
#' sedimentation-rate grid, and expanded to a coherent nine-term sine target.
#' Every validation window is evaluated with the nearest lower-depth anchor
#' (forward) and the nearest higher-depth anchor (backward) with frozen
#' weights. The consensus correlation at a rate is the smaller of the two
#' directional correlations; at the record edges the available direction is
#' used (\code{consensus_policy = "supported"}) or the window is left empty
#' (\code{"strict"}). Stationary Gaussian AR(1) Monte Carlo nulls, one per
#' unique window interval, repeat the complete anchor training and frozen
#' validation. Global p-values use the null maximum over the complete rate
#' grid within each window, local p-values the same rate.
#' pCOCO = rho * |log10(p_global)| is the ridge score; the sedimentation-rate
#' ridge is tracked as in Acycle (\code{trackEcocoRidge}).
#' Translated from Acycle v3.1 ecocoCrossfitCore.m, ecoco.m (Blocked branch),
#' ecocoWindowSpectra.m, ecocoEstimateAr1.m, ecocoPhysicalCenters.m,
#' ecocoPhysicalWindowValues.m, cocoNonnegativeLeakageSolve.m and zeropad2.m
#' (Mingsong Li). Uses the COCO helper functions of this package
#' (helpers_coco.R and the ecoco_* functions).
#'
#' @param data Record to analyse (matrix or data frame, two columns): depth in metres in column 1 and the proxy in
#' column 2. Non-finite rows are dropped, duplicated depths averaged, and unevenly spaced data are linearly
#' interpolated to the median sample spacing (see \code{\link{coco_prepare_regular_data}}). Required
#' (\code{Default=NULL}).
#' @param orbit9 The nine astronomical target periods in kyr (numeric vector of nine distinct positive values: long
#' eccentricity, four short eccentricity terms, obliquity, three precession terms), e.g. from
#' \code{\link{orbit9_at_age}} for the age of the record.
#' \code{Default=c(405.6912, 130.6979, 123.8532, 98.8517, 94.8856, 40.9897, 23.6820, 22.3758, 18.9519)}.
#' @param window Length of the sliding validation and anchor windows, in metres (numeric, > 0; every window must
#' contain at least 8 samples). \code{NULL} uses 25\% of the record length (Acycle GUI default).
#' \code{Default=NULL}.
#' @param sr Sedimentation rates to test, in cm/kyr (numeric vector, > 0, strictly increasing). \code{NULL} uses the
#' Acycle GUI default: from \code{2 * dt * 0.065 * 100} up to 100 cm/kyr (reduced so that the longest period fits in
#' the record), in steps of 0.1 cm/kyr or at most 300 steps. \code{Default=NULL}.
#' @param step_depth Spacing of the output (validation) window centres, in metres (numeric, > 0); only used with
#' \code{window_mode="physical-depth"}. \code{NULL} uses the sample interval \code{dt}, or
#' \code{dt * ceiling(n / 300)} for more than 300 samples. \code{Default=NULL}.
#' @param step Step between successive validation windows, in samples (integer, >= 1); only used with
#' \code{window_mode="legacy-count"}. \code{Default=1}.
#' @param pad Length (number of points) to which every window is zero-padded for the periodogram (integer, at least
#' the number of points per window). \code{NULL} uses the Acycle default: 5000 for up to 2500 samples, 10000 for up
#' to 5000 samples, and \code{floor(n / 5000) * 5000 + 5000} above that. \code{Default=NULL}.
#' @param red Background (red-noise) removal applied to every window spectrum (integer): \code{0} none, \code{1}
#' classical AR(1) background, \code{2} robust AR(1) background, \code{3} smoothed-window-average (SWA) background.
#' Negative power after subtraction is set to zero. \code{Default=0}.
#' @param n_sim Number of AR(1) Monte Carlo simulations used for the global and local p-values (integer, >= 0); each
#' simulation repeats the anchor training and validation. \code{Default=2000}.
#' @param method Correlation coefficient between the window spectrum and the target spectrum (character):
#' \code{"Pearson"} or \code{"Spearman"}. \code{Default="Pearson"}.
#' @param max_frequency Highest frequency of the target spectrum that enters the correlation, in cycles/kyr
#' (numeric, at least \code{max(1/orbit9)}); \code{NULL} uses \code{1.2 * max(1/orbit9)}. \code{Default=NULL}.
#' @param seed Seed of the random number generator (integer), so that the Monte Carlo results can be reproduced;
#' \code{NULL} leaves the current random state unchanged. \code{Default=1}.
#' @param anchor_fraction Spacing of the anchor windows (where the target amplitudes are trained) as a fraction of
#' \code{window} (numeric, > 0 and <= 2); smaller values give more, overlapping anchors. \code{Default=0.5}.
#' @param window_mode How the windows are defined (character): \code{"physical-depth"} (windows of fixed length in
#' metres centred every \code{step_depth}, as in the Acycle wrapper and GUI) or \code{"legacy-count"} (windows of a
#' fixed number of samples moved by \code{step} samples; needs evenly sampled data).
#' \code{Default="physical-depth"}.
#' @param edge_pad Zero-pad the record by half a window at both ends and restrict the output centres to the original
#' record, so that windows reach the record ends (logical, TRUE or FALSE; Acycle GUI setting for Blocked eCOCO).
#' Only used with \code{window_mode="physical-depth"}. \code{Default=TRUE}.
#' @param center_limits Depths of the first and last output window centre, in metres (numeric vector
#' \code{c(first, last)}, or \code{NULL}); overrides the limits set by \code{edge_pad}. \code{NULL} uses the record
#' ends (with \code{edge_pad=TRUE}) or \code{c(first depth + window/2, last depth - window/2)}.
#' \code{Default=NULL}.
#' @param compute_local_p Also compute the local p-values, comparing each rate with the null correlations at the same
#' rate (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param consensus_policy How windows at the record edges, which have an anchor on only one side, are treated
#' (character): \code{"supported"} uses the available direction, \code{"strict"} leaves windows without both a
#' forward and a backward anchor empty. \code{Default="supported"}.
#' @param batch_size Number of Monte Carlo simulations processed together (integer, >= 1); larger batches are faster
#' but use more memory, the results do not depend on it. \code{Default=20}.
#' @param genplot Plot the result with \code{\link{plot_ecoco_cv}} (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param verbose Print progress messages (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' @references Li, M., Kump, L.R., Hinnov, L.A., Mann, M.E., 2018. Tracking
#' variable sedimentation rates and astronomical forcing in Phanerozoic
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
#' set.seed(2)
#' dt <- 0.05
#' depth <- seq(0, 30, by = dt)
#' tt <- depth * 100 / 5          # 5 cm/kyr
#' x <- sin(2 * pi * tt / 405.7) + 0.8 * sin(2 * pi * tt / 98.9) +
#'   0.5 * sin(2 * pi * tt / 40.99) + 0.4 * sin(2 * pi * tt / 23.68) +
#'   0.3 * rnorm(length(tt))
#' res <- ecoco_blocked(cbind(depth, x), window = 10, sr = seq(3, 8, 0.5),
#'                      step_depth = 1, pad = 400, n_sim = 20)
#' res$track[, c("depth", "sr")]
#' }
#'
#' @return A list (class "ecoco_cv") with elements \code{method},
#' \code{sr} (rate grid), \code{depth} (window centres), rate-by-window
#' matrices \code{rho}, \code{p_global}, \code{p_local}, \code{norbit},
#' \code{pcoco} and \code{score} (the published consensus surface), the lists
#' \code{forward}, \code{backward}, \code{consensus}, \code{strict}
#' (each with rho, p_global, p_local, norbit, pcoco, best_rate, ...),
#' \code{support} (per-window "bidirectional", "forward-only",
#' "backward-only" or "none"), \code{windows} and \code{anchors}
#' (data frames with bounds, anchor links and trained rates/weights),
#' \code{geometry} (rate masks), \code{rho_m} (per-window AR(1)),
#' \code{best} (best rate per window) and \code{track} (tracked ridge:
#' depth, sr, rho, p_global, norbit, score, sr_low, sr_high).
#' @param centre_limits British spelling of \code{center_limits}; if given (not \code{NULL}) it overrides
#' \code{center_limits}. \code{Default=NULL}.
#' @export
#' @importFrom stats rnorm median
#' @importFrom utils txtProgressBar setTxtProgressBar
ecoco_blocked <- function(data = NULL,
    orbit9 = c(405.6912, 130.6979, 123.8532, 98.8517,
        94.8856, 40.9897, 23.6820, 22.3758,
        18.9519),
    window = NULL, sr = NULL, step_depth = NULL,
    step = 1, pad = NULL, red = 0, n_sim = 2000,
    method = "Pearson", max_frequency = NULL, seed = 1,
    anchor_fraction = 0.5,
    window_mode = "physical-depth", edge_pad = TRUE,
    center_limits = NULL, compute_local_p = TRUE,
    consensus_policy = "supported", batch_size = 20,
    genplot = FALSE, verbose = FALSE,
    centre_limits = NULL){
  if (!is.null(centre_limits)) {
    center_limits <- centre_limits
  }  # British spelling synonym
  # Internal helper functions (local to this function)

  # MATLAB-compatible scalar helpers

  ecocox_round <- function(x){
    # MATLAB round: halves away from zero (R round() rounds half to even)
    sign(x) * floor(abs(x) + 0.5)
  }

  ecocox_eps <- function(x = 1){
    # MATLAB eps(x): spacing of doubles at |x|
    x <- abs(x)
    if (!is.finite(x)){
      return(NaN)
    }
    if (x < 2.2250738585072014e-308){
      return(2^-1074)
    }
    2^(floor(log2(x)) - 52)
  }

  ecocox_max_finite <- function(x){
    # MATLAB max() ignoring NaN; NaN when nothing is finite
    x <- x[is.finite(x)]
    if (length(x) < 1){
      return(NaN)
    }
    max(x)
  }

  ecocox_col_max_finite <- function(M){
    M <- as.matrix(M)
    out <- rep(NaN, ncol(M))
    for (k in seq_len(ncol(M))){
      out[k] <- ecocox_max_finite(M[, k])
    }
    out
  }

  ecocox_first_max <- function(v){
    # index of the first finite maximum (MATLAB max over finite entries)
    ok <- which(is.finite(v))
    if (length(ok) < 1){
      return(NA_integer_)
    }
    ok[which.max(v[ok])]
  }

  ecocox_sampling_tolerance <- function(spacing){
    # cocoSamplingTolerance.m
    1e-5 * abs(spacing)
  }

  ecocox_all_period_rate_range <- function(periods, dz, n){
    # cocoAllPeriodRateRange.m
    c(max(200 * dz / periods), min(100 * n * dz / periods))
  }

  ecocox_default_pad <- function(npts){
    # defaultPad in eCOCOGUI.m
    if (npts <= 2500){
      return(5000)
    } else if (npts <= 5000){
      return(10000)
    }
    floor(npts / 5000) * 5000 + 5000
  }

  ecocox_default_sr <- function(n, dt, orbit9){
    # defaultSedRange in eCOCOGUI.m
    sedmin <- 0
    sedmax <- 100
    sedstep <- 0.1
    fh <- 0.065
    fnyq <- sedmin / (2 * dt)
    if (fh > fnyq){
      sedmin <- 2 * dt * fh * 100
    }
    fray <- sedmax / (n * dt)
    flow <- 1 / max(orbit9)
    if (fray > flow){
      sedmax <- n * dt * flow * 100
    }
    if ((sedmax - sedmin) / sedstep > 300){
      sedstep <- (sedmax - sedmin) / 300
    }
    seq(sedmin, sedmax, by = sedstep)
  }

  # Variance / AR(1) helpers

  ecocox_resolved_variance <- function(original, detrended){
    # cocoResolvedDetrendedVariance.m
    original <- as.numeric(original)
    detrended <- as.numeric(detrended)
    n <- length(detrended)
    out <- list(resolved = FALSE, sd = 0)
    ref <- max(abs(original))
    rs <- max(abs(detrended))
    if (n < 2 || ref <= 0 || rs <= 0){
      return(out)
    }
    scaled <- detrended / rs
    energy <- sum(scaled^2)
    if (!is.finite(energy) || energy <= 0){
      return(out)
    }
    rel <- (rs / ref) * sqrt(energy / n)
    s <- rs * sqrt(energy / (n - 1))
    out$sd <- s
    out$resolved <- is.finite(rel) && rel > 256 * .Machine$double.eps &&
        is.finite(s) && s > 0
    out
  }

  ecocox_estimate_ar1 <- function(values){
    # ecocoEstimateAr1.m: column-wise conditional LS AR(1) of detrended data
    values <- as.matrix(values)
    ns <- ncol(values)
    rho <- rep(NaN, ns)
    sdv <- rep(NaN, ns)
    valid <- rep(FALSE, ns)
    for (k in seq_len(ns)){
      orig <- values[, k]
      if (length(orig) < 4 || any(!is.finite(orig))){
        next
      }
      x <- detrend_linear(orig)
      rv <- ecocox_resolved_variance(orig, x)
      if (!rv$resolved || !is.finite(rv$sd) || rv$sd <= 0){
        next
      }
      sc <- max(abs(x))
      if (!is.finite(sc) || sc <= 0){
        next
      }
      z <- x / sc
      a <- z[-length(z)]
      b <- z[-1]
      den <- sum(a * a)
      if (!is.finite(den) || den <= 0){
        next
      }
      est <- sum(a * b) / den
      if (!is.finite(est)){
        next
      }
      rho[k] <- min(max(est, -0.999), 0.999)
      sdv[k] <- rv$sd
      valid[k] <- TRUE
    }
    list(rho = rho, sd = sdv, valid = valid)
  }

  # Spectra with background removal (ecocoWindowSpectra.m, spectrumBatch)

  ecocox_window_spectra <- function(values, dt, nfft, red = 0){
    values <- as.matrix(values)
    if (any(!is.finite(values))){
      stop("Every window value must be finite before spectral analysis.")
    }
    if (nfft < nrow(values)){
      stop("pad must be at least the number of samples in a window.")
    }
    if (red %in% c(0, 1, 2)){
      sp <- coco_window_spectra(values, dt, nfft, red)
      return(list(power = as.matrix(sp$power), f = sp$f))
    }
    if (red != 3){
      stop("red must be 0, 1, 2 or 3")
    }
    if (!exists("swa_background", mode = "function")){
      stop("red = 3 requires swa_background() (R/swa_background.R)")
    }
    sp <- coco_window_spectra(values, dt, nfft, 0)
    P <- as.matrix(sp$power)
    f <- sp$f
    if (nrow(P) < 33){
      stop("Smoothed-window removal requires at least 33 one-sided bins.")
    }
    for (k in seq_len(ncol(P))){
      p <- P[, k]
      pos <- p[is.finite(p) & p > 0]
      if (length(pos) < 1){
        stop("SWA removal requires positive spectral power.")
      }
      fl <- max(2.2250738585072014e-308, max(pos) * .Machine$double.eps)
      bg <- swa_background(f, log10(pmax(p, fl)), nrow(values),
          verbose = FALSE)$swa
      bg <- as.numeric(bg)
      if (length(bg) != length(p) || any(!is.finite(bg) | bg <= 0)){
        stop("The selected background model returned invalid values.")
      }
      p <- p - bg
      p[p < 0] <- 0
      P[, k] <- p
    }
    list(power = P, f = f)
  }

  # Exact four-variable nonnegative least squares (cocoNonnegativeLeakageSolve.m)

  ecocox_nnls4 <- function(M, E){
    E <- as.matrix(E)
    ns <- ncol(E)
    cs <- rep(0, ns)
    for (k in seq_len(ns)){
      cs[k] <- max(E[, k])
    }
    cs[!is.finite(cs) | cs <= 0] <- 1
    S <- E / matrix(cs, 4, ns, byrow = TRUE)
    G <- matrix(0, 4, ns)
    best <- colSums(S^2)
    for (code in 1:15){
      act <- bitwAnd(code, c(1L, 2L, 4L, 8L)) > 0
      D <- M[, act, drop = FALSE]
      sv <- svd(D, nu = 0, nv = 0)$d
      tol <- max(dim(D)) * ecocox_eps(max(sv))
      if (sum(sv > tol) < sum(act)){
        next
      }
      if (ncol(D) == 4){
        cand <- solve(D, S)
      } else {
        cand <- qr.coef(qr(D, tol = 1e-300), S)
      }
      cand <- matrix(cand, nrow = sum(act))
      feasible <- rep(FALSE, ns)
      for (k in seq_len(ns)){
        ft <- 256 * .Machine$double.eps * max(1, max(abs(cand[, k])))
        feasible[k] <- all(cand[, k] >= -ft)
      }
      if (!any(feasible)){
        next
      }
      cand <- pmax(cand, 0)
      res <- colSums((D %*% cand - S)^2)
      better <- feasible & is.finite(res) & res < best
      if (any(better)){
        G[, better] <- 0
        G[act, better] <- cand[, better]
        best[better] <- res[better]
      }
    }
    G * matrix(cs, 4, ns, byrow = TRUE)
  }

  # Shared four-group coherent-nine geometry
  # (prepareGeometry in ecocoCrossfitCore.m and in cvcoco.m, four-group branch)

  ecocox_group_overlap <- function(center, rayleigh, upper, resolved, gindex){
    active <- which(resolved)
    if (length(active) < 2){
      return(FALSE)
    }
    lo <- pmax(0, center - rayleigh)
    hi <- pmin(upper, center + rayleigh)
    tol <- 128 * ecocox_eps(max(c(1, abs(lo), abs(hi), abs(upper))))
    for (a in 1:(length(active) - 1)){
      i1 <- active[a]
      for (b in (a + 1):length(active)){
        i2 <- active[b]
        if (gindex[i1] == gindex[i2]){
          next
        }
        if (min(hi[i1], hi[i2]) - max(lo[i1], lo[i2]) > tol){
          return(TRUE)
        }
      }
    }
    FALSE
  }

  ecocox_geometry <- function(n, dz, periods, sr_grid, nfft, max_frequency,
      engine = c("blocked", "cvcoco")){
    engine <- match.arg(engine)
    gindex <- coco_orbit_groups(periods)
    fs <- 1 / dz
    nf <- floor(nfft / 2) + 1
    frequency <- (0:(nf - 1)) * (fs / nfft)
    rayleigh <- fs / n
    nyquist <- fs / 2
    n_rate <- length(sr_grid)
    n_orbit <- length(periods)
    w1 <- rep(2, nf)
    w1[1] <- 1
    if (nfft %% 2 == 0){
      w1[nf] <- 1
    }
    geom <- list(n = n, dz = dz, fs = fs, nfft = nfft, frequency = frequency,
        rayleigh = rayleigh, nyquist = nyquist, sr_grid = sr_grid,
        group_index = gindex, engine = engine,
        one_sided = w1, spectrum_scale = 1 / (fs * n),
        basis = vector("list", n_rate),
        group_basis = vector("list", n_rate),
        band_index = vector("list", n_rate),
        solve_matrix = vector("list", n_rate),
        leakage = vector("list", n_rate),
        corr_index = vector("list", n_rate),
        group_unit_energy = matrix(NaN, n_rate, 4),
        full_rcond = rep(NaN, n_rate),
        active_rcond = rep(NaN, n_rate),
        valid_rate = rep(FALSE, n_rate),
        strict_training = rep(FALSE, n_rate),
        partial_training = rep(FALSE, n_rate),
        active_group = matrix(FALSE, n_rate, 4),
        active_orbit = matrix(FALSE, n_rate, n_orbit),
        resolved = matrix(FALSE, n_rate, n_orbit),
        orbit_count = rep(0, n_rate),
        resolved_group_count = matrix(0, n_rate, 4),
        overlap = rep(FALSE, n_rate),
        all_nine_range = ecocox_all_period_rate_range(periods, dz, n))
    pad_rows <- nfft - n
    for (m in seq_len(n_rate)){
      sr <- sr_grid[m]
      spatial <- 100 / (periods * sr)
      res <- spatial >= rayleigh & spatial < nyquist
      geom$orbit_count[m] <- sum(res)
      geom$resolved[m, ] <- res
      for (g in 1:4){
        geom$resolved_group_count[m, g] <- sum(res & gindex == g)
      }
      tfreq <- frequency * sr / 100
      cidx <- which(tfreq <= max_frequency + 16 * ecocox_eps(max(1, max_frequency)))
      geom$corr_index[[m]] <- cidx
      if (length(cidx) < 3 || !any(res)){
        next
      }
      lo <- rep(0, n_orbit)
      hi <- rep(0, n_orbit)
      for (j in seq_len(n_orbit)){
        inb <- which(abs(frequency - spatial[j]) <= rayleigh)
        if (length(inb) < 1){
          inb <- which.min(abs(frequency - spatial[j]))
        }
        lo[j] <- inb[1]
        hi[j] <- inb[length(inb)]
      }
      tt <- (0:(n - 1)) * (100 * dz / sr)
      ph <- matrix(2 * pi * tt, n, n_orbit) /
          matrix(periods, n, n_orbit, byrow = TRUE)
      us <- detrend_linear(sin(ph))
      uc <- detrend_linear(cos(ph))
      if (pad_rows > 0){
        us <- rbind(us, matrix(0, pad_rows, n_orbit))
        uc <- rbind(uc, matrix(0, pad_rows, n_orbit))
      }
      H <- mvfft(us)[1:nf, , drop = FALSE]
      Hc <- mvfft(uc)[1:nf, , drop = FALSE]
      if (any(!is.finite(Re(H))) || any(!is.finite(Im(H))) ||
          any(!is.finite(Re(Hc))) || any(!is.finite(Im(Hc)))){
        next
      }
      unit_power <- 0.5 * (Mod(H)^2 + Mod(Hc)^2) * w1 * (1 / (fs * n))
      gt <- matrix(0, nf, 4)
      upper <- min(nyquist, max_frequency * 100 / sr)
      for (j in seq_len(n_orbit)){
        if (res[j]){
          gt[, gindex[j]] <- gt[, gindex[j]] + unit_power[, j]
        }
      }
      bands <- vector("list", 4)
      for (g in 1:4){
        mem <- which(gindex == g & res)
        bm <- rep(FALSE, nf)
        for (mm in mem){
          bm[lo[mm]:hi[mm]] <- TRUE
        }
        bm <- bm & frequency <= upper + 64 * ecocox_eps(max(1, upper))
        bands[[g]] <- which(bm)
        if (length(bands[[g]]) > 0){
          geom$group_unit_energy[m, g] <- sum(gt[bands[[g]], g]) * (fs / nfft)
        }
      }
      geom$band_index[[m]] <- bands
      geom$overlap[m] <- ecocox_group_overlap(spatial, rayleigh, upper, res,
          gindex)
      L <- matrix(0, 4, 4)
      for (g in 1:4){
        if (length(bands[[g]]) > 0){
          L[g, ] <- colSums(gt[bands[[g]], , drop = FALSE]) * (fs / nfft)
        }
      }
      geom$leakage[[m]] <- L
      if (all(is.finite(L))){
        geom$full_rcond[m] <- rcond(L)
      }
      if (engine == "blocked"){
        GB <- matrix(0 + 0i, nf, 4)
        for (g in 1:4){
          mem <- gindex == g & res
          if (any(mem)){
            GB[, g] <- rowSums(H[, mem, drop = FALSE])
          }
        }
        geom$group_basis[[m]] <- GB
      } else {
        geom$basis[[m]] <- H
      }
      geom$valid_rate[m] <- any(gt > 0)
      ag <- geom$resolved_group_count[m, ] > 0 &
          is.finite(geom$group_unit_energy[m, ]) & geom$group_unit_energy[m, ] > 0
      ag[is.na(ag)] <- FALSE
      geom$active_group[m, ] <- ag
      geom$active_orbit[m, ] <- res & ag[gindex]
      Ms <- L
      for (g in which(!ag)){
        Ms[g, ] <- 0
        Ms[, g] <- 0
        Ms[g, g] <- 1
      }
      geom$solve_matrix[[m]] <- Ms
      if (any(ag)){
        A <- L[ag, ag, drop = FALSE]
        if (all(is.finite(A)) && all(diag(A) > 0)){
          geom$active_rcond[m] <- rcond(A)
        }
      }
      geom$partial_training[m] <- geom$valid_rate[m] && any(ag) &&
          is.finite(geom$active_rcond[m]) && geom$active_rcond[m] >= 1e-10 &&
          !geom$overlap[m]
      ue <- geom$group_unit_energy[m, ]
      geom$strict_training[m] <- all(res) && all(is.finite(ue) & ue > 0) &&
          is.finite(geom$full_rcond[m]) && geom$full_rcond[m] >= 1e-10 &&
          !geom$overlap[m]
    }
    geom$training_rate <- geom$strict_training
    geom
  }

  ecocox_target_power <- function(geom, m, amp){
    # fourGroupTargetPower (blocked: coherent group basis; cvcoco: expanded)
    amp <- as.matrix(amp)
    if (geom$engine == "blocked"){
      tf <- geom$group_basis[[m]] %*% amp
    } else {
      ex <- amp[geom$group_index, , drop = FALSE] *
          as.numeric(geom$resolved[m, ])
      tf <- geom$basis[[m]] %*% ex
    }
    Mod(tf)^2 * geom$one_sided * geom$spectrum_scale
  }

  ecocox_train <- function(power, geom, method, keep_curve = TRUE){
    # trainFourGroup (ecocoCrossfitCore.m) / fourGroupAdaptiveTrain (cvcoco.m)
    power <- as.matrix(power)
    n_rate <- length(geom$sr_grid)
    ns <- ncol(power)
    curve <- NULL
    if (keep_curve){
      curve <- matrix(NaN, n_rate, ns)
    }
    best_score <- rep(-Inf, ns)
    best_rate <- rep(NaN, ns)
    amp_best <- matrix(NaN, 4, ns)
    df <- geom$fs / geom$nfft
    tol2 <- (1e-8)^2
    for (m in seq_len(n_rate)){
      if (!geom$training_rate[m]){
        next
      }
      ag <- geom$active_group[m, ]
      E <- matrix(0, 4, ns)
      for (g in 1:4){
        if (!ag[g]){
          next
        }
        idx <- geom$band_index[[m]][[g]]
        if (length(idx) < 1){
          stop("An eligible partial-orbit training group has no frequency band.")
        }
        E[g, ] <- colSums(power[idx, , drop = FALSE]) * df
      }
      if (any(!is.finite(E)) || any(E < 0)){
        stop("A four-group integration band has invalid energy.")
      }
      E[!ag, ] <- 0
      GP <- ecocox_nnls4(geom$solve_matrix[[m]], E)
      GP[!ag, ] <- 0
      gmax <- rep(0, ns)
      for (k in seq_len(ns)){
        gmax[k] <- max(GP[, k])
      }
      GP[GP <= matrix(gmax * tol2, 4, ns, byrow = TRUE)] <- 0
      GA <- sqrt(GP)
      tp <- ecocox_target_power(geom, m, GA)
      idx <- geom$corr_index[[m]]
      r <- .column_correlation(power[idx, , drop = FALSE],
          tp[idx, , drop = FALSE], method)
      if (keep_curve){
        curve[m, ] <- r
      }
      better <- is.finite(r) & r > best_score
      if (any(better)){
        best_score[better] <- r[better]
        best_rate[better] <- geom$sr_grid[m]
        amp_best[, better] <- GA[, better]
      }
    }
    best_score[!is.finite(best_rate)] <- NaN
    list(curve = curve, best_rate = best_rate, best_score = best_score,
        amplitude = amp_best)
  }

  ecocox_normalize <- function(amp){
    # normalizeGroupAmplitudes (crossfit) / normalized part of groupAmplitudes
    amp <- as.matrix(amp)
    sc <- ecocox_col_max_finite(amp)
    W <- amp / matrix(sc, nrow(amp), ncol(amp), byrow = TRUE)
    for (k in seq_len(ncol(amp))){
      if (is.finite(sc[k]) && sc[k] == 0 && all(is.finite(amp[, k]))){
        W[, k] <- 0
      }
      if (!is.finite(sc[k]) || sc[k] < 0){
        W[, k] <- NaN
      }
    }
    small <- !is.na(W) & W >= 0 & W <= 1e-8
    W[small] <- 0
    W
  }

  ecocox_validate <- function(power, geom, weights, method){
    # validateFrozen (crossfit) / fixedValidate curve (cvcoco)
    power <- as.matrix(power)
    weights <- as.matrix(weights)
    n_rate <- length(geom$sr_grid)
    ns <- ncol(power)
    if (ncol(weights) == 1 && ns > 1){
      weights <- matrix(weights, 4, ns)
    }
    curve <- matrix(NaN, n_rate, ns)
    for (m in seq_len(n_rate)){
      if (!geom$valid_rate[m]){
        next
      }
      tp <- ecocox_target_power(geom, m, weights)
      idx <- geom$corr_index[[m]]
      curve[m, ] <- .column_correlation(power[idx, , drop = FALSE],
          tp[idx, , drop = FALSE], method)
    }
    curve
  }

  ecocox_active_orbit_count <- function(geom, weights){
    # participation counts per rate for one frozen weight vector
    ag <- as.numeric(weights) > 1e-8
    ag[is.na(ag)] <- FALSE
    ao <- ag[geom$group_index]
    as.numeric(geom$resolved %*% as.numeric(ao))
  }

  ecocox_min_na <- function(a, b){
    # MATLAB min(a,b): NaN ignored unless both are NaN
    out <- pmin(a, b, na.rm = TRUE)
    out[is.na(a) & is.na(b)] <- NaN
    out
  }

  ecocox_max_na <- function(a, b){
    out <- pmax(a, b, na.rm = TRUE)
    out[is.na(a) & is.na(b)] <- NaN
    out
  }

  ecocox_p_floor <- function(n_sim){
    if (n_sim > 0){
      return(1 / (n_sim + 1))
    }
    NaN
  }

  ecocox_pcoco <- function(rho, p){
    ps <- p
    ps[!is.finite(ps) | ps <= 0] <- NaN
    ps[!is.na(ps) & ps > 1] <- 1
    rho * abs(log10(ps))
  }

  # Physical-depth window sampling (ecocoPhysicalWindowValues.m)

  ecocox_physical_window_values <- function(data, window, centers, nominal_dt,
      min_source = 8){
    depth <- data[, 1]
    proxy <- data[, 2]
    nd <- length(depth)
    h <- ecocox_round(window / (2 * nominal_dt))
    if (h < 1){
      stop("WINDOW is too short relative to the sampling interval.")
    }
    npt <- 2 * h + 1
    adt <- window / (2 * h)
    rel <- (-h:h) * adt
    rel[1] <- -window / 2
    rel[h + 1] <- 0
    rel[npt] <- window / 2
    nw <- length(centers)
    V <- matrix(NaN, npt, nw)
    first <- rep(NA_integer_, nw)
    last <- rep(NA_integer_, nw)
    count <- rep(0, nw)
    obs_center <- rep(NaN, nw)
    obs_span <- rep(NaN, nw)
    for (w in seq_len(nw)){
      lo <- centers[w] - window / 2
      hi <- centers[w] + window / 2
      tol <- 64 * ecocox_eps(max(c(1, abs(depth[1]), abs(depth[nd]), abs(lo),
          abs(hi), abs(window)))) * max(1, nd)
      if (lo < depth[1] - tol || hi > depth[nd] + tol){
        stop(sprintf(paste0("Physical window %d requests bounds %.12g to %.12g, ",
            "outside the data support %.12g to %.12g."),
            w, lo, hi, depth[1], depth[nd]))
      }
      ii <- which(depth >= lo - tol & depth <= hi + tol)
      count[w] <- length(ii)
      if (length(ii) < min_source){
        stop(sprintf("Physical window %d contains %d observations; at least %d are required.",
            w, length(ii), min_source))
      }
      first[w] <- ii[1]
      last[w] <- ii[length(ii)]
      ld <- depth[ii]
      lp <- proxy[ii]
      obs_center[w] <- 0.5 * (ld[1] + ld[length(ld)])
      obs_span[w] <- ld[length(ld)] - ld[1]
      q <- centers[w] + rel
      q[1] <- lo
      q[h + 1] <- centers[w]
      q[npt] <- hi
      # aligned fast path: every query coincides with a source point
      aligned <- FALSE
      if (length(ld) >= npt){
        si <- 1
        ai <- rep(0L, npt)
        aligned <- TRUE
        for (qi in seq_len(npt)){
          while (si <= length(ld) && ld[si] < q[qi] - tol){
            si <- si + 1
          }
          if (si > length(ld) || abs(ld[si] - q[qi]) > tol){
            aligned <- FALSE
            break
          }
          ai[qi] <- si
          si <- si + 1
        }
      }
      if (aligned){
        V[, w] <- lp[ai]
      } else {
        V[, w] <- interp1_extrap(ld, lp, q)
      }
    }
    list(values = V, n = npt, dt = adt, first = first, last = last,
        count = count, observed_center = obs_center, observed_span = obs_span)
  }

  ecocox_zeropad <- function(data, window, padding = 1){
    # zeropad2.m with ensureHalfWindowCoverage = TRUE (padding 1 zero,
    # 2 mirror, 3 mean)
    data <- data[order(data[, 1]), , drop = FALSE]
    x <- data[, 1]
    y <- data[, 2] - mean(data[, 2])
    dt <- mean(diff(x))
    n <- ceiling(window / 2 / dt)
    x1 <- x[1] - (n:1) * dt
    x2 <- x[length(x)] + (1:n) * dt
    if (padding == 1){
      y1 <- rep(0, n)
      y2 <- rep(0, n)
    } else if (padding == 2){
      if (n > length(y)){
        stop("Mirror padding cannot exceed the input data length.")
      }
      y1 <- y[n:1]
      y2 <- y[length(y):(length(y) - n + 1)]
    } else if (padding == 3){
      y1 <- rep(mean(y[1:n]), n)
      y2 <- rep(mean(y[(length(y) - n + 1):length(y)]), n)
    } else {
      stop("padding must be 1 (zero), 2 (mirror) or 3 (mean)")
    }
    rbind(cbind(x1, y1), cbind(x, y), cbind(x2, y2))
  }

  # Result helpers shared by Blocked and Interleaved eCOCO

  ecocox_values_at <- function(M, index){
    out <- rep(NaN, length(index))
    for (j in seq_along(index)){
      if (!is.na(index[j]) && index[j] >= 1 && index[j] <= nrow(M)){
        out[j] <- M[index[j], j]
      }
    }
    out
  }

  ecocox_direction <- function(rho, p_global, p_local, norbit, sr, supported){
    nw <- ncol(rho)
    best_rate <- rep(NaN, nw)
    best_index <- rep(NA_integer_, nw)
    best_rho <- rep(NaN, nw)
    for (j in seq_len(nw)){
      b <- ecocox_first_max(rho[, j])
      if (!is.na(b)){
        best_rate[j] <- sr[b]
        best_index[j] <- b
        best_rho[j] <- rho[b, j]
      }
    }
    pc <- ecocox_pcoco(rho, p_global)
    list(rho = rho, p_global = p_global, p_local = p_local, norbit = norbit,
        pcoco = pc, score = pc, best_rate = best_rate,
        best_index = best_index, best_rho = best_rho,
        best_p_global = ecocox_values_at(p_global, best_index),
        supported = supported)
  }

  ecocox_finish <- function(out, sr, depth, verbose){
    # bestPerWindowSummary + trackEcocoRidge (ecoco.m) via ecoco_track_ridge()
    best <- ecoco_best_per_window(sr, depth, out$rho, out$p_global, out$norbit,
        out$score)
    track <- ecoco_track_ridge(out$score, sr, depth, out$rho, out$p_global,
        out$norbit, best)
    out$best <- best
    out$track <- track
    if (verbose){
      ok <- is.finite(track$sr)
      med <- NaN
      if (any(ok)){
        med <- stats::median(track$sr[ok])
      }
      message(sprintf("%s: %d windows, %d with a tracked rate (median %.4g cm/kyr)",
          out$method, length(depth), sum(ok), med))
    }
    out
  }

  ecocox_blocked_core <- function(dat, orbit9, window, dt, step, red, pad, sr,
      n_sim, method, max_frequency, seed,
      anchor_fraction, compute_local_p,
      window_mode, step_depth, center_limits,
      batch_size, verbose){
    depth <- dat[, 1]
    nd <- length(depth)
    if (any(diff(depth) <= 0)){
      stop("Depths must be strictly increasing")
    }
    if (window_mode == "legacy-count" &&
        any(abs(diff(depth) - dt) > ecocox_sampling_tolerance(dt))){
      stop("legacy-count mode requires evenly sampled data")
    }
    n_win <- 2 * ecocox_round(window / (2 * dt)) + 1
    if (n_win < 8){
      stop("window must contain at least eight samples")
    }
    if (window_mode == "legacy-count" && n_win > nd){
      stop("window is longer than the data")
    }
    if (pad < n_win){
      stop("pad must be at least the window point count")
    }
    dtol <- max(ecocox_sampling_tolerance(min(window, dt)),
        128 * ecocox_eps(max(1, max(abs(depth)))))
    if (window_mode == "physical-depth"){
      vcent <- ecoco_physical_centers(depth, window, step_depth, center_limits)
      anchor_step_depth <- anchor_fraction * window
      acent <- ecoco_physical_centers(depth, window, anchor_step_depth,
          center_limits)
      comb <- c(vcent, acent)
      ord <- order(comb)
      sc <- comb[ord]
      mtol <- max(ecocox_sampling_tolerance(min(max(sc[length(sc)] - sc[1], dt),
          dt)),
          128 * ecocox_eps(max(1, max(abs(depth)))))
      uc <- rep(0, length(sc))
      sidx <- rep(0L, length(sc))
      nu <- 0
      for (k in seq_along(sc)){
        if (nu == 0 || abs(sc[k] - uc[nu]) > mtol){
          nu <- nu + 1
          uc[nu] <- sc[k]
        }
        sidx[k] <- nu
      }
      block_center <- uc[seq_len(nu)]
      cidx <- rep(0L, length(comb))
      cidx[ord] <- sidx
      vblock <- cidx[seq_along(vcent)]
      ablock <- cidx[length(vcent) + seq_along(acent)]
      wv <- ecocox_physical_window_values(dat, window, block_center, dt, 8)
      values <- wv$values
      n_win <- wv$n
      adt <- wv$dt
      vbounds <- cbind(vcent - window / 2, vcent + window / 2)
      abounds <- cbind(acent - window / 2, acent + window / 2)
      block_bounds <- cbind(block_center - window / 2, block_center + window / 2)
      vstart <- wv$first[vblock]
      vend <- wv$last[vblock]
      block_count <- wv$count
      ratio <- anchor_step_depth / dt
      anchor_step <- NaN
      if (abs(ratio - ecocox_round(ratio)) <= 1e-10 * max(1, abs(ratio))){
        anchor_step <- ecocox_round(ratio)
      }
    } else {
      last_start <- nd - n_win + 1
      vstart <- seq(1, last_start, by = step)
      if (vstart[length(vstart)] != last_start){
        vstart <- c(vstart, last_start)
      }
      anchor_step <- max(1, ecocox_round(anchor_fraction * (n_win - 1)))
      astart <- seq(1, last_start, by = anchor_step)
      if (astart[length(astart)] != last_start){
        astart <- c(astart, last_start)
      }
      vend <- vstart + n_win - 1
      aend <- astart + n_win - 1
      vcent <- 0.5 * (depth[vstart] + depth[vend])
      acent <- 0.5 * (depth[astart] + depth[aend])
      vbounds <- cbind(depth[vstart], depth[vend])
      abounds <- cbind(depth[astart], depth[aend])
      bstart <- sort(unique(c(vstart, astart)))
      vblock <- match(vstart, bstart)
      ablock <- match(astart, bstart)
      values <- matrix(0, n_win, length(bstart))
      for (k in seq_along(bstart)){
        values[, k] <- dat[bstart[k]:(bstart[k] + n_win - 1), 2]
      }
      block_center <- 0.5 * (depth[bstart] + depth[bstart + n_win - 1])
      block_bounds <- cbind(depth[bstart], depth[bstart + n_win - 1])
      block_count <- rep(n_win, length(bstart))
      adt <- dt
      anchor_step_depth <- anchor_step * dt
    }
    n_pos <- length(vcent)
    n_anchor <- length(acent)
    n_block <- ncol(values)
    n_rate <- length(sr)
    fwd <- rep(NA_integer_, n_pos)
    bwd <- rep(NA_integer_, n_pos)
    for (j in seq_len(n_pos)){
      if (window_mode == "physical-depth"){
        cand <- which(abounds[, 2] < vbounds[j, 1] - dtol)
      } else {
        cand <- which(aend < vstart[j])
      }
      if (length(cand) > 0){
        fwd[j] <- cand[length(cand)]
      }
      if (window_mode == "physical-depth"){
        cand <- which(abounds[, 1] > vbounds[j, 2] + dtol)
      } else {
        cand <- which(astart > vend[j])
      }
      if (length(cand) > 0){
        bwd[j] <- cand[1]
      }
    }
    geom <- ecocox_geometry(n_win, adt, orbit9, sr, pad, max_frequency,
        "blocked")
    partial <- FALSE
    if (!any(geom$strict_training) && any(geom$partial_training)){
      partial <- TRUE
      geom$training_rate <- geom$partial_training
      warning(sprintf(paste0("Blocked eCOCO has no tested rate with complete ",
          "all-nine training geometry; continuing with the ",
          "partial-orbit target at %d of %d rates ",
          "(degraded/exploratory result)."),
          sum(geom$training_rate), n_rate))
    }
    if (!any(geom$training_rate)){
      stop("No tested sedimentation rate has a resolvable, well-conditioned training geometry.")
    }
    if (!any(geom$valid_rate)){
      stop("No tested sedimentation rate has usable validation geometry.")
    }
    if (verbose){
      message(sprintf("Blocked eCOCO: %d windows, %d anchors, %d points per window, %d rates",
          n_pos, n_anchor, n_win, n_rate))
    }
    block_power <- ecocox_window_spectra(values, adt, pad, red)$power
    ar <- ecocox_estimate_ar1(values)
    apow <- block_power[, ablock, drop = FALSE]
    atr <- ecocox_train(apow, geom, method, TRUE)
    aw <- ecocox_normalize(atr$amplitude)
    avalid <- is.finite(atr$best_rate) & colSums(!is.finite(aw)) == 0 &
        ar$valid[ablock]
    a_best_index <- match(atr$best_rate, sr)
    a_partial <- rep(FALSE, n_anchor)
    ok <- avalid & !is.na(a_best_index)
    a_partial[ok] <- geom$partial_training[a_best_index[ok]] &
        !geom$strict_training[a_best_index[ok]]
    has_f <- !is.na(fwd)
    has_b <- !is.na(bwd)
    for (j in seq_len(n_pos)){
      vv <- ar$valid[vblock[j]]
      if (has_f[j]){
        has_f[j] <- vv && avalid[fwd[j]]
      }
      if (has_b[j]){
        has_b[j] <- vv && avalid[bwd[j]]
      }
    }
    vpow <- block_power[, vblock, drop = FALSE]
    rho_f <- matrix(NaN, n_rate, n_pos)
    rho_b <- matrix(NaN, n_rate, n_pos)
    if (any(has_f)){
      iw <- which(has_f)
      rho_f[, iw] <- ecocox_validate(vpow[, iw, drop = FALSE], geom,
          aw[, fwd[iw], drop = FALSE], method)
    }
    if (any(has_b)){
      iw <- which(has_b)
      rho_b[, iw] <- ecocox_validate(vpow[, iw, drop = FALSE], geom,
          aw[, bwd[iw], drop = FALSE], method)
    }
    both <- has_f & has_b
    fonly <- has_f & !has_b
    bonly <- !has_f & has_b
    rho_c <- matrix(NaN, n_rate, n_pos)
    rho_c[, both] <- ecocox_min_na(rho_f[, both], rho_b[, both])
    rho_c[, fonly] <- rho_f[, fonly]
    rho_c[, bonly] <- rho_b[, bonly]
    support <- rep("none", n_pos)
    support[both] <- "bidirectional"
    support[fonly] <- "forward-only"
    support[bonly] <- "backward-only"
    rho_s <- rho_c
    rho_s[, !both] <- NaN

    # participation (active orbit) counts
    nor_f <- matrix(NaN, n_rate, n_pos)
    nor_b <- matrix(NaN, n_rate, n_pos)
    nor_c <- matrix(NaN, n_rate, n_pos)
    nor_s <- matrix(NaN, n_rate, n_pos)
    for (j in seq_len(n_pos)){
      if (has_f[j]){
        nor_f[, j] <- ecocox_active_orbit_count(geom, aw[, fwd[j]])
      }
      if (has_b[j]){
        nor_b[, j] <- ecocox_active_orbit_count(geom, aw[, bwd[j]])
      }
      if (has_f[j] && has_b[j]){
        w2 <- as.numeric(aw[, fwd[j]] > 1e-8 & aw[, bwd[j]] > 1e-8)
        nor_c[, j] <- ecocox_active_orbit_count(geom, w2)
        nor_s[, j] <- nor_c[, j]
      } else if (has_f[j]){
        nor_c[, j] <- nor_f[, j]
      } else if (has_b[j]){
        nor_c[, j] <- nor_b[, j]
      }
    }

    cg_f <- matrix(0, n_rate, n_pos)
    cg_b <- matrix(0, n_rate, n_pos)
    cg_c <- matrix(0, n_rate, n_pos)
    cl_f <- matrix(0, n_rate, n_pos)
    cl_b <- matrix(0, n_rate, n_pos)
    cl_c <- matrix(0, n_rate, n_pos)
    if (n_sim > 0){
      if (!is.null(seed)){
        set.seed(seed)
      }
      s_rho <- ar$rho
      s_sd <- ar$sd
      s_rho[!ar$valid] <- 0
      s_sd[!ar$valid] <- 1
      first <- 1
      pb <- NULL
      if (verbose){
        pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
      }
      while (first <= n_sim){
        last <- min(first + batch_size - 1, n_sim)
        nb <- last - first + 1
        rc <- rep(s_rho, nb)
        scol <- rep(s_sd, nb)
        Z <- matrix(stats::rnorm(n_win * n_block * nb), n_win, n_block * nb)
        Z[1, ] <- scol * Z[1, ]
        isc <- scol * sqrt(pmax(0, 1 - rc^2))
        for (r in 2:n_win){
          Z[r, ] <- rc * Z[r - 1, ] + isc * Z[r, ]
        }
        npow <- ecocox_window_spectra(Z, adt, pad, red)$power
        acol <- as.numeric(matrix(ablock, n_anchor, nb) +
            matrix((0:(nb - 1)) * n_block, n_anchor, nb,
                byrow = TRUE))
        ntr <- ecocox_train(npow[, acol, drop = FALSE], geom, method, FALSE)
        nw <- ecocox_normalize(ntr$amplitude)
        if (any(!is.finite(ntr$best_rate)) || any(!is.finite(nw))){
          stop("A null anchor failed the Blocked eCOCO training pipeline.")
        }
        for (j in seq_len(n_pos)){
          if (!(has_f[j] || has_b[j])){
            next
          }
          vcol <- vblock[j] + (0:(nb - 1)) * n_block
          nvp <- npow[, vcol, drop = FALSE]
          nf_c <- matrix(NaN, n_rate, nb)
          nb_c <- matrix(NaN, n_rate, nb)
          if (has_f[j]){
            wcol <- fwd[j] + (0:(nb - 1)) * n_anchor
            nf_c <- ecocox_validate(nvp, geom, nw[, wcol, drop = FALSE], method)
            cg_f[, j] <- cg_f[, j] + ecocox_global_exceed(nf_c, rho_f[, j],
                geom$valid_rate)
            if (compute_local_p){
              cl_f[, j] <- cl_f[, j] + ecocox_local_exceed(nf_c, rho_f[, j])
            }
          }
          if (has_b[j]){
            wcol <- bwd[j] + (0:(nb - 1)) * n_anchor
            nb_c <- ecocox_validate(nvp, geom, nw[, wcol, drop = FALSE], method)
            cg_b[, j] <- cg_b[, j] + ecocox_global_exceed(nb_c, rho_b[, j],
                geom$valid_rate)
            if (compute_local_p){
              cl_b[, j] <- cl_b[, j] + ecocox_local_exceed(nb_c, rho_b[, j])
            }
          }
          if (has_f[j] && has_b[j]){
            nc_c <- ecocox_min_na(nf_c, nb_c)
          } else if (has_f[j]){
            nc_c <- nf_c
          } else {
            nc_c <- nb_c
          }
          cg_c[, j] <- cg_c[, j] + ecocox_global_exceed(nc_c, rho_c[, j],
              geom$valid_rate)
          if (compute_local_p){
            cl_c[, j] <- cl_c[, j] + ecocox_local_exceed(nc_c, rho_c[, j])
          }
        }
        if (!is.null(pb)){
          utils::setTxtProgressBar(pb, last)
        }
        first <- last + 1
      }
      if (!is.null(pb)){
        close(pb)
      }
    }
    pg_f <- ecocox_p_counts(cg_f, rho_f, n_sim)
    pg_b <- ecocox_p_counts(cg_b, rho_b, n_sim)
    pg_c <- ecocox_p_counts(cg_c, rho_c, n_sim)
    pg_s <- pg_c
    pg_s[, !both] <- NaN
    if (compute_local_p){
      pl_f <- ecocox_p_counts(cl_f, rho_f, n_sim)
      pl_b <- ecocox_p_counts(cl_b, rho_b, n_sim)
      pl_c <- ecocox_p_counts(cl_c, rho_c, n_sim)
    } else {
      pl_f <- rho_f * NaN
      pl_b <- rho_b * NaN
      pl_c <- rho_c * NaN
    }
    pl_s <- pl_c
    pl_s[, !both] <- NaN
    p_dirmax <- matrix(NaN, n_rate, n_pos)
    p_dirmax[, both] <- ecocox_max_na(pg_f[, both], pg_b[, both])
    p_dirmax[, fonly] <- pg_f[, fonly]
    p_dirmax[, bonly] <- pg_b[, bonly]
    p_both <- matrix(NaN, n_rate, n_pos)
    p_both[, both] <- ecocox_max_na(pg_f[, both], pg_b[, both])

    forward <- ecocox_direction(rho_f, pg_f, pl_f, nor_f, sr, has_f)
    backward <- ecocox_direction(rho_b, pg_b, pl_b, nor_b, sr, has_b)
    consensus <- ecocox_direction(rho_c, pg_c, pl_c, nor_c, sr, has_f | has_b)
    consensus$p_directional_max <- p_dirmax
    consensus$p_both <- p_both
    strict <- ecocox_direction(rho_s, pg_s, pl_s, nor_s, sr, both)
    strict$p_directional_max <- p_both
    windows <- data.frame(center = vcent, lower = vbounds[, 1],
        upper = vbounds[, 2], start_index = vstart,
        end_index = vend, block = vblock,
        forward_anchor = fwd, backward_anchor = bwd,
        has_forward = has_f, has_backward = has_b,
        support = support, rho_m = ar$rho[vblock],
        valid_null = ar$valid[vblock],
        stringsAsFactors = FALSE)
    anchors <- data.frame(center = acent, lower = abounds[, 1],
        upper = abounds[, 2], block = ablock,
        best_rate = atr$best_rate,
        best_rho = atr$best_score,
        w_long_ecc = aw[1, ], w_short_ecc = aw[2, ],
        w_obliquity = aw[3, ], w_precession = aw[4, ],
        valid = avalid, partial_orbit = a_partial,
        rho_m = ar$rho[ablock])
    list(method = "Blocked eCOCO", sr = sr, depth = vcent,
        rho = rho_c, p_global = pg_c, p_local = pl_c, norbit = nor_c,
        pcoco = consensus$pcoco, score = consensus$score,
        forward = forward, backward = backward, consensus = consensus,
        strict = strict, support = support, windows = windows,
        anchors = anchors, anchor_train_curve = atr$curve,
        block_center = block_center, block_count = block_count,
        rho_m = ar$rho[vblock], window = window, window_points = n_win,
        analysis_dt = adt, sampling_dt = dt, step_depth = step_depth,
        anchor_step_depth = anchor_step_depth, anchor_step = anchor_step,
        anchor_fraction = anchor_fraction, window_mode = window_mode,
        geometry = list(valid_rate = geom$valid_rate,
            training_rate = geom$training_rate,
            strict_training = geom$strict_training,
            partial_training = geom$partial_training,
            orbit_count = geom$orbit_count,
            active_group = geom$active_group,
            leakage_rcond = geom$active_rcond,
            full_leakage_rcond = geom$full_rcond,
            group_overlap = geom$overlap,
            rayleigh = geom$rayleigh, nyquist = geom$nyquist),
        degraded = partial, n_sim = n_sim,
        p_floor = ecocox_p_floor(n_sim), seed = seed,
        red = red, pad = pad, max_frequency = max_frequency, orbit9 = orbit9,
        correlation = method)
  }

  ecocox_global_exceed <- function(null_curve, observed, valid_rate){
    null_curve <- as.matrix(null_curve)
    rows <- valid_rate & rowSums(is.finite(null_curve)) > 0
    if (!any(rows)){
      stop("A null validation curve has no finite sedimentation-rate statistic.")
    }
    nmax <- ecocox_col_max_finite(null_curve[rows, , drop = FALSE])
    if (any(!is.finite(nmax))){
      stop("A Monte Carlo realization has no finite SR-global maximum.")
    }
    cnt <- rowSums(matrix(nmax, length(observed), length(nmax), byrow = TRUE) >=
        matrix(observed, length(observed), length(nmax)))
    cnt[!is.finite(observed)] <- 0
    cnt
  }

  ecocox_local_exceed <- function(null_curve, observed){
    null_curve <- as.matrix(null_curve)
    ge <- null_curve >= matrix(observed, nrow(null_curve), ncol(null_curve))
    ge[is.na(ge)] <- FALSE
    cnt <- rowSums(ge)
    cnt[!is.finite(observed)] <- 0
    cnt
  }

  ecocox_p_counts <- function(count, observed, n_sim){
    p <- observed * NaN
    if (n_sim <= 0){
      return(p)
    }
    ok <- is.finite(observed)
    p[ok] <- (count[ok] + 1) / (n_sim + 1)
    p
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

  # MATLAB periodogram (rectangular window, one-sided PSD)
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

  # Interpolation
  interp1_na <- function(x, y, xi){
    # Linear interpolation, NA outside range (MATLAB interp1 default).
    o <- order(x)
    x <- x[o]
    y <- y[o]
    keep <- !duplicated(x)
    stats::approx(x[keep], y[keep], xout = xi, method = "linear", rule = 1,
        ties = "ordered")$y
  }

  interp1_extrap <- function(x, y, xi){
    # Linear interpolation with linear extrapolation (MATLAB 'linear','extrap').
    out <- interp1_na(x, y, xi)
    n <- length(x)
    if (n >= 2){
      lo <- xi < x[1]
      hi <- xi > x[n]
      if (any(lo)){
        out[lo] <- y[1] + (xi[lo] - x[1]) * (y[2] - y[1]) / (x[2] - x[1])
      }
      if (any(hi)){
        out[hi] <- y[n] + (xi[hi] - x[n]) * (y[n] - y[n - 1]) / (x[n] - x[n - 1])
      }
    }
    out
  }

  # Moving median with shrinking endpoints
  move_median <- function(x, w){
    x <- as.vector(x)
    n <- length(x)
    halfw <- floor(w / 2)
    if (w %% 2 == 1){
      back <- halfw
      ahead <- halfw
    } else {
      back <- max(0, halfw - 1)
      ahead <- halfw
    }
    out <- numeric(n)
    for (i in seq_len(n)){
      out[i] <- stats::median(x[max(1, i - back):min(n, i + ahead)])
    }
    out
  }

  # AR(1) red-noise backgrounds
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

  minirhos0 <- function(s0, fn, ft, pxxsmooth, linlog = 2){
    cospara <- cos(pi * as.vector(ft) / fn)
    power <- as.vector(pxxsmooth)
    dist_grid <- function(rho_g, s0_g){
      # AR(1) shape for every (frequency, rho) pair: rows = frequencies
      cmat <- matrix(cospara, length(cospara), length(rho_g))
      rmat <- matrix(rho_g, length(cospara), length(rho_g), byrow = TRUE)
      shape <- (1 - rmat^2) / (1 - 2 * cmat * rmat + rmat^2)
      if (linlog == 1){
        a <- colSums(shape^2)
        b <- as.vector(crossprod(shape, power))
        cc <- sum(power^2)
        outer(a, s0_g^2) - 2 * outer(b, s0_g) + cc
      } else {
        base <- log(shape) - log(power)
        rs <- colSums(base)
        rss <- colSums(base^2)
        ls <- log(s0_g)
        matrix(rss, length(rho_g), length(s0_g)) + 2 * outer(rs, ls) +
            length(power) * matrix(ls^2, length(rho_g), length(s0_g), byrow = TRUE)
      }
    }
    nn <- 50
    rhoi <- seq(0.001, 0.999, length.out = nn)
    s0i <- seq(0.2 * s0, 5 * s0, length.out = nn)
    d <- dist_grid(rhoi, s0i)
    k <- which.min(d)
    x <- (k - 1) %% nn + 1
    y <- (k - 1) %/% nn + 1
    mm <- nn / 2
    for (it in 1:3){
      ds <- s0i[2] - s0i[1]
      dr <- rhoi[2] - rhoi[1]
      s0max <- s0i[y] + ds
      s0min <- s0i[y] - ds
      rmax <- min(rhoi[x] + dr, 0.9999)
      rmin <- rhoi[x] - dr
      if (rmin <= 0){
        rmin <- 1e-4
      }
      rhoi <- seq(rmin, rmax, length.out = mm)
      s0i <- seq(s0min, s0max, length.out = mm)
      d <- dist_grid(rhoi, s0i)
      k <- which.min(d)
      x <- (k - 1) %% mm + 1
      y <- (k - 1) %/% mm + 1
    }
    list(rho = rhoi[x], s0 = s0i[y])
  }

  redconf_any <- function(f_angular, pxx, dt, smoothwin = 0.25, linlog = 2){
    # Robust AR(1) background (Mann & Lees 1996 style): median-smoothed spectrum,
    # grid-search fit of rho and S0. f_angular is radians/sample (2*pi*f*dt).
    fn <- 1 / (2 * dt)
    ft <- f_angular / pi * fn
    smoothn <- max(1, round(smoothwin * length(pxx)))
    pxxsmooth <- move_median(pxx, smoothn)
    pos <- pxxsmooth[pxxsmooth > 0]
    if (length(pos)){
      floor_val <- max(.Machine$double.xmin, max(pos) * .Machine$double.eps)
    } else {
      floor_val <- 1e-300
    }
    pxxsmooth <- pmax(pxxsmooth, floor_val)
    s0 <- mean(pxxsmooth)
    fit <- minirhos0(s0, fn, ft, pxxsmooth, linlog)
    fit$s0 * (1 - fit$rho^2) / (1 - 2 * fit$rho * cos(pi * ft / fn) + fit$rho^2)
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

  # Batched spectra with background removal (ecocoWindowSpectra.m)
  coco_window_spectra <- function(values, dt, pad, red = 1, already_detrended = FALSE){
    values <- as.matrix(values)
    if (!already_detrended){
      values <- detrend_linear(values)
    }
    sp <- ml_periodogram(values, pad, 1 / dt)
    P <- as.matrix(sp$p)
    f <- sp$f
    P[P < 0] <- 0
    if (red == 0){
      return(list(power = P, f = f))
    }
    if (red == 3){
      return(coco_swa_subtract(P, f, nrow(values)))
    }
    for (j in seq_len(ncol(P))){
      p <- P[, j]
      if (red == 1){
        bg <- theored_ar1_ml(values[, j], f, mean(p), dt)
      } else {
        bg <- redconf_any(2 * pi * f * dt, p, dt, 0.25, 2)
      }
      p <- p - bg
      p[p < 0] <- 0
      P[, j] <- p
    }
    list(power = P, f = f)
  }

  # @title Remove the smoothed-window-average background (red = 3)
  # @description Helper of \code{coco_window_spectra}. Acycle red = 3: every
  # one-sided periodogram (at least 33 bins) is floored at
  # \code{max(realmin, max(positive power) * eps)}, converted to log10, passed
  # to specswa (\code{swa_background}), and the SWA background is subtracted
  # with negative values set to zero. Translated from Acycle v3.1
  # corrcoefslices_rankNew.m, redNoisePeriodogramMC.m and ecocoWindowSpectra.m
  # (Mingsong Li); verified in verification/ecoco_ext (red = 3 cases).
  # @keywords internal
  # @noRd
  coco_swa_subtract <- function(P, f, n_data){
    if (nrow(P) < 33){
      stop(paste0("red = 3 (smoothed-window average) requires at least 33 ",
          "one-sided spectral bins; increase pad"))
    }
    for (j in seq_len(ncol(P))){
      p <- P[, j]
      pos <- p[is.finite(p) & p > 0]
      if (length(pos) < 1){
        sc <- 1
      } else {
        sc <- max(pos)
      }
      fl <- max(2.2250738585072014e-308, sc * .Machine$double.eps)
      bg <- as.numeric(swa_background(f, log10(pmax(p, fl)), n_data,
          verbose = FALSE)$swa)
      if (length(bg) != length(p) || any(!is.finite(bg) | bg <= 0)){
        stop("The smoothed-window background is nonpositive or nonfinite.")
      }
      p <- p - bg
      p[!is.finite(p) | p < 0] <- 0
      P[, j] <- p
    }
    list(power = P, f = f)
  }

  .column_correlation <- function(x, y, method){
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (method == "Spearman"){
      xr <- x
      for (j in seq_len(ncol(x))){
        xr[, j] <- rank(x[, j])
      }
      yr <- y
      for (j in seq_len(ncol(y))){
        yr[, j] <- rank(y[, j])
      }
      x <- as.matrix(xr)
      y <- as.matrix(yr)
    }
    n <- nrow(x)
    xc <- sweep(x, 2, colMeans(x))
    yc <- sweep(y, 2, colMeans(y))
    ssx <- colSums(xc^2)
    ssy <- colSums(yc^2)
    maxx <- numeric(ncol(x))
    for (j in seq_len(ncol(x))){
      maxx[j] <- max(abs(x[, j]))
    }
    maxy <- numeric(ncol(y))
    for (j in seq_len(ncol(y))){
      maxy[j] <- max(abs(y[, j]))
    }
    tolx <- n * (256 * .eps * pmax(maxx, .Machine$double.xmin))^2
    toly <- n * (256 * .eps * pmax(maxy, .Machine$double.xmin))^2
    r <- colSums(xc * yc) / sqrt(ssx * ssy)
    r[ssx <= tolx | ssy <= toly] <- 0
    r[!is.finite(r)] <- NA
    r
  }

  # Physical-depth windows (ecocoPhysicalCenters.m / ecocoPhysicalWindowValues.m)
  ecoco_physical_centers <- function(depth, window, step_depth, center_limits = NULL){
    if (is.null(center_limits)){
      lim <- c(depth[1] + window / 2, depth[length(depth)] - window / 2)
    } else {
      lim <- center_limits
    }
    if (lim[2] < lim[1]){
      stop("Data span is shorter than the window")
    }
    ratio <- (lim[2] - lim[1]) / step_depth
    nint <- max(0, floor(ratio + 64 * .eps * max(1, abs(ratio))))
    if (nint + 1 > 10000){
      stop("More than 10000 windows; increase the step")
    }
    cen <- lim[1] + (0:nint) * step_depth
    cen[cen <= lim[2] + 64 * .eps * max(1, abs(lim), abs(cen))]
  }

  # Small helpers
  .eps <- .Machine$double.eps

  # Main code
  method <- match.arg(method, c("Pearson", "Spearman"))
  window_mode <- match.arg(window_mode, c("physical-depth", "legacy-count"))
  consensus_policy <- match.arg(consensus_policy, c("supported", "strict"))
  orbit9 <- as.numeric(orbit9)
  if (length(orbit9) != 9 || length(unique(orbit9)) != 9 || any(orbit9 <= 0)){
    stop("orbit9 must contain nine distinct positive periods")
  }
  if (is.null(max_frequency)){
    max_frequency <- 1.2 * max(1 / orbit9)
  }
  if (max_frequency < max(1 / orbit9) - 64 * ecocox_eps(max(1, max(1 / orbit9)))){
    stop("max_frequency must include the highest nominal orbital frequency")
  }
  if (!(red %in% 0:3)){
    stop("red must be 0, 1, 2 or 3")
  }
  if (!is.numeric(anchor_fraction) || anchor_fraction <= 0 ||
      anchor_fraction > 2){
    stop("anchor_fraction must be in (0, 2]")
  }
  prep <- coco_prepare_regular_data(data, verbose = verbose)
  dat <- unname(prep$data)
  dt <- prep$dt
  n_raw <- nrow(dat)
  if (is.null(window)){
    window <- max(1e-6, 0.25 * abs(dat[n_raw, 1] - dat[1, 1]))
  }
  if (is.null(step_depth)){
    step_depth <- dt
    if (n_raw > 300){
      step_depth <- dt * ceiling(n_raw / 300)
    }
  }
  if (is.null(pad)){
    pad <- ecocox_default_pad(n_raw)
  }
  if (is.null(sr)){
    sr <- ecocox_default_sr(n_raw, dt, orbit9)
  }
  sr <- as.numeric(sr)
  if (any(diff(sr) <= 0) || any(sr <= 0)){
    stop("sr must be positive and strictly increasing")
  }
  if (edge_pad && window_mode == "physical-depth"){
    if (is.null(center_limits)){
      center_limits <- dat[c(1, n_raw), 1]
    }
    dat <- ecocox_zeropad(dat, window, 1)
    prep <- coco_prepare_regular_data(dat, verbose = verbose)
    dat <- unname(prep$data)
    dt <- prep$dt
  }
  res <- ecocox_blocked_core(dat, orbit9, window, dt, step, red, pad, sr,
      n_sim, method, max_frequency, seed,
      anchor_fraction, compute_local_p, window_mode,
      step_depth, center_limits, batch_size, verbose)
  if (consensus_policy == "strict"){
    res$rho <- res$strict$rho
    res$p_global <- res$strict$p_global
    res$p_local <- res$strict$p_local
    res$norbit <- res$strict$norbit
    res$pcoco <- res$strict$pcoco
    res$score <- res$strict$score
  }
  res$consensus_policy <- consensus_policy
  res <- ecocox_finish(res, sr, res$depth, verbose)
  class(res) <- "ecoco_cv"
  if (genplot){
    plot_ecoco_cv(res)
  }
  res
}
