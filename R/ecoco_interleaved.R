#' @title Interleaved eCOCO (odd/even cross-validated evolutionary COCO)
#'
#' @description Interleaved eCOCO evaluates complete sliding windows of the
#' cleaned (finite, sorted, duplicate-averaged) record. Within each window
#' the globally odd observations train a frozen four-group coherent-nine
#' astronomical target that is validated on the globally even observations,
#' and vice versa (Interleaved cvCOCO applied per window). Each fold is
#' linearly regularised at its own median spacing after the split. The
#' consensus correlation is the same-rate minimum of the two directional
#' curves; its global p-value is a plus-one Monte Carlo p against the null
#' maximum over the complete rate grid. Each null realisation is one
#' stationary Gaussian AR(1) process on the complete window order, split and
#' interpolated like the data. Windows that cannot be resolved for the rate
#' grid are returned empty. pCOCO = rho * |log10(p_global)| is the ridge
#' score. Translated from Acycle v3.1 ecocoInterleavedCore.m,
#' interleavedcvcoco.m and cvcoco.m (Mingsong Li). Requires the cocoR helpers
#' (part of this package).
#'
#' @param data Record to analyse (matrix or data frame, two columns): depth in metres in column 1 and the proxy in
#' column 2 (at least 16 unique depths). Give the raw, possibly irregular data; do not interpolate beforehand,
#' because the odd and even observations are split before each fold is regularised. Non-finite rows are dropped,
#' rows sorted and duplicated depths averaged. Required (\code{Default=NULL}).
#' @param orbit9 The nine astronomical target periods in kyr (numeric vector of nine distinct positive values: long
#' eccentricity, four short eccentricity terms, obliquity, three precession terms), e.g. from
#' \code{\link{orbit9_at_age}} for the age of the record.
#' \code{Default=c(405.6912, 130.6979, 123.8532, 98.8517, 94.8856, 40.9897, 23.6820, 22.3758, 18.9519)}.
#' @param window Length of the sliding window, in metres (numeric, > 0, at most the record span; windows with fewer
#' than 16 observations are left empty). \code{NULL} uses 25\% of the record length. \code{Default=NULL}.
#' @param sr Sedimentation rates to test, in cm/kyr (numeric vector, > 0, strictly increasing with a constant step).
#' \code{NULL} uses the Acycle GUI default: from \code{2 * dt * 0.065 * 100} up to 100 cm/kyr (reduced so that the
#' longest period fits in the record), in steps of 0.1 cm/kyr or at most 300 steps. \code{Default=NULL}.
#' @param step_depth Spacing of the window centres, in metres (numeric, > 0); only used with
#' \code{window_mode="physical-depth"}. \code{NULL} uses \code{dt}, or \code{dt * ceiling(n / 300)} for more than
#' 300 samples. \code{Default=NULL}.
#' @param step Step between successive windows, in samples (integer, >= 1); only used with
#' \code{window_mode="legacy-count"}. \code{Default=1}.
#' @param pad Length (number of points) to which the periodogram of each fold is zero-padded; each fold uses
#' \code{max(pad, fold length)} (integer, > 0). \code{NULL} uses the Acycle default: 5000 for up to 2500 samples,
#' 10000 for up to 5000 samples, and \code{floor(n / 5000) * 5000 + 5000} above that. \code{Default=NULL}.
#' @param red Background (red-noise) removal applied to every spectrum (integer): \code{0} none, \code{1} classical
#' AR(1) background, \code{2} robust AR(1) background, \code{3} smoothed-window-average (SWA) background.
#' \code{Default=0}.
#' @param n_sim Number of AR(1) Monte Carlo simulations per window used for the global and local p-values (integer,
#' >= 0). \code{Default=2000}.
#' @param method Correlation coefficient between the data and target spectra (character): \code{"Pearson"} or
#' \code{"Spearman"}. \code{Default="Pearson"}.
#' @param max_frequency Highest frequency of the target spectrum that enters the correlation, in cycles/kyr
#' (numeric, at least \code{max(1/orbit9)}); \code{NULL} uses \code{1.2 * max(1/orbit9)}. \code{Default=NULL}.
#' @param seed Parent seed of the random number generator (integer). Each window uses a seed derived from it and
#' from the window position (Acycle deriveWindowSeed), so results are reproducible and independent of the other
#' windows; \code{NULL} leaves the current random state unchanged. \code{Default=1}.
#' @param window_mode How the windows are defined (character): \code{"physical-depth"} (windows of fixed length in
#' metres centred every \code{step_depth}) or \code{"legacy-count"} (windows of a fixed, even number of points,
#' \code{window / dt + 1} rounded, at least 16, moved by \code{step} samples). \code{Default="physical-depth"}.
#' @param dt Nominal sampling interval, in metres (numeric, > 0), used for the number of points per window in
#' \code{"legacy-count"} mode and for the defaults of \code{step_depth} and \code{sr}. \code{NULL} uses the median
#' spacing of the cleaned data. \code{Default=NULL}.
#' @param batch_size Number of Monte Carlo simulations processed together (integer, >= 1); larger batches are faster
#' but use more memory. \code{Default=20}.
#' @param genplot Plot the result with \code{\link{plot_ecoco_cv}} (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param verbose Print messages and a progress bar over the windows (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' @references Li, M., Kump, L.R., Hinnov, L.A., Mann, M.E., 2018. Tracking
#' variable sedimentation rates and astronomical forcing in Phanerozoic
#' paleoclimate proxy series with evolutionary correlation coefficients and
#' hypothesis testing. Earth and Planetary Science Letters 501, 165-179.
#' \doi{10.1016/j.epsl.2018.08.041}
#'
#' @examples
#' \donttest{
#' set.seed(2)
#' dt <- 0.05
#' depth <- seq(0, 30, by = dt)
#' tt <- depth * 100 / 5
#' x <- sin(2 * pi * tt / 405.7) + 0.8 * sin(2 * pi * tt / 98.9) +
#'   0.5 * sin(2 * pi * tt / 40.99) + 0.4 * sin(2 * pi * tt / 23.68) +
#'   0.3 * rnorm(length(tt))
#' res <- ecoco_interleaved(cbind(depth, x), window = 15,
#'                          sr = seq(3, 8, 0.5), step_depth = 2, pad = 400,
#'                          n_sim = 20)
#' res$track[, c("depth", "sr")]
#' }
#'
#' @return A list (class "ecoco_cv") with \code{method}, \code{sr},
#' \code{depth}, rate-by-window matrices \code{rho}, \code{p_global},
#' \code{p_local}, \code{norbit}, \code{pcoco}, \code{score}, lists
#' \code{forward} (odd -> even), \code{backward} (even -> odd),
#' \code{consensus} and \code{strict} (consensus curve with the larger
#' directional p), \code{windows} (data frame: bounds, point counts, parity
#' phase, derived seed, success, window-level p-values, fold spacings,
#' all-nine rate ranges, AR(1) rho), \code{best} and \code{track}.
#' @export
#' @importFrom stats rnorm median
#' @importFrom utils txtProgressBar setTxtProgressBar
ecoco_interleaved <- function(data = NULL,
    orbit9 = c(405.6912, 130.6979, 123.8532,
        98.8517, 94.8856, 40.9897, 23.6820,
        22.3758, 18.9519),
    window = NULL, sr = NULL, step_depth = NULL,
    step = 1, pad = NULL, red = 0, n_sim = 2000,
    method = "Pearson", max_frequency = NULL,
    seed = 1, window_mode = "physical-depth",
    dt = NULL, batch_size = 20, genplot = FALSE,
    verbose = FALSE){
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

  ecocox_seed <- function(seed){
    # map a MATLAB uint32 seed to an R integer seed
    as.integer(seed %% 2147483647)
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

  ecocox_estimate_rho_m <- function(values){
    # estimateRhoM in cvcoco.m
    x <- as.numeric(values)
    x <- x - mean(x)
    sc <- max(abs(x))
    if (length(x) < 4 || any(!is.finite(x)) || !is.finite(sc) || sc <= 0){
      stop("Each segment needs at least four finite, variable detrended values.")
    }
    z <- x / sc
    a <- z[-length(z)]
    b <- z[-1]
    den <- sum(a * a)
    if (!is.finite(den) || den <= 0){
      stop("The conditional AR(1) denominator is zero or nonfinite.")
    }
    r <- sum(a * b) / den
    min(max(r, -0.999), 0.999)
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

  ecocox_curve_best <- function(curve, sr_grid){
    # fixedValidate best score / rate (strictly greater, first rate wins)
    curve <- as.matrix(curve)
    ns <- ncol(curve)
    score <- rep(NaN, ns)
    rate <- rep(NaN, ns)
    index <- rep(NA_integer_, ns)
    for (k in seq_len(ns)){
      b <- ecocox_first_max(curve[, k])
      if (!is.na(b)){
        score[k] <- curve[b, k]
        rate[k] <- sr_grid[b]
        index[k] <- b
      }
    }
    list(score = score, rate = rate, index = index)
  }

  ecocox_min_finite_pair <- function(a, b){
    # minFinitePair (cvcoco.m): NaN unless both finite
    out <- a * NaN
    ok <- is.finite(a) & is.finite(b)
    out[ok] <- pmin(a[ok], b[ok])
    out
  }

  ecocox_max_finite_pair <- function(a, b){
    out <- a * NaN
    ok <- is.finite(a) & is.finite(b)
    out[ok] <- pmax(a[ok], b[ok])
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

  ecocox_interp_columns <- function(x, Y, xi){
    # interp1(x, Y, xi, 'linear') for every column of Y (NaN outside range)
    Y <- as.matrix(Y)
    nx <- length(x)
    out <- matrix(NaN, length(xi), ncol(Y))
    inside <- xi >= x[1] & xi <= x[nx]
    if (!any(inside)){
      return(out)
    }
    xin <- xi[inside]
    k <- findInterval(xin, x, rightmost.closed = TRUE, all.inside = TRUE)
    s <- (xin - x[k]) / (x[k + 1] - x[k])
    val <- Y[k, , drop = FALSE] + (Y[k + 1, , drop = FALSE] -
        Y[k, , drop = FALSE]) * s
    exact <- xin == x[k + 1]
    if (any(exact)){
      val[exact, ] <- Y[k[exact] + 1, , drop = FALSE]
    }
    out[inside, ] <- val
    out
  }

  ecocox_clean_unique <- function(data){
    # finite rows, sorted, duplicate depths averaged
    data <- as.matrix(data)[, 1:2, drop = FALSE]
    n_in <- nrow(data)
    data <- data[is.finite(data[, 1]) & is.finite(data[, 2]), , drop = FALSE]
    n_fin <- nrow(data)
    data <- data[order(data[, 1]), , drop = FALSE]
    ud <- unique(data[, 1])
    if (length(ud) < nrow(data)){
      grp <- match(data[, 1], ud)
      s <- rowsum(data[, 2], grp, reorder = TRUE)
      cnt <- rowsum(rep(1, nrow(data)), grp, reorder = TRUE)
      data <- cbind(ud, as.numeric(s / cnt))
    }
    list(data = data, n_input = n_in, n_finite = n_fin,
        n_duplicate = n_fin - nrow(data))
  }

  # One held-out interleaved window (cvcoco.m, SplitMode 'interleaved',
  # TargetModel 'four-group-coherent-nine')

  ecocox_regularize_fold <- function(raw){
    sp <- diff(raw[, 1])
    sp <- sp[is.finite(sp) & sp > 0]
    dz <- stats::median(sp)
    if (length(sp) < 1 || !is.finite(dz) || dz <= 0){
      stop("A fold has no positive finite depth spacing.")
    }
    tol <- ecocox_sampling_tolerance(dz)
    if (max(abs(sp - dz)) <= tol){
      return(raw)
    }
    ex <- (raw[nrow(raw), 1] - raw[1, 1]) / dz
    rex <- ecocox_round(ex)
    if (abs(ex - rex) <= 1e-10 * max(1, abs(ex))){
      ni <- rex
    } else {
      ni <- floor(ex)
    }
    if (ni + 1 < 4){
      stop("A fold contains fewer than four points on the median-spacing grid.")
    }
    grid <- raw[1, 1] + (0:ni) * dz
    etol <- 16 * ecocox_eps(max(1, max(abs(raw[c(1, nrow(raw)), 1])))) *
        max(1, length(grid))
    if (abs(grid[length(grid)] - raw[nrow(raw), 1]) <= etol){
      grid[length(grid)] <- raw[nrow(raw), 1]
    }
    vals <- interp1_na(raw[, 1], raw[, 2], grid)
    ok <- is.finite(grid) & is.finite(vals)
    out <- cbind(grid[ok], vals[ok])
    if (nrow(out) < 4){
      stop("A fold contains fewer than four finite interpolated points.")
    }
    out
  }

  ecocox_unresolvable <- function(msg){
    structure(class = c("ecocox_unresolvable", "error", "condition"),
        list(message = msg, call = NULL))
  }

  ecocox_interleaved_window <- function(data, orbit9, pad, sr_grid, red, n_sim,
      method, max_frequency, phase,
      batch_size = 20){
    clean <- data
    nc <- nrow(clean)
    if (nc < 8){
      stop("At least eight finite observations are required.")
    }
    odd <- seq(1, nc, by = 2)
    even <- seq(2, nc, by = 2)
    if (phase == 0){
      ia <- odd
      ib <- even
    } else {
      ia <- even
      ib <- odd
    }
    raw_a <- clean[ia, , drop = FALSE]
    raw_b <- clean[ib, , drop = FALSE]
    if (nrow(raw_a) < 4 || nrow(raw_b) < 4){
      stop("Each odd/even fold must contain at least four unique observations.")
    }
    A <- ecocox_regularize_fold(raw_a)
    B <- ecocox_regularize_fold(raw_b)
    dza <- stats::median(diff(A[, 1]))
    dzb <- stats::median(diff(B[, 1]))
    det_a <- detrend_linear(A[, 2])
    det_b <- detrend_linear(B[, 2])
    rva <- ecocox_resolved_variance(A[, 2], det_a)
    rvb <- ecocox_resolved_variance(B[, 2], det_b)
    if (!rva$resolved || !rvb$resolved){
      stop("Each held-out fold must retain resolved variance after detrending.")
    }
    rho_ma <- ecocox_estimate_rho_m(det_a)
    rho_mb <- ecocox_estimate_rho_m(det_b)
    det_full <- detrend_linear(clean[, 2])
    rvf <- ecocox_resolved_variance(clean[, 2], det_full)
    if (!rvf$resolved){
      stop("The complete window must retain resolved variance for the null.")
    }
    sd_full <- rvf$sd
    rho_m <- ecocox_estimate_rho_m(det_full)

    ga <- ecocox_geometry(nrow(A), dza, orbit9, sr_grid, max(pad, nrow(A)),
        max_frequency, "cvcoco")
    gb <- ecocox_geometry(nrow(B), dzb, orbit9, sr_grid, max(pad, nrow(B)),
        max_frequency, "cvcoco")
    if (!any(ga$valid_rate) || !any(gb$valid_rate)){
      stop(ecocox_unresolvable(
          "No tested sedimentation rate has usable geometry in both folds."))
    }
    degraded <- FALSE
    comp_a <- "complete-nine"
    comp_b <- "complete-nine"
    if (!any(ga$strict_training) || !any(gb$strict_training)){
      if (any(ga$partial_training) && any(gb$partial_training)){
        degraded <- TRUE
        if (!any(ga$strict_training)){
          ga$training_rate <- ga$partial_training
          comp_a <- "partial-orbit"
        }
        if (!any(gb$strict_training)){
          gb$training_rate <- gb$partial_training
          comp_b <- "partial-orbit"
        }
      } else {
        stop(ecocox_unresolvable(
            "No tested rate supports four-group training in both folds."))
      }
    }
    spa <- ecocox_window_spectra(A[, 2], dza, ga$nfft, red)$power
    spb <- ecocox_window_spectra(B[, 2], dzb, gb$nfft, red)$power
    tra <- ecocox_train(spa, ga, method, TRUE)
    trb <- ecocox_train(spb, gb, method, TRUE)
    amp9a <- tra$amplitude[ga$group_index, , drop = FALSE]
    amp9b <- trb$amplitude[gb$group_index, , drop = FALSE]
    fa <- ecocox_normalize(ecocox_group_rms(amp9a, ga$group_index))
    fb <- ecocox_normalize(ecocox_group_rms(amp9b, gb$group_index))
    if (any(!is.finite(fa)) || any(!is.finite(fb))){
      stop("At least one training fold produced nonfinite frozen target weights.")
    }
    cab <- ecocox_validate(spb, gb, fa, method)[, 1]
    cba <- ecocox_validate(spa, ga, fb, method)[, 1]
    bab <- ecocox_curve_best(cab, sr_grid)
    bba <- ecocox_curve_best(cba, sr_grid)
    ccons <- ecocox_min_finite_pair(cab, cba)
    bcons <- ecocox_curve_best(ccons, sr_grid)
    score_sym <- ecocox_min_finite_pair(bab$score, bba$score)

    n_rate <- length(sr_grid)
    null_sym <- rep(NaN, n_sim)
    null_cons <- rep(NaN, n_sim)
    null_ab <- rep(NaN, n_sim)
    null_ba <- rep(NaN, n_sim)
    lex_ab <- rep(0, n_rate)
    lex_ba <- rep(0, n_rate)
    lex_c <- rep(0, n_rate)
    lva_ab <- rep(0, n_rate)
    lva_ba <- rep(0, n_rate)
    lva_c <- rep(0, n_rate)
    if (n_sim > 0){
      ntot <- nc
      isd <- sqrt(max(0, 1 - rho_m^2))
      first <- 1
      while (first <= n_sim){
        last <- min(first + batch_size - 1, n_sim)
        nb <- last - first + 1
        Z <- matrix(stats::rnorm(ntot * nb), ntot, nb)
        if (ntot > 1){
          Z[2:ntot, ] <- isd * Z[2:ntot, , drop = FALSE]
        }
        for (r in seq_len(ntot)[-1]){
          Z[r, ] <- Z[r, ] + rho_m * Z[r - 1, ]
        }
        Z <- sd_full * Z
        sa <- ecocox_interp_columns(raw_a[, 1], Z[ia, , drop = FALSE], A[, 1])
        sb <- ecocox_interp_columns(raw_b[, 1], Z[ib, , drop = FALSE], B[, 1])
        if (any(!is.finite(sa)) || any(!is.finite(sb))){
          stop("The fixed odd/even interpolation map produced nonfinite values.")
        }
        pa <- ecocox_window_spectra(sa, dza, ga$nfft, red)$power
        pb <- ecocox_window_spectra(sb, dzb, gb$nfft, red)$power
        ta <- ecocox_train(pa, ga, method, FALSE)
        wa <- ecocox_normalize(ecocox_group_rms(
            ta$amplitude[ga$group_index, , drop = FALSE], ga$group_index))
        sab <- ecocox_validate(pb, gb, wa, method)
        tb <- ecocox_train(pb, gb, method, FALSE)
        wb <- ecocox_normalize(ecocox_group_rms(
            tb$amplitude[gb$group_index, , drop = FALSE], gb$group_index))
        sba <- ecocox_validate(pa, ga, wb, method)
        scons <- ecocox_min_finite_pair(sab, sba)
        obs_ab <- matrix(cab, n_rate, nb)
        obs_ba <- matrix(cba, n_rate, nb)
        obs_c <- matrix(ccons, n_rate, nb)
        fin <- is.finite(sab) & is.finite(obs_ab)
        lva_ab <- lva_ab + rowSums(fin)
        lex_ab <- lex_ab + rowSums(fin & sab >= obs_ab)
        fin <- is.finite(sba) & is.finite(obs_ba)
        lva_ba <- lva_ba + rowSums(fin)
        lex_ba <- lex_ba + rowSums(fin & sba >= obs_ba)
        fin <- is.finite(scons) & is.finite(obs_c)
        lva_c <- lva_c + rowSums(fin)
        lex_c <- lex_c + rowSums(fin & scons >= obs_c)
        null_cons[first:last] <- ecocox_col_max_finite(scons)
        s_ab <- ecocox_curve_best(sab, sr_grid)$score
        s_ba <- ecocox_curve_best(sba, sr_grid)$score
        null_ab[first:last] <- s_ab
        null_ba[first:last] <- s_ba
        null_sym[first:last] <- ecocox_min_finite_pair(s_ab, s_ba)
        first <- last + 1
      }
      if (any(!is.finite(null_sym)) || any(!is.finite(null_cons)) ||
          any(!is.finite(null_ab)) || any(!is.finite(null_ba))){
        stop("A Monte Carlo realization produced a nonfinite statistic.")
      }
    }
    p_ab <- ecocox_mc_p(null_ab, bab$score)
    p_ba <- ecocox_mc_p(null_ba, bba$score)
    p_sym <- ecocox_mc_p(null_sym, score_sym)
    p_cons <- ecocox_mc_p(null_cons, bcons$score)
    pc_ab <- ecocox_max_p_curve(null_ab, cab)
    pc_ba <- ecocox_max_p_curve(null_ba, cba)
    pc_c <- ecocox_max_p_curve(null_cons, ccons)
    pl_ab <- ecocox_point_p_curve(lex_ab, lva_ab, cab)
    pl_ba <- ecocox_point_p_curve(lex_ba, lva_ba, cba)
    pl_c <- ecocox_point_p_curve(lex_c, lva_c, ccons)
    ma <- gb$resolved & matrix(ecocox_active_flags(gb, fa), n_rate, 9,
        byrow = TRUE)
    mb <- ga$resolved & matrix(ecocox_active_flags(ga, fb), n_rate, 9,
        byrow = TRUE)
    shared <- c(max(ga$all_nine_range[1], gb$all_nine_range[1]),
        min(ga$all_nine_range[2], gb$all_nine_range[2]))
    list(rho_forward = cab, rho_backward = cba, rho_consensus = ccons,
        p_global_forward = pc_ab, p_global_backward = pc_ba,
        p_global_consensus = pc_c,
        p_local_forward = pl_ab, p_local_backward = pl_ba,
        p_local_consensus = pl_c,
        p_directional_max = ecocox_max_finite_pair(pc_ab, pc_ba),
        norbit_forward = rowSums(ma), norbit_backward = rowSums(mb),
        norbit_consensus = rowSums(ma & mb),
        best_forward = bab, best_backward = bba, best_consensus = bcons,
        p_forward = p_ab, p_backward = p_ba, p_sym = p_sym,
        p_consensus = p_cons, p_robust = ecocox_max_finite_pair(p_ab, p_ba),
        spacing_odd = dza, spacing_even = dzb,
        raw_count_odd = nrow(raw_a), raw_count_even = nrow(raw_b),
        regular_count_odd = nrow(A), regular_count_even = nrow(B),
        all_nine_odd = ga$all_nine_range, all_nine_even = gb$all_nine_range,
        all_nine_shared = shared, rho_m = rho_m, rho_m_odd = rho_ma,
        rho_m_even = rho_mb, degraded = degraded, completeness_odd = comp_a,
        completeness_even = comp_b, weights_odd = fa[, 1],
        weights_even = fb[, 1], best_rate_train_odd = tra$best_rate,
        best_rate_train_even = trb$best_rate)
  }

  ecocox_group_rms <- function(amp9, gindex){
    # raw group amplitudes of groupAmplitudes (cvcoco.m): RMS within group
    amp9 <- as.matrix(amp9)
    out <- matrix(NaN, 4, ncol(amp9))
    for (g in 1:4){
      a <- amp9[gindex == g, , drop = FALSE]
      out[g, ] <- sqrt(colMeans(a^2))
    }
    out
  }

  ecocox_active_flags <- function(geom, weights){
    ag <- as.numeric(weights) > 1e-8
    ag[is.na(ag)] <- FALSE
    ag[geom$group_index]
  }

  ecocox_mc_p <- function(null, observed){
    null <- null[is.finite(null)]
    if (length(null) < 1 || !is.finite(observed)){
      return(NaN)
    }
    (sum(null >= observed) + 1) / (length(null) + 1)
  }

  ecocox_max_p_curve <- function(null, curve){
    out <- rep(NaN, length(curve))
    null <- null[is.finite(null)]
    if (length(null) < 1){
      return(out)
    }
    for (i in which(is.finite(curve))){
      out[i] <- (1 + sum(null >= curve[i])) / (length(null) + 1)
    }
    out
  }

  ecocox_point_p_curve <- function(exceed, valid, curve){
    out <- rep(NaN, length(curve))
    ok <- is.finite(curve) & valid > 0
    out[ok] <- (exceed[ok] + 1) / (valid[ok] + 1)
    out
  }

  ecocox_physical_windows_interleaved <- function(depth, window, step_depth){
    # buildPhysicalDepthWindows in ecocoInterleavedCore.m
    nd <- length(depth)
    rs <- depth[1]
    re <- depth[nd]
    tol <- ecocox_sampling_tolerance(min(window, step_depth))
    fc <- rs + window / 2
    lc <- re - window / 2
    avail <- max(0, lc - fc)
    nw <- floor((avail + tol) / step_depth) + 1
    if (nw > 10000){
      stop("The requested physical-depth step would create more than 10000 windows.")
    }
    centers <- fc + (0:(nw - 1)) * step_depth
    lower <- centers - window / 2
    upper <- centers + window / 2
    lower[abs(lower - rs) <= tol] <- rs
    upper[abs(upper - re) <= tol] <- re
    starts <- rep(NA_real_, nw)
    ends <- rep(NA_real_, nw)
    count <- rep(0, nw)
    oc <- rep(NaN, nw)
    os <- rep(NaN, nw)
    li <- 1
    ri <- 0
    for (w in seq_len(nw)){
      while (li <= nd && depth[li] < lower[w] - tol){
        li <- li + 1
      }
      ri <- max(ri, li - 1)
      while (ri < nd && depth[ri + 1] <= upper[w] + tol){
        ri <- ri + 1
      }
      if (li <= nd && li <= ri && depth[li] <= upper[w] + tol){
        starts[w] <- li
        ends[w] <- ri
        count[w] <- ri - li + 1
        oc[w] <- 0.5 * (depth[li] + depth[ri])
        os[w] <- depth[ri] - depth[li]
      }
    }
    list(starts = starts, ends = ends, centers = centers, lower = lower,
        upper = upper, count = count, observed_center = oc, observed_span = os)
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

  # Small helpers
  .eps <- .Machine$double.eps

  # Main code
  method <- match.arg(method, c("Pearson", "Spearman"))
  window_mode <- match.arg(window_mode, c("physical-depth", "legacy-count"))
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
  cl <- ecocox_clean_unique(data)
  dat <- cl$data
  nd <- nrow(dat)
  if (nd < 16){
    stop("At least 16 unique depths are required.")
  }
  depth <- dat[, 1]
  med_dt <- stats::median(diff(depth))
  if (is.null(dt)){
    dt <- med_dt
  }
  if (is.null(window)){
    window <- max(1e-6, 0.25 * abs(depth[nd] - depth[1]))
  }
  if (is.null(step_depth)){
    step_depth <- dt
    if (nd > 300){
      step_depth <- dt * ceiling(nd / 300)
    }
  }
  if (is.null(pad)){
    pad <- ecocox_default_pad(nd)
  }
  if (is.null(sr)){
    sr <- ecocox_default_sr(nd, dt, orbit9)
  }
  sr <- as.numeric(sr)
  if (any(diff(sr) <= 0) || any(sr <= 0)){
    stop("sr must be positive and strictly increasing")
  }
  if (length(sr) > 2){
    spc <- diff(sr)
    if (any(abs(spc - stats::median(spc)) > 256 * ecocox_eps(max(1, max(abs(sr)))))){
      stop("sr must have a constant step")
    }
  }
  if (window_mode == "physical-depth"){
    if (window > depth[nd] - depth[1] + ecocox_sampling_tolerance(window)){
      stop("window exceeds the cleaned record span")
    }
    pw <- ecocox_physical_windows_interleaved(depth, window, step_depth)
    starts <- pw$starts
    ends <- pw$ends
    centers <- pw$centers
    count <- pw$count
    lower <- pw$lower
    upper <- pw$upper
  } else {
    nwp <- 2 * ecocox_round((window / dt + 1) / 2)
    if (nwp < 16){
      stop("window must contain at least 16 points")
    }
    if (nwp > nd){
      stop("window is longer than the cleaned data")
    }
    starts <- seq(1, nd - nwp + 1, by = step)
    ends <- starts + nwp - 1
    count <- ends - starts + 1
    lower <- depth[starts]
    upper <- depth[ends]
    centers <- 0.5 * (lower + upper)
  }
  n_w <- length(starts)
  n_rate <- length(sr)
  phase <- (starts - 1) %% 2
  wseed <- rep(NaN, n_w)
  for (w in seq_len(n_w)){
    if (is.finite(starts[w]) && is.finite(ends[w])){
      ident <- 104729 * (starts[w] - 1) + 130363 * (ends[w] - 1)
    } else {
      ident <- 32452843 * (w - 1)
    }
    parent <- 0
    if (!is.null(seed)){
      parent <- seed
    }
    wseed[w] <- parent + ident
    wseed[w] <- wseed[w] %% 2^32
  }
  M <- list()
  keys <- c("rho_forward", "rho_backward", "rho_consensus",
      "p_global_forward", "p_global_backward", "p_global_consensus",
      "p_local_forward", "p_local_backward", "p_local_consensus",
      "p_directional_max", "norbit_forward", "norbit_backward",
      "norbit_consensus")
  for (k in keys){
    M[[k]] <- matrix(NaN, n_rate, n_w)
  }
  wtab <- data.frame(center = centers, lower = lower, upper = upper,
      start_index = starts, end_index = ends,
      point_count = count, phase = phase, seed = wseed,
      success = FALSE, failure = "",
      best_rate_forward = NaN, best_rate_backward = NaN,
      best_rate_consensus = NaN, best_rho_forward = NaN,
      best_rho_backward = NaN, best_rho_consensus = NaN,
      best_index_consensus = NA_integer_,
      p_odd_to_even = NaN, p_even_to_odd = NaN,
      p_robust = NaN, p_sym = NaN, p_consensus = NaN,
      spacing_odd = NaN, spacing_even = NaN,
      raw_count_odd = NaN, raw_count_even = NaN,
      regular_count_odd = NaN, regular_count_even = NaN,
      all_nine_low = NaN, all_nine_high = NaN,
      rho_m = NaN, partial_orbit = FALSE,
      stringsAsFactors = FALSE)
  pb <- NULL
  if (verbose){
    message(sprintf("Interleaved eCOCO: %d windows, %d rates", n_w, n_rate))
    pb <- utils::txtProgressBar(min = 0, max = n_w, style = 3)
  }
  for (w in seq_len(n_w)){
    if (!is.finite(starts[w]) || !is.finite(ends[w]) || count[w] < 16){
      wtab$failure[w] <- "fewer than 16 raw observations"
      next
    }
    if (n_sim > 0 && !is.null(seed)){
      set.seed(ecocox_seed(wseed[w]))
    }
    wr <- tryCatch(
        ecocox_interleaved_window(dat[starts[w]:ends[w], , drop = FALSE],
            orbit9, pad, sr, red, n_sim, method,
            max_frequency, phase[w], batch_size),
        ecocox_unresolvable = function(e){
          e
        })
    if (inherits(wr, "ecocox_unresolvable")){
      wtab$failure[w] <- conditionMessage(wr)
      next
    }
    M$rho_forward[, w] <- wr$rho_forward
    M$rho_backward[, w] <- wr$rho_backward
    M$rho_consensus[, w] <- wr$rho_consensus
    M$p_global_forward[, w] <- wr$p_global_forward
    M$p_global_backward[, w] <- wr$p_global_backward
    M$p_global_consensus[, w] <- wr$p_global_consensus
    M$p_local_forward[, w] <- wr$p_local_forward
    M$p_local_backward[, w] <- wr$p_local_backward
    M$p_local_consensus[, w] <- wr$p_local_consensus
    M$p_directional_max[, w] <- wr$p_directional_max
    M$norbit_forward[, w] <- wr$norbit_forward
    M$norbit_backward[, w] <- wr$norbit_backward
    M$norbit_consensus[, w] <- wr$norbit_consensus
    wtab$success[w] <- TRUE
    wtab$best_rate_forward[w] <- wr$best_forward$rate
    wtab$best_rate_backward[w] <- wr$best_backward$rate
    wtab$best_rate_consensus[w] <- wr$best_consensus$rate
    wtab$best_rho_forward[w] <- wr$best_forward$score
    wtab$best_rho_backward[w] <- wr$best_backward$score
    wtab$best_rho_consensus[w] <- wr$best_consensus$score
    wtab$best_index_consensus[w] <- wr$best_consensus$index
    wtab$p_odd_to_even[w] <- wr$p_forward
    wtab$p_even_to_odd[w] <- wr$p_backward
    wtab$p_robust[w] <- wr$p_robust
    wtab$p_sym[w] <- wr$p_sym
    wtab$p_consensus[w] <- wr$p_consensus
    wtab$spacing_odd[w] <- wr$spacing_odd
    wtab$spacing_even[w] <- wr$spacing_even
    wtab$raw_count_odd[w] <- wr$raw_count_odd
    wtab$raw_count_even[w] <- wr$raw_count_even
    wtab$regular_count_odd[w] <- wr$regular_count_odd
    wtab$regular_count_even[w] <- wr$regular_count_even
    wtab$all_nine_low[w] <- wr$all_nine_shared[1]
    wtab$all_nine_high[w] <- wr$all_nine_shared[2]
    wtab$rho_m[w] <- wr$rho_m
    wtab$partial_orbit[w] <- wr$degraded
    if (!is.null(pb)){
      utils::setTxtProgressBar(pb, w)
    }
  }
  if (!is.null(pb)){
    close(pb)
  }
  if (!any(wtab$success)){
    stop("Every sliding window was frequency-unresolvable for the rate grid.")
  }
  if (any(wtab$partial_orbit)){
    warning(sprintf(paste0("Interleaved eCOCO used partial-orbit training in ",
        "%d of %d successful windows (exploratory)."),
        sum(wtab$partial_orbit), sum(wtab$success)))
  }
  wtab$best_rate_all_nine_resolved <- wtab$success &
      is.finite(wtab$best_rate_consensus) &
      wtab$best_rate_consensus > wtab$all_nine_low &
      wtab$best_rate_consensus <= wtab$all_nine_high
  succ <- wtab$success
  forward <- ecocox_direction(M$rho_forward, M$p_global_forward,
      M$p_local_forward, M$norbit_forward, sr, succ)
  backward <- ecocox_direction(M$rho_backward, M$p_global_backward,
      M$p_local_backward, M$norbit_backward, sr, succ)
  consensus <- ecocox_direction(M$rho_consensus, M$p_global_consensus,
      M$p_local_consensus, M$norbit_consensus, sr,
      succ)
  consensus$best_rate <- wtab$best_rate_consensus
  consensus$best_index <- wtab$best_index_consensus
  consensus$best_rho <- wtab$best_rho_consensus
  consensus$best_p_global <- ecocox_values_at(M$p_global_consensus,
      wtab$best_index_consensus)
  forward$best_rate <- wtab$best_rate_forward
  forward$best_rho <- wtab$best_rho_forward
  backward$best_rate <- wtab$best_rate_backward
  backward$best_rho <- wtab$best_rho_backward
  consensus$p_directional_max <- M$p_directional_max
  pl_strict <- ecocox_max_finite_pair(M$p_local_forward, M$p_local_backward)
  strict <- ecocox_direction(M$rho_consensus, M$p_directional_max, pl_strict,
      M$norbit_consensus, sr, succ)
  strict$best_rate <- consensus$best_rate
  strict$best_index <- consensus$best_index
  strict$best_rho <- consensus$best_rho
  strict$best_p_global <- ecocox_values_at(M$p_directional_max,
      consensus$best_index)
  res <- list(method = "Interleaved eCOCO", sr = sr, depth = centers,
      rho = M$rho_consensus, p_global = M$p_global_consensus,
      p_local = M$p_local_consensus, norbit = M$norbit_consensus,
      pcoco = consensus$pcoco, score = consensus$score,
      forward = forward, backward = backward, consensus = consensus,
      strict = strict, windows = wtab, rho_m = wtab$rho_m,
      window = window, step_depth = step_depth, step = step,
      window_mode = window_mode, sampling_dt = med_dt,
      n_input = cl$n_input, n_finite = cl$n_finite,
      n_duplicate = cl$n_duplicate, degraded = any(wtab$partial_orbit),
      n_sim = n_sim, p_floor = ecocox_p_floor(n_sim),
      seed = seed, red = red, pad = pad,
      max_frequency = max_frequency, orbit9 = orbit9,
      correlation = method)
  res <- ecocox_finish(res, sr, centers, verbose)
  class(res) <- "ecoco_cv"
  if (genplot){
    plot_ecoco_cv(res)
  }
  res
}
