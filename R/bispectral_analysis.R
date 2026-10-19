#' @title Bispectral analysis: bispectrum, biphase and squared bicoherence
#'
#' @description The \code{\link{bispectral_analysis}} function estimates the
#' complex bispectrum \eqn{B(f_1,f_2) = mean[X(f_1) X(f_2) conj(X(f_1+f_2))]},
#' the biphase \eqn{arg B} and the Kim and Powers (1979) magnitude-squared
#' bicoherence (bounded in [0,1]) of an evenly or unevenly sampled series.
#' Two estimators are available: a direct FFT bispectrum averaged over
#' tapered (WOSA/Welch) segments, and a two-dimensional hexagonal
#' frequency-smoothed direct estimate on the native DFT grid. The numerical
#' domain is the non-redundant positive principal domain
#' (\eqn{k_2 \le k_1}, \eqn{k_1+k_2} below Nyquist).
#' Before estimation the series is cleaned (non-finite rows removed, sorted,
#' duplicate coordinates averaged), regularised by interpolation when the
#' sampling is irregular, detrended and optionally standardised.
#' Optional inference: an analytical Beta(1,R-1) pointwise reference, pointwise
#' Monte Carlo surrogate thresholds, or a global maximum-statistic surrogate
#' test (family-wise error control over the full computed principal domain)
#' using Fourier-transform (FT) phase-randomised or IAAFT surrogates
#' (Schreiber and Schmitz, 1996) with plus-one Monte Carlo p-values.
#' High bicoherence indicates repeatable quadratic phase coupling at
#' \eqn{(f_1,f_2,f_1+f_2)}; it is not by itself evidence of causality or energy
#' transfer.
#' Translated from Acycle v3.1 bispectralAnalyze.m, bispectralDefaults.m,
#' bispectralPreprocess.m, bispectralEstimate.m and bispectralSignificance.m
#' (Mingsong Li).
#'
#' @param data Input data set: a matrix or data frame with 2 columns, the first
#' being depth/time and the second the proxy value.
#' @param estimator Estimator: "wosa" (segment averaged) or
#' "frequency-smoothed". \code{NULL} chooses as Acycle does: "wosa", unless
#' three 32-point segments do not fit, then "frequency-smoothed".
#' \code{Default=NULL}
#' @param n_segments Number of WOSA segments (>= 3). \code{NULL} uses
#' min(8, maximum number of 32-point segments) as in Acycle. \code{Default=NULL}
#' @param overlap Segment overlap in percent, in [0,90). \code{Default=50}
#' @param window Taper: "hann", "hamming", "blackman" or "rectangular".
#' \code{Default="hann"}
#' @param segment_detrend Detrending applied to each segment (or the whole
#' record for the frequency-smoothed estimator): "none", "mean" or "linear".
#' \code{Default="mean"}
#' @param smoothing_span Odd span (>= 3) of the hexagonal frequency-smoothing
#' kernel (3 gives 7 triads). \code{Default=3}
#' @param smoothing_kernel "daniell" (uniform) or "cosine" weights.
#' \code{Default="daniell"}
#' @param detrend Whole-record detrending: "none", "mean", "linear" or
#' "polynomial". \code{Default="linear"}
#' @param poly_order Polynomial order (1-10) for \code{detrend="polynomial"}.
#' \code{Default=2}
#' @param standardize Standardise the series to zero mean and unit standard
#' deviation after detrending. \code{Default=TRUE}
#' @param interpolate "auto" (interpolate only irregular series), "always" or
#' "never". \code{Default="auto"}
#' @param interp_method Interpolation method: "linear", "pchip" or "makima".
#' \code{Default="linear"}
#' @param sample_interval Sample interval of the regular grid; \code{NULL} uses
#' the median original spacing. May not be coarser than the median spacing.
#' \code{Default=NULL}
#' @param irregular_tol Maximum relative spacing departure treated as regular
#' (input_policy "prepare"). \code{Default=0.01}
#' @param gap_warning_factor Largest gap (in sample intervals) above which a
#' warning is stored. \code{Default=5}
#' @param input_policy "prepare" (clean the input) or "strict" (refuse
#' non-finite, unsorted or duplicated rows). \code{Default="prepare"}
#' @param nfft FFT length for the WOSA estimator; \code{NULL} uses
#' round(segment length x zero_pad). Ignored by the frequency-smoothed
#' estimator. \code{Default=NULL}
#' @param zero_pad Zero-padding factor for the WOSA estimator. \code{Default=1}
#' @param max_freq_bins Maximum number of computed frequency-axis bins (>= 16);
#' larger grids are decimated with a constant stride. \code{Default=512}
#' @param significance Inference: "none", "analytical", "pointwise"
#' (pointwise surrogate thresholds) or "global" (maximum-statistic surrogate
#' test). \code{Default="none"}
#' @param conf_level Confidence level in (0.5,1). \code{Default=0.95}
#' @param n_sim Number of accepted surrogates (19-99999). \code{Default=999}
#' @param surrogate_type "iaaft" or "ft" (phase randomisation).
#' \code{Default="iaaft"}
#' @param iaaft_iter Maximum number of IAAFT iterations. \code{Default=200}
#' @param iaaft_tol Maximum accepted relative Fourier-amplitude error of an
#' IAAFT surrogate (non-DC, non-Nyquist bins); worse surrogates are rejected
#' and replaced. \code{Default=0.02}
#' @param max_attempts Maximum number of surrogate attempts; \code{NULL} allows
#' n_sim + max(10, ceiling(0.1 n_sim)). \code{Default=NULL}
#' @param seed Seed for the random number generator (\code{NULL} leaves the
#' generator untouched). \code{Default=NULL}
#' @param genplot Plot the result with \code{\link{plot_bispectral}}.
#' \code{Default=FALSE}
#' @param verbose Print progress (progress bar for the surrogates).
#' \code{Default=FALSE}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' (Acycle bispectral module).
#'
#' @references
#' Kim, Y.C., Powers, E.J., 1979. Digital bispectral analysis and its
#' applications to nonlinear wave interactions. IEEE Transactions on Plasma
#' Science 7, 120-131.
#'
#' Choudhury, M.A.A.S., Shah, S.L., Thornhill, N.F., 2008. Diagnosis of
#' Process Nonlinearities and Valve Stiction. Springer.
#'
#' Schreiber, T., Schmitz, A., 1996. Improved surrogate data for nonlinearity
#' tests. Physical Review Letters 77, 635-638.
#'
#' Li, M., Hinnov, L., Kump, L., 2019. Acycle: Time-series analysis software
#' for paleoclimate research and education. Computers & Geosciences 127,
#' 12-22.
#'
#' @examples
#' \donttest{
#' # quadratic phase coupling: f3 = f1 + f2 with phase3 = phase1 + phase2 + 0.4
#' set.seed(1)
#' seg_len <- 256
#' n_seg <- 16
#' t <- 0:(seg_len * n_seg - 1)
#' y <- numeric(length(t))
#' f1 <- 30 / seg_len
#' f2 <- 17 / seg_len
#' for (s in seq_len(n_seg)){
#'   idx <- (s - 1) * seg_len + seq_len(seg_len)
#'   p1 <- 2 * pi * runif(1)
#'   p2 <- 2 * pi * runif(1)
#'   y[idx] <- cos(2 * pi * f1 * t[idx] + p1) +
#'     0.85 * cos(2 * pi * f2 * t[idx] + p2) +
#'     0.65 * cos(2 * pi * (f1 + f2) * t[idx] + p1 + p2 + 0.4) +
#'     0.18 * rnorm(seg_len)
#' }
#' res <- bispectral_analysis(data = cbind(t, y), n_segments = 16,
#'   overlap = 0, detrend = "none", max_freq_bins = 256,
#'   significance = "global", n_sim = 19, seed = 1, genplot = TRUE)
#' }
#' @return
#' A list with elements
#' \code{frequency} (computed frequency axis, cycles per unit),
#' \code{bispectrum} (complex matrix; rows = f2, columns = f1; NA outside the
#' principal domain), \code{bispectrum_magnitude}, \code{biphase} (radians,
#' NA where undefined), \code{bicoherence_squared}, \code{denominator},
#' \code{valid_mask}, \code{principal_domain_mask}, \code{power} (segment-
#' averaged power on the frequency axis), \code{pair_index},
#' \code{pair_bicoherence_squared}, \code{meta} (estimator settings and
#' diagnostics), \code{processed_data} (2-column processed series),
#' \code{preprocessing} (list), \code{significance} (list with
#' \code{method}, \code{threshold}, \code{threshold_matrix}, \code{p_value},
#' \code{significant_mask}, \code{surrogate_maxima},
#' \code{surrogate_spectral_errors}, attempt/rejection counts and warnings),
#' \code{significant_mask} and \code{warnings}.
#'
#' @param standardise British spelling of \code{standardize}; if given it overrides \code{standardize}.
#' @export
#' @importFrom stats fft
#' @importFrom stats median
#' @importFrom stats sd
#' @importFrom stats approx
#' @importFrom stats runif
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar

bispectral_analysis <- function(data = NULL,
                                estimator = NULL,
                                n_segments = NULL,
                                overlap = 50,
                                window = "hann",
                                segment_detrend = "mean",
                                smoothing_span = 3,
                                smoothing_kernel = "daniell",
                                detrend = "linear",
                                poly_order = 2,
                                standardize = TRUE,
                                interpolate = "auto",
                                interp_method = "linear",
                                sample_interval = NULL,
                                irregular_tol = 0.01,
                                gap_warning_factor = 5,
                                input_policy = "prepare",
                                nfft = NULL,
                                zero_pad = 1,
                                max_freq_bins = 512,
                                significance = "none",
                                conf_level = 0.95,
                                n_sim = 999,
                                surrogate_type = "iaaft",
                                iaaft_iter = 200,
                                iaaft_tol = 0.02,
                                max_attempts = NULL,
                                seed = NULL,
                                genplot = FALSE,
                                verbose = FALSE,
  standardise = NULL){
  if (!is.null(standardise)) standardize <- standardise  # British spelling synonym
  # ---- internal helper functions (local to this function) ----

  acy_detrend_linear <- function(x, t = NULL){
    x <- as.numeric(x)
    n <- length(x)
    if (n < 2){
      return(x - mean(x))
    }
    if (is.null(t)){
      t <- seq_len(n)
    }
    X <- cbind(1, t)
    beta <- qr.solve(X, x)
    x - as.numeric(X %*% beta)
  }

  # ---------------------------------------------------------------------------
  # MATLAB interp1(...,'pchip')
  # ---------------------------------------------------------------------------

  acy_pchip_slopes <- function(x, y){
    n <- length(x)
    h <- diff(x)
    del <- diff(y) / h
    d <- numeric(n)
    if (n == 2){
      d[] <- del[1]
      return(d)
    }
    for (k in 2:(n - 1)){
      if (del[k - 1] * del[k] > 0){
        w1 <- 2 * h[k] + h[k - 1]
        w2 <- h[k] + 2 * h[k - 1]
        d[k] <- (w1 + w2) / (w1 / del[k - 1] + w2 / del[k])
      } else {
        d[k] <- 0
      }
    }
    d[1] <- acy_pchip_end(h[1], h[2], del[1], del[2])
    d[n] <- acy_pchip_end(h[n - 1], h[n - 2], del[n - 1], del[n - 2])
    d
  }

  acy_pchip_end <- function(h1, h2, del1, del2){
    d <- ((2 * h1 + h2) * del1 - h1 * del2) / (h1 + h2)
    if (sign(d) != sign(del1)){
      d <- 0
    } else if ((sign(del1) != sign(del2)) && (abs(d) > abs(3 * del1))){
      d <- 3 * del1
    }
    d
  }

  acy_pchip <- function(x, y, xi){
    x <- as.numeric(x)
    y <- as.numeric(y)
    xi <- as.numeric(xi)
    n <- length(x)
    if (n < 2){
      stop("pchip needs at least two points")
    }
    d <- acy_pchip_slopes(x, y)
    h <- diff(x)
    del <- diff(y) / h
    idx <- findInterval(xi, x, all.inside = TRUE)
    s <- xi - x[idx]
    hk <- h[idx]
    c3 <- (d[idx] - 2 * del[idx] + d[idx + 1]) / hk^2
    c2 <- (3 * del[idx] - 2 * d[idx] - d[idx + 1]) / hk
    yi <- y[idx] + s * (d[idx] + s * (c2 + s * c3))
    yi[xi < x[1] | xi > x[n]] <- NA
    yi
  }

  # ===========================================================================
  # internal helpers
  # ===========================================================================

  # @keywords internal
  # @noRd
  bisp_eps <- function(v){
    # MATLAB eps(v) for finite v
    v <- abs(v)
    if (v < 2^-1022){
      return(2^-1074)
    }
    2^(floor(log2(v)) - 52)
  }

  # @keywords internal
  # @noRd
  bisp_round <- function(x){
    # MATLAB round (half away from zero)
    sign(x) * floor(abs(x) + 0.5)
  }

  # @keywords internal
  # @noRd
  bisp_defaults <- function(data, overlap = 50){
    # bispectralDefaults.m: feasible segment count / estimator for short records
    x <- data[, 1]
    keep <- is.finite(data[, 1]) & is.finite(data[, 2])
    n <- length(unique(x[keep]))
    max_seg <- floor(1 + (n / 32 - 1) / (1 - overlap / 100))
    if (max_seg >= 3){
      list(estimator = "wosa", n_segments = min(8, max_seg))
    } else {
      list(estimator = "frequency-smoothed", n_segments = 8)
    }
  }

  # @keywords internal
  # @noRd
  bisp_canonical_estimator <- function(value){
    value <- tolower(trimws(value))
    if (value %in% c("wosa", "welch", "welch/wosa", "segmented")){
      return("wosa")
    }
    if (value %in% c("frequency-smoothed", "frequency smoothed", "smooth",
                     "direct-smoothed")){
      return("frequency-smoothed")
    }
    stop("estimator must be wosa or frequency-smoothed.")
  }

  # @keywords internal
  # @noRd
  bisp_canonical_window <- function(value){
    value <- tolower(trimws(value))
    if (value %in% c("hann", "hanning")){
      return("hann")
    }
    if (value %in% c("hamming", "blackman")){
      return(value)
    }
    if (value %in% c("rectangular", "rectwin", "none")){
      return("rectangular")
    }
    stop("window must be hann, hamming, blackman, or rectangular.")
  }

  # @keywords internal
  # @noRd
  bisp_canonical_segment_detrend <- function(value){
    value <- tolower(trimws(value))
    if (!(value %in% c("none", "mean", "linear"))){
      stop("segment_detrend must be none, mean, or linear.")
    }
    value
  }

  # @keywords internal
  # @noRd
  bisp_canonical_kernel <- function(value){
    value <- tolower(trimws(value))
    if (value %in% c("daniell", "uniform", "boxcar")){
      return("daniell")
    }
    if (value %in% c("cosine", "raised-cosine")){
      return("cosine")
    }
    stop("smoothing_kernel must be daniell or cosine.")
  }

  # @keywords internal
  # @noRd
  bisp_canonical_surrogate <- function(value){
    value <- tolower(trimws(value))
    if (value %in% c("phase", "ft", "fourier")){
      return("ft")
    }
    if (value == "iaaft"){
      return("iaaft")
    }
    stop("surrogate type must be ft (phase) or iaaft.")
  }

  # @keywords internal
  # @noRd
  bisp_canonical_significance <- function(value){
    value <- tolower(trimws(value))
    if (value %in% c("none", "off", "")){
      return("none")
    }
    if (value %in% c("analytical", "analytic", "beta")){
      return("analytical")
    }
    if (value %in% c("pointwise", "surrogate-pointwise", "pointwise-surrogate")){
      return("pointwise")
    }
    if (value %in% c("global", "surrogate-global", "global-surrogate",
                     "max-surrogate")){
      return("global")
    }
    stop("significance must be none, analytical, pointwise, or global.")
  }

  # @keywords internal
  # @noRd
  bisp_is_int <- function(v){
    length(v) == 1 && is.numeric(v) && is.finite(v) && v == round(v)
  }

  # @keywords internal
  # @noRd
  bisp_validate_options <- function(opts){
    # bispectralAnalyze.m validateOptions + bispectralEstimate.m checks
    if (!(opts$input_policy %in% c("prepare", "strict"))){
      stop("input_policy must be prepare or strict.")
    }
    if (!(bisp_is_int(opts$poly_order) && opts$poly_order >= 0)){
      stop("poly_order must be a finite nonnegative integer.")
    }
    int_pos <- list(n_segments = opts$n_segments, max_freq_bins = opts$max_freq_bins,
                    smoothing_span = opts$smoothing_span, n_sim = opts$n_sim,
                    iaaft_iter = opts$iaaft_iter)
    for (i in seq_len(length(int_pos))){
      if (!(bisp_is_int(int_pos[[i]]) && int_pos[[i]] > 0)){
        stop(paste0(names(int_pos)[i], " must be a finite positive integer."))
      }
    }
    if (!(length(opts$zero_pad) == 1 && is.finite(opts$zero_pad) && opts$zero_pad > 0)){
      stop("zero_pad must be a finite positive scalar.")
    }
    if (!(length(opts$overlap) == 1 && is.finite(opts$overlap) &&
          opts$overlap >= 0 && opts$overlap < 90)){
      stop("overlap must be in the interval [0,90).")
    }
    if (!is.null(opts$sample_interval)){
      if (!(length(opts$sample_interval) == 1 && is.finite(opts$sample_interval) &&
            opts$sample_interval > 0)){
        stop("sample_interval must be NULL or positive.")
      }
    }
    if (!is.null(opts$nfft)){
      if (!(bisp_is_int(opts$nfft) && opts$nfft > 0)){
        stop("nfft must be NULL or a positive integer.")
      }
    }
    if (!(length(opts$iaaft_tol) == 1 && is.finite(opts$iaaft_tol) &&
          opts$iaaft_tol > 0 && opts$iaaft_tol <= 1)){
      stop("iaaft_tol must lie in (0,1].")
    }
    if (!is.null(opts$max_attempts)){
      if (!(bisp_is_int(opts$max_attempts) && opts$max_attempts >= opts$n_sim)){
        stop("max_attempts must be NULL or an integer no smaller than n_sim.")
      }
    }
    if (!(length(opts$irregular_tol) == 1 && is.finite(opts$irregular_tol) &&
          opts$irregular_tol >= 0)){
      stop("irregular_tol must be a finite nonnegative scalar.")
    }
    if (!(length(opts$gap_warning_factor) == 1 && is.finite(opts$gap_warning_factor) &&
          opts$gap_warning_factor > 0)){
      stop("gap_warning_factor must be a finite positive scalar.")
    }
    invisible(TRUE)
  }

  # ---------------------------------------------------------------------------
  # preprocessing (bispectralPreprocess.m)
  # ---------------------------------------------------------------------------

  # @keywords internal
  # @noRd
  bisp_preprocess <- function(data, opts){
    data <- matrix(as.numeric(data[, 1:2]), ncol = 2)
    meta <- list(original_count = nrow(data), warnings = character(0),
                 input_policy = opts$input_policy)
    if (opts$input_policy == "strict"){
      if (any(!is.finite(data))){
        stop("Strict input policy does not remove nonfinite rows.")
      }
      if (nrow(data) < 32){
        stop("At least 32 finite observations are required.")
      }
      if (is.unsorted(data[, 1])){
        stop("Strict input policy does not sort coordinates.")
      }
      x <- data[, 1]
      y <- data[, 2]
      if (any(diff(x) == 0)){
        stop("Strict input policy does not combine duplicate coordinates.")
      }
      meta$removed_nonfinite_rows <- 0
      meta$was_sorted <- TRUE
      meta$duplicate_rows_combined <- 0
    } else {
      keep <- is.finite(data[, 1]) & is.finite(data[, 2])
      data <- data[keep, , drop = FALSE]
      meta$removed_nonfinite_rows <- meta$original_count - nrow(data)
      if (nrow(data) < 32){
        stop("At least 32 finite observations are required.")
      }
      was_sorted <- !is.unsorted(data[, 1])
      data <- data[order(data[, 1], method = "radix"), , drop = FALSE]
      x <- data[, 1]
      y <- data[, 2]
      meta$was_sorted <- was_sorted
      xu <- unique(x)
      if (length(xu) < length(x)){
        grp <- match(x, xu)
        sums <- rowsum(y, grp, reorder = TRUE)
        cnt <- rowsum(rep(1, length(y)), grp, reorder = TRUE)
        y <- as.numeric(sums / cnt)
        x <- xu
      }
      meta$duplicate_rows_combined <- nrow(data) - length(x)
      if (length(x) < 32){
        stop("Fewer than 32 distinct coordinates remain after combining duplicates.")
      }
    }

    dx <- diff(x)
    if (any(!is.finite(dx)) || any(dx <= 0)){
      stop("Coordinates must be finite and contain at least 32 distinct values.")
    }
    med_dx <- stats::median(dx)
    if (is.null(opts$sample_interval)){
      dt <- med_dx
      source <- "median spacing"
    } else {
      dt <- as.numeric(opts$sample_interval)
      source <- "user input"
    }
    coarse_tol <- max(32 * bisp_eps(max(1, abs(med_dx))), 1e-10 * med_dx)
    if (dt > med_dx + coarse_tol){
      stop(paste0("sample_interval is coarser than the original median spacing. ",
                  "Interpolation is not an anti-aliasing filter: low-pass filter ",
                  "and resample the series first."))
    }
    max_err <- max(abs(dx - dt))
    rel_dep <- max_err / dt
    gap <- max(dx) / dt
    mode <- opts$interpolate
    if (!(mode %in% c("auto", "always", "never"))){
      stop("interpolate must be auto, always, or never.")
    }
    strict_rel <- 1e-5
    strict_round <- 128 * bisp_eps(max(1, abs(dt)))
    strict_abs <- strict_rel * dt
    if (opts$input_policy == "strict"){
      active_tol <- strict_rel
    } else {
      active_tol <- opts$irregular_tol
    }
    irregular <- rel_dep > active_tol
    if (mode == "auto"){
      do_interp <- irregular
    } else if (mode == "always"){
      do_interp <- TRUE
    } else {
      do_interp <- FALSE
    }
    if (irregular && !do_interp){
      stop(paste0("Bispectral FFT estimators require regular sampling. Use ",
                  "interpolate = 'auto' or 'always', or a suitable sample_interval."))
    }

    if (do_interp){
      n_int <- floor((x[length(x)] - x[1]) / dt +
                       100 * bisp_eps(max(1, abs(x[length(x)] / dt))))
      if (n_int + 1 < 32){
        stop("The requested interpolation interval leaves fewer than 32 samples.")
      } else if (n_int + 1 > 1e7){
        stop("The requested interpolation grid exceeds 10 million samples.")
      }
      xr <- x[1] + (0:n_int) * dt
      im <- opts$interp_method
      if (im == "linear"){
        yr <- stats::approx(x, y, xout = xr, method = "linear", rule = 1)$y
      } else if (im == "pchip"){
        yr <- acy_pchip(x, y, xr)
      } else if (im == "makima"){
        yr <- bisp_makima(x, y, xr)
      } else {
        stop("interp_method must be linear, pchip, or makima.")
      }
      keep <- is.finite(yr)
      x <- xr[keep]
      y <- yr[keep]
      if (irregular){
        meta$warnings <- c(meta$warnings, sprintf(paste0(
          "Coordinate-spacing departure %.6g ppm exceeded the %.6g ppm ",
          "threshold. The data were interpolated onto a regular grid with ",
          "sample interval %.9g using %s."), 1e6 * rel_dep, 1e6 * active_tol,
          dt, toupper(im)))
      }
    }
    if (gap > opts$gap_warning_factor){
      meta$warnings <- c(meta$warnings, sprintf(paste0(
        "Largest original gap is %.2f sampling intervals; interpolated values ",
        "across large gaps should be interpreted cautiously."), gap))
    }

    trend <- numeric(length(y))
    method <- opts$detrend
    if (method == "none"){
      trend <- numeric(length(y))
    } else if (method == "mean"){
      trend[] <- mean(y)
      y <- y - trend
    } else if (method == "linear"){
      trend <- y - acy_detrend_linear(y)
      y <- y - trend
    } else if (method == "polynomial"){
      ord <- opts$poly_order
      if (!(bisp_is_int(ord) && ord >= 1 && ord <= 10 && ord < length(y) - 1)){
        stop("poly_order must be an integer from 1 to 10 and below N-1.")
      }
      xs <- (x - mean(x)) / stats::sd(x)
      V <- matrix(1, length(xs), ord + 1)
      for (j in seq_len(ord)){
        V[, ord + 1 - j] <- xs^j
      }
      coef <- qr.solve(V, y)
      trend <- as.numeric(V %*% coef)
      y <- y - trend
    } else {
      stop("detrend must be none, mean, linear, or polynomial.")
    }

    scale <- 1
    center <- mean(y)
    if (opts$standardize == TRUE){
      scale <- stats::sd(y)
      if (!(is.finite(scale) && scale > 10 * bisp_eps(max(1, max(abs(y)))))){
        stop("The series is constant, or nearly constant, after detrending.")
      }
      y <- (y - center) / scale
    }
    processed <- cbind(x, y)
    colnames(processed) <- NULL
    meta$final_count <- nrow(processed)
    meta$original_start <- data[1, 1]
    meta$original_end <- data[nrow(data), 1]
    meta$sample_interval <- dt
    meta$sample_interval_source <- source
    meta$nyquist <- 1 / (2 * dt)
    meta$original_median_spacing <- med_dx
    meta$original_minimum_spacing <- min(dx)
    meta$original_maximum_spacing <- max(dx)
    meta$maximum_spacing_error <- max_err
    meta$relative_spacing_departure <- rel_dep
    meta$largest_gap_factor <- gap
    meta$was_irregular <- irregular
    meta$was_interpolated <- do_interp
    meta$strict_spacing_absolute_tolerance <- strict_abs
    meta$strict_spacing_roundoff_tolerance <- strict_round
    meta$interp_method <- opts$interp_method
    meta$detrend <- method
    meta$trend <- trend
    meta$was_standardized <- opts$standardize
    meta$standardization_center <- center
    meta$standardization_scale <- scale
    list(processed = processed, meta = meta)
  }

  # @keywords internal
  # @noRd
  bisp_makima <- function(x, y, xi){
    # MATLAB makima (modified Akima) interpolation, NA outside [x1, xn]
    n <- length(x)
    h <- diff(x)
    del <- diff(y) / h
    if (n == 2){
      s <- c(del[1], del[1])
    } else {
      d0 <- 2 * del[1] - del[2]
      dm1 <- 2 * d0 - del[1]
      dn <- 2 * del[n - 1] - del[n - 2]
      dn1 <- 2 * dn - del[n - 1]
      de <- c(dm1, d0, del, dn, dn1)
      w <- abs(diff(de)) + abs((de[-length(de)] + de[-1]) / 2)
      w1 <- w[1:n]
      w2 <- w[3:(n + 2)]
      dl1 <- de[2:(n + 1)]
      dl2 <- de[3:(n + 2)]
      w12 <- w1 + w2
      s <- (w2 / w12) * dl1 + (w1 / w12) * dl2
      zero <- w12 == 0
      s[zero] <- (dl1[zero] + dl2[zero]) / 2
    }
    idx <- findInterval(xi, x, all.inside = TRUE)
    t <- xi - x[idx]
    hk <- h[idx]
    c3 <- (s[idx] - 2 * del[idx] + s[idx + 1]) / hk^2
    c2 <- (3 * del[idx] - 2 * s[idx] - s[idx + 1]) / hk
    yi <- y[idx] + t * (s[idx] + t * (c2 + t * c3))
    yi[xi < x[1] | xi > x[n]] <- NA
    yi
  }

  # ---------------------------------------------------------------------------
  # estimator (bispectralEstimate.m)
  # ---------------------------------------------------------------------------

  # @keywords internal
  # @noRd
  bisp_window <- function(n, name){
    k <- 0:(n - 1)
    den <- max(1, n - 1)
    if (name == "hann"){
      w <- 0.5 - 0.5 * cos(2 * pi * k / den)
    } else if (name == "hamming"){
      w <- 0.54 - 0.46 * cos(2 * pi * k / den)
    } else if (name == "blackman"){
      w <- 0.42 - 0.5 * cos(2 * pi * k / den) + 0.08 * cos(4 * pi * k / den)
    } else {
      w <- rep(1, n)
    }
    w
  }

  # @keywords internal
  # @noRd
  bisp_hex_kernel <- function(h, kernel){
    # MATLAB meshgrid(-h:h,-h:h) visited in column-major order
    v <- -h:h
    m <- length(v)
    a <- matrix(v, m, m, byrow = TRUE)
    b <- matrix(v, m, m)
    keep <- abs(a + b) <= h
    offsets <- cbind(a[keep], b[keep])
    if (kernel == "daniell"){
      w <- rep(1, nrow(offsets))
    } else {
      g1 <- 0.5 * (1 + cos(pi * offsets[, 1] / (h + 1)))
      g2 <- 0.5 * (1 + cos(pi * offsets[, 2] / (h + 1)))
      g3 <- 0.5 * (1 + cos(pi * rowSums(offsets) / (h + 1)))
      w <- g1 * g2 * g3
    }
    list(offsets = offsets, weights = w / sum(w))
  }

  # @keywords internal
  # @noRd
  bisp_build_state <- function(n, dt, opts){
    est <- opts$estimator
    if (n < 32){
      stop("At least 32 samples are required.")
    }
    if (est == "wosa"){
      nseg <- opts$n_segments
      if (nseg < 3){
        stop("n_segments must be an integer of at least 3.")
      }
      ov <- opts$overlap / 100
      L <- floor(n / (1 + (nseg - 1) * (1 - ov)))
      if (L < 32){
        stop("This segment count/overlap leaves fewer than 32 samples per segment.")
      }
      starts <- bisp_round(seq(1, n - L + 1, length.out = nseg))
      if (length(unique(starts)) != nseg){
        stop("The requested settings do not produce distinct segments.")
      }
      if (is.null(opts$nfft)){
        nfft <- bisp_round(L * opts$zero_pad)
      } else {
        nfft <- opts$nfft
      }
      if (!(is.finite(nfft) && nfft >= L)){
        stop("nfft must be at least the segment length.")
      }
      win <- bisp_window(L, opts$window)
      offsets <- matrix(0, 0, 2)
      weights <- numeric(0)
      eff <- nseg
      actual_overlap <- 100 * (1 - stats::median(diff(starts)) / L)
      rayleigh <- 1 / (L * dt)
    } else {
      nfft <- n
      L <- n
      starts <- 1
      win <- bisp_window(n, opts$window)
      span <- opts$smoothing_span
      if (!(span >= 3 && span %% 2 == 1)){
        stop("smoothing_span must be an odd integer of at least 3.")
      }
      hk <- bisp_hex_kernel((span - 1) / 2, opts$smoothing_kernel)
      offsets <- hk$offsets
      weights <- hk$weights
      eff <- 1 / sum(weights^2)
      actual_overlap <- NA
      rayleigh <- 1 / (n * dt)
    }
    pos_max <- floor((nfft - 1) / 2)
    min_bin <- 1
    max_bin <- pos_max
    hw <- 0
    if (est == "frequency-smoothed"){
      hw <- max(abs(offsets))
      min_bin <- max(min_bin, hw + 1)
    }
    if (max_bin - min_bin + 1 < 2){
      stop("Too few positive FFT bins remain for bispectral estimation.")
    }
    all_bins <- min_bin:max_bin
    if (opts$max_freq_bins < 16){
      stop("max_freq_bins must be an integer of at least 16.")
    }
    stride <- max(1, ceiling(length(all_bins) / opts$max_freq_bins))
    axis_bins <- all_bins[seq(1, length(all_bins), by = stride)]
    frequency <- axis_bins / (nfft * dt)
    na <- length(axis_bins)
    bin1 <- matrix(axis_bins, na, na, byrow = TRUE)
    bin2 <- matrix(axis_bins, na, na)
    valid <- bin2 <= bin1 & (bin1 + bin2 <= max_bin)
    if (est == "frequency-smoothed"){
      valid <- valid & (bin1 - hw >= 1) & (bin2 - hw >= 1) &
        (bin1 + bin2 + hw <= max_bin)
    }
    if (!any(valid)){
      stop("No triads remain in the selected nonredundant principal domain.")
    }
    pair_index <- which(valid)
    pb1 <- bin1[pair_index]
    pb2 <- bin2[pair_index]

    win <- win / sqrt(mean(win^2))
    enbw <- L * sum(win^2) / (sum(win)^2)
    coverage <- numeric(n)
    for (i in seq_len(length(starts))){
      idx <- starts[i]:(starts[i] + L - 1)
      coverage[idx] <- coverage[idx] + 1
    }
    warn <- character(0)
    if (est == "wosa" && opts$n_segments < 8){
      warn <- c(warn, "Fewer than 8 segments gives a high-variance bicoherence estimate.")
    }
    if (est == "wosa" && opts$overlap > 0 && opts$significance == "analytical"){
      warn <- c(warn, paste0("Overlapped segments are correlated; analytical Beta ",
                             "significance is approximate and retained only for ",
                             "compatibility. Formal inference uses IAAFT ",
                             "maximum-statistic FWER control."))
    } else if (est == "frequency-smoothed"){
      warn <- c(warn, paste0("Frequency-smoothed triads share FFT coefficients; ",
                             "the kernel effective count is diagnostic, not ",
                             "strict degrees of freedom."))
    }
    meta <- list(estimator = est, record_length = n, sample_interval = dt,
                 nyquist = 1 / (2 * dt), nfft = nfft,
                 zero_padding_factor_actual = nfft / L,
                 segment_length = L, segment_count = length(starts),
                 segment_starts = starts, requested_overlap = opts$overlap,
                 actual_median_overlap = actual_overlap, window = opts$window,
                 window_enbw_bins = enbw, rayleigh_resolution = rayleigh,
                 frequency_bin_stride = stride, computed_axis_bin_count = na,
                 native_positive_frequency_maximum = pos_max / (nfft * dt),
                 principal_sum_frequency_maximum = max(pb1 + pb2) / (nfft * dt),
                 triad_count = length(pair_index),
                 effective_realization_count = eff,
                 smoothing_offsets = offsets, smoothing_weights = weights,
                 smoothing_bandwidth = hw * rayleigh,
                 unused_sample_count = sum(coverage == 0),
                 minimum_coverage = min(coverage),
                 maximum_coverage = max(coverage),
                 warnings = warn)
    list(estimator = est, n = n, dt = dt, nfft = nfft, seg_len = L,
         starts = starts, window = win, axis_bins = axis_bins,
         frequency = frequency, valid = valid, pair_index = pair_index,
         pb1 = pb1, pb2 = pb2, offsets = offsets, weights = weights, meta = meta)
  }

  # @keywords internal
  # @noRd
  bisp_segment_detrend <- function(y, method){
    if (method == "mean"){
      y <- y - mean(y)
    } else if (method == "linear"){
      y <- acy_detrend_linear(y)
    }
    y
  }

  # @keywords internal
  # @noRd
  bisp_pad_fft <- function(x, nfft){
    if (nfft > length(x)){
      x <- c(x, numeric(nfft - length(x)))
    }
    stats::fft(x)
  }

  # @keywords internal
  # @noRd
  bisp_estimate_pairs <- function(y, state, opts){
    # returns the pair (triad) vectors only; used for surrogates
    L <- state$seg_len
    if (state$estimator == "wosa"){
      nseg <- length(state$starts)
      spectra <- matrix(0 + 0i, state$nfft, nseg)
      for (i in seq_len(nseg)){
        first <- state$starts[i]
        seg <- y[first:(first + L - 1)]
        seg <- bisp_segment_detrend(seg, opts$segment_detrend)
        spectra[, i] <- bisp_pad_fft(seg * state$window, state$nfft) / L
      }
      a <- spectra[state$pb1 + 1, , drop = FALSE] * spectra[state$pb2 + 1, , drop = FALSE]
      cc <- spectra[state$pb1 + state$pb2 + 1, , drop = FALSE]
      sum_triple <- rowSums(a * Conj(cc))
      sum_a2 <- rowSums(Mod(a)^2)
      sum_c2 <- rowSums(Mod(cc)^2)
      den <- sum_a2 * sum_c2
      bisp <- sum_triple / nseg
      power <- rowMeans(Mod(spectra[state$axis_bins + 1, , drop = FALSE])^2)
    } else {
      series <- bisp_segment_detrend(y, opts$segment_detrend)
      spec <- bisp_pad_fft(series * state$window, state$nfft) / L
      np <- length(state$pb1)
      sum_triple <- complex(np)
      sum_a2 <- numeric(np)
      sum_c2 <- numeric(np)
      for (i in seq_len(nrow(state$offsets))){
        k1 <- state$pb1 + state$offsets[i, 1]
        k2 <- state$pb2 + state$offsets[i, 2]
        k3 <- k1 + k2
        w <- state$weights[i]
        a <- spec[k1 + 1] * spec[k2 + 1]
        cc <- spec[k3 + 1]
        sum_triple <- sum_triple + w * a * Conj(cc)
        sum_a2 <- sum_a2 + w * Mod(a)^2
        sum_c2 <- sum_c2 + w * Mod(cc)^2
      }
      bisp <- sum_triple
      den <- sum_a2 * sum_c2
      power <- Mod(spec[state$axis_bins + 1])^2
    }
    b2 <- rep(NA_real_, length(den))
    ok <- is.finite(den) & den > 0
    b2[ok] <- Mod(sum_triple[ok])^2 / den[ok]
    b2[ok] <- pmin(1, pmax(0, b2[ok]))
    list(bispectrum = bisp, b2 = b2, denominator = den, power = power)
  }

  # @keywords internal
  # @noRd
  bisp_estimate <- function(y, state, opts){
    y <- as.numeric(y)
    if (any(!is.finite(y))){
      stop("The estimator requires a finite preprocessed series.")
    }
    p <- bisp_estimate_pairs(y, state, opts)
    na <- length(state$frequency)
    B <- matrix(complex(real = NA_real_, imaginary = NA_real_), na, na)
    b2 <- matrix(NA_real_, na, na)
    den <- matrix(NA_real_, na, na)
    B[state$pair_index] <- p$bispectrum
    b2[state$pair_index] <- p$b2
    den[state$pair_index] <- p$denominator
    mag <- Mod(B)
    biphase <- Arg(B)
    defined <- state$valid & is.finite(b2) & is.finite(Re(B)) &
      is.finite(Im(B)) & mag > 0
    defined[is.na(defined)] <- FALSE
    biphase[!defined] <- NA
    valid <- state$valid & is.finite(b2)
    meta <- state$meta
    meta$invalid_denominator_triad_count <- sum(state$valid & !is.finite(b2))
    list(frequency = state$frequency,
         bispectrum = B,
         bispectrum_magnitude = mag,
         biphase = biphase,
         bicoherence_squared = b2,
         denominator = den,
         principal_domain_mask = state$valid,
         valid_mask = valid,
         invalid_denominator_mask = state$valid & !is.finite(b2),
         power = p$power,
         pair_index = state$pair_index,
         pair_bicoherence_squared = p$b2,
         meta = meta)
  }

  # ---------------------------------------------------------------------------
  # surrogates (bispectralSurrogate.m)
  # ---------------------------------------------------------------------------

  # @keywords internal
  # @noRd
  bisp_draw <- function(n, type){
    # random numbers consumed by one surrogate (same order as MATLAB code)
    if (type == "ft"){
      stats::runif(max(0, floor((n + 1) / 2) - 1))
    } else {
      sample.int(n)
    }
  }

  # @keywords internal
  # @noRd
  bisp_spectral_error <- function(y, target){
    n <- length(y)
    if (n < 3){
      return(0)
    }
    pos <- 2:floor((n + 1) / 2)
    cand <- Mod(stats::fft(y))
    d <- cand[pos] - target[pos]
    sqrt(sum(d^2)) / max(sqrt(sum(target[pos]^2)), .Machine$double.xmin)
  }

  # @keywords internal
  # @noRd
  bisp_surrogate_core <- function(y, type, iterations, draw){
    n <- length(y)
    if (type == "ft"){
      spec <- stats::fft(y)
      target <- Mod(spec)
      npos <- max(0, floor((n + 1) / 2) - 1)
      pos <- 1 + seq_len(npos)
      if (length(pos) > 0){
        spec[pos] <- Mod(spec[pos]) * exp(1i * 2 * pi * draw)
        spec[n - pos + 2] <- Conj(spec[pos])
      }
      sur <- Re(stats::fft(spec, inverse = TRUE)) / n
      return(list(surrogate = sur, iterations = 1,
                  relative_spectral_error = bisp_spectral_error(sur, target)))
    }
    target <- Mod(stats::fft(y))
    sorted_target <- sort(y)
    sur <- y[draw]
    prev <- Inf
    completed <- 0
    err <- NA_real_
    for (i in seq_len(iterations)){
      ph <- Arg(stats::fft(sur))
      proj <- Re(stats::fft(target * exp(1i * ph), inverse = TRUE)) / n
      rk <- order(proj, method = "radix")
      nxt <- numeric(n)
      nxt[rk] <- sorted_target
      sur <- nxt
      err <- bisp_spectral_error(sur, target)
      completed <- i
      if (i > 1 && abs(prev - err) <= 1e-10 * max(1, prev)){
        break
      }
      prev <- err
    }
    list(surrogate = sur, iterations = completed, relative_spectral_error = err)
  }

  # ---------------------------------------------------------------------------
  # significance (bispectralSignificance.m)
  # ---------------------------------------------------------------------------

  # @keywords internal
  # @noRd
  bisp_critical_rank <- function(m, alpha){
    if (!(bisp_is_int(m) && m >= 1)){
      return(NA_real_)
    }
    allowed <- floor(alpha * (m + 1)) - 1
    allowed <- max(-1, min(m, allowed))
    while (allowed < m && (allowed + 2) / (m + 1) <= alpha){
      allowed <- allowed + 1
    }
    while (allowed >= 0 && (allowed + 1) / (m + 1) > alpha){
      allowed <- allowed - 1
    }
    if (allowed < 0){
      return(NA_real_)
    }
    m - allowed
  }

  # @keywords internal
  # @noRd
  bisp_count_ge <- function(sorted_values, query){
    # number of sorted_values >= query (NA for non-finite queries)
    out <- rep(NA_real_, length(query))
    ok <- is.finite(query)
    out[ok] <- length(sorted_values) -
      findInterval(query[ok], sorted_values, left.open = TRUE)
    out
  }

  # @keywords internal
  # @noRd
  bisp_significance <- function(y, opts, observed, state, verbose = FALSE){
    method <- opts$significance
    na <- length(observed$frequency)
    empty <- matrix(NA_real_, na, na)
    sig <- list(method = method, conf_level = opts$conf_level, threshold = NA_real_,
                threshold_matrix = empty, p_value = empty,
                significant_mask = matrix(FALSE, na, na),
                surrogate_maxima = numeric(0),
                surrogate_spectral_errors = numeric(0),
                attempt_spectral_errors = numeric(0),
                surrogate_type = "none", n_sim = 0, attempt_count = 0,
                rejected_count = 0, rejected_iaaft = 0, rejected_nonfinite = 0,
                max_attempts = 0, family_triad_count = 0,
                warnings = character(0))
    if (method == "none"){
      return(sig)
    }
    conf <- opts$conf_level
    if (!(length(conf) == 1 && is.finite(conf) && conf > 0.5 && conf < 1)){
      stop("conf_level must lie strictly between 0.5 and 1.")
    }
    alpha <- 1 - conf
    obs <- observed$pair_bicoherence_squared
    pidx <- observed$pair_index
    family <- is.finite(obs)
    sig$family_triad_count <- sum(family)

    if (method == "analytical"){
      R <- state$meta$effective_realization_count
      if (!(is.finite(R) && R > 1)){
        stop("At least two effective realizations are needed for the Beta approximation.")
      }
      thr <- 1 - alpha^(1 / (R - 1))
      pp <- (1 - obs)^(R - 1)
      pp[!is.finite(obs)] <- NA
      tm <- empty
      tm[pidx] <- thr
      pv <- empty
      pv[pidx] <- pp
      sm <- matrix(FALSE, na, na)
      s <- obs >= thr
      s[is.na(s)] <- FALSE
      sm[pidx] <- s
      sig$threshold <- thr
      sig$threshold_matrix <- tm
      sig$p_value <- pv
      sig$significant_mask <- sm & observed$valid_mask
      sig$warnings <- paste0("Analytical Beta inference is retained for ",
                             "compatibility and method validation only.")
      return(sig)
    }

    is_global <- method == "global"
    if (is_global && !any(family)){
      stop("Maximum-statistic inference requires at least one finite observed b^2.")
    }
    nsim <- opts$n_sim
    if (!(bisp_is_int(nsim) && nsim >= 19 && nsim <= 99999)){
      stop("Use an integer between 19 and 99999 surrogates (n_sim).")
    }
    if (1 / (nsim + 1) > alpha){
      stop(sprintf("n_sim = %d cannot reach alpha = %.6g; use at least %d surrogates.",
                   nsim, alpha, ceiling(1 / alpha - 1)))
    }
    stype <- opts$surrogate_type
    if (is.null(opts$max_attempts)){
      max_att <- nsim + max(10, ceiling(0.1 * nsim))
    } else {
      max_att <- opts$max_attempts
    }
    np <- length(obs)
    if (!is_global){
      if (8 * np * nsim > 1024^3){
        stop("Pointwise surrogate thresholds would need more than 1 GiB; reduce max_freq_bins or n_sim.")
      }
      sur_values <- matrix(NA_real_, np, nsim)
    }
    sur_max <- rep(NA_real_, nsim)
    sur_err <- rep(NA_real_, nsim)
    att_err <- rep(NA_real_, max_att)
    exceed <- numeric(np)
    finite_count <- numeric(np)
    att <- 0
    acc <- 0
    rej_iaaft <- 0
    rej_nonfinite <- 0
    n <- length(y)
    if (verbose == TRUE){
      pb <- utils::txtProgressBar(min = 0, max = nsim, style = 3)
    }
    while (acc < nsim && att < max_att){
      att <- att + 1
      draw <- bisp_draw(n, stype)
      s <- bisp_surrogate_core(y, stype, opts$iaaft_iter, draw)
      err <- s$relative_spectral_error
      att_err[att] <- err
      if (stype == "iaaft" && !(is.finite(err) && err <= opts$iaaft_tol)){
        rej_iaaft <- rej_iaaft + 1
        next
      }
      vals <- bisp_estimate_pairs(s$surrogate, state, opts)$b2
      if (is_global){
        fv <- vals[family]
        has <- length(fv) > 0 && all(is.finite(fv))
      } else {
        fv <- vals[is.finite(vals)]
        has <- length(fv) > 0
      }
      if (!has){
        rej_nonfinite <- rej_nonfinite + 1
        next
      }
      acc <- acc + 1
      sur_err[acc] <- err
      ge <- vals >= obs
      ge[is.na(ge)] <- FALSE
      exceed <- exceed + ge
      finite_count <- finite_count + is.finite(vals)
      sur_max[acc] <- max(fv)
      if (!is_global){
        sur_values[, acc] <- vals
      }
      if (verbose == TRUE){
        utils::setTxtProgressBar(pb, acc)
      }
    }
    if (verbose == TRUE){
      close(pb)
    }
    if (acc < nsim){
      if (rej_iaaft > 0){
        stop(sprintf(paste0("Only %d of %d surrogates passed after %d attempts; %d IAAFT ",
                            "surrogate(s) exceeded spectral-error tolerance %.4g. ",
                            "Increase iaaft_iter or adjust iaaft_tol/max_attempts."),
                     acc, nsim, att, rej_iaaft, opts$iaaft_tol))
      }
      stop(sprintf("Only %d of %d finite null surrogates were obtained after %d attempts.",
                   acc, nsim, att))
    }

    pp <- (1 + exceed) / (1 + finite_count)
    pp[!is.finite(obs)] <- NA
    pv <- empty
    pv[pidx] <- pp
    sm <- matrix(FALSE, na, na)
    tm <- empty
    if (!is_global){
      sur_values[!is.finite(sur_values)] <- Inf
      o <- order(row(sur_values), sur_values, method = "radix")
      sorted_vals <- matrix(sur_values[o], np, nsim, byrow = TRUE)
      thr_pairs <- rep(NA_real_, np)
      counts <- sort(unique(finite_count[finite_count > 0]))
      for (ci in seq_len(length(counts))){
        rank <- bisp_critical_rank(counts[ci], alpha)
        if (is.na(rank)){
          next
        }
        rows <- which(finite_count == counts[ci])
        thr_pairs[rows] <- sorted_vals[cbind(rows, rep(rank, length(rows)))]
      }
      thr_pairs[!is.finite(obs)] <- NA
      thr_pairs[!is.finite(thr_pairs)] <- NA
      tm[pidx] <- thr_pairs
      s <- obs > thr_pairs & pp <= alpha
      s[is.na(s)] <- FALSE
      sm[pidx] <- s
      sig$threshold_comparison <- paste0("strictly greater than the plus-one Monte ",
                                         "Carlo critical order statistic, with p <= alpha")
      sig$warnings <- paste0("Pointwise surrogate inference does not control the ",
                             "family-wise error rate across the full map.")
    } else {
      fmax <- sort(sur_max[is.finite(sur_max)])
      if (length(fmax) < ceiling(0.9 * nsim)){
        stop("Too many surrogate estimates lacked finite bicoherence values.")
      }
      rank <- bisp_critical_rank(length(fmax), alpha)
      if (is.na(rank)){
        stop("The accepted surrogate count cannot attain the requested significance level.")
      }
      thr <- fmax[rank]
      tm[pidx[family]] <- thr
      gp <- (1 + bisp_count_ge(fmax, obs)) / (length(fmax) + 1)
      gp[!is.finite(obs)] <- NA
      pv[pidx] <- gp
      s <- obs > thr & gp <= alpha
      s[is.na(s)] <- FALSE
      sm[pidx] <- s
      sig$threshold <- thr
      sig$threshold_comparison <- paste0("strictly greater than the plus-one Monte ",
                                         "Carlo critical order statistic, with ",
                                         "global p <= alpha")
      if (stype == "ft"){
        sig$warnings <- paste0("FT phase maximum-statistic inference is retained for ",
                               "compatibility and validation; formal maps use IAAFT.")
      }
    }
    if (rej_iaaft > 0){
      sig$warnings <- c(sig$warnings, sprintf(paste0(
        "Rejected and replaced %d IAAFT surrogate(s) whose relative ",
        "Fourier-amplitude error exceeded %.4g."), rej_iaaft, opts$iaaft_tol))
    }
    if (rej_nonfinite > 0){
      sig$warnings <- c(sig$warnings, sprintf(paste0(
        "Rejected and replaced %d surrogate estimate(s) lacking b^2 values ",
        "required by the fixed inference family."), rej_nonfinite))
    }
    sig$threshold_matrix <- tm
    sig$p_value <- pv
    sig$significant_mask <- sm & observed$valid_mask
    sig$surrogate_maxima <- sur_max
    sig$surrogate_spectral_errors <- sur_err
    sig$attempt_spectral_errors <- att_err[seq_len(att)]
    sig$surrogate_type <- stype
    sig$n_sim <- nsim
    sig$attempt_count <- att
    sig$rejected_count <- rej_iaaft + rej_nonfinite
    sig$rejected_iaaft <- rej_iaaft
    sig$rejected_nonfinite <- rej_nonfinite
    sig$max_attempts <- max_att
    sig
  }

  # ---- main code ----
  if (is.null(data)){
    stop("data must be a 2-column matrix [depth/time, value]")
  }
  data <- as.matrix(data)
  if (!is.numeric(data) || ncol(data) < 2){
    stop("data must be a numeric matrix with at least 2 columns")
  }
  auto <- bisp_defaults(data, overlap)
  if (is.null(estimator)){
    estimator <- auto$estimator
  }
  if (is.null(n_segments)){
    n_segments <- auto$n_segments
  }
  opts <- list(estimator = bisp_canonical_estimator(estimator),
               n_segments = n_segments,
               overlap = overlap,
               window = bisp_canonical_window(window),
               segment_detrend = bisp_canonical_segment_detrend(segment_detrend),
               smoothing_span = smoothing_span,
               smoothing_kernel = bisp_canonical_kernel(smoothing_kernel),
               detrend = tolower(trimws(detrend)),
               poly_order = poly_order,
               standardize = standardize,
               interpolate = tolower(trimws(interpolate)),
               interp_method = tolower(trimws(interp_method)),
               sample_interval = sample_interval,
               irregular_tol = irregular_tol,
               gap_warning_factor = gap_warning_factor,
               input_policy = tolower(trimws(input_policy)),
               nfft = nfft,
               zero_pad = zero_pad,
               max_freq_bins = max_freq_bins,
               significance = bisp_canonical_significance(significance),
               conf_level = conf_level,
               n_sim = n_sim,
               surrogate_type = bisp_canonical_surrogate(surrogate_type),
               iaaft_iter = iaaft_iter,
               iaaft_tol = iaaft_tol,
               max_attempts = max_attempts)
  bisp_validate_options(opts)

  if (verbose == TRUE){
    cat("Preprocessing series\n")
  }
  pre <- bisp_preprocess(data, opts)
  processed <- pre$processed
  dt <- pre$meta$sample_interval

  if (verbose == TRUE){
    cat("Estimating bispectrum and bicoherence\n")
  }
  state <- bisp_build_state(nrow(processed), dt, opts)
  est <- bisp_estimate(processed[, 2], state, opts)

  if (!is.null(seed)){
    set.seed(seed)
  }
  sig <- bisp_significance(processed[, 2], opts, est, state, verbose)

  warn <- unique(c(state$meta$warnings, pre$meta$warnings, sig$warnings))
  if (verbose == TRUE && length(warn) > 0){
    for (i in seq_len(length(warn))){
      cat("Note:", warn[i], "\n")
    }
  }
  res <- est
  res$processed_data <- processed
  res$input_data <- data[, 1:2, drop = FALSE]
  res$preprocessing <- pre$meta
  res$significance <- sig
  res$significant_mask <- sig$significant_mask
  res$warnings <- warn
  res$options <- opts
  if (genplot == TRUE){
    plot_bispectral(res)
  }
  res
}
