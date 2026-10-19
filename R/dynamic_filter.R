#' @title Dynamic (depth-varying) Gaussian band-pass filter
#'
#' @description The \code{\link{dynamic_filter}} function extracts a cycle
#' whose frequency changes along the depth/time axis (e.g. because of a
#' changing sedimentation rate) using a pass band that follows user-defined
#' control points. The method is the dynamic filter of Nicolas Thibault and
#' Giovanni Rizzi (2019), adapted by Mingsong Li (2020) for Acycle's
#' "Frequency Stabilization" tool, in the deterministic form of Acycle's
#' current core: the series is divided into sliding windows; in each window
#' the pass band is obtained by linear interpolation of the control points
#' at the window centre (held constant beyond the first/last control point),
#' the window is mean-centred, tapered by the square root of a
#' half-sample-shifted Hann window, Gaussian band-pass filtered (half-power
#' edges at the lower/upper frequency, zero padded to
#' \code{2^nextpow2(2 * window_samples - 1)}), tapered again and
#' overlap-added. The result is divided by the summed Hann weights.
#' In Acycle the control points are clicked on an evolutionary spectrum;
#' here they are supplied as arguments, e.g. built with
#' \code{\link{dynamic_filter_controls}} from a period curve tracked with
#' WaverideR (\code{track_period_wavelet}, \code{completed_series}).
#' Translated from Acycle v3.1 acycleDynamicFilter.m (strict core; legacy
#' flow dynamic_filter.m / dynamic_filter_lang.m, GUI DynamicFilter.m)
#' (Mingsong Li).
#'
#' @param data Input data set, a matrix or data frame with 2 columns: the
#' first column being depth/time (strictly increasing, regularly sampled to
#' a relative spacing tolerance of 1e-5; interpolate beforehand otherwise)
#' and the second column being the proxy value.
#' @param controls Pass band control points. Either a matrix/data frame with
#' columns \code{depth}, \code{f_low}, \code{f_high} and optionally
#' \code{f_centre} (unnamed columns are taken in that order), or a list with
#' elements \code{lower} and \code{upper}, each a 2-column matrix
#' [depth, frequency] (Acycle's separate lower/upper boundary picks, which
#' may be at different depths). Depths must be strictly increasing and lie
#' within the data range. A single row gives a constant pass band. Because
#' the Acycle Gaussian response is symmetric, a supplied \code{f_centre}
#' must equal \code{(f_low + f_high) / 2} (it is only checked).
#' @param window_length Window length in depth units (first-to-last sample
#' span). \code{Default=NULL}: 20 percent of the record length (Acycle GUI
#' default). A window of about two times the longest period of the target
#' cycle is recommended.
#' @param step_length Step between successive windows in depth units.
#' \code{Default=NULL}: 4 times the sample interval, or
#' \code{(record length - window_length) / 500} when that gives more than
#' 500 steps (Acycle GUI default).
#' @param genplot Generate a plot using \code{\link{plot_dynamic_filter}}
#' \code{Default=FALSE}.
#' @param x_lab label for the x-axis \code{Default="depth"}.
#' @param y_lab label for the y-axis \code{Default="proxy"}.
#' @param verbose Show a progress bar \code{Default=FALSE}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Nicolas
#' Thibault and Giovanni Rizzi (2019), adapted by Mingsong Li (Acycle).
#'
#' @references
#' Li, M., Hinnov, L., Kump, L., 2019. Acycle: Time-series analysis software
#' for paleoclimate research and education. Computers & Geosciences 127,
#' 12-22. <doi:10.1016/j.cageo.2019.02.011>
#'
#' Thibault, N., Rizzi, G. (2019). Dynamic filtering routine, distributed
#' with Acycle (Frequency Stabilization).
#'
#' @examples
#' \donttest{
#' # cycle whose period grows linearly from 5 to 15 m (changing sed. rate)
#' depth <- seq(0, 300, by = 0.25)
#' period <- 5 + 10 * depth / 300
#' phase <- 2 * pi * cumsum(1 / period) * 0.25
#' set.seed(1)
#' proxy <- sin(phase) + 0.8 * sin(2 * pi * depth / 2.1) + rnorm(length(depth), sd = 0.3)
#' ctrl <- dynamic_filter_controls(tracked = cbind(depth, period),
#'   rel_bandwidth = 0.3, n_points = 30)
#' res <- dynamic_filter(data = cbind(depth, proxy), controls = ctrl,
#'   window_length = 40, step_length = 1, genplot = FALSE)
#' plot(depth, sin(phase), type = "l")
#' lines(res$filtered[, 1], res$filtered[, 2], col = "red")
#' }
#'
#' @return A list with elements
#' \describe{
#' \item{filtered}{matrix with columns \code{depth}, \code{filtered},
#' \code{coverage_count} (number of windows covering each sample) and
#' \code{coverage_weight} (sum of Hann weights).}
#' \item{windows}{matrix with one row per window and columns
#' \code{centre}, \code{f_low}, \code{f_high}, \code{f_centre},
#' \code{start_depth}, \code{end_depth}, \code{start_index},
#' \code{end_index}.}
#' \item{lower_controls, upper_controls}{the 2-column [depth, frequency]
#' control points used.}
#' \item{sample_interval}{median sample spacing.}
#' \item{nyquist}{Nyquist frequency.}
#' \item{window_samples, step_samples}{window and step in samples.}
#' \item{window_length, step_length}{effective (sample-rounded) window and
#' step in depth units.}
#' \item{fft_length}{zero-padded FFT length per window.}
#' \item{terminal_window_appended}{whether a final window ending on the last
#' sample was added.}
#' }
#'
#' @export
#' @importFrom stats fft
#' @importFrom stats approx
#' @importFrom stats median
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar

dynamic_filter <- function(data = NULL,
                           controls = NULL,
                           window_length = NULL,
                           step_length = NULL,
                           genplot = FALSE,
                           x_lab = "depth",
                           y_lab = "proxy",
                           verbose = FALSE){
  # ---- internal helper functions (local to this function) ----

  # ---------------------------------------------------------------------------
  # internal helpers
  # ---------------------------------------------------------------------------

  # @title Strict series check (Acycle validateSeries)
  # @keywords internal
  # @noRd
  dynf_check_series <- function(data = NULL, fun = "dynamic_filter"){
    if (is.null(data)){
      stop(paste0(fun, ": data must be a 2-column matrix [depth, value]"))
    }
    data <- as.matrix(data)
    if (ncol(data) != 2 || !is.numeric(data)){
      stop(paste0(fun, ": data must be a numeric 2-column matrix [depth, value]"))
    }
    storage.mode(data) <- "double"
    if (nrow(data) < 8){
      stop(paste0(fun, ": at least eight rows are required"))
    }
    if (any(!is.finite(data))){
      stop(paste0(fun, ": every depth and value must be finite"))
    }
    spacing <- diff(data[, 1])
    if (any(spacing <= 0)){
      stop(paste0(fun, ": depths must be sorted and strictly increasing"))
    }
    dt <- median(spacing)
    rel_dep <- max(abs(spacing - dt)) / dt
    if (!(is.finite(rel_dep)) || rel_dep > 1e-5){
      stop(sprintf(paste0("%s: relative spacing departure %.3g exceeds 1e-5; ",
                          "interpolate the data to a regular grid first"),
                   fun, rel_dep))
    }
    list(data = data, sample_interval = dt, nyquist = 1 / (2 * dt))
  }

  # @title Positive finite scalar check
  # @keywords internal
  # @noRd
  dynf_positive_scalar <- function(value = NULL, name = "value"){
    if (is.null(value) || length(value) != 1 || !is.numeric(value) ||
        !is.finite(value) || value <= 0){
      stop(paste0(name, " must be a finite positive number"))
    }
    as.numeric(value)
  }

  # @title Exact next power of two (2^nextpow2(n))
  # @keywords internal
  # @noRd
  dynf_pow2 <- function(n = 1){
    p <- 1
    while (p < n){
      p <- p * 2
    }
    p
  }

  # @title Gaussian response with half-power edges at f_low and f_high
  # @keywords internal
  # @noRd
  dynf_gauss_response <- function(abs_freq = NULL, f_low = NULL, f_high = NULL){
    f_centre <- f_low + (f_high - f_low) / 2
    offset <- (abs_freq - f_centre) / ((f_high - f_low) / 2)
    exp(-log(sqrt(2)) * offset^2)
  }

  # @title Half-sample-shifted Hann weight (strictly positive)
  # @keywords internal
  # @noRd
  dynf_hann_shifted <- function(n = 8){
    phase <- ((0:(n - 1)) + 0.5) / n
    0.5 - 0.5 * cos(2 * pi * phase)
  }

  # @title Linear interpolation with constant end hold (Acycle interpolateClamped)
  # @keywords internal
  # @noRd
  dynf_interp_clamped <- function(points = NULL, query = NULL){
    if (nrow(points) == 1){
      return(rep(points[1, 2], length(query)))
    }
    approx(points[, 1], points[, 2], xout = query, method = "linear",
           rule = 2, ties = "ordered")$y
  }

  # @title Split controls into lower/upper [depth, frequency] matrices
  # @keywords internal
  # @noRd
  dynf_split_controls <- function(controls = NULL){
    if (is.null(controls)){
      stop("controls must be given: [depth, f_low, f_high] or list(lower, upper)")
    }
    if (is.list(controls) && !is.data.frame(controls)){
      if (is.null(controls$lower) || is.null(controls$upper)){
        stop("a controls list must contain elements 'lower' and 'upper'")
      }
      return(list(lower = as.matrix(controls$lower),
                  upper = as.matrix(controls$upper)))
    }
    cm <- as.matrix(controls)
    if (ncol(cm) < 3){
      stop("controls must have columns depth, f_low, f_high (and optionally f_centre)")
    }
    nm <- colnames(cm)
    if (!is.null(nm) && all(c("depth", "f_low", "f_high") %in% nm)){
      d <- as.numeric(cm[, "depth"])
      fl <- as.numeric(cm[, "f_low"])
      fh <- as.numeric(cm[, "f_high"])
      if ("f_centre" %in% nm){
        fc <- as.numeric(cm[, "f_centre"])
      } else {
        fc <- NULL
      }
    } else {
      d <- as.numeric(cm[, 1])
      fl <- as.numeric(cm[, 2])
      fh <- as.numeric(cm[, 3])
      if (ncol(cm) >= 4){
        fc <- as.numeric(cm[, 4])
      } else {
        fc <- NULL
      }
    }
    if (!is.null(fc)){
      mid <- fl + (fh - fl) / 2
      if (any(!is.finite(fc)) || any(abs(fc - mid) > 1e-8 * abs(mid))){
        stop(paste0("f_centre must equal (f_low + f_high) / 2: the Acycle ",
                    "Gaussian pass band is symmetric about its centre"))
      }
    }
    list(lower = cbind(d, fl), upper = cbind(d, fh))
  }

  # @title Validate [depth, frequency] control points (Acycle validateControlPoints)
  # @keywords internal
  # @noRd
  dynf_check_points <- function(points = NULL, name = "lower", d_start = 0,
                                d_end = 1){
    points <- as.matrix(points)
    if (ncol(points) != 2 || nrow(points) < 1){
      stop(paste0(name, " control points must be a K x 2 matrix [depth, frequency]"))
    }
    points <- matrix(as.numeric(points), ncol = 2)
    if (nrow(points) > 1e4){
      stop(paste0(name, " control points exceed the maximum of 10000 rows"))
    }
    if (any(!is.finite(points))){
      stop(paste0("every ", name, " control depth and frequency must be finite"))
    }
    if (nrow(points) > 1 && any(diff(points[, 1]) <= 0)){
      stop(paste0(name, " control depths must be strictly increasing"))
    }
    if (any(points[, 1] < d_start) || any(points[, 1] > d_end)){
      stop(sprintf("every %s control depth must lie inside the data range [%.10g, %.10g]",
                   name, d_start, d_end))
    }
    if (any(points[, 2] <= 0)){
      stop(paste0("every ", name, " control frequency must be positive"))
    }
    points
  }

  # ---- main code ----
  series <- dynf_check_series(data, fun = "dynamic_filter")
  depth <- series$data[, 1]
  values <- series$data[, 2]
  n <- length(values)
  dt <- series$sample_interval
  nyquist <- series$nyquist
  d_start <- depth[1]
  d_end <- depth[n]
  span <- d_end - d_start

  ctrl <- dynf_split_controls(controls)
  lower <- dynf_check_points(ctrl$lower, "lower", d_start, d_end)
  upper <- dynf_check_points(ctrl$upper, "upper", d_start, d_end)
  check_depth <- sort(unique(c(d_start, lower[, 1], upper[, 1], d_end)))
  check_low <- dynf_interp_clamped(lower, check_depth)
  check_high <- dynf_interp_clamped(upper, check_depth)
  if (any(check_high >= nyquist)){
    first <- which(check_high >= nyquist)[1]
    stop(sprintf(paste0("upper frequency %.10g at depth %.10g must be ",
                        "strictly below the Nyquist frequency %.10g"),
                 check_high[first], check_depth[first], nyquist))
  }
  if (any(check_low >= check_high)){
    first <- which(check_low >= check_high)[1]
    stop(sprintf(paste0("pass band must satisfy f_low < f_high; at depth ",
                        "%.10g they are %.10g and %.10g"),
                 check_depth[first], check_low[first], check_high[first]))
  }

  # Acycle GUI defaults
  if (is.null(window_length)){
    window_length <- 0.2 * span
  }
  window_length <- dynf_positive_scalar(window_length, "window_length")
  if (is.null(step_length)){
    step_length <- 4 * dt
    if ((span - window_length) / dt > 500){
      step_length <- abs(span - window_length) / 500
    }
  }
  step_length <- dynf_positive_scalar(step_length, "step_length")

  if (window_length > span){
    stop(sprintf("window_length %.10g exceeds the data span %.10g",
                 window_length, span))
  }
  window_intervals <- floor(window_length / dt + 0.5)
  window_samples <- window_intervals + 1
  if (window_samples < 8){
    stop(sprintf(paste0("window_length resolves to %d samples; at least ",
                        "eight are required"), window_samples))
  }
  if (window_samples > n){
    stop("window_length resolves to more samples than the data contain")
  }
  step_samples <- floor(step_length / dt + 0.5)
  if (step_samples < 1){
    stop("step_length resolves to zero samples; use at least half the sample interval")
  }
  if (step_samples > window_samples){
    stop("step_length is longer than the window, which would leave gaps")
  }

  last_start <- n - window_samples + 1
  start_idx <- seq(1, last_start, by = step_samples)
  terminal <- start_idx[length(start_idx)] != last_start
  if (terminal == TRUE){
    start_idx <- c(start_idx, last_start)
  }
  end_idx <- start_idx + window_samples - 1
  centres <- depth[start_idx] + (depth[end_idx] - depth[start_idx]) / 2
  n_win <- length(start_idx)
  if (n_win > 1e4){
    stop(sprintf(paste0("the request produces %d windows (maximum 10000); ",
                        "increase step_length"), n_win))
  }

  win_low <- dynf_interp_clamped(lower, centres)
  win_high <- dynf_interp_clamped(upper, centres)
  win_centre <- win_low + (win_high - win_low) / 2

  fft_length <- dynf_pow2(2 * window_samples - 1)
  if (fft_length > 2^20){
    stop("the window requires an FFT longer than 2^20; reduce window_length")
  }
  df <- 1 / (fft_length * dt)
  width <- win_high - win_low
  tol <- 64 * .Machine$double.eps * 2^floor(log2(max(1, df)))
  if (any(width < df - tol)){
    first <- which(width < df - tol)[1]
    stop(sprintf(paste0("window %d has pass band width %.10g, below the ",
                        "FFT resolution %.10g; increase window_length or ",
                        "widen the pass band"), first, width[first], df))
  }

  value_scale <- max(abs(values))
  if (value_scale == 0){
    scaled <- rep(0, n)
  } else {
    scaled <- values / value_scale
  }
  weight <- dynf_hann_shifted(window_samples)
  analysis_weight <- sqrt(weight)
  cov_weight <- rep(0, n)
  cov_count <- rep(0, n)
  accum <- rep(0, n)
  bins <- 0:(fft_length - 1)
  abs_freq <- pmin(bins, fft_length - bins) * df
  pad <- rep(0, fft_length - window_samples)

  if (verbose == TRUE){
    pb <- txtProgressBar(min = 0, max = n_win, style = 3)
  }
  for (w in seq_len(n_win)){
    ind <- start_idx[w]:end_idx[w]
    segment <- scaled[ind]
    centred <- segment - mean(segment)
    response <- dynf_gauss_response(abs_freq, win_low[w], win_high[w])
    coefs <- fft(c(centred * analysis_weight, pad))
    time_domain <- fft(coefs * response, inverse = TRUE) / fft_length
    local <- Re(time_domain[1:window_samples]) * analysis_weight
    accum[ind] <- accum[ind] + local
    cov_weight[ind] <- cov_weight[ind] + weight
    cov_count[ind] <- cov_count[ind] + 1
    if (verbose == TRUE){
      setTxtProgressBar(pb, w)
    }
  }
  if (verbose == TRUE){
    close(pb)
  }
  if (any(cov_weight <= 0) || any(cov_count < 1)){
    stop("the window schedule failed to cover every sample")
  }
  if (value_scale == 0){
    filtered <- rep(0, n)
  } else {
    filtered <- (accum / cov_weight) * value_scale
  }

  out_filtered <- cbind(depth = depth, filtered = filtered,
                        coverage_count = cov_count,
                        coverage_weight = cov_weight)
  out_windows <- cbind(centre = centres, f_low = win_low, f_high = win_high,
                       f_centre = win_centre, start_depth = depth[start_idx],
                       end_depth = depth[end_idx], start_index = start_idx,
                       end_index = end_idx)
  colnames(lower) <- c("depth", "frequency")
  colnames(upper) <- c("depth", "frequency")
  res <- list(filtered = out_filtered, windows = out_windows,
              lower_controls = lower, upper_controls = upper,
              sample_interval = dt, nyquist = nyquist,
              window_samples = window_samples, step_samples = step_samples,
              window_length = window_intervals * dt,
              step_length = step_samples * dt, fft_length = fft_length,
              terminal_window_appended = terminal)

  if (genplot == TRUE){
    plot_dynamic_filter(data = series$data, dyn_filter = res,
                        x_lab = x_lab, y_lab = y_lab)
  }
  res
}
