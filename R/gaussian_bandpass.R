#' @title Gaussian band-pass filter (Acycle)
#'
#' @description The \code{\link{gaussian_bandpass}} function applies the
#' Acycle Gaussian band-pass filter to a regularly sampled series. The value
#' mean is removed, the series is zero padded to \code{10 * 2^nextpow2(n)}
#' samples and multiplied in the Fourier domain by a symmetric Gaussian
#' response centred on \code{(f_low + f_high) / 2} whose lower and upper
#' frequencies are exact half-power (-3 dB, gain \code{1/sqrt(2)}) points.
#' The amplitude envelope is the modulus of the analytic signal (MATLAB
#' \code{hilbert}) of the filtered series. The Gaussian parameterisation
#' mimics the Gaussian filter of Analyseries 2.0 (Paillard et al., 1996),
#' as implemented by K. Ramamurthy (2013) and corrected by L. Hinnov (2016).
#' Translated from Acycle v3.1 acycleBandpassFilter.m (method "gaussian")
#' and gaussfilter.m (Mingsong Li).
#'
#' @param data Numeric two-column matrix or data frame: depth/time in column 1 (strictly increasing and regularly
#' sampled, relative spacing tolerance 1e-5; interpolate first if needed) and the proxy value in column 2; at least
#' eight rows, all values finite. Required (\code{Default=NULL}).
#' @param f_low Lower half-power (-3 dB) frequency of the pass band, in cycles per unit of the depth/time column
#' (numeric, > 0; e.g. 1/25 to keep periods shorter than 25 m). Required (\code{Default=NULL}).
#' @param f_high Upper half-power (-3 dB) frequency of the pass band, in cycles per unit of the depth/time column
#' (numeric, > \code{f_low} and strictly below the Nyquist frequency \code{1 / (2 * sampling interval)}). Required
#' (\code{Default=NULL}).
#' @param genplot Plot the demeaned data with the filtered signal and its envelope, and the filter response
#' (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param x_lab Label of the x-axis of the data plot (character). \code{Default="depth"}.
#' @param y_lab Label of the y-axis of the data plot (character). \code{Default="proxy"}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' (Acycle), Linda Hinnov and K. Ramamurthy (gaussfilter).
#'
#' @references
#' Li, M., Hinnov, L., Kump, L., 2019. Acycle: Time-series analysis software
#' for paleoclimate research and education. Computers & Geosciences 127,
#' 12-22. <doi:10.1016/j.cageo.2019.02.011>
#'
#' Paillard, D., Labeyrie, L., Yiou, P., 1996. Macintosh program performs
#' time-series analysis. Eos 77, 379. <doi:10.1029/96EO00259>
#'
#' @examples
#' \donttest{
#' depth <- seq(0, 200, by = 0.5)
#' proxy <- sin(2 * pi * depth / 20) + 0.5 * sin(2 * pi * depth / 5)
#' res <- gaussian_bandpass(data = cbind(depth, proxy), f_low = 0.04,
#'   f_high = 0.06, genplot = FALSE)
#' plot(res$filtered[, 1], res$filtered[, 2], type = "l")
#' }
#'
#' @return A list with elements
#' \describe{
#' \item{filtered}{matrix with columns \code{depth}, \code{filtered} and
#' \code{envelope}.}
#' \item{response}{matrix with columns \code{frequency} and \code{gain}: the
#' one-sided Gaussian response (0 to Nyquist of the padded FFT grid).}
#' \item{f_low, f_high, f_centre}{the pass band used.}
#' \item{sample_interval}{median sample spacing.}
#' \item{nyquist}{Nyquist frequency \code{1 / (2 * sample_interval)}.}
#' \item{nfft}{length of the zero-padded FFT.}
#' \item{input_mean}{mean removed before filtering.}
#' }
#'
#' @export
#' @importFrom stats fft
#' @importFrom stats median
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics abline
#' @importFrom graphics legend

gaussian_bandpass <- function(data = NULL,
    f_low = NULL,
    f_high = NULL,
    genplot = FALSE,
    x_lab = "depth",
    y_lab = "proxy"){
  # Internal helper functions (local to this function)

  # internal helpers

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

  # @title Analytic signal (MATLAB hilbert)
  # @keywords internal
  # @noRd
  dynf_hilbert <- function(x = NULL){
    n <- length(x)
    h <- rep(0, n)
    h[1] <- 1
    if (n %% 2 == 0){
      h[n / 2 + 1] <- 1
      if (n > 2){
        h[2:(n / 2)] <- 2
      }
    } else {
      if (n > 1){
        h[2:((n + 1) / 2)] <- 2
      }
    }
    fft(fft(x) * h, inverse = TRUE) / n
  }

  # @title Scaled mean (Acycle stableMean)
  # @keywords internal
  # @noRd
  dynf_stable_mean <- function(x = NULL){
    s <- max(abs(x))
    if (s == 0){
      return(0)
    }
    s * mean(x / s)
  }

  # @title Scaled sample variance (Acycle stableVariance)
  # @keywords internal
  # @noRd
  dynf_stable_variance <- function(x = NULL){
    r <- x - dynf_stable_mean(x)
    s <- max(abs(r))
    if (s == 0){
      return(0)
    }
    (s * sqrt(sum((r / s)^2) / (length(x) - 1)))^2
  }

  # Main code
  series <- dynf_check_series(data, fun = "gaussian_bandpass")
  f_low <- dynf_positive_scalar(f_low, "f_low")
  f_high <- dynf_positive_scalar(f_high, "f_high")
  if (f_high <= f_low){
    stop("f_high must be greater than f_low")
  }
  if (f_high >= series$nyquist){
    stop(sprintf("f_high must be strictly below the Nyquist frequency %.10g",
        series$nyquist))
  }
  x <- series$data[, 2]
  n <- length(x)
  input_mean <- dynf_stable_mean(x)
  centred <- x - input_mean
  if (!(dynf_stable_variance(x) > 0)){
    stop("the value column must have a finite positive variance")
  }
  f_centre <- (f_low + f_high) / 2
  dt <- series$sample_interval

  nfft <- 10 * dynf_pow2(n)
  df <- 1 / (nfft * dt)
  idx <- 0:(nfft - 1)
  abs_freq <- pmin(idx, nfft - idx) * df
  response <- dynf_gauss_response(abs_freq, f_low, f_high)
  padded <- c(centred, rep(0, nfft - n))
  time_domain <- fft(fft(padded) * response, inverse = TRUE) / nfft
  filtered <- Re(time_domain[1:n])
  envelope <- Mod(dynf_hilbert(filtered))

  out_filtered <- cbind(depth = series$data[, 1], filtered = filtered,
      envelope = envelope)
  out_response <- cbind(frequency = (0:(nfft / 2)) * df,
      gain = response[1:(nfft / 2 + 1)])
  res <- list(filtered = out_filtered, response = out_response,
      f_low = f_low, f_high = f_high, f_centre = f_centre,
      sample_interval = dt, nyquist = series$nyquist, nfft = nfft,
      input_mean = input_mean)

  if (genplot == TRUE){
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
    plot(series$data[, 1], centred, type = "l", col = "grey50",
        xlab = x_lab, ylab = y_lab)
    lines(series$data[, 1], filtered, col = "red", lwd = 1.5)
    lines(series$data[, 1], envelope, col = "blue")
    lines(series$data[, 1], -envelope, col = "blue")
    legend("topright", legend = c("data (demeaned)", "filtered", "envelope"),
        col = c("grey50", "red", "blue"), lty = 1, bty = "n", cex = 0.8)
    xmax <- min(series$nyquist, f_high + 4 * (f_high - f_low))
    plot(out_response[, 1], out_response[, 2], type = "l",
        xlim = c(0, xmax), xlab = "frequency", ylab = "gain")
    abline(v = c(f_low, f_high), lty = 2, col = "grey40")
    abline(h = 1 / sqrt(2), lty = 3, col = "grey40")
  }
  res
}
