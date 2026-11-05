#' @title Power decomposition analysis (strict evolutionary MTM band-power ratio)
#'
#' @description The \code{\link{pda}} function calculates, in moving windows,
#' the ratio between the multitaper (MTM) power inside a set of target
#' frequency bands and the power inside a total frequency band (power
#' decomposition analysis, PDA; Li et al., 2016). Each window is linearly
#' detrended, its adaptive-weighted MTM power spectral density is computed on
#' an explicit FFT grid, and band powers are obtained by trapezoidal
#' integration of the PSD with piecewise-linear interpolation at the exact
#' band limits. The input is divided by one global scale (maximum absolute
#' value) so that powers are comparable between windows.
#' The data must already be regularly sampled (relative spacing tolerance
#' 1e-5); no sorting, interpolation or padding is performed.
#' Defaults for \code{step}, \code{nfft} and \code{total_band} follow the
#' Acycle compatibility wrapper \code{pda.m} (one sample interval, next power
#' of two of the window length, DC to Nyquist).
#' Translated from Acycle v3.1 acyclePowerDecomposition.m and pda.m
#' (Mingsong Li).
#'
#' @param data Numeric two-column matrix or data frame: depth/time in column 1 (strictly increasing and evenly
#' spaced, relative spacing tolerance 1e-5; interpolate first if needed) and the proxy value in column 2; at least
#' eight rows, all values finite. Required (\code{Default=NULL}).
#' @param target_bands Target frequency bands whose summed power is divided by the power of \code{total_band}, in
#' cycles per unit of the depth/time column: either a B x 2 matrix with the lower limit in column 1 and the upper
#' limit in column 2 (rows ordered by increasing lower limit, non-overlapping), or a numeric vector of paired limits
#' \code{c(lower1, upper1, lower2, upper2, ...)}. Every band must lie inside \code{total_band}, below the Nyquist
#' frequency and be at least as wide as the FFT resolution. \code{Default=c(1/45, 1/25)}.
#' @param window Length of the moving window (first-to-last coordinate span), in the unit of the depth/time column
#' (numeric, > 0 and not longer than the record; rounded to a whole number of sample intervals, at least 8 samples).
#' \code{Default=500}.
#' @param step Distance between the starts of successive windows, in the unit of the depth/time column (numeric,
#' > 0; rounded to whole samples). \code{NULL} uses one sample interval. \code{Default=NULL}.
#' @param nw Time-bandwidth product of the DPSS (Slepian) tapers of the multitaper spectrum (numeric, >= 0.5 and less
#' than half the number of samples per window); larger values smooth the spectrum more. \code{Default=2}.
#' @param nfft Length of the FFT (even integer, at least the number of samples per window); larger values give a
#' finer frequency grid. \code{NULL} uses the next power of two of the number of samples per window.
#' \code{Default=NULL}.
#' @param total_band Total frequency band \code{c(lower, upper)} in cycles per unit of the depth/time column (numeric
#' vector of length 2, 0 <= lower < upper <= Nyquist frequency), the denominator of the power ratio. \code{NULL} uses
#' \code{c(0, Nyquist frequency)}. \code{Default=NULL}.
#' @param genplot Plot the data and the power ratio against depth/time (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param verbose Show a progress bar over the windows (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' and Linda Hinnov (Acycle).
#'
#' @references
#' Li, M., Huang, C., Hinnov, L., Ogg, J., Chen, Z.-Q., Zhang, Y., 2016.
#' Obliquity-forced climate during the Early Triassic hothouse in China.
#' Geology 44, 623-626. \doi{10.1130/G37970.1}
#'
#' Li, M., Hinnov, L., Kump, L., 2019. Acycle: Time-series analysis software
#' for paleoclimate research and education. Computers & Geosciences 127,
#' 12-22. \doi{10.1016/j.cageo.2019.02.011}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 2000, by = 2)
#' x <- sin(2 * pi * t / 40) + 0.5 * sin(2 * pi * t / 100) +
#'   0.5 * as.numeric(arima.sim(list(ar = 0.7), n = length(t)))
#' res <- pda(data = cbind(t, x), target_bands = c(1/45, 1/35),
#'   window = 300, step = 10, nw = 2, genplot = FALSE)
#' head(res$result)
#' }
#' @return A list with elements
#' \code{result}: matrix with columns window center, ratio (target power /
#' total power), target power (union of all target bands) and total power;
#' \code{band_powers}: matrix with the window center followed by one power
#' column per target band;
#' \code{meta}: list with the resolved sampling interval, window and step
#' samples, effective window/step lengths, nfft, frequency resolution,
#' Nyquist frequency, value scale and number of unused trailing samples.
#'
#' @export
#' @importFrom stats median
#' @importFrom stats approx
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom graphics par
#' @importFrom graphics plot

pda <- function(data = NULL,
    target_bands = c(1/45, 1/25),
    window = 500,
    step = NULL,
    nw = 2,
    nfft = NULL,
    total_band = NULL,
    genplot = FALSE,
    verbose = FALSE){
  # Internal helper functions (local to this function)

  # small utilities

  acy_nextpow2 <- function(n){
    if (n <= 1){
      return(0)
    }
    ceiling(log2(n) - 1e-12)
  }

  acy_dpss <- function(n, nw, k = NULL){
    if (is.null(k)){
      k <- round(2 * nw) - 1
    }
    key <- paste(n, format(nw, digits = 17), k, sep = "_")
    if (exists(key, envir = acy_dpss_cache, inherits = FALSE)){
      return(get(key, envir = acy_dpss_cache))
    }
    w <- nw / n
    i <- 0:(n - 1)
    dg <- ((n - 1 - 2 * i) / 2)^2 * cos(2 * pi * w)
    od <- (1:(n - 1)) * (n - (1:(n - 1))) / 2
    if (n <= 2048){
      Tm <- diag(dg, n, n)
      for (j in 1:(n - 1)){
        Tm[j, j + 1] <- od[j]
        Tm[j + 1, j] <- od[j]
      }
      ev <- eigen(Tm, symmetric = TRUE)
      E <- ev$vectors[, 1:k, drop = FALSE]
    } else {
      E <- acy_dpss_tridiag(dg, od, k)
    }
    tc <- i - (n - 1) / 2
    for (j in 1:k){
      E[, j] <- E[, j] / sqrt(sum(E[, j]^2))
      if (j %% 2 == 1){
        if (sum(E[, j]) < 0){
          E[, j] <- -E[, j]
        }
      } else {
        if (sum(tc * E[, j]) < 0){
          E[, j] <- -E[, j]
        }
      }
    }
    V <- numeric(k)
    for (j in 1:k){
      V[j] <- acy_dpss_concentration(E[, j], w)
    }
    out <- list(E = E, V = V)
    assign(key, out, envir = acy_dpss_cache)
    out
  }

  acy_dpss_concentration <- function(e, w){
    n <- length(e)
    m <- 2^acy_nextpow2(2 * n)
    fe <- stats::fft(c(e, rep(0, m - n)))
    r <- Re(stats::fft(Mod(fe)^2, inverse = TRUE)) / m
    lag <- 1:(n - 1)
    a <- sin(2 * pi * w * lag) / (pi * lag)
    2 * w * r[1] + 2 * sum(a * r[lag + 1])
  }

  acy_dpss_tridiag <- function(dg, od, k){
    n <- length(dg)
    sturm_count <- function(lambda){
      cnt <- 0
      q <- dg[1] - lambda
      if (q < 0){
        cnt <- cnt + 1
      }
      for (j in 2:n){
        if (q == 0){
          q <- 1e-300
        }
        q <- dg[j] - lambda - od[j - 1]^2 / q
        if (q < 0){
          cnt <- cnt + 1
        }
      }
      cnt
    }
    bound <- max(abs(dg)) + 2 * max(abs(od))
    E <- matrix(0, n, k)
    for (j in 1:k){
      target <- n - j
      lo <- -bound
      hi <- bound
      for (it in 1:200){
        mid <- (lo + hi) / 2
        if (sturm_count(mid) > target){
          hi <- mid
        } else {
          lo <- mid
        }
        if (hi - lo <= 1e-14 * max(1, abs(mid))){
          break
        }
      }
      lambda <- (lo + hi) / 2 + 1e-10 * max(1, abs(hi))
      v <- rep(1, n) / sqrt(n)
      for (it in 1:6){
        v <- acy_tridiag_solve(od, dg - lambda, od, v)
        v <- v / sqrt(sum(v^2))
      }
      E[, j] <- v
    }
    E
  }

  acy_tridiag_solve <- function(a, b, c, d){
    n <- length(b)
    cp <- numeric(n)
    dp <- numeric(n)
    cp[1] <- c[1] / b[1]
    dp[1] <- d[1] / b[1]
    for (i in 2:n){
      den <- b[i] - a[i - 1] * cp[i - 1]
      if (den == 0){
        den <- 1e-300
      }
      if (i < n){
        cp[i] <- c[i] / den
      }
      dp[i] <- (d[i] - a[i - 1] * dp[i - 1]) / den
    }
    x <- numeric(n)
    x[n] <- dp[n]
    for (i in (n - 1):1){
      x[i] <- dp[i] - cp[i] * x[i + 1]
    }
    x
  }

  # MATLAB periodogram / pmtm (one-sided, adaptive weights)

  acy_datawrap <- function(x, nfft){
    n <- length(x)
    if (n <= nfft){
      return(c(x, rep(0, nfft - n)))
    }
    out <- rep(0, nfft)
    nseg <- ceiling(n / nfft)
    xx <- c(x, rep(0, nseg * nfft - n))
    for (s in 1:nseg){
      out <- out + xx[((s - 1) * nfft + 1):(s * nfft)]
    }
    out
  }

  acy_onesided <- function(S, nfft){
    if (nfft %% 2 == 0){
      sel <- 1:(nfft / 2 + 1)
      P <- S[sel]
      if (length(P) > 2){
        P[2:(length(P) - 1)] <- 2 * P[2:(length(P) - 1)]
      }
    } else {
      sel <- 1:((nfft + 1) / 2)
      P <- S[sel]
      if (length(P) > 1){
        P[2:length(P)] <- 2 * P[2:length(P)]
      }
    }
    P
  }

  acy_freq_vector <- function(nfft, fs = NULL){
    if (nfft %% 2 == 0){
      m <- nfft / 2 + 1
    } else {
      m <- (nfft + 1) / 2
    }
    if (is.null(fs)){
      return(2 * pi * (0:(m - 1)) / nfft)
    }
    fs * (0:(m - 1)) / nfft
  }

  acy_pmtm <- function(x, nw = 4, nfft = NULL, fs = NULL, drop_last = TRUE,
      method = "adapt"){
    x <- as.numeric(x)
    n <- length(x)
    if (is.null(nfft)){
      nfft <- max(256, 2^acy_nextpow2(n))
    }
    if (drop_last == TRUE){
      k <- round(2 * nw) - 1
    } else {
      k <- round(2 * nw)
    }
    if (k < 1){
      k <- 1
    }
    tp <- acy_dpss(n, nw, k)
    E <- tp$E
    V <- tp$V
    Sk <- matrix(0, nfft, k)
    for (j in 1:k){
      Sk[, j] <- Mod(stats::fft(acy_datawrap(E[, j] * x, nfft)))^2
    }
    if (method == "adapt" && k > 1){
      sig2 <- sum(x^2) / n
      S <- (Sk[, 1] + Sk[, 2]) / 2
      S1 <- rep(0, nfft)
      tol <- 0.0005 * sig2 / nfft
      a <- sig2 * (1 - V)
      it <- 0
      while (sum(abs(S - S1) / nfft) > tol && it < 100){
        it <- it + 1
        num <- 0
        den <- 0
        for (j in 1:k){
          b <- S / (S * V[j] + a[j])
          wk <- b^2 * V[j]
          num <- num + wk * Sk[, j]
          den <- den + wk
        }
        S1 <- S
        S <- num / den
      }
    } else if (method == "eigen"){
      S <- rep(0, nfft)
      for (j in 1:k){
        S <- S + V[j] * Sk[, j]
      }
      S <- S / sum(V)
    } else {
      S <- rowMeans(Sk)
    }
    P <- acy_onesided(S, nfft)
    if (is.null(fs)){
      P <- P / (2 * pi)
    } else {
      P <- P / fs
    }
    list(pxx = P, f = acy_freq_vector(nfft, fs), tapers = E, eigenvalues = V)
  }

  # internal helpers

  # @keywords internal
  # @noRd
  pda_eps <- function(x){
    x <- abs(x)
    if (x == 0){
      return(2^-1074)
    }
    2^(floor(log2(x)) - 52)
  }

  # @keywords internal
  # @noRd
  pda_linspace <- function(a, b, n){
    if (n < 1){
      return(numeric(0))
    }
    if (n == 1){
      return(b)
    }
    out <- a + (0:(n - 1)) * (b - a) / (n - 1)
    out[1] <- a
    out[n] <- b
    out
  }

  # @keywords internal
  # @noRd
  pda_band_matrix <- function(target_bands){
    if (is.matrix(target_bands)){
      bands <- target_bands
    } else {
      tb <- as.numeric(target_bands)
      if (length(tb) %% 2 != 0 || length(tb) < 2){
        stop("target_bands must contain paired band limits")
      }
      bands <- matrix(tb, ncol = 2, byrow = TRUE)
    }
    bands <- matrix(as.numeric(bands), ncol = 2)
    if (ncol(bands) != 2 || nrow(bands) < 1 || any(!is.finite(bands))){
      stop("target_bands must be a finite B x 2 matrix [lower, upper]")
    }
    bands
  }

  # @keywords internal
  # @noRd
  pda_detrend_segment <- function(segment){
    # linearDetrendSegment of acyclePowerDecomposition.m
    offset <- min(segment) / 2 + max(segment) / 2
    centered <- segment - offset
    local_scale <- max(abs(centered))
    if (!(is.finite(local_scale) && local_scale > 0)){
      stop("A window has no finite value variation")
    }
    n <- length(segment)
    ordinal <- pda_linspace(-1, 1, n)
    design <- cbind(1, ordinal)
    scaled <- centered / local_scale
    beta <- qr.solve(design, scaled)
    residual <- (scaled - as.numeric(design %*% beta)) * local_scale
    threshold <- 16 * .Machine$double.eps * max(1, max(abs(segment)))
    if (!(max(abs(residual)) > threshold)){
      stop("A window has no detrended variation above floating-point resolution")
    }
    residual
  }

  # @keywords internal
  # @noRd
  pda_integrate_band <- function(freq, power, lower, upper){
    interior <- freq > lower & freq < upper
    fi <- c(lower, freq[interior], upper)
    pi_ <- stats::approx(freq, power, xout = fi, method = "linear", rule = 1)$y
    if (any(!is.finite(pi_)) || any(pi_ < 0)){
      stop("Piecewise-linear PSD interpolation failed")
    }
    n <- length(fi)
    sum(diff(fi) * (pi_[-1] + pi_[-n]) / 2)
  }

  # Slepian (DPSS) tapers, MATLAB dpss() conventions (unit energy)

  acy_dpss_cache <- new.env(parent = emptyenv())

  # Main code

  data <- as.matrix(data)
  if (!is.numeric(data) || ncol(data) != 2){
    stop("data must be a numeric matrix with two columns [depth/time, value]")
  }
  n_rows <- nrow(data)
  if (n_rows < 8){
    stop("Power decomposition requires at least eight input rows")
  }
  if (any(!is.finite(data))){
    stop("Every data coordinate and value must be finite")
  }
  coordinates <- as.numeric(data[, 1])
  values <- as.numeric(data[, 2])
  spacing <- diff(coordinates)
  if (any(spacing <= 0)){
    stop("data coordinates must be sorted, unique and strictly increasing")
  }
  dt <- stats::median(spacing)
  spacing_error <- max(abs(spacing - dt))
  if (!(is.finite(dt) && dt > 0) || spacing_error / dt > 1e-5){
    stop("Relative spacing departure exceeds 1e-5; interpolate before power decomposition")
  }
  nyquist <- 1 / (2 * dt)
  span <- coordinates[n_rows] - coordinates[1]

  bands <- pda_band_matrix(target_bands)
  if (any(bands[, 1] < 0) || any(bands[, 2] <= bands[, 1])){
    stop("Every target band must satisfy 0 <= lower < upper")
  }
  if (nrow(bands) > 1){
    if (any(diff(bands[, 1]) <= 0)){
      stop("target_bands rows must be ordered by strictly increasing lower limit")
    }
    if (any(bands[-1, 1] < bands[-nrow(bands), 2])){
      stop("Target bands must not overlap")
    }
  }

  if (!(is.numeric(window) && length(window) == 1 && is.finite(window) && window > 0)){
    stop("window must be a finite positive scalar")
  }
  if (is.null(step)){
    step <- dt
  }
  if (!(is.numeric(step) && length(step) == 1 && is.finite(step) && step > 0)){
    stop("step must be a finite positive scalar")
  }
  if (window > span){
    stop("window exceeds the data coordinate span")
  }
  window_intervals <- floor(window / dt + 0.5)
  window_samples <- window_intervals + 1
  if (window_samples < 8){
    stop("window resolves to fewer than eight samples")
  }
  if (window_samples > n_rows){
    stop("window resolves to more samples than input rows")
  }
  step_samples <- floor(step / dt + 0.5)
  if (step_samples < 1){
    stop("step resolves to zero samples")
  }
  if (!(is.numeric(nw) && length(nw) == 1 && is.finite(nw) && nw > 0)){
    stop("nw must be a finite positive scalar")
  }
  if (nw < 0.5 || nw >= window_samples / 2){
    stop("nw must be at least 0.5 and strictly less than window_samples/2")
  }
  if (is.null(nfft)){
    nfft <- 2^acy_nextpow2(window_samples)
  }
  if (!(is.numeric(nfft) && length(nfft) == 1 && is.finite(nfft) &&
      nfft == trunc(nfft) && nfft >= window_samples && nfft %% 2 == 0)){
    stop("nfft must be an even integer at least as large as the window sample count")
  }
  if (is.null(total_band)){
    total_band <- c(0, nyquist)
  }
  total_band <- as.numeric(total_band)
  if (length(total_band) != 2 || any(!is.finite(total_band))){
    stop("total_band must be a finite two-element vector")
  }
  if (total_band[1] < 0 || total_band[2] <= total_band[1]){
    stop("total_band must satisfy 0 <= lower < upper")
  }
  if (total_band[2] > nyquist){
    stop("total_band upper frequency exceeds the Nyquist frequency")
  }
  if (any(bands[, 2] > nyquist)){
    stop("A target band upper frequency exceeds the Nyquist frequency")
  }
  if (any(bands[, 1] < total_band[1]) || any(bands[, 2] > total_band[2])){
    stop("Every target band must lie completely inside total_band")
  }

  fs <- 1 / dt
  resolution <- fs / nfft
  res_tol <- 64 * pda_eps(max(1, resolution))
  band_widths <- bands[, 2] - bands[, 1]
  if (any(band_widths < resolution - res_tol)){
    stop("A target band is narrower than the FFT resolution; increase nfft or widen the band")
  }
  if ((total_band[2] - total_band[1]) < resolution - res_tol){
    stop("total_band is narrower than the FFT resolution")
  }

  start_idx <- seq(1, n_rows - window_samples + 1, by = step_samples)
  end_idx <- start_idx + window_samples - 1
  centers <- coordinates[start_idx] + (coordinates[end_idx] - coordinates[start_idx]) / 2
  n_win <- length(start_idx)
  n_band <- nrow(bands)

  value_scale <- max(abs(values))
  if (!(is.finite(value_scale) && value_scale > 0)){
    stop("The value column must contain nonzero finite variation")
  }
  scaled <- values / value_scale

  ratio <- numeric(n_win)
  target_power <- numeric(n_win)
  total_power <- numeric(n_win)
  band_power <- matrix(0, n_win, n_band)

  if (verbose == TRUE){
    pb <- utils::txtProgressBar(min = 0, max = n_win, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (w in seq_len(n_win)){
    segment <- scaled[start_idx[w]:end_idx[w]]
    residual <- pda_detrend_segment(segment)
    spec <- acy_pmtm(residual, nw = nw, nfft = nfft, fs = fs)
    power <- spec$pxx
    freq <- spec$f
    freq[1] <- 0
    freq[length(freq)] <- nyquist
    total_power[w] <- pda_integrate_band(freq, power, total_band[1], total_band[2])
    for (b in seq_len(n_band)){
      band_power[w, b] <- pda_integrate_band(freq, power, bands[b, 1], bands[b, 2])
    }
    target_power[w] <- sum(band_power[w, ])
    if (!(is.finite(total_power[w]) && total_power[w] > 0)){
      stop(paste0("Window ", w, " has zero or nonfinite power inside total_band"))
    }
    raw_ratio <- target_power[w] / total_power[w]
    if (raw_ratio < 0 || raw_ratio > 1 + 4096 * .Machine$double.eps || !is.finite(raw_ratio)){
      stop(paste0("Window ", w, " produced a power ratio outside [0, 1]"))
    }
    if (raw_ratio > 1){
      raw_ratio <- 1
    }
    ratio[w] <- raw_ratio
    if (verbose == TRUE){
      utils::setTxtProgressBar(pb, w)
    }
  }

  result <- cbind(centers, ratio, target_power, total_power)
  colnames(result) <- c("center", "ratio", "target_power", "total_power")
  band_powers <- cbind(centers, band_power)
  colnames(band_powers) <- c("center", paste0("band_", seq_len(n_band)))

  meta <- list(sample_interval = dt,
      nyquist = nyquist,
      window_samples = window_samples,
      effective_window_length = window_intervals * dt,
      step_samples = step_samples,
      effective_step_length = step_samples * dt,
      nw = nw,
      nfft = nfft,
      frequency_resolution = resolution,
      target_bands = bands,
      total_band = total_band,
      value_scale = value_scale,
      window_count = n_win,
      unused_trailing_samples = n_rows - end_idx[n_win])

  if (genplot == TRUE){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
    graphics::plot(coordinates, values, type = "l", xlab = "depth/time", ylab = "value")
    graphics::plot(centers, ratio, type = "l", xlab = "depth/time", ylab = "Power ratio",
        xlim = range(coordinates))
  }

  list(result = result, band_powers = band_powers, meta = meta)
}
