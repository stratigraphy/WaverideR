#' @title Legacy power decomposition analysis (Acycle PDA menu / DYNOT core)
#'
#' @description The \code{\link{pda_legacy}} function reproduces the power
#' decomposition analysis as run from the Acycle menu
#' "Power Decomposition Analysis" and as used inside DYNOT (\code{pdan.m};
#' Li et al., 2016). Preprocessing as in the Acycle menu: rows with NA are
#' removed, data are sorted, duplicated depths are replaced by their mean,
#' unevenly sampled data (relative tolerance 1e-5) are linearly interpolated
#' at the median sampling interval, and the series is optionally padded by
#' half a window at both ends (\code{zeropad2.m}). In each window of
#' \code{trunc(window/dt)} samples the data are linearly detrended and the
#' adaptive MTM spectrum (\code{pmtm(x, nw, pad)}) is computed. Target and
#' total powers are sums of spectral bins selected with the historical
#' proportional bin mapping (lower = ceiling(nbin*f/fNyq),
#' upper = trunc(nbin*f/fNyq)); overlapping target bands are counted once.
#' Window centers are evenly spaced between first coordinate + window/2 and
#' last coordinate - window/2.
#' Translated from Acycle v3.1 pdan.m, zeropad2.m, interpolate.m,
#' findduplicate.m and AC.m (menu_pda_Callback) (Mingsong Li).
#'
#' @param data Two-column numeric matrix or data frame: depth/time in column 1 and the proxy value in column 2 (further
#' columns are ignored). Rows with \code{NA} are removed, the data are sorted, duplicated depths are averaged and
#' unevenly sampled data are interpolated at the median sampling interval. Required (\code{Default=NULL}).
#' @param target_bands Target frequency bands whose summed power is divided by the power of the total band, in cycles
#' per unit of the depth/time column: a numeric vector of paired limits \code{c(lower1, upper1, lower2, upper2, ...)}
#' (any order within a pair) or a B x 2 matrix with one band per row; overlapping bands are counted once.
#' \code{Default=c(1/45, 1/25)}.
#' @param window Length of the moving window in the unit of the depth/time column (numeric, > 0); each window holds
#' \code{trunc(window / dt)} samples. \code{Default=500}.
#' @param nw Time-bandwidth product of the DPSS (Slepian) tapers of the multitaper spectrum (numeric, > 0); larger
#' values smooth the spectrum more. \code{Default=2}.
#' @param f_min Lower cutoff frequency of the total band, in cycles per unit of the depth/time column (numeric, >= 0).
#' \code{Default=0}.
#' @param f_max Upper cutoff frequency of the total band, in cycles per unit of the depth/time column (numeric, >
#' \code{f_min}); values above the Nyquist frequency are clipped. \code{NULL} uses the Nyquist frequency.
#' \code{Default=NULL}.
#' @param step Step between successive windows, in number of samples (positive integer). \code{Default=1}.
#' @param pad Length of the FFT used by the multitaper spectrum, i.e. the zero-padding number of \code{pmtm} (positive
#' integer); larger values give a finer frequency grid. \code{Default=5000}.
#' @param pad_type Padding of both ends of the (demeaned) series by half a window before the analysis, so that the
#' windows cover the whole record (integer): \code{0} none, \code{1} zeros, \code{2} mirror image of the record,
#' \code{3} mean of the first/last half window, \code{4} random values with the mean and standard deviation of the
#' first/last half window. \code{Default=1}.
#' @param genplot Plot the power ratio against depth/time (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param verbose Show a progress bar over the windows (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' and Linda Hinnov (Acycle); zeropad2 padding options by Nicolas Thibault
#' and Giovanni Rizzi.
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
#' res <- pda_legacy(data = cbind(t, x), target_bands = c(1/45, 1/35),
#'   window = 300, step = 5, pad = 1000, pad_type = 1)
#' head(res)
#' }
#' @return A matrix with columns window center, ratio (target / total),
#' target power and total power.
#'
#' @export
#' @importFrom stats median
#' @importFrom stats approx
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom graphics par
#' @importFrom graphics plot

pda_legacy <- function(data = NULL,
    target_bands = c(1/45, 1/25),
    window = 500,
    nw = 2,
    f_min = 0,
    f_max = NULL,
    step = 1,
    pad = 5000,
    pad_type = 1,
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

  # @keywords internal
  # @noRd
  pda_round <- function(x){
    sign(x) * floor(abs(x) + 0.5)
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
  pda_colon <- function(a, b, d){
    # MATLAB/Octave colon operator a:d:b for d > 0
    q <- (b - a) / d
    n <- floor(q + 3 * .Machine$double.eps * max(1, abs(q))) + 1
    if (n < 1){
      return(numeric(0))
    }
    out <- a + (0:(n - 1)) * d
    if (out[n] > b){
      out[n] <- b
    }
    out
  }

  # @keywords internal
  # @noRd
  pda_is_uneven <- function(x){
    # acycleSamplingIsUneven.m with relative tolerance 1e-5
    if (length(x) < 2){
      return(FALSE)
    }
    d <- diff(x)
    if (any(!is.finite(d)) || any(d <= 0)){
      return(TRUE)
    }
    md <- stats::median(d)
    max(abs(d - md)) / md > 1e-5
  }

  # @keywords internal
  # @noRd
  pda_findduplicate <- function(data){
    # findduplicate.m: duplicated depths replaced by one row with mean values
    x <- data[, 1]
    b <- sort(unique(x))
    if (length(b) == length(x)){
      return(data)
    }
    y <- data[, -1, drop = FALSE]
    d <- matrix(0, length(b), ncol(y))
    for (i in seq_along(b)){
      sel <- which(x == b[i])
      for (j in seq_len(ncol(y))){
        d[i, j] <- sum(y[sel, j]) / length(sel)
      }
    }
    cbind(b, d)
  }

  # @keywords internal
  # @noRd
  pda_interpolate <- function(data, rate){
    # interpolate.m: linear interpolation on data(1,1):rate:data(end,1)
    n <- nrow(data)
    xe <- pda_colon(data[1, 1], data[n, 1], rate)
    ye <- stats::approx(data[, 1], data[, 2], xout = xe, method = "linear",
        rule = 1, ties = "ordered")$y
    cbind(xe, ye)
  }

  # @keywords internal
  # @noRd
  pda_zeropad <- function(data, win, padding = 1){
    # zeropad2.m (ensureHalfWindowCoverage = false)
    data <- data[order(data[, 1], data[, 2]), 1:2, drop = FALSE]
    x <- data[, 1]
    y <- data[, 2]
    ny <- length(y)
    y <- y - sum(y) / ny
    dx <- diff(x)
    dt <- sum(dx) / length(dx)
    n <- pda_round(win / 2 / dt)
    if (n > 0){
      x1x <- x[1] - (n:1) * dt
      x2x <- x[ny] + (1:n) * dt
    } else {
      x1x <- numeric(0)
      x2x <- numeric(0)
    }
    if (padding == 1){
      x1y <- rep(0, n)
      x2y <- rep(0, n)
    } else if (padding == 2){
      if (n > ny){
        stop("Mirror padding cannot exceed the input data length")
      }
      if (n > 0){
        x1y <- y[n:1]
        x2y <- y[ny:(ny - n + 1)]
      } else {
        x1y <- numeric(0)
        x2y <- numeric(0)
      }
    } else if (padding == 3){
      m1 <- sum(y[seq_len(n)]) / n
      m2 <- sum(y[(ny - n + 1):ny]) / n
      x1y <- rep(m1, n)
      x2y <- rep(m2, n)
    } else if (padding == 4){
      m1 <- sum(y[seq_len(n)]) / n
      m2 <- sum(y[(ny - n + 1):ny]) / n
      s1 <- stats::sd(y[seq_len(n)])
      s2 <- stats::sd(y[(ny - n + 1):ny])
      x1y <- stats::rnorm(n) * s1 + m1
      x2y <- stats::rnorm(n) * s2 + m2
    } else {
      stop("padding must be either 1, 2, 3, or 4")
    }
    rbind(cbind(x1x, x1y), cbind(x, y), cbind(x2x, x2y), deparse.level = 0)
  }

  # @keywords internal
  # @noRd
  pda_band_indices <- function(bands, n_bin, nyquist){
    # legacyBandIndices of pdan.m
    idx <- integer(0)
    for (b in seq_len(nrow(bands))){
      lower <- ceiling(n_bin * bands[b, 1] / nyquist)
      upper <- trunc(n_bin * bands[b, 2] / nyquist)
      if (lower == 0){
        lower <- 1
      }
      lower <- max(1, lower)
      upper <- min(n_bin, upper)
      if (upper >= lower){
        idx <- c(idx, lower:upper)
      }
    }
    sort(unique(idx))
  }

  # @keywords internal
  # @noRd
  pda_pdan <- function(data, f3, window, nw = 2, ftmin = 0, fterm = NULL,
      step = 1, pad = 1000, verbose = FALSE){
    # pdan.m
    data <- as.matrix(data)
    if (ncol(data) != 2 || nrow(data) < 2 || any(!is.finite(data))){
      stop("data must be a finite N x 2 matrix")
    }
    coordinate <- as.numeric(data[, 1])
    dt <- coordinate[2] - coordinate[1]
    spacing <- diff(coordinate)
    if (any(spacing <= 0) || max(abs(spacing - dt)) / dt > 1e-5){
      stop("data coordinates must be strictly increasing and regularly sampled (relative tolerance 1e-5)")
    }
    nyquist <- 1 / (2 * dt)
    if (is.null(fterm)){
      fterm <- nyquist
    }
    if (!(is.finite(window) && window > 0)){
      stop("window must be positive")
    }
    if (!(is.finite(nw) && nw > 0)){
      stop("nw must be positive")
    }
    if (!(is.finite(step) && step > 0 && step == trunc(step))){
      stop("step must be a positive integer (number of samples)")
    }
    if (!(is.finite(pad) && pad > 0 && pad == trunc(pad))){
      stop("pad must be a positive integer")
    }
    f3 <- as.numeric(f3)
    if (length(f3) < 2 || length(f3) %% 2 != 0 || any(!is.finite(f3))){
      stop("target bands must be a finite vector of paired limits")
    }
    window_samples <- trunc(window / dt)
    n_rows <- nrow(data)
    if (window_samples < 2 || window_samples > n_rows){
      stop("window must resolve to between 2 and the number of rows samples")
    }
    if (nw < 0.5 || nw >= window_samples / 2){
      stop("nw must be at least 0.5 and less than window_samples/2")
    }
    fterm <- min(fterm, nyquist)
    if (ftmin < 0 || fterm <= ftmin){
      stop("Frequency limits must satisfy 0 <= f_min < f_max after clipping")
    }
    bands <- matrix(f3, ncol = 2, byrow = TRUE)
    lo <- pmin(bands[, 1], bands[, 2])
    hi <- pmax(bands[, 1], bands[, 2])
    bands <- cbind(lo, hi)
    bands[bands > nyquist] <- nyquist
    if (any(bands[, 1] < 0)){
      stop("Target-frequency limits must be nonnegative")
    }
    start_idx <- seq(1, n_rows - window_samples + 1, by = step)
    n_win <- length(start_idx)
    n_bin <- floor(pad / 2) + 1
    total_idx <- pda_band_indices(matrix(c(ftmin, fterm), 1, 2), n_bin, nyquist)
    target_idx <- pda_band_indices(bands, n_bin, nyquist)
    if (length(total_idx) == 0){
      stop("The total frequency range contains no FFT bins")
    }
    values <- as.numeric(data[, 2])
    target_power <- numeric(n_win)
    total_power <- numeric(n_win)
    if (verbose == TRUE){
      pb <- utils::txtProgressBar(min = 0, max = n_win, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    for (w in seq_len(n_win)){
      first <- start_idx[w]
      segment <- acy_detrend_linear(values[first:(first + window_samples - 1)])
      spectrum <- acy_pmtm(segment, nw = nw, nfft = pad)$pxx
      total_power[w] <- sum(spectrum[total_idx])
      if (length(target_idx) > 0){
        target_power[w] <- sum(spectrum[target_idx])
      }
      if (verbose == TRUE){
        utils::setTxtProgressBar(pb, w)
      }
    }
    if (any(!is.finite(total_power)) || any(total_power <= 0)){
      stop("A window has zero or nonfinite total-band power")
    }
    ratio <- target_power / total_power
    centers <- pda_linspace(coordinate[1] + window / 2,
        coordinate[n_rows] - window / 2, n_win)
    pow <- cbind(centers, ratio, target_power, total_power)
    colnames(pow) <- c("center", "ratio", "target_power", "total_power")
    pow
  }

  # Slepian (DPSS) tapers, MATLAB dpss() conventions (unit energy)

  acy_dpss_cache <- new.env(parent = emptyenv())

  # Main code

  data <- as.matrix(data)
  if (!is.numeric(data) || ncol(data) < 2){
    stop("data must be a numeric matrix with two columns [depth/time, value]")
  }
  data <- data[, 1:2, drop = FALSE]
  data <- data[!(is.na(data[, 1]) | is.na(data[, 2])), , drop = FALSE]
  if (any(diff(data[, 1]) < 0)){
    data <- data[order(data[, 1], data[, 2]), , drop = FALSE]
  }
  if (any(diff(data[, 1]) == 0)){
    data <- pda_findduplicate(data)
  }
  if (pda_is_uneven(data[, 1])){
    data <- pda_interpolate(data, stats::median(diff(data[, 1])))
  }
  if (!(pad_type %in% c(0, 1, 2, 3, 4))){
    stop("pad_type must be 0, 1, 2, 3 or 4")
  }
  if (pad_type > 0){
    data <- pda_zeropad(data, window, pad_type)
  }
  if (is.null(f_max)){
    f_max <- 1 / (2 * (data[2, 1] - data[1, 1]))
  }
  pow <- pda_pdan(data, target_bands, window, nw, f_min, f_max, step, pad,
      verbose = verbose)

  if (genplot == TRUE){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::plot(pow[, 1], pow[, 2], type = "l", xlab = "depth/time",
        ylab = "Power ratio")
  }
  pow
}
