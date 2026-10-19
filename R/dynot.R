#' @title Dynamic noise after orbital tuning (DYNOT)
#'
#' @description The \code{\link{dynot}} function runs the DYNOT sedimentary
#' noise model of Li et al. (2018a,b) on an (orbitally tuned, usually
#' unevenly sampled) time series. In each Monte Carlo realisation a sampling
#' interval is drawn from a Weibull distribution fitted (maximum likelihood)
#' to the sampling intervals of the data (values outside
#' [\code{sample_min}, \code{sample_max}] are replaced by one shared uniform
#' draw for all low values and one for all high values, as in Acycle), a
#' window length is drawn uniformly from [\code{window_min}, \code{window_max}],
#' a time-bandwidth product from \code{nw_min + k/2} and, for each target
#' orbital frequency f, the lower and upper band limits are
#' f -/+ fz * nw / window with fz uniform in [\code{fz_min}, \code{fz_max}].
#' The series is pchip-interpolated at the drawn sampling interval, linearly
#' detrended, padded by half a window (\code{pad_type}) and analysed with the
#' legacy power decomposition (\code{pdan.m}, see \code{\link{pda_legacy}}).
#' The power-ratio curves are interpolated onto a regular grid, rescaled by
#' (grand mean / realisation mean) of the power ratio, normalised to a
#' maximum of 1 if needed, and summarized by percentiles (median plus the
#' selected confidence intervals). The noise is 1 - power ratio.
#' Translated from Acycle v3.1 DYNOS.m (runDynot, loadDataReady,
#' sanitizeData), dynotIterationWorker.m, dynotPercentilePlan.m, pdan.m and
#' zeropad2.m (Mingsong Li).
#'
#' @param data Matrix with two columns: time (kyr) and proxy value. Data are
#' sorted, duplicated times are averaged and non-finite rows removed.
#' \code{Default=NULL}
#' @param target_periods Target orbital periods (same unit as time). Default
#' are the nine Acycle default periods (405.6912, 130.6979, 123.8532,
#' 98.8517, 94.8856, 40.9897, 23.6820, 22.3758, 18.9519 kyr).
#' @param window_min Minimum window length. \code{Default=300}
#' @param window_max Maximum window length. \code{Default=500}
#' @param nw_min Minimum time-bandwidth product. \code{Default=2}
#' @param nw_max Maximum time-bandwidth product; 2*(nw_max-nw_min) must be an
#' integer. \code{Default=2}
#' @param pad FFT length of the MTM (zero-padding number). \code{Default=1000}
#' @param step Step between PDA windows in samples. \code{Default=5}
#' @param n_out Number of regular output grid points. \code{Default=1000}
#' @param pad_type Half-window padding: 0 = none, 1 = zero, 2 = mirror,
#' 3 = mean, 4 = random. \code{Default=1}
#' @param fz_min Minimum bandwidth multiplier. \code{Default=0.9}
#' @param fz_max Maximum bandwidth multiplier. \code{Default=1.2}
#' @param f_min Lower cutoff frequency (clips target limits and total band).
#' \code{Default=0.001}
#' @param f_max Upper cutoff frequency. \code{Default=1}
#' @param sample_min Lower limit of simulated sampling intervals. If NULL the
#' 5th percentile of the data sampling intervals. \code{Default=NULL}
#' @param sample_max Upper limit of simulated sampling intervals. If NULL the
#' 95th percentile of the data sampling intervals. \code{Default=NULL}
#' @param intervals Confidence intervals (in percent) to report, any of
#' 50, 68, 80, 90, 95; the median is always reported.
#' \code{Default=c(50, 68, 80, 90, 95)}
#' @param n_sim Number of Monte Carlo simulations. \code{Default=1000}
#' @param seed Seed for the random number generator (NULL = not set).
#' \code{Default=NULL}
#' @param verbose Show a progress bar. \code{Default=FALSE}
#' @param genplot Plot the result with \code{\link{plot_dynot}}.
#' \code{Default=FALSE}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' (Acycle).
#'
#' @references
#' Li, M., Hinnov, L.A., Huang, C., Ogg, J.G., 2018a. Sedimentary noise and
#' sea levels linked to land-ocean water exchange and obliquity forcing.
#' Nature Communications 9, 1004. \doi{10.1038/s41467-018-03454-y}
#'
#' Li, M., Huang, C., Hinnov, L., Chen, W., Ogg, J., Tian, W., 2018b.
#' Astrochronology of the Anisian stage (Middle Triassic) at the Guandao
#' reference section, South China. Earth and Planetary Science Letters 482,
#' 591-606. \doi{10.1016/j.epsl.2017.11.042}
#'
#' Li, M., Huang, C., Hinnov, L., Ogg, J., Chen, Z.-Q., Zhang, Y., 2016.
#' Obliquity-forced climate during the Early Triassic hothouse in China.
#' Geology 44, 623-626. \doi{10.1130/G37970.1}
#'
#' @examples
#' \donttest{
#' set.seed(2)
#' t <- cumsum(runif(900, 1, 3))
#' x <- sin(2 * pi * t / 405) + sin(2 * pi * t / 100) +
#'   sin(2 * pi * t / 41) + rnorm(900) * (0.2 + t / max(t))
#' res <- dynot(data = cbind(t, x), n_sim = 20, n_out = 200, seed = 1)
#' plot_dynot(res)
#' }
#' @return A list with elements
#' \code{percentile}: matrix [time, noise percentiles] on the grid points
#' where the median is defined (column names give the percentile levels);
#' \code{median}: matrix [time, median noise];
#' \code{percentiles}: percentile levels; \code{median_index};
#' \code{interval_count}: number of confidence intervals;
#' \code{grid}: full output grid; \code{noise}: matrix (n_out x n_sim) of
#' adjusted noise (1 - adjusted power ratio) of every realisation;
#' \code{power_mean}: mean power ratio of every realisation;
#' \code{draws}: list with the simulated sampling intervals, windows, nw,
#' band limits (matrix f3m) and the underlying uniform draws;
#' \code{weibull}: fitted Weibull scale and shape; \code{sample_range};
#' \code{data}: the cleaned input data.
#'
#' @export
#' @importFrom stats runif
#' @importFrom stats quantile
#' @importFrom stats approx
#' @importFrom stats uniroot
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar

dynot <- function(data = NULL,
                  target_periods = c(405.6912, 130.6979, 123.8532, 98.8517,
                                     94.8856, 40.9897, 23.6820, 22.3758,
                                     18.9519),
                  window_min = 300,
                  window_max = 500,
                  nw_min = 2,
                  nw_max = 2,
                  pad = 1000,
                  step = 5,
                  n_out = 1000,
                  pad_type = 1,
                  fz_min = 0.9,
                  fz_max = 1.2,
                  f_min = 0.001,
                  f_max = 1,
                  sample_min = NULL,
                  sample_max = NULL,
                  intervals = c(50, 68, 80, 90, 95),
                  n_sim = 1000,
                  seed = NULL,
                  verbose = FALSE,
                  genplot = FALSE){
  # ---- internal helper functions (local to this function) ----

  # ---------------------------------------------------------------------------
  # small utilities
  # ---------------------------------------------------------------------------

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

  # ---------------------------------------------------------------------------
  # MATLAB periodogram / pmtm (one-sided, adaptive weights)
  # ---------------------------------------------------------------------------

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

  # ---------------------------------------------------------------------------
  # internal helpers
  # ---------------------------------------------------------------------------

  # @keywords internal
  # @noRd
  dynot_sanitize <- function(data){
    # DYNOS.m sanitizeData: sort, average duplicates, drop non-finite rows
    data <- as.matrix(data)
    if (!is.numeric(data) || ncol(data) < 2){
      stop("data must be a numeric matrix with two columns [time, value]")
    }
    data <- data[, 1:2, drop = FALSE]
    data <- data[order(data[, 1], na.last = TRUE), , drop = FALSE]
    fin_x <- !is.na(data[, 1])
    good <- data[fin_x, , drop = FALSE]
    if (nrow(good) > 1 && any(diff(good[, 1]) == 0)){
      good <- pda_findduplicate(good)
    }
    data <- rbind(good, data[!fin_x, , drop = FALSE])
    keep <- is.finite(data[, 1]) & is.finite(data[, 2])
    out <- data[keep, , drop = FALSE]
    dimnames(out) <- NULL
    out
  }

  # @keywords internal
  # @noRd
  dynot_prctile <- function(x, p){
    # MATLAB prctile (midpoint rule, NaN ignored) == quantile type 5
    x <- x[!is.na(x)]
    if (length(x) == 0){
      return(rep(NA_real_, length(p)))
    }
    as.numeric(stats::quantile(x, probs = p / 100, type = 5, names = FALSE))
  }

  # @keywords internal
  # @noRd
  dynot_row_prctile <- function(mat, p){
    out <- matrix(NA_real_, nrow(mat), length(p))
    for (r in seq_len(nrow(mat))){
      out[r, ] <- dynot_prctile(mat[r, ], p)
    }
    out
  }

  # @keywords internal
  # @noRd
  dynot_wblfit <- function(x){
    # Weibull maximum likelihood (MATLAB wblfit via the extreme-value fit of
    # log(x)); solved to near machine precision instead of TolX = 1e-6
    if (any(x <= 0)){
      stop("Sampling intervals must be positive for the Weibull fit")
    }
    lx <- log(x)
    x_max <- max(lx)
    x_range <- x_max - min(lx)
    if (x_range == 0){
      return(c(exp(lx[1]), Inf))
    }
    x0 <- (lx - x_max) / x_range
    mean_x0 <- sum(x0) / length(x0)
    lkeq <- function(s){
      w <- exp(x0 / s)
      s + mean_x0 - sum(x0 * w) / sum(w)
    }
    upper <- 1
    while (lkeq(upper) <= 0){
      upper <- upper * 2
    }
    lower <- upper / 2
    while (lkeq(lower) > 0){
      lower <- lower / 2
    }
    s <- stats::uniroot(lkeq, c(lower, upper), tol = 1e-15, maxiter = 2000)$root
    mu <- s * log(sum(exp(x0 / s)) / length(x0))
    c(exp(x_range * mu + x_max), 1 / (x_range * s))
  }

  # @keywords internal
  # @noRd
  dynot_wblinv <- function(u, scale, shape){
    # MATLAB wblinv: scale * (-log(1-u))^(1/shape)
    scale * (-log(1 - u))^(1 / shape)
  }

  # @keywords internal
  # @noRd
  dynot_percentile_plan <- function(intervals){
    # dynotPercentilePlan.m
    endpoints <- list(c(25, 75), c(15.865, 84.135), c(10, 90), c(5, 95),
                      c(2.5, 97.5))
    levels <- c(50, 68, 80, 90, 95)
    pp <- 50
    count <- 0
    for (i in seq_along(levels)){
      if (levels[i] %in% intervals){
        pp <- c(pp, endpoints[[i]])
        count <- count + 1
      }
    }
    pp <- sort(pp)
    list(percentiles = pp, interval_count = count,
         median_index = which(pp == 50)[1])
  }

  # @keywords internal
  # @noRd
  dynot_interp1 <- function(x, y, xi){
    # MATLAB interp1 linear, NaN outside the range
    if (length(x) < 2){
      stop("At least two points are needed for interpolation")
    }
    stats::approx(x, y, xout = xi, method = "linear", rule = 1,
                  ties = "ordered")$y
  }

  # @keywords internal
  # @noRd
  dynot_iteration <- function(data, y_grid, sample, window, nw, f3, pad_type,
                              f_min, f_max, step, pad){
    # dynotIterationWorker.m
    n <- nrow(data)
    coords <- pda_colon(data[1, 1], data[n, 1], sample)
    dat <- cbind(coords, acy_pchip(data[, 1], data[, 2], coords))
    beta <- qr.solve(cbind(coords, 1), dat[, 2])
    dat[, 2] <- dat[, 2] - (beta[1] * coords + beta[2])
    if (pad_type > 0){
      dat <- pda_zeropad(dat, window, pad_type)
    }
    power <- pda_pdan(dat, f3, window, nw, f_min, f_max, step, pad)
    column <- dynot_interp1(power[, 1], power[, 2], y_grid)
    vals <- power[, 2]
    list(column = column, mean = mean(vals[!is.na(vals)]))
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

  # ---------------------------------------------------------------------------
  # Slepian (DPSS) tapers, MATLAB dpss() conventions (unit energy)
  # ---------------------------------------------------------------------------

  acy_dpss_cache <- new.env(parent = emptyenv())

  # ---- main code ----

  data <- dynot_sanitize(data)
  if (nrow(data) < 10){
    stop("DYNOT needs at least 10 data points")
  }
  sample_rate <- diff(data[, 1])
  parmhat <- dynot_wblfit(sample_rate)
  sample_range <- dynot_prctile(sample_rate, c(5, 95))
  if (is.null(sample_min)){
    sample_min <- sample_range[1]
  }
  if (is.null(sample_max)){
    sample_max <- sample_range[2]
  }
  if (!(is.finite(n_sim) && n_sim >= 1 && n_sim == trunc(n_sim))){
    stop("n_sim must be a positive integer")
  }
  if (!(is.finite(n_out) && n_out >= 1 && n_out == trunc(n_out))){
    stop("n_out must be a positive integer")
  }
  if (!(pad_type %in% c(0, 1, 2, 3, 4))){
    stop("pad_type must be 0, 1, 2, 3 or 4")
  }
  if (any(!is.finite(target_periods)) || any(target_periods <= 0)){
    stop("Target orbital cycles must be finite positive values")
  }
  if (sample_min <= 0 || sample_max <= 0 || sample_max < sample_min){
    stop("Sample rate range is invalid")
  }
  if (!all(intervals %in% c(50, 68, 80, 90, 95))){
    stop("intervals must be a subset of c(50, 68, 80, 90, 95)")
  }
  f <- 1 / target_periods

  # random draws (order as in DYNOS.m runDynot)
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (nw_max == nw_min){
    randnw <- rep(0, n_sim)
  } else {
    k_nw <- 2 * (nw_max - nw_min)
    if (!(k_nw >= 1 && abs(k_nw - round(k_nw)) < 1e-12)){
      stop("2*(nw_max - nw_min) must be a positive integer")
    }
    randnw <- sample.int(round(k_nw), n_sim, replace = TRUE)
  }
  u_weibull <- stats::runif(n_sim)
  r_lower <- stats::runif(1)
  r_upper <- stats::runif(1)
  u_window <- stats::runif(n_sim)
  u_fz <- matrix(0, n_sim, 2 * length(f))
  for (i in seq_along(f)){
    u_fz[, 2 * i - 1] <- stats::runif(n_sim)
    u_fz[, 2 * i] <- stats::runif(n_sim)
  }

  samplez <- dynot_wblinv(u_weibull, parmhat[1], parmhat[2])
  lower_replacement <- sample_min + (sample_max - sample_min) * r_lower
  samplez[samplez < sample_min] <- lower_replacement
  upper_replacement <- sample_min + (sample_max - sample_min) * r_upper
  samplez[samplez > sample_max] <- upper_replacement
  windowz <- window_min + (window_max - window_min) * u_window
  nwz <- nw_min + randnw / 2
  bw <- nwz / windowz
  f3m <- matrix(0, n_sim, 2 * length(f))
  for (i in seq_along(f)){
    fz <- fz_min + (fz_max - fz_min) * u_fz[, 2 * i - 1]
    f3m[, 2 * i - 1] <- f[i] - fz * bw
    fz <- fz_min + (fz_max - fz_min) * u_fz[, 2 * i]
    f3m[, 2 * i] <- f[i] + fz * bw
  }
  f3m[f3m < f_min] <- f_min
  f3m[f3m > f_max] <- f_max

  n_data <- nrow(data)
  y_grid <- pda_linspace(data[1, 1], data[n_data, 1], n_out)
  powy <- matrix(NA_real_, n_out, n_sim)
  powmean <- rep(NA_real_, n_sim)

  if (verbose == TRUE){
    pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  for (i in seq_len(n_sim)){
    it <- dynot_iteration(data, y_grid, samplez[i], windowz[i], nwz[i],
                          f3m[i, ], pad_type, f_min, f_max, step, pad)
    powy[, i] <- it$column
    powmean[i] <- it$mean
    if (verbose == TRUE){
      utils::setTxtProgressBar(pb, i)
    }
  }

  pm <- mean(powmean[!is.na(powmean)])
  adj <- pm / powmean
  powyadjust <- powy * matrix(adj, n_out, n_sim, byrow = TRUE)
  if (any(!is.na(powyadjust))){
    maxp <- max(powyadjust, na.rm = TRUE)
    if (is.finite(maxp) && maxp > 1){
      powyadjust <- powyadjust / maxp
    }
  }

  plan <- dynot_percentile_plan(intervals)
  pp <- plan$percentiles
  powyadjustp <- dynot_row_prctile(powyadjust, pp)
  mask <- !is.na(powyadjustp[, 1])
  y_grid_nan <- y_grid[mask]
  powyad_p_nan <- 1 - powyadjustp[mask, , drop = FALSE]

  percentile <- cbind(y_grid_nan, powyad_p_nan)
  colnames(percentile) <- c("x", paste0("p", pp))
  med <- cbind(y_grid_nan, powyad_p_nan[, plan$median_index])
  colnames(med) <- c("x", "median")

  res <- list(percentile = percentile,
              median = med,
              percentiles = pp,
              median_index = plan$median_index,
              interval_count = plan$interval_count,
              grid = y_grid,
              noise = 1 - powyadjust,
              power_mean = powmean,
              draws = list(samplez = samplez, windowz = windowz, nwz = nwz,
                           f3m = f3m, randnw = randnw, u_weibull = u_weibull,
                           r_lower = r_lower, r_upper = r_upper,
                           u_window = u_window, u_fz = u_fz),
              weibull = c(scale = parmhat[1], shape = parmhat[2]),
              sample_range = c(sample_min, sample_max),
              data = data)
  if (genplot == TRUE){
    plot_dynot(res)
  }
  res
}
