#' @title Sliding-window spectral moments (mean frequency and bandwidth)
#'
#' @description The \code{\link{spectral_moments}} function calculates the
#' first two spectral moments, the mean frequency (\eqn{\mu_f}) and the
#' bandwidth (\eqn{B}), of a uniformly sampled series in a sliding window
#' (Sinnesael et al., 2018). Because both moments scale with the
#' sedimentation rate, their evolution along a stratigraphic series can be
#' used to detect changes in sedimentation rate.
#'
#' Within every window the (globally mean-centred) segment is linearly
#' detrended, its one-sided rectangular-window periodogram is computed on
#' an even NFFT grid from DC up to the Nyquist frequency and
#' \deqn{\mu_f = \sum f P(f) / \sum P(f), \quad
#' B = \sqrt{\sum (f-\mu_f)^2 P(f) / \sum P(f)}}
#' are calculated.
#'
#' Two computational paths are available, mirroring the Acycle GUI
#' (SpectralMomentsGUI.m):
#' \itemize{
#' \item the strict deterministic core (default; used whenever
#' \code{sr_mean=NULL} and \code{edge_padding} is \code{"none"} or
#' \code{"zero"}). Translated from Acycle v3.1 acycleSpectralMoments.m
#' (Mingsong Li). The window length is rounded (half away from zero) to an
#' integer number of sampling intervals, the step to an integer number of
#' samples, NFFT = \code{pad} x window samples (+1 if odd), and window
#' centres are the midpoints of the actual window coordinates.
#' \item the legacy path (used when \code{sr_mean} is given or when
#' \code{edge_padding} is \code{"mirror"} or \code{"mean"}). Translated from
#' Acycle v3.1 spectralmoments.m and zeropad2.m (Mingsong Li, based on
#' SpectralMoments_with_SR_estimation.m and DetrendSignalPP.m by Matthias
#' Sinnesael). When \code{sr_mean} is given the long-term bandwidth trend is
#' estimated with the overlap-add polynomial baseline model of
#' DetrendSignalPP.m (Acycle's 'poly' option) and converted into a relative
#' sedimentation-rate curve scaled to \code{sr_mean}
#' (\eqn{sr = sr_{mean} B_1 / B_{trend}}, re-centred on \code{sr_mean}). As in
#' Acycle, the trend window is \code{round(window/step_samples)} bandwidth
#' values, i.e. \code{window} expressed in the unit of the depth column is
#' used as a number of samples, so the smoothing depends on the depth unit.
#' }
#'
#' @param data Uniformly sampled record (required): matrix or data frame with depth or time in column 1 (strictly
#' increasing, all values finite) and the proxy in column 2; at least three rows. \code{Default=NULL}.
#' @param window Length of the running window, in the unit of the first column of \code{data} (numeric, > 0). In the
#' strict core it is rounded to an integer number of sampling intervals. \code{NULL} uses 25\% of the record length
#' (the Acycle GUI default). \code{Default=NULL}.
#' @param step Step of the running window, in the unit of the first column (numeric, > 0); rounded to an integer number
#' of samples in the strict core. \code{NULL} uses the median sampling interval. \code{Default=NULL}.
#' @param pad Zero-padding factor of the periodogram (numeric, >= 1): NFFT is \code{pad} times the number of samples in
#' a window (+1 if odd). Larger values give a finer frequency grid. \code{Default=5}.
#' @param edge_padding Padding of the record edges by half a window (character, case-insensitive): \code{"none"} or
#' \code{"zero"} (zeros; strict core), \code{"mirror"} (mirrored data) or \code{"mean"} (the mean value) (both use the
#' legacy path, as in the Acycle GUI). \code{Default="none"}.
#' @param sr_mean Mean sedimentation rate (numeric, e.g. in cm/kyr). When given the legacy path is used and the
#' bandwidth trend and an Acycle-style sedimentation rate curve scaled to \code{sr_mean} are returned; \code{NULL}
#' computes the moments only. \code{Default=NULL}.
#' @param freq_range Frequency band over which the moments are taken, in cycles per depth unit (numeric vector
#' \code{c(lower, upper)}); \code{NULL} uses the whole spectrum up to the Nyquist frequency, as in the paper. The
#' moments of a record with white noise are biased towards the cutoff frequency (section 2.3.3 of the paper), so a band
#' around the astronomical lines makes them follow the sedimentation rate rather than the noise. \code{Default=NULL}.
#' @param noise_floor Subtract the median of every window's periodogram (an estimate of the white-noise floor) before
#' the moments are taken; negative values are set to zero (logical). In the synthetic test of the examples (precession,
#' obliquity, eccentricity, white noise 0.3 times the signal, rate varying 1-3 cm/kyr, Nyquist 25 cycles/m) the
#' anchored bandwidth trend correlates 0.00 with the true rate over the whole spectrum (the noise plateau up to the
#' Nyquist frequency dominates the second moment), 0.97 over \code{freq_range = c(0.1, 8)}, and 0.99 with the floor
#' removed as well (91 \% of the record within 20 \%); the mean frequency is less sensitive (0.98 over the whole
#' spectrum). \code{Default=FALSE}.
#' @param anchor Anchor point that turns the smoothed trend of the moment into a sedimentation rate curve, as in
#' Sinnesael et al. (2018): every spatial frequency, and so every moment, scales with 1 / sedimentation rate (their
#' equation 12). Give either a numeric vector \code{c(depth, period_kyr, frequency)} (the component with period
#' \code{period_kyr} sits at \code{frequency} cycles per depth unit at \code{depth}, so the rate there is
#' \code{1 / (period_kyr * frequency)}), or \code{list(depth = , rate = )} with the rate in depth units per kyr. The
#' result then holds \code{moment_trend} and \code{sedrate_anchored}. This is the procedure of the paper (U1437: the 0.45
#' cycle/m precession line read as 21 kyr; ODP 982 and the synthetic insolation: the obliquity line at the top of the
#' record), as opposed to the Acycle \code{sr_mean} path, which rescales by the first raw bandwidth value and
#' re-centres the curve additively on a mean rate. \code{NULL} gives no anchored rate. \code{Default=NULL}.
#' @param moment Which moment carries the sedimentation rate when \code{anchor} is given (character):
#' \code{"bandwidth"} (the paper's choice, less sensitive to the noise bias than the mean frequency) or
#' \code{"mean_frequency"}. \code{Default=c("bandwidth", "mean_frequency")}, i.e. \code{"bandwidth"}.
#' @param trend Smoothing of the moment before anchoring (character): \code{"poly"} (the overlap-add polynomial
#' baseline of Verhelst 2000 used in the paper; falls back to lowess when the trend window is not shorter than the
#' series) or \code{"lowess"} (\code{\link[stats]{lowess}} with span \code{trend_window} / record length).
#' \code{Default=c("poly", "lowess")}, i.e. \code{"poly"}.
#' @param trend_window Length of the smoothing window of the trend, in depth units (numeric, > 0); only used with
#' \code{anchor}. The paper takes a long-term trend: the result is the first-order change of the rate, not its detail.
#' \code{NULL} uses three analysis windows (\code{3 * window}). \code{Default=NULL}.
#' @param trend_order Polynomial order of the \code{"poly"} baseline (integer, >= 0). \code{Default=2}.
#' @param genplot Plot the result with \code{\link{plot_spectral_moments}} (logical). \code{Default=FALSE}.
#' @param verbose Show a progress bar (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' (Acycle), Matthias Sinnesael (spectral moments and DetrendSignalPP), Milan
#' Zivanovic.
#'
#' @references
#' Sinnesael, M., Zivanovic, M., De Vleeschouwer, D., Claeys, P. (2018).
#' Spectral moments in cyclostratigraphy: advantages and disadvantages
#' compared to more classic approaches. Paleoceanography and
#' Paleoclimatology, 33, 493-510. <doi:10.1029/2017PA003293>
#'
#' Li, M., Hinnov, L., Kump, L. (2019). Acycle: Time-series analysis software
#' for paleoclimate research and education. Computers & Geosciences, 127,
#' 12-22. <doi:10.1016/j.cageo.2019.02.011>
#'
#' @examples
#' \donttest{
#' # Two cycles (40 and 10 cm period in "age" units) recorded with a
#' # sedimentation rate that doubles halfway the record (depth in cm): the
#' # mean frequency halves in the upper part.
#' dt <- 0.5
#' depth <- seq(0, 1000, by = dt)
#' age <- ifelse(depth <= 500, depth, 500 + (depth - 500) / 2)
#' set.seed(1)
#' proxy <- sin(2 * pi * age / 40) + 0.5 * sin(2 * pi * age / 10) +
#'   stats::rnorm(length(depth), sd = 0.2)
#' sm <- spectral_moments(data = cbind(depth, proxy), window = 200,
#'                        step = 5, pad = 5, genplot = FALSE)
#' plot_spectral_moments(sm, x_lab = "depth (cm)")
#'
#' # with a sedimentation-rate estimate (legacy path, 'poly' trend)
#' sm_sr <- spectral_moments(data = cbind(depth, proxy), window = 200,
#'                           step = 5, pad = 5, sr_mean = 1.5)
#' plot_spectral_moments(sm_sr, x_lab = "depth (cm)")
#' # the sedimentation rate of Sinnesael et al. (2018): the obliquity line is identified at
#' # 5 m (41 kyr at f_obl cycles/m), the bandwidth trend is anchored there; moments over a band
#' set.seed(3)
#' t <- seq(0, 3000, 0.5)
#' x <- 0.8 * sin(2 * pi * t / 41) + sin(2 * pi * t / 23.7) + 0.5 * sin(2 * pi * t / 405) +
#'   0.3 * rnorm(length(t))
#' rate <- 0.02 + 0.01 * sin(2 * pi * t / 1500)
#' z <- cumsum(c(0, rate[-1] * diff(t)))
#' zg <- seq(0, max(z), 0.02)
#' d <- data.frame(zg, approx(z, x, zg)$y)
#' f_obl <- 1 / (41 * approx(z, rate, 5)$y)     # the obliquity line at 5 m depth
#' sm <- spectral_moments(d, window = 5, step = 0.1, anchor = c(5, 41, f_obl),
#'                        freq_range = c(0.1, 8), noise_floor = TRUE)
#' plot_spectral_moments(sm, freq_max = 8)
#' }
#'
#' @return A list with elements
#' \itemize{
#' \item \code{result}: matrix with columns \code{depth} (window centre),
#' \code{mean_frequency}, \code{bandwidth} and, when \code{sr_mean} was
#' given, \code{bandwidth_trend} and \code{sedrate}.
#' \item \code{method}: \code{"strict"} or \code{"legacy"}.
#' \item \code{meta}: list with the effective settings (sampling interval,
#' window/step in samples and physical units, NFFT, padding, Nyquist, number
#' of windows).
#' \item with \code{anchor}: \code{result} also holds \code{moment_trend} and \code{sedrate_anchored}, and \code{anchor}
#' records the anchor used.
#' \item \code{spectra}: the periodogram of every window (\code{frequency}, \code{depth} of the window centres,
#' \code{power} matrix frequency x window), the evolutive spectrum that \code{\link{plot_spectral_moments}} draws
#' as an image with the two moments on top.
#' }
#'
#' @export
#' @importFrom stats median
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
spectral_moments <- function(data = NULL,
    window = NULL,
    step = NULL,
    pad = 5,
    edge_padding = "none",
    sr_mean = NULL,
    freq_range = NULL,
    noise_floor = FALSE,
    anchor = NULL,
    moment = c("bandwidth", "mean_frequency"),
    trend = c("poly", "lowess"),
    trend_window = NULL,
    trend_order = 2,
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

  acy_periodogram <- function(x, nfft = NULL, fs = NULL){
    x <- as.numeric(x)
    n <- length(x)
    if (is.null(nfft)){
      nfft <- max(256, 2^acy_nextpow2(n))
    }
    X <- stats::fft(acy_datawrap(x, nfft))
    S <- Mod(X)^2 / n
    P <- acy_onesided(S, nfft)
    if (is.null(fs)){
      P <- P / (2 * pi)
    } else {
      P <- P / fs
    }
    list(pxx = P, f = acy_freq_vector(nfft, fs))
  }

  # internal helpers

  # @keywords internal
  # @noRd
  smom_round <- function(x){
    # MATLAB round(): half away from zero
    sign(x) * floor(abs(x) + 0.5)
  }

  # @keywords internal
  # @noRd
  smom_eps <- function(x){
    # MATLAB eps(x) for finite double x
    x <- abs(x)
    if (x == 0){
      return(2^-1074)
    }
    e <- floor(log2(x))
    if (2^e > x){
      e <- e - 1
    } else if (2^(e + 1) <= x){
      e <- e + 1
    }
    2^(max(e, -1022) - 52)
  }

  # @keywords internal
  # @noRd
  smom_core <- function(data, window, step, pad, edge_padding, verbose = FALSE){
    # strict core, acycleSpectralMoments.m
    x <- as.numeric(data[, 1])
    y <- as.numeric(data[, 2])
    n_rows <- length(x)
    spacing <- diff(x)
    if (any(spacing <= 0)){
      stop("Coordinates must be finite, unique, and strictly increasing.")
    }
    dt <- stats::median(spacing)
    rel_dep <- max(abs(spacing - dt)) / dt
    if (rel_dep > 1e-5){
      stop(sprintf(paste0("Relative spacing departure %.17g exceeds the ",
          "near-uniform tolerance 1e-5; interpolate first."),
          rel_dep))
    }
    nyq <- 1 / (2 * dt)
    x_start <- x[1]
    x_end <- x[n_rows]
    span <- x_end - x_start

    if (!(is.numeric(window) && length(window) == 1 && is.finite(window) &&
        window > 0)){
      stop("window must be a finite positive numeric scalar.")
    }
    if (!(is.numeric(step) && length(step) == 1 && is.finite(step) &&
        step > 0)){
      stop("step must be a finite positive numeric scalar.")
    }
    if (!(is.numeric(pad) && length(pad) == 1 && is.finite(pad) && pad >= 1 &&
        pad == trunc(pad))){
      stop("pad must be a finite positive integer scalar.")
    }

    # window
    if (window > span){
      stop(sprintf("window %.17g exceeds the input coordinate span %.17g.",
          window, span))
    }
    win_ratio <- window / dt
    win_int <- smom_round(win_ratio)
    if (win_int < 2){
      stop(paste0("window/dt rounds to fewer than two intervals; at least ",
          "three samples are required."))
    }
    win_samp <- win_int + 1
    if (win_samp > n_rows){
      stop("window rounds to more samples than the input rows.")
    }
    # step
    if (step > span){
      stop(sprintf("step %.17g exceeds the input coordinate span %.17g.",
          step, span))
    }
    step_ratio <- step / dt
    step_samp <- smom_round(step_ratio)
    if (step_samp < 1){
      stop("step/dt rounds to zero samples. Increase step.")
    }
    # nfft
    nominal_nfft <- pad * win_samp
    nfft <- nominal_nfft + nominal_nfft %% 2

    # prepare series: global scaling and mean centring
    vscale <- max(abs(y))
    if (vscale == 0){
      yn <- rep(0, n_rows)
    } else {
      yn <- y / vscale
      yn <- yn - mean(yn)
    }
    left_n <- 0
    right_n <- 0
    xp <- x
    yp <- yn
    if (edge_padding == "zero"){
      left_n <- floor(win_int / 2)
      right_n <- ceiling(win_int / 2)
      lx <- numeric(0)
      rx <- numeric(0)
      if (left_n > 0){
        lx <- x[1] - (left_n:1) * dt
      }
      if (right_n > 0){
        rx <- x[n_rows] + (1:right_n) * dt
      }
      xp <- c(lx, x, rx)
      yp <- c(rep(0, left_n), yn, rep(0, right_n))
    }

    last_start <- length(yp) - win_samp + 1
    cand_starts <- seq(1, last_start, by = step_samp)
    cand_ends <- cand_starts + win_samp - 1
    cand_centers <- xp[cand_starts] + (xp[cand_ends] - xp[cand_starts]) / 2
    tol <- max(64 * smom_eps(max(1, max(abs(c(x_start, x_end))))),
        64 * smom_eps(max(1, abs(dt))))
    inside <- cand_centers >= x_start - tol & cand_centers <= x_end + tol
    starts <- cand_starts[inside]
    ends <- cand_ends[inside]
    centers <- cand_centers[inside]
    centers[abs(centers - x_start) <= tol] <- x_start
    centers[abs(centers - x_end) <= tol] <- x_end
    n_win <- length(starts)
    if (n_win < 1){
      stop("The window, step and edge policy retain no window centre.")
    }

    fs <- 2 * nyq
    mean_f <- numeric(n_win)
    bw <- numeric(n_win)
    ordinal <- (0:(win_samp - 1)) - (win_samp - 1) / 2
    denom <- sum(ordinal^2)
    freq <- NULL
    spec_mat <- NULL
    if (verbose == TRUE){
      pb <- utils::txtProgressBar(min = 0, max = n_win, style = 3)
    }
    for (w in seq_len(n_win)){
      seg <- yp[starts[w]:ends[w]]
      seg_scale <- max(abs(seg))
      if (!(is.finite(seg_scale) && seg_scale > 0)){
        stop(sprintf(paste0("Window %d is constant after global mean ",
            "centring and has no resolvable spectral power."),
            w))
      }
      scaled <- seg / seg_scale
      centered <- scaled - mean(scaled)
      slope <- sum(ordinal * centered) / denom
      resid <- centered - slope * ordinal
      peak <- max(abs(resid))
      resolution <- 16 * .Machine$double.eps * max(1, max(abs(scaled)))
      if (!(is.finite(peak) && peak > resolution)){
        stop(sprintf(paste0("Window %d has zero or numerically unresolvable ",
            "residual power after linear detrending."), w))
      }
      resid <- resid / peak
      pg <- acy_periodogram(resid, nfft = nfft, fs = fs)
      pw <- pg$pxx
      if (is.null(spec_mat)) {
        spec_mat <- matrix(NA_real_, length(pw), n_win)
      }
      spec_mat[, w] <- pw
      if (is.null(freq)){
        freq <- pg$f
        freq[length(freq)] <- nyq
      }
      tot <- sum(pw)
      if (!(is.finite(tot) && tot > 0)){
        stop(sprintf("Window %d has zero or nonfinite periodogram power.", w))
      }
      wts <- pw / tot
      fnorm <- freq / nyq
      m_norm <- sum(fnorm * wts)
      v_norm <- sum((fnorm - m_norm)^2 * wts)
      vtol <- 128 * smom_eps(max(1, abs(m_norm)))
      if (v_norm < 0 && v_norm >= -vtol){
        v_norm <- 0
      }
      m_norm <- min(m_norm, 1)
      mean_f[w] <- m_norm * nyq
      bw[w] <- sqrt(v_norm) * nyq
      if (verbose == TRUE){
        utils::setTxtProgressBar(pb, w)
      }
    }
    if (verbose == TRUE){
      close(pb)
    }

    res <- cbind(depth = centers, mean_frequency = mean_f, bandwidth = bw)
    spectra <- list(frequency = freq, depth = centers, power = spec_mat)
    meta <- list(sample_interval = dt,
        physical_nyquist = nyq,
        requested_window = window,
        window_samples = win_samp,
        effective_window = win_int * dt,
        requested_step = step,
        step_samples = step_samp,
        effective_step = step_samp * dt,
        pad = pad,
        nfft = nfft,
        frequency_resolution = freq[2] - freq[1],
        edge_padding = edge_padding,
        left_padding_samples = left_n,
        right_padding_samples = right_n,
        window_count = n_win,
        window_start_index = starts)
    list(result = res, method = "strict", meta = meta, spectra = spectra)
  }

  # @keywords internal
  # @noRd
  smom_zeropad2 <- function(data, win, padding = 1){
    # zeropad2.m (padding 1 = zero, 2 = mirror, 3 = mean), without resorting
    data <- data[order(data[, 1]), , drop = FALSE]
    x <- data[, 1]
    y <- data[, 2]
    y <- y - mean(y)
    data[, 2] <- y
    dt <- mean(diff(x))
    n <- smom_round(win / 2 / dt)
    ny <- length(y)
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
        stop("Mirror padding cannot exceed the input data length.")
      }
      if (n > 0){
        x1y <- y[n:1]
        x2y <- y[ny:(ny - n + 1)]
      } else {
        x1y <- numeric(0)
        x2y <- numeric(0)
      }
    } else if (padding == 3){
      if (n > 0){
        x1y <- rep(mean(y[1:n]), n)
        x2y <- rep(mean(y[(ny - n + 1):ny]), n)
      } else {
        x1y <- numeric(0)
        x2y <- numeric(0)
      }
    } else {
      stop("padding must be 1 (zero), 2 (mirror) or 3 (mean)")
    }
    rbind(cbind(x1x, x1y), data, cbind(x2x, x2y), deparse.level = 0)
  }

  # @keywords internal
  # @noRd
  smom_linspace <- function(d1, d2, n){
    if (n < 2){
      return(d2)
    }
    out <- d1 + (0:(n - 1)) * ((d2 - d1) / (n - 1))
    out[n] <- d2
    out
  }

  # @keywords internal
  # @noRd
  smom_hanning <- function(n){
    # MATLAB hanning(n): symmetric Hann window without the zero end points
    0.5 * (1 - cos(2 * pi * (1:n) / (n + 1)))
  }

  # @keywords internal
  # @noRd
  smom_detrend_pp <- function(signal, windowsize, fs, porder){
    # DetrendSignalPP.m: overlap-add polynomial baseline, returns the trend
    duration <- length(signal)
    sigo <- numeric(duration)
    normo <- numeric(duration)
    if (windowsize < 1){
      return(sigo / normo)
    }
    win <- smom_hanning(windowsize)
    n_pos <- max(0, smom_round(duration - windowsize + 1))
    tt <- (0:(windowsize - 1)) / fs
    tmat <- matrix(0, windowsize, porder)
    for (i in seq_len(porder)){
      tmat[, i] <- tt^(i - 1)
    }
    norms <- sqrt(colSums(tmat^2))
    tn <- tmat
    for (i in seq_len(porder)){
      tn[, i] <- tmat[, i] / norms[i]
    }
    ooff <- 0
    for (i in seq_len(n_pos)){
      seg <- signal[i:(i + windowsize - 1)]
      theta <- qr.solve(tn, seg)
      seg_est <- as.numeric(tn %*% theta)
      idx <- (1:windowsize) + ooff
      sigo[idx] <- sigo[idx] + seg_est * win
      normo[idx] <- normo[idx] + win
      ooff <- ooff + 1
    }
    sigo / normo
  }

  # @keywords internal
  # @noRd
  smom_legacy <- function(data, window, step, pad, sr_mean = 0, want_sr = FALSE,
      verbose = FALSE){
    # spectralmoments.m with padedge = 0 and smoothmodel = 'poly'
    if (sr_mean < 0){
      stop("Mean sedimentation rate must be > 0")
    }
    if (pad < 1){
      stop("padding must be > 0")
    }
    if (step <= 0){
      stop("Step must be larger than 0")
    }
    x <- as.numeric(data[, 1])
    y <- as.numeric(data[, 2])
    dt <- stats::median(diff(x))
    if (step < dt){
      step <- dt
    }
    npts <- length(x)
    if (window >= npts * dt){
      stop("Window must be smaller than data length")
    }
    n_win <- floor(window / dt)
    n_step <- floor(step / dt)
    if (n_win < 2){
      stop("window must contain at least two samples")
    }
    npts_new1 <- npts - n_win + 1
    npts_new <- ceiling(npts_new1 / n_step)
    nfft <- pad * n_win
    uf <- numeric(npts_new)
    bw <- numeric(npts_new)
    if (verbose == TRUE){
      pb <- utils::txtProgressBar(min = 0, max = npts_new, style = 3)
    }
    for (j in seq_len(npts_new)){
      i <- 1 + (j - 1) * n_step
      zz <- y[i:(i + n_win - 1)]
      cc <- acy_detrend_linear(zz)
      pg <- acy_periodogram(cc, nfft = nfft, fs = 1 / dt)
      pxx <- pg$pxx
      f <- pg$f
      if (j == 1) {
        spec_mat <- matrix(NA_real_, length(pxx), npts_new)
      }
      spec_mat[, j] <- pxx
      uf[j] <- sum(f * pxx / sum(pxx))
      bw[j] <- sqrt(sum((f - uf[j])^2 * pxx / sum(pxx)))
      if (verbose == TRUE){
        utils::setTxtProgressBar(pb, j)
      }
    }
    if (verbose == TRUE){
      close(pb)
    }
    depth <- smom_linspace(x[1] + window / 2, x[npts] - window / 2, npts_new)
    res <- cbind(depth = depth, mean_frequency = uf, bandwidth = bw)
    meta <- list(sample_interval = dt,
        requested_window = window,
        window_samples = n_win,
        effective_window = n_win * dt,
        requested_step = step,
        step_samples = n_step,
        effective_step = n_step * dt,
        pad = pad,
        nfft = nfft,
        window_count = npts_new)
    if (want_sr == TRUE){
      windowsize <- smom_round(window / n_step)
      if (windowsize < 2 || windowsize > length(bw)){
        warning(paste0("The bandwidth-trend window round(window/step_samples) = ",
            windowsize, " is not between 2 and the number of ",
            "windows (", length(bw), "); the trend and ",
            "sedimentation rate are NaN (Acycle behaviour)."))
      }
      bw_trend <- smom_detrend_pp(bw, windowsize, n_step, 2)
      sr <- sr_mean * bw[1] / bw_trend
      sr <- sr - mean(sr)
      sr <- sr + sr_mean
      res <- cbind(res, bandwidth_trend = bw_trend, sedrate = sr)
      meta$sr_mean <- sr_mean
      meta$trend_model <- "poly"
      meta$trend_window_samples <- windowsize
    }
    list(result = res, method = "legacy", meta = meta, spectra = list(frequency = f, depth = depth, power = spec_mat))
  }

  # Main code
  data <- as.matrix(data)
  if (ncol(data) < 2){
    stop("data must be a matrix with two columns: [depth/time, value]")
  }
  data <- data[, 1:2, drop = FALSE]
  if (any(!is.finite(data))){
    stop("All coordinates and values of data must be finite.")
  }
  if (nrow(data) < 3){
    stop("Spectral moments require at least three input rows.")
  }
  edge_padding <- tolower(edge_padding)
  if (!(edge_padding %in% c("none", "zero", "mirror", "mean"))){
    stop("edge_padding must be one of 'none', 'zero', 'mirror', 'mean'.")
  }
  x <- data[, 1]
  dx <- diff(x)
  if (is.null(window)){
    window <- 0.25 * abs(x[length(x)] - x[1])
  }
  if (is.null(step)){
    step <- stats::median(dx)
  }

  use_core <- is.null(sr_mean) && (edge_padding %in% c("none", "zero"))
  if (use_core == TRUE){
    out <- smom_core(data = data, window = window, step = step, pad = pad,
        edge_padding = edge_padding, verbose = verbose)
  } else {
    if (any(dx <= 0)){
      stop("The depth/time column must be strictly increasing.")
    }
    dat <- data
    if (edge_padding != "none"){
      pad_type <- 1
      if (edge_padding == "mirror"){
        pad_type <- 2
      } else if (edge_padding == "mean"){
        pad_type <- 3
      }
      dat <- smom_zeropad2(dat, window, pad_type)
    }
    if (is.null(sr_mean)){
      out <- smom_legacy(data = dat, window = window, step = step, pad = pad,
          sr_mean = 0, want_sr = FALSE, verbose = verbose)
    } else {
      out <- smom_legacy(data = dat, window = window, step = step, pad = pad,
          sr_mean = sr_mean, want_sr = TRUE,
          verbose = verbose)
    }
    out$meta$edge_padding <- edge_padding
  }

  # Moments over a frequency band and/or above the white-noise floor (the bias of section 2.3.3 / Appendix B
  # of the paper: white noise adds a flat spectrum whose moments depend on the cutoff frequency, not on the rate)
  if ((!is.null(freq_range) || isTRUE(noise_floor)) && !is.null(out$spectra)) {
    f <- out$spectra$frequency
    P <- out$spectra$power
    keep <- rep(TRUE, length(f))
    if (!is.null(freq_range)) {
      keep <- f >= freq_range[1] & f <= freq_range[2]
    }
    mf <- bw <- rep(NA_real_, ncol(P))
    for (j in seq_len(ncol(P))) {
      p <- P[keep, j]
      fk <- f[keep]
      if (isTRUE(noise_floor)) {
        p <- pmax(p - stats::median(p, na.rm = TRUE), 0)
      }
      tot <- sum(p, na.rm = TRUE)
      if (!is.finite(tot) || tot <= 0) {
        next
      }
      w <- p / tot
      mf[j] <- sum(fk * w, na.rm = TRUE)
      bw[j] <- sqrt(max(0, sum((fk - mf[j])^2 * w, na.rm = TRUE)))
    }
    out$result[, "mean_frequency"] <- mf
    out$result[, "bandwidth"] <- bw
    out$meta$freq_range <- freq_range
    out$meta$noise_floor <- noise_floor
  }
  # Sedimentation rate as in Sinnesael et al. (2018): the smoothed trend of a moment, anchored at one depth
  moment <- match.arg(moment)
  trend <- match.arg(trend)
  if (!is.null(anchor)) {
    res <- out$result
    z <- res[, "depth"]
    m <- res[, if (moment == "bandwidth") { "bandwidth" } else { "mean_frequency" }]
    if (is.null(trend_window)) {
      trend_window <- 3 * window
    }
    dz <- stats::median(diff(z))
    wsamp <- max(3, smom_round(trend_window / dz))
    tr <- if (trend == "poly" && wsamp < length(m)) {
      smom_detrend_pp(m, wsamp, 1 / dz, trend_order + 1)
    } else {
      stats::lowess(z, m, f = min(1, trend_window / diff(range(z))))$y
    }
    an <- as.list(anchor)
    z0 <- as.numeric(an[[1]])
    z0 <- min(max(z0, min(z)), max(z))
    r0 <- if (!is.null(an$rate)) { as.numeric(an$rate) } else if (length(an) >= 3) { 1 / (as.numeric(an[[2]]) * as.numeric(an[[3]])) } else { as.numeric(an[[2]]) }
    t0 <- stats::approx(z, tr, z0, rule = 2)$y
    sr <- r0 * t0 / tr                       # every spatial frequency scales with 1 / rate (equation 12 of the paper)
    res <- cbind(res, moment_trend = tr, sedrate_anchored = sr)
    out$result <- res
    out$anchor <- list(depth = z0, rate = r0, moment = moment, trend = trend, trend_window = trend_window, trend_order = trend_order)
  }
  if (genplot == TRUE){
    plot_spectral_moments(out)
  }
  out
}
