# Shared internal helpers of the MESA functions

# Burg recursion: coefficients of the prediction error filter A(z) = 1 + a1 z^-1 + ...
# for every order up to m_max and the prediction error power of every order
.mesa_burg <- function(x, m_max){
  n <- length(x)
  f <- x
  b <- x
  a <- 1
  e <- numeric(m_max + 1)
  e[1] <- sum(x^2) / n
  coefs <- vector("list", m_max)
  for (m in seq_len(m_max)){
    len <- length(f)
    ff <- f[2:len]
    bb <- b[1:(len - 1)]
    den <- sum(ff^2) + sum(bb^2)
    if (den <= 0){
      k <- 0
    } else {
      k <- -2 * sum(ff * bb) / den
    }
    a <- c(a, 0) + k * c(0, rev(a))
    f <- ff + k * bb
    b <- bb + k * ff
    e[m + 1] <- e[m] * (1 - k^2)
    coefs[[m]] <- a
  }
  list(coefs = coefs, e = e)
}

# Fit the AR model: fixed order or order chosen by an information criterion
.mesa_fit <- function(x, order = NULL, order_max = NULL, criterion = "fpe"){
  n <- length(x)
  if (is.null(order_max)){
    order_max <- min(floor(n / 3), 300)
  }
  order_max <- min(order_max, n - 3)
  if (!is.null(order)){
    if (order < 1 || order > n - 3){
      stop("order must be between 1 and N - 3")
    }
    order_max <- max(order_max, order)
  }
  if (order_max < 1){
    stop("the series is too short for an AR model")
  }
  rec <- .mesa_burg(x, order_max)
  e <- pmax(rec$e, .Machine$double.xmin)
  m <- 0:order_max
  fpe <- e * (n + m + 1) / (n - m - 1)
  aic <- n * log(e) + 2 * m
  mdl <- n * log(e) + m * log(n)
  cat_c <- numeric(order_max + 1)
  inv_e <- (n - m) / (n * e)
  for (j in m){
    if (j == 0){
      cat_c[j + 1] <- -(1 + 1 / n) / e[1]
    } else {
      cat_c[j + 1] <- mean((n - seq_len(j)) / (n * e[seq_len(j) + 1])) - inv_e[j + 1]
    }
  }
  criteria <- data.frame(order = m, fpe = fpe, aic = aic, mdl = mdl, cat = cat_c)
  chosen <- "fixed"
  if (is.null(order)){
    order <- which.min(criteria[[criterion]][-1])
    chosen <- criterion
  }
  list(order = order, a = rec$coefs[[order]], noise_variance = e[order + 1],
      criteria = criteria, criterion = chosen)
}

# One-sided power spectral density of the AR model at the given frequencies
.mesa_psd <- function(a, noise_variance, dt, freq){
  lag <- 0:(length(a) - 1)
  z <- exp(-2i * pi * outer(freq * dt, lag))
  resp <- Mod(z %*% a)^2
  as.numeric(2 * dt * noise_variance / resp)
}

# Local maxima of a spectrum with prominence (power over the higher neighbouring minimum)
.mesa_peaks <- function(freq, power, n_peaks = 10){
  nf <- length(power)
  empty <- data.frame(frequency = numeric(0), period = numeric(0),
      power = numeric(0), prominence = numeric(0))
  if (nf < 3){
    return(empty)
  }
  mid <- 2:(nf - 1)
  idx <- mid[power[mid] > power[mid - 1] & power[mid] >= power[mid + 1]]
  if (!length(idx)){
    return(empty)
  }
  prom <- numeric(length(idx))
  for (k in seq_along(idx)){
    i <- idx[k]
    l <- i
    while (l > 1 && power[l - 1] < power[l]){
      l <- l - 1
    }
    r <- i
    while (r < nf && power[r + 1] < power[r]){
      r <- r + 1
    }
    prom[k] <- power[i] / max(power[l], power[r])
  }
  pk <- data.frame(frequency = freq[idx], period = 1 / freq[idx],
      power = power[idx], prominence = prom)
  pk <- pk[order(-pk$power), ]
  pk <- utils::head(pk, n_peaks)
  rownames(pk) <- NULL
  pk
}

# Samples of the evenly sampled record on the grid between lo and hi, with the part outside
# the record filled according to pad ("reflect", "mean", "zero" or "none")
.mesa_segment <- function(x, t0, dt, lo, hi, pad = "reflect"){
  n <- length(x)
  i_lo <- ceiling((lo - t0) / dt - 1e-9) + 1
  i_hi <- floor((hi - t0) / dt + 1e-9) + 1
  idx <- i_lo:i_hi
  inside <- idx >= 1 & idx <= n
  if (pad == "none"){
    return(x[idx[inside]])
  }
  out <- numeric(length(idx))
  if (pad == "reflect" && n > 1){
    per <- 2 * (n - 1)
    j <- (idx - 1) %% per
    j <- ifelse(j > n - 1, per - j, j)
    out <- x[j + 1]
  } else {
    out[inside] <- x[idx[inside]]
    if (pad == "mean" && any(inside)){
      out[!inside] <- mean(x[idx[inside]])
    }
  }
  out
}

# Read, check, regularise and clean a record: returns t, x, dt (x demeaned or detrended)
.mesa_prepare <- function(data, regularise = TRUE, demean = TRUE, detrend = FALSE){
  data <- as.data.frame(data)
  if (ncol(data) < 2){
    stop("data must have two columns: time (or depth) and value")
  }
  data <- data[stats::complete.cases(data[, 1:2]), 1:2]
  data <- data[order(data[, 1]), ]
  t <- data[, 1]
  x <- data[, 2]
  if (length(t) < 8){
    stop("at least 8 samples are needed")
  }
  dts <- diff(t)
  if (any(dts <= 0)){
    stop("the first column (time or depth) must be strictly increasing")
  }
  dt <- stats::median(dts)
  if (max(abs(dts - dt)) > 1e-3 * dt){
    if (!regularise){
      stop("the series is not evenly sampled; use regularise = TRUE or interpolate it first")
    }
    message("MESA: series is unevenly sampled, interpolated to dt = ", signif(dt, 6))
    tn <- seq(t[1], t[length(t)], by = dt)
    x <- stats::approx(t, x, xout = tn)$y
    t <- tn
  }
  n <- length(x)
  if (detrend){
    x <- x - stats::lm.fit(cbind(1, seq_len(n)), x)$fitted.values
  } else if (demean){
    x <- x - mean(x)
  }
  list(t = t, x = x, dt = dt)
}

#' @title Maximum entropy spectral analysis (MESA, Burg)
#'
#' @description Estimates the power spectrum of a time series with the maximum entropy
#' method of Burg (1967, 1975). An autoregressive (AR) model is fitted to the series with the
#' Burg recursion, which minimises the forward and backward prediction error, and the power
#' spectrum is the spectrum of that AR model. Compared with the periodogram or the multi-taper
#' method (\code{astrochron::mtm}), MESA gives very sharp spectral lines for short records and
#' needs no window, which is useful for the few cycles typical of cyclostratigraphic data.
#' The price is that the number of peaks and their width depend on the model order, so the order
#' is chosen with an information criterion (default the final prediction error, FPE) or set by hand.
#' Spurious peak splitting at high orders is a known property of the method; compare several
#' orders or criteria (see \code{criteria} in the output) before interpreting a peak.
#'
#' @param data Record to analyse: a two-column matrix or data frame with depth or time in column 1 (strictly
#' increasing after sorting) and the proxy value in column 2; rows with \code{NA} are removed and at least 8 samples
#' are needed. It should be evenly sampled, see \code{regularise}. Required (\code{Default=NULL}).
#' @param order Order of the AR model, i.e. the number of AR coefficients (integer, between 1 and N - 3 for N
#' samples). Higher orders give sharper and more peaks. \code{NULL} selects the order with \code{criterion}.
#' \code{Default=NULL}.
#' @param order_max Largest AR order considered when the order is selected (integer, >= 1; it is capped at N - 3).
#' \code{NULL} uses \code{min(floor(N / 3), 300)}. \code{Default=NULL}.
#' @param criterion Order-selection criterion (character), one of \code{"fpe"} (final prediction error, Akaike
#' 1969), \code{"aic"} (Akaike information criterion), \code{"mdl"} (minimum description length of Rissanen, the
#' most conservative) or \code{"cat"} (criterion autoregressive transfer of Parzen); ignored when \code{order} is
#' given. \code{Default=c("fpe", "aic", "mdl", "cat")}, i.e. \code{"fpe"}.
#' @param demean Remove the mean of the series before the analysis (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param detrend Remove a linear trend before the analysis instead of only the mean (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param regularise Linearly interpolate an unevenly sampled series to its median sampling interval, with a message
#' (logical, TRUE or FALSE); with \code{FALSE} an unevenly sampled series gives an error. \code{Default=TRUE}.
#' @param nfreq Number of frequencies, evenly spaced between \code{fmin} and \code{fmax}, at which the spectrum is
#' evaluated (integer, > 1). \code{Default=2048}.
#' @param fmin Lowest frequency of the spectrum, in cycles per unit of the first column of \code{data} (numeric,
#' >= 0 and < \code{fmax}). \code{Default=0}.
#' @param fmax Highest frequency of the spectrum, in cycles per unit of the first column of \code{data} (numeric,
#' > \code{fmin}; values above the Nyquist frequency are reduced to it). \code{NULL} uses the Nyquist frequency.
#' \code{Default=NULL}.
#' @param n_peaks Maximum number of spectral peaks (local maxima) returned, the strongest ones (integer, > 0).
#' \code{Default=10}.
#' @param window Part of the record to analyse, \code{c(from, to)} in the unit of the first column of \code{data}
#' (numeric vector of length 2, from < to). The window may extend beyond the record; the missing part is filled
#' according to \code{pad}, so that a spectrum of a window centred on the very top or bottom can be made.
#' \code{NULL} uses the entire record. \code{Default=NULL}.
#' @param pad How samples of \code{window} that lie outside the record are filled (character): \code{"reflect"}
#' (mirror the record at its edge), \code{"mean"} (the mean of the part of the window inside the record),
#' \code{"zero"} (zero, i.e. the mean of the record after demeaning) or \code{"none"} (the window is cut at the edge
#' of the record). Padding only adds samples outside the record; it never changes the record. Only used with
#' \code{window}. \code{Default=c("reflect", "mean", "zero", "none")}, i.e. \code{"reflect"}.
#'
#' @details The power spectral density is one-sided and scaled so that its integral over frequency
#' equals the variance of the (demeaned) series. The AR coefficients are returned in the convention
#' \code{x[t] = sum(ar[j] * x[t - j]) + e[t]}, the same as \code{stats::ar.burg}.
#'
#' @return A list of class \code{mesa_spectrum} with elements
#' \describe{
#' \item{spectrum}{data frame with \code{frequency}, \code{period} and \code{power}.}
#' \item{peaks}{data frame with the strongest local maxima: \code{frequency}, \code{period},
#' \code{power}, and \code{prominence} (peak power divided by the lower of the two neighbouring minima).}
#' \item{order}{the AR order used.}
#' \item{ar}{the AR coefficients.}
#' \item{noise_variance}{the prediction error variance of the final model.}
#' \item{criteria}{data frame with the criteria (fpe, aic, mdl, cat) for every order 0 to the maximum order.}
#' \item{dt}{sampling interval.}
#' \item{n}{number of samples.}
#' }
#'
#' @references
#' Burg, J.P., 1967. Maximum entropy spectral analysis. 37th Annual International Meeting,
#' Society of Exploration Geophysicists, Oklahoma City.
#'
#' Burg, J.P., 1975. Maximum entropy spectral analysis. PhD thesis, Stanford University.
#'
#' Akaike, H., 1969. Fitting autoregressive models for prediction. Annals of the Institute of
#' Statistical Mathematics 21, 243-247.
#'
#' Parzen, E., 1974. Some recent advances in time series modeling. IEEE Transactions on Automatic
#' Control 19, 723-730.
#'
#' Rissanen, J., 1978. Modeling by shortest data description. Automatica 14, 465-471.
#'
#' Ulrych, T.J., Bishop, T.N., 1975. Maximum entropy spectral analysis and autoregressive
#' decomposition. Reviews of Geophysics 13, 183-200.
#'
#' @examples
#' # two cycles (100 and 41 units) plus noise in a short record
#' set.seed(1)
#' t <- seq(0, 600, by = 2)
#' y <- sin(2 * pi * t / 100) + 0.6 * sin(2 * pi * t / 41) + rnorm(length(t), sd = 0.5)
#' m <- mesa_spectrum(data.frame(t, y))
#' m$order
#' m$peaks
#' plot_mesa_spectrum(m)
#'
#' # a fixed, low order gives a smoother spectrum
#' m8 <- mesa_spectrum(data.frame(t, y), order = 8)
#' plot_mesa_spectrum(m8, peaks = FALSE)
#'
#' @importFrom stats approx median complete.cases lm.fit
#' @export
mesa_spectrum <- function(data = NULL,
    order = NULL,
    order_max = NULL,
    criterion = c("fpe", "aic", "mdl", "cat"),
    demean = TRUE,
    detrend = FALSE,
    regularise = TRUE,
    nfreq = 2048,
    fmin = 0,
    fmax = NULL,
    n_peaks = 10,
    window = NULL,
    pad = c("reflect", "mean", "zero", "none")) {
  # Main code
  criterion <- match.arg(criterion)
  pad <- match.arg(pad)
  prep <- .mesa_prepare(data, regularise = regularise, demean = demean, detrend = detrend)
  t <- prep$t
  x <- prep$x
  dt <- prep$dt
  if (!is.null(window)){
    if (length(window) != 2 || window[2] <= window[1]){
      stop("window must be c(from, to) with from < to")
    }
    x <- .mesa_segment(x, t[1], dt, window[1], window[2], pad)
    if (length(x) < 8){
      stop("the window holds fewer than 8 samples")
    }
    if (pad == "mean" || demean){
      x <- x - mean(x)
    }
  }
  n <- length(x)
  fit <- .mesa_fit(x, order = order, order_max = order_max, criterion = criterion)
  nyq <- 1 / (2 * dt)
  if (is.null(fmax)){
    fmax <- nyq
  }
  fmax <- min(fmax, nyq)
  if (fmin < 0 || fmin >= fmax){
    stop("fmin must be >= 0 and smaller than fmax")
  }
  freq <- seq(fmin, fmax, length.out = nfreq)
  power <- .mesa_psd(fit$a, fit$noise_variance, dt, freq)
  spec <- data.frame(frequency = freq,
      period = ifelse(freq > 0, 1 / freq, Inf),
      power = power)
  pk <- .mesa_peaks(freq, power, n_peaks)
  out <- list(spectrum = spec, peaks = pk, order = fit$order, ar = -fit$a[-1],
      noise_variance = fit$noise_variance, criteria = fit$criteria,
      criterion = fit$criterion, dt = dt, n = n, window = window)
  class(out) <- "mesa_spectrum"
  out
}

#' @export
print.mesa_spectrum <- function(x, ...){
  cat("Maximum entropy (Burg) spectrum: N =", x$n, ", dt =", signif(x$dt, 5),
      ", AR order =", x$order, "(", x$criterion, ")\n")
  if (nrow(x$peaks)){
    cat("Strongest peaks:\n")
    print(format(x$peaks, digits = 4), row.names = FALSE)
  }
  invisible(x)
}
