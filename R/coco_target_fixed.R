#' @title Fixed-weight astronomical target spectrum (Fixed-target COCO)
#'
#' @description Builds the fixed-relative-weight astronomical target used by
#' Acycle's Fixed-target COCO: nine (or other) orbital sine terms with the
#' historical amplitudes 1.0 (eccentricity), 0.8 (obliquity) and 0.6
#' (precession) (cocoFixedTargetWeights.m), sampled on the data grid mapped
#' to time at sedimentation rate \code{sr}, linearly detrended, and turned
#' into a rectangular-window periodogram of length max(pad, n) that is
#' linearly interpolated at the requested temporal frequencies (0 outside
#' the target band). Unlike the adaptive target, the amplitudes never depend
#' on the data spectrum. Translated from Acycle v3.1 freq2targetFixed.m and
#' cocoFixedTargetWeights.m (Mingsong Li).
#'
#' @param data Two-column matrix (depth, value) defining the sampling grid;
#' only the depth column is used. NULL uses 5 x max(orbit9) unit-spaced
#' samples. \code{Default=NULL}
#' @param pad Periodogram length. \code{Default=5000}
#' @param orbit9 Orbital periods (kyr).
#' \code{Default=c(405.6912, 130.6979, 123.8532, 98.8517, 94.8856, 40.9897, 23.6820, 22.3758, 18.9519)}
#' @param freq Temporal frequencies (cycles/kyr) at which the target power
#' is returned. \code{Default=NULL}
#' @param sr Sedimentation rate (cm/kyr). \code{Default=NULL}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' @references Li, M., Kump, L.R., Hinnov, L.A., Mann, M.E., 2018. Tracking
#' variable sedimentation rates and astronomical forcing in Phanerozoic
#' paleoclimate proxy series with evolutionary correlation coefficients and
#' hypothesis testing. Earth and Planetary Science Letters 501, 165-179.
#'
#' @examples
#' \donttest{
#' depth <- seq(0, 20, by = 0.05)
#' f <- seq(0, 0.06, by = 0.0005)
#' la <- coco_target_fixed(cbind(depth, 0), pad = 5000, freq = f, sr = 5)
#' plot(f, la, type = "l", xlab = "Frequency (cycles/kyr)", ylab = "Power")
#' }
#'
#' @return Numeric vector of target power at \code{freq}.
#' @export
#' @importFrom stats sd approx median
coco_target_fixed <- function(data = NULL, pad = 5000,
                              orbit9 = c(405.6912, 130.6979, 123.8532,
                                         98.8517, 94.8856, 40.9897, 23.6820,
                                         22.3758, 18.9519),
                              freq = NULL, sr = NULL){
  # ---- internal helper functions (local to this function) ----

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

  # --- MATLAB periodogram (rectangular window, one-sided PSD) -------------------
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

  # ---- main code ----
  freq <- as.numeric(freq)
  la <- rep(0, length(freq))
  if (length(freq) < 1){
    return(la)
  }
  orbit9 <- as.numeric(orbit9)
  w <- coco_fixed_target_weights(orbit9)
  w[!is.finite(orbit9) | orbit9 <= 0] <- 0
  ok <- is.finite(orbit9) & orbit9 > 0 & is.finite(w) & w > 0
  orbit9 <- orbit9[ok]
  w <- w[ok]
  if (length(orbit9) < 1){
    return(la)
  }
  if (is.null(pad) || !is.finite(pad) || pad <= 0){
    pad <- 5000
  }
  fs <- 1
  nsmp <- max(2, ceiling(5 * max(orbit9)))
  if (!is.null(data) && ncol(as.matrix(data)) >= 1 && !is.null(sr) &&
      is.finite(sr) && sr > 0){
    d <- as.matrix(data)[, 1]
    d <- d[is.finite(d)]
    if (length(d) >= 2){
      d <- sort(d)
      dz <- diff(d)
      dz <- dz[is.finite(dz) & dz > 0]
      if (length(dz) > 0){
        nsmp <- length(d)
        fs <- sr / (100 * stats::median(dz))
      }
    }
  }
  if (!is.finite(fs) || fs <= 0){
    fs <- 1
  }
  nsmp <- max(2, nsmp)
  nfft <- max(ceiling(pad), nsmp)
  tt <- (0:(nsmp - 1)) / fs
  y0 <- rep(0, nsmp)
  for (i in seq_along(orbit9)){
    y0 <- y0 + w[i] * sin(2 * pi / orbit9[i] * tt)
  }
  if (stats::sd(y0) == 0){
    return(la)
  }
  y0 <- detrend_linear(y0)
  sp <- ml_periodogram(y0, nfft, fs)
  tf <- sp$f
  tp <- sp$p
  okf <- is.finite(tf) & is.finite(tp)
  tf <- tf[okf]
  tp <- tp[okf]
  if (length(tf) < 1){
    return(la)
  }
  keep <- !duplicated(tf)
  o <- order(tf[keep])
  tf <- tf[keep][o]
  tp <- tp[keep][o]
  la <- stats::approx(tf, tp, xout = freq, method = "linear", rule = 1)$y
  la[!is.finite(la)] <- 0
  la
}
