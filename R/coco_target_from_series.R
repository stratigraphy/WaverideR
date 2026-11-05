#' @title Target periodogram from an astronomical time series (e.g. La2004 ETP)
#' @description Periodogram of an astronomical series, or of an ETP series built
#'   from eccentricity, obliquity and precession as in gentarget.m.
#' @param time_kyr Time axis of the astronomical series, in kyr (numeric vector, evenly spaced; the series is
#' sorted by it and its first step is used as the sampling interval). Required, no default.
#' @param value Astronomical series to use as target (numeric vector, same length as \code{time_kyr}), e.g. an ETP
#' series; \code{NULL} builds an ETP series from \code{ecc}, \code{obl} and \code{prec}. \code{Default=NULL}.
#' @param ecc Eccentricity series (numeric vector, same length and order as \code{time_kyr}), used with \code{obl}
#' and \code{prec} instead of \code{value}: the ETP series is \code{ne * z(ecc) + nt * z(obl) - np * z(prec)},
#' with \code{z} the standardisation to zero mean and unit standard deviation (as in Acycle gentarget.m).
#' \code{Default=NULL}.
#' @param obl Obliquity series (numeric vector, same length and order as \code{time_kyr}), used to build the ETP
#' series when \code{value} is \code{NULL}. \code{Default=NULL}.
#' @param prec Precession (index) series (numeric vector, same length and order as \code{time_kyr}), used with a
#' negative sign to build the ETP series when \code{value} is \code{NULL}. \code{Default=NULL}.
#' @param ne Weight of the standardised eccentricity in the ETP series (numeric). \code{Default=1}.
#' @param nt Weight of the standardised obliquity in the ETP series (numeric). \code{Default=0.6}.
#' @param np Weight of the standardised precession in the ETP series (numeric). \code{Default=0.5}.
#' @param f1 Lowest frequency of the target that is returned, in cycles/kyr (numeric, >= 0). \code{Default=0}.
#' @param f2 Highest frequency of the target that is returned, in cycles/kyr (numeric, > \code{f1}).
#' \code{Default=0.06}.
#' @param pad Length (number of points) to which the series is zero-padded for the periodogram (integer; at least
#' the series length is used); a larger value gives a finer frequency grid. \code{Default=10000}.
#' @return two-column matrix: frequency (cycles/kyr), power
#' @export
coco_target_from_series <- function(time_kyr, value = NULL, ecc = NULL, obl = NULL,
    prec = NULL, ne = 1, nt = 0.6, np = 0.5,
    f1 = 0, f2 = 0.06, pad = 10000){
  # Internal helper functions (local to this function)

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

  select_interval <- function(x, lo, hi){
    x[x[, 1] >= lo & x[, 1] <= hi, , drop = FALSE]
  }

  # Main code
  if (is.null(value)){
    if (is.null(ecc) || is.null(obl) || is.null(prec)){
      stop("Supply either 'value' or all of 'ecc', 'obl', 'prec'")
    }
    z <- function(v){
      (v - mean(v)) / stats::sd(v)
    }
    value <- ne * z(ecc) + nt * z(obl) - np * z(prec)
  }
  o <- order(time_kyr)
  time_kyr <- time_kyr[o]
  value <- value[o]
  dtt <- time_kyr[2] - time_kyr[1]
  sp <- ml_periodogram(value, max(pad, length(value)), 1 / dtt)
  select_interval(cbind(sp$f, sp$p), f1, f2)
}
