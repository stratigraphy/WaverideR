#' @title Target periodogram from a list of orbital periods (period2spectrum.m)
#' @description Sum of sinusoids at the orbital periods, sampled from t1 to t2
#'   (kyr), and its zero-padded periodogram between f1 and f2 (cycles/kyr).
#' @param orbit Astronomical periods of the target sinusoids, in kyr (numeric vector, > 0; sorted longest first
#' internally), e.g. \code{c(405, 125, 95, 41, 23, 22, 19)}. \code{NULL} uses the nine periods at \code{age_ma}.
#' \code{Default=NULL}.
#' @param age_ma Age in Ma (single numeric value, 0 to 4429) at which the nine periods are taken from
#' \code{\link{orbit9_at_age}} (La2004 g and s frequencies with the AstroGeo22 precession constant); only used
#' when \code{orbit} is \code{NULL}, and one of the two must be given. \code{Default=NULL}.
#' @param t1,t2 Start and end of the time axis on which the synthetic sum of sinusoids is sampled, in kyr
#' (numeric, \code{t2 > t1}); the length \code{t2 - t1} sets the width of the target peaks.
#' \code{Default=300000} and \code{Default=308000}.
#' @param step Time step of the synthetic series, in kyr (numeric, > 0). \code{Default=1}.
#' @param f1 Lowest frequency of the target that is returned, in cycles/kyr (numeric, >= 0). \code{Default=0}.
#' @param f2 Highest frequency of the target that is returned, in cycles/kyr (numeric, > \code{f1}); \code{NULL} uses
#' 1.2 times the highest orbital frequency, as in the Acycle GUI. \code{Default=NULL}.
#' @param pad Length (number of points) to which the synthetic series is zero-padded for the periodogram (integer;
#' at least the series length is used); a larger value gives a finer frequency grid. \code{Default=10000}.
#' @param weights Amplitudes of the sinusoids, one per period in \code{orbit} after sorting longest first (numeric
#' vector). \code{NULL} uses 1.0 for eccentricity, 0.8 for obliquity and 0.6 for precession (for nine periods via
#' \code{\link{coco_fixed_target_weights}}; otherwise the three longest get 1.0, the fourth 0.8 and the rest 0.6, as
#' in Acycle period2spectrum.m). \code{Default=NULL}.
#' @details The returned target carries attribute "orbit" so coco_classic() uses
#'   the same periods to count contributing parameters.
#' @return two-column matrix: frequency (cycles/kyr), power
#' @export
coco_target_from_periods <- function(orbit = NULL, age_ma = NULL, t1 = 300000, t2 = 308000,
    step = 1, f1 = 0, f2 = NULL, pad = 10000,
    weights = NULL){
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
  if (is.null(orbit)){
    if (is.null(age_ma)){
      stop("Supply 'orbit' (kyr) or 'age_ma'")
    }
    orbit <- orbit9_at_age(age_ma)
  }
  orbit <- sort(orbit, decreasing = TRUE)
  if (is.null(f2)){
    f2 <- 1.2 * max(1 / orbit)
  }
  if (is.null(weights)){
    weights <- rep(0.6, length(orbit))
    weights[seq_len(min(3, length(orbit)))] <- 1
    if (length(orbit) >= 4){
      weights[4] <- 0.8
    }
    if (length(orbit) == 9){
      weights <- coco_fixed_target_weights(orbit)
    }
  }
  x <- seq(t1, t2, by = step)
  sines <- matrix(0, length(x), length(orbit))
  for (i in seq_along(orbit)){
    sines[, i] <- weights[i] * sin(2 * pi / orbit[i] * x)
  }
  y <- rowSums(sines)
  sp <- ml_periodogram(y, max(pad, length(y)), 1 / step)
  tg <- select_interval(cbind(sp$f, sp$p), f1, f2)
  attr(tg, "orbit") <- orbit
  tg
}
