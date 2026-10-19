#' @title Target periodogram from a list of orbital periods (period2spectrum.m)
#' @description Sum of sinusoids at the orbital periods, sampled from t1 to t2
#'   (kyr), and its zero-padded periodogram between f1 and f2 (cycles/kyr).
#' @param orbit periods in kyr, longest first (e.g. c(405,125,95,41,23,22,19))
#'   \code{Default=NULL}
#' @param age_ma alternative to `orbit`: use the nine periods at this age
#'   (calculate_orbit9(); La2004 g/s + AstroGeo22 precession constant)
#'   \code{Default=NULL}
#' @param t1,t2,step time axis of the synthetic series (kyr)
#' @param f1,f2 frequency interval kept (cycles/kyr); f2 defaults to 1.2 x the
#'   highest orbital frequency (Acycle GUI)
#' @param pad zero padding of the periodogram \code{Default=10000}
#' @param weights sinusoid amplitudes; default 1.0 eccentricity, 0.8 obliquity,
#'   0.6 precession (as in period2spectrum.m)
#' @details The returned target carries attribute "orbit" so coco_classic() uses
#'   the same periods to count contributing parameters.
#' @return two-column matrix: frequency (cycles/kyr), power
#' @export
coco_target_from_periods <- function(orbit = NULL, age_ma = NULL, t1 = 300000, t2 = 308000,
                                     step = 1, f1 = 0, f2 = NULL, pad = 10000,
                                     weights = NULL){
  # ---- internal helper functions (local to this function) ----

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

  select_interval <- function(x, lo, hi){
    x[x[, 1] >= lo & x[, 1] <= hi, , drop = FALSE]
  }

  # ---- main code ----
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
