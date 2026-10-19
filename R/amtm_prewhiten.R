# amtm_spectrum.R
#
# Adaptive multitaper (aMTM) PSD background estimation and periodic signal
# detection (gamma-test + harmonic F-test), R translation of the Python package
# amtm_signaldetect v1.0.6 (https://github.com/hsalinasGIT/amtm_signaldetect)
# and of the parts of Prieto's 'multitaper' Python package (v1.2.0,
# https://github.com/gaprieto/multitaper) that amtm_signaldetect calls
# (MTSpec with iadapt = 0: dpss normalisation, eigenspec, adaptspec, rspec, ftest).
#
# ---------------------------------------------------------------------------
# Original licences (both MIT):
#
# amtm_signaldetect -- MIT License, Copyright (c) 2024 Hector A Salinas
# multitaper        -- MIT License, Copyright (c) 2022 German A. Prieto
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ---------------------------------------------------------------------------
# The AR(1) background model, the AIC / MERIT model-selection criteria and the
# AR(1) starting value follow the IDL package SPD_MTM (c) Simone Di Matteo 2020,
# Apache License 2.0 (SPEDAS: general/science/spd_mtm/models/ar1/*.pro,
# spd_mtm_modelgof.pro, spd_mtm_fitmodel.pro). They are not part of the Python
# package. The prewhitening option follows the "Python PreWhiten Filter Example"
# distributed with amtm_signaldetect.
#
# Uses from helpers_matlab_equivalents.R: acy_dpss(), acy_check_two_column().
# ---------------------------------------------------------------------------



#' @title Autoregressive prewhitening of a series for aMTM analysis
#'
#' @description The \code{\link{amtm_prewhiten}} function prewhitens an evenly
#' spaced series with an AR(p) filter whose coefficients are obtained from the
#' autocovariance sequence of a low-variance aMTM spectrum (NW = 4, K = 7, zero
#' padding to \code{2^(ceiling(log2(N)) + nfft_add)}) through the Yule-Walker
#' equations. The one-step predictions minus the observations form the prewhitened
#' series; the first \code{2*order} values are dropped and the mean is removed.
#' The spectrum of the prewhitened series can be "recoloured" by dividing it by
#' \eqn{|1 - \sum_j \phi_j e^{-i 2\pi f j \Delta t}|^2}.
#' Translated from the "Python PreWhiten Filter Example" (prewhiten(),
#' correct_spec_for_pw()) distributed with amtm_signaldetect (MIT licence), itself
#' based on an R prewhiten() routine (Wesley Burr / David Riegert).
#'
#' @param data Input data set: 2-column matrix [depth/time, value], evenly spaced
#' @param order Order of the AR filter \code{Default=2}
#' @param nw Time-bandwidth product of the spectrum used for the ACVS \code{Default=4}
#' @param k Number of tapers of that spectrum \code{Default=7}
#' @param nfft_add Zero-padding exponent increment \code{Default=2}
#'
#' @author Michiel Arts (R translation); original Python code: Hector A. Salinas et al.
#' @references Percival, D. B., & Walden, A. T. (2020). Spectral Analysis for
#' Univariate Time Series. Cambridge University Press. \doi{10.1017/9781139235723}
#' @examples
#' \donttest{
#' set.seed(3)
#' x <- as.numeric(stats::filter(rnorm(400), 0.9, method = "recursive"))
#' pw <- amtm_prewhiten(cbind(seq_along(x), x), order = 1)
#' pw$ar_coef
#' }
#' @return A list with \code{data} (2-column matrix of the prewhitened series, with
#' the depth/time values of the retained samples), \code{ar_coef} (AR coefficients
#' phi_1..phi_p) and \code{acvs} (autocovariance sequence used).
#' @export
#' @importFrom stats fft
amtm_prewhiten <- function(data = NULL, order = 2, nw = 4, k = 7, nfft_add = 2){
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

  acy_check_two_column <- function(data, name = "data"){
    data <- as.matrix(data)
    if (ncol(data) < 2){
      stop(paste0(name, " must have at least two columns: [depth/time, value]"))
    }
    data <- data[, 1:2, drop = FALSE]
    keep <- is.finite(data[, 1]) & is.finite(data[, 2])
    data <- data[keep, , drop = FALSE]
    data <- data[order(data[, 1]), , drop = FALSE]
    data
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

  # @title Prieto-normalised DPSS tapers
  # @description Unit-energy DPSS (acy_dpss, equal to scipy.signal.windows.dpss to
  # ~1e-13) with Prieto's "positive standard" sign convention (taper positive at
  # sample floor(N/2) (0-based)) and concentration eigenvalues.
  # @keywords internal
  # @noRd
  amtm_dpss <- function(n, nw, k){
    tp <- acy_dpss(n, nw, k)
    vn <- tp$E
    for (j in seq_len(k)){
      vn[, j] <- vn[, j] / sqrt(sum(vn[, j]^2))
    }
    if (n %% 2 == 1){
      lh <- (n + 1) / 2
    } else {
      lh <- n / 2
    }
    for (j in seq_len(k)){
      if (vn[lh + 1, j] < 0){
        vn[, j] <- -vn[, j]
      }
    }
    list(vn = vn, lamb = tp$V)
  }

  # @title Adaptive multitaper spectrum (Prieto MTSpec, iadapt = 0)
  # @description Returns the full (two-sided, nfft) adaptive spectrum scaled so that
  # sum(spec)*df = var(x), the adaptive weights, eigencoefficients and tapers.
  # @keywords internal
  # @noRd
  amtm_mtspec <- function(x, nw = 3, k = NULL, dt = 1, nfft = NULL){
    npts <- length(x)
    if (is.null(k)){
      k <- as.integer(round(2 * nw - 1))
    }
    if (is.null(nfft)){
      nfft <- 2 * npts + 1
    }
    nfft <- as.integer(nfft)
    if (nfft < npts){
      nfft <- as.integer(2 * npts + 1)
    }
    x <- x - mean(x)
    xvar <- mean((x - mean(x))^2)
    if (nfft %% 2 == 0){
      nf <- nfft / 2 + 1
    } else {
      nf <- (nfft + 1) / 2
    }
    val <- 1 / (nfft * dt)
    freq <- (0:(nf - 1)) * val
    df <- 2 * val - 1 * val
    tp <- amtm_dpss(npts, nw, k)
    vn <- tp$vn
    lamb <- tp$lamb
    yk <- matrix(0 + 0i, nfft, k)
    pad <- numeric(nfft - npts)
    for (j in seq_len(k)){
      yk[, j] <- stats::fft(c(vn[, j] * x, pad))
    }
    sk <- Mod(yk)^2
    ad <- amtm_adaptspec(sk, lamb)
    sscal <- xvar / (sum(ad$spec) * df)
    list(freq = freq, spec = sscal * ad$spec, wt = ad$wt, yk = yk, vn = vn,
         lamb = lamb, nf = nf, nfft = nfft, df = df, xvar = xvar)
  }

  # @title Adaptive weighting (Prieto utils.adaptspec, iadapt = 0)
  # @description Thomson adaptive weights iterated (max 1000) until
  # max|(S-S_last)/(S+S_last)| <= 9.5e-7. The returned weights are those of the last
  # iteration (computed from the previous estimate), as in Python.
  # @keywords internal
  # @noRd
  amtm_adaptspec <- function(sk, lamb){
    nfft <- nrow(sk)
    k <- ncol(sk)
    lamb1 <- 1 - lamb
    dfa <- 1 / (nfft - 1)
    varsk <- colSums(sk) * dfa
    dvar <- mean(varsk)
    bk <- dvar * lamb1
    sqlamb <- sqrt(lamb)
    rerr <- 9.5e-7
    sbar <- (sk[, 1] + sk[, 2]) / 2
    bkm <- matrix(bk, nfft, k, byrow = TRUE)
    wt <- NULL
    for (i in 1:1000){
      slast <- sbar
      wt1 <- outer(sbar, sqlamb)
      wt2 <- outer(sbar, lamb) + bkm
      wt <- pmin(wt1 / wt2, 1)
      wt2s <- wt^2
      sbar <- rowSums(wt2s * sk) / rowSums(wt2s)
      oerr <- max(abs((sbar - slast) / (sbar + slast)))
      if (!(oerr > rerr)){
        break
      }
    }
    list(spec = sbar, wt = wt)
  }

  # ---------------------------------------------------------------------------
  # Slepian (DPSS) tapers, MATLAB dpss() conventions (unit energy)
  # ---------------------------------------------------------------------------

  acy_dpss_cache <- new.env(parent = emptyenv())

  # ---- main code ----
  data <- acy_check_two_column(data)
  tt <- as.numeric(data[, 1])
  xt <- as.numeric(data[, 2])
  n <- length(xt)
  if (n <= 2 * order + 2){
    stop("series too short for the requested AR order")
  }
  nfft <- as.integer(2^(ceiling(log2(n)) + nfft_add))
  xm <- xt - mean(xt)
  mt <- amtm_mtspec(x = xm, nw = nw, k = k, dt = 1, nfft = nfft)
  acv <- Re(stats::fft(mt$spec, inverse = TRUE)) / nfft
  tm <- matrix(0, order, order)
  for (i in seq_len(order)){
    for (j in seq_len(order)){
      tm[i, j] <- acv[abs(i - j) + 1]
    }
  }
  phi <- as.numeric(solve(tm, acv[2:(order + 1)]))
  nfit <- n - order
  fit <- numeric(nfit)
  for (j in seq_len(order)){
    fit <- fit + phi[j] * xm[(order + 1 - j):(n - j)]
  }
  pwx <- fit - xm[(order + 1):n]
  pwx <- pwx[(order + 1):nfit]
  pwx <- pwx - mean(pwx)
  list(data = cbind(tt[(2 * order + 1):n], pwx), ar_coef = phi, acvs = acv)
}
