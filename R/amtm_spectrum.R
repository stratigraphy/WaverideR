# amtm_spectrum.R
#
# Adaptive multitaper (aMTM) PSD background estimation and periodic signal
# detection (gamma-test + harmonic F-test), R translation of the Python package
# amtm_signaldetect v1.0.6 (https://github.com/hsalinasGIT/amtm_signaldetect)
# and of the parts of Prieto's 'multitaper' Python package (v1.2.0,
# https://github.com/gaprieto/multitaper) that amtm_signaldetect calls
# (MTSpec with iadapt = 0: dpss normalisation, eigenspec, adaptspec, rspec, ftest).
#
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
# The AR(1) background model, the AIC / MERIT model-selection criteria and the
# AR(1) starting value follow the IDL package SPD_MTM (c) Simone Di Matteo 2020,
# Apache License 2.0 (SPEDAS: general/science/spd_mtm/models/ar1/*.pro,
# spd_mtm_modelgof.pro, spd_mtm_fitmodel.pro). They are not part of the Python
# package. The prewhitening option follows the "Python PreWhiten Filter Example"
# distributed with amtm_signaldetect.
#
# Uses from helpers_matlab_equivalents.R: acy_dpss(), acy_check_two_column().

#' @title Adaptive multitaper spectrum with background fit and periodic signal detection
#'
#' @description The \code{\link{amtm_spectrum}} function computes the adaptive
#' multitaper (aMTM) power spectral density (PSD) of an evenly spaced series, fits
#' the PSD background by maximum likelihood with several noise models, selects the
#' best model and detects periodic signals with the dual gamma-test / harmonic
#' F-test criterion of Di Matteo et al. (2021) and Salinas et al. (2025).
#'
#' Method: (1) Thomson's adaptively weighted multitaper PSD (Prieto's 'multitaper'
#' implementation: unit-energy DPSS tapers, adaptive weights iterated until the
#' relative change is below 9.5e-7, PSD rescaled so that its integral equals the
#' variance, one-sided by doubling all bins except 0 and Nyquist). The half degrees
#' of freedom are \eqn{\alpha_k = (\sum_j w_{jk}^2)^2 / \sum_j w_{jk}^4}.
#' (2) Background models, fitted on the frequency range [2NW/T, f_Nyq - 2NW/T]
#' (or a user range) by minimising the aMTM log-likelihood
#' \eqn{M = 2\sum_k [\alpha_k S_k/B_k + \ln(\Gamma(\alpha_k) S_k) - \alpha_k\ln(\alpha_k S_k/B_k)]}:
#' white noise WHT \eqn{B = c} (analytical), power law PL \eqn{B = c f^{-\beta}},
#' bending power law BPL \eqn{B = c f^{-\beta}/(1 + (f/f_b)^{\gamma-\beta})}, and
#' first-order autoregressive AR1 \eqn{B = c/(1 - 2\rho\cos(\pi f/f_{Nyq}) + \rho^2)}.
#' (3) Model selection by the MERIT function (default in SPD_MTM) or AIC
#' (\eqn{M + 2 n_{free}}, with 0/1/1/3 free parameters besides c for WHT/PL/AR1/BPL).
#' Adjacent aMTM frequencies are correlated over the 2W band, so these criteria tend
#' to favour flexible models (BPL is sometimes chosen for pure white noise): inspect
#' \code{model_stats} and use \code{bkg_model} to force a model when needed.
#' The Python code minimises M with scipy SLSQP (or Powell) on the raw parameters;
#' here \code{stats::nlminb} is used on (log c, beta[, gamma, log fb]) with the same
#' bounds and starting values. Both reach the same optimum for PL (parameters within
#' ~1e-5 relative); for BPL R reaches the same or a lower (better) M, because SLSQP
#' can stop in a local minimum.
#' (4) Gamma-test \eqn{\gamma_k = S_k/B_k}, whose threshold is taken from the
#' mixture-of-gamma CDF built from the histogram (bin width 0.1) of \eqn{\alpha_k};
#' harmonic F-test with threshold \code{qf(conf, 2, 2K-2)}. A frequency is a
#' detected periodic signal when an F-test peak above its threshold falls inside a
#' band where \eqn{\gamma_k} exceeds its threshold and that band is at least
#' \eqn{NW/(2T)} wide.
#'
#' Translated from the Python package amtm_signaldetect v1.0.6
#' (Hector A. Salinas; functions get_amtm_specs, get_spectra_dofA,
#' get_background_psdfit, get_gamtest_confs, get_ftest_confs, get_gftest_confpeaks)
#' and from Prieto's Python 'multitaper' v1.2.0 (MTSpec, adaptspec, ftest), both MIT
#' licensed. The AR1 model and the AIC/MERIT criteria follow the IDL code SPD_MTM
#' (Di Matteo et al., 2021). This is not an Acycle function.
#'
#' Cyclostratigraphic interpretation: for a depth series the frequency axis is in
#' cycles per depth unit (e.g. cycles/m if depth is in m); the period (thickness) of a
#' detected cycle is 1/frequency. The fitted background is the null hypothesis of
#' "no cycle": use PL/BPL/AR1 for the red spectra typical of sedimentary proxies, and
#' treat the gamma-test as the test for excess power above the red background and
#' the F-test as the test for a phase-coherent sinusoid. Because sedimentation-rate
#' changes smear orbital lines, a real Milankovitch signal in depth often passes the
#' gamma-test over a band but gives only a moderate F-test; inspect both, use
#' \code{freq_range} to exclude trends/very high frequencies, and compare the ratios
#' of detected periods with expected orbital ratios (e.g. 405:100:~20 kyr) rather
#' than relying on single detections. The frequency resolution is \eqn{2NW/T} with
#' \eqn{T} the series length in depth units.
#'
#' @param data Input record: a matrix or data frame with depth/time in column 1 (evenly spaced; rows are sorted
#' by increasing depth/time) and the proxy value in column 2 (no \code{NA}/\code{Inf} values). Required
#' (\code{Default=NULL} gives an error).
#' @param nw Time-bandwidth product NW of the multitaper spectrum (numeric, > 0, typically 2 to 4). Larger values
#' give a smoother spectrum with lower variance but a coarser frequency resolution (2NW/T). \code{Default=3}.
#' @param k Number of DPSS tapers (integer, >= 2, at most about \code{2*nw - 1}); \code{NULL} uses
#' \code{round(2*nw - 1)}, as in the Python code. \code{Default=NULL}.
#' @param nfft Number of FFT points, i.e. zero padding (integer); \code{NULL} uses the series length N (the Python
#' default); values smaller than N are replaced by 2N+1, as in 'multitaper'. \code{Default=NULL}.
#' @param models Background models to fit (character vector, any combination of \code{"WHT"} (white noise),
#' \code{"PL"} (power law), \code{"AR1"} (first-order autoregressive, red noise) and \code{"BPL"} (bending power
#' law); case-insensitive). \code{Default=c("WHT", "PL", "AR1", "BPL")}.
#' @param selection Criterion used to choose the background model among \code{models} (character, case-insensitive):
#' \code{"MERIT"} (MERIT goodness-of-fit function of SPD_MTM) or \code{"AIC"} (Akaike information criterion);
#' the model with the lowest value is selected. \code{Default="MERIT"}.
#' @param bkg_model Background model to use for the detection instead of the automatically selected one
#' (character, one of the names in \code{models}); \code{NULL} uses the model chosen by \code{selection}.
#' \code{Default=NULL}.
#' @param freq_range Frequency range used for the background fit and the detection: numeric vector
#' \code{c(f_low, f_high)} with \code{f_low < f_high}, in cycles per unit of column 1 of \code{data} (e.g.
#' cycles/m); \code{NULL} uses [2NW/T, f_Nyq - 2NW/T]. At least 12 frequencies must fall in the range.
#' \code{Default=NULL}.
#' @param conf_level Confidence level of the gamma-test and F-test thresholds used for the detection (numeric,
#' strictly between 0 and 1, e.g. 0.95 or 0.99). \code{Default=0.99}.
#' @param prewhiten Prewhiten the series with an AR(\code{pw_order}) filter (via \code{\link{amtm_prewhiten}} with
#' its default \code{nw=4}, \code{k=7}) before the analysis (logical, TRUE or FALSE); the spectrum, fits and
#' detection are then made on the prewhitened series and a recoloured PSD is returned in \code{prewhiten}.
#' \code{Default=FALSE}.
#' @param pw_order Order of the AR prewhitening filter (integer, >= 1); only used when \code{prewhiten=TRUE}.
#' \code{Default=2}.
#' @param genplot Plot the results with \code{\link{plot_amtm_spectrum}} (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param verbose Print progress messages, the fitted parameters of each model and the selected model (logical,
#' TRUE or FALSE). \code{Default=FALSE}.
#'
#' @author Michiel Arts (R translation); original Python code: Hector A. Salinas
#' (amtm_signaldetect) and German A. Prieto (multitaper); original method and IDL code
#' SPD_MTM: Simone Di Matteo
#'
#' @references
#' Di Matteo, S., Viall, N. M., & Kepko, L. (2021). Power spectral density background
#' estimate and signal detection via the multitaper method. Journal of Geophysical
#' Research: Space Physics, 126, e2020JA028748. \doi{10.1029/2020JA028748}
#'
#' Salinas, H. A., et al. (2025). Journal of Geophysical Research: Space Physics.
#' \doi{10.1029/2025JA034006} (open-source Python aMTM signal detection)
#'
#' Prieto, G. A. (2022). The multitaper spectrum analysis package in Python.
#' Seismological Research Letters, 93(3), 1922-1929. \doi{10.1785/0220210332}
#'
#' Thomson, D. J. (1982). Spectrum estimation and harmonic analysis. Proceedings of
#' the IEEE, 70(9), 1055-1096. \doi{10.1109/PROC.1982.12433}
#'
#' @examples
#' \donttest{
#' # synthetic depth series: AR(1) red noise + 2 cycles (periods 5 m and 1.25 m)
#' set.seed(1)
#' n <- 600
#' depth <- seq(0, by = 0.1, length.out = n)
#' red <- as.numeric(stats::filter(rnorm(n), 0.8, method = "recursive"))
#' x <- red + 1.2 * sin(2 * pi * depth / 5) + 0.8 * sin(2 * pi * depth / 1.25)
#' res <- amtm_spectrum(data = cbind(depth, x), nw = 3, genplot = FALSE)
#' res$selected_model
#' res$detected           # frequency in cycles/m, period in m
#' plot_amtm_spectrum(res)
#' }
#'
#' @return A list with elements
#' \itemize{
#'   \item \code{freq}, \code{psd}, \code{alpha}, \code{ftest}: aMTM frequency
#'   (cycles per unit), one-sided PSD, half degrees of freedom and harmonic F-test over
#'   [0, f_Nyq].
#'   \item \code{fit_index}: indices of \code{freq} used for the fit/detection;
#'   \code{fit_freq}: those frequencies.
#'   \item \code{background}: matrix (fit range x models) of fitted backgrounds.
#'   \item \code{params}: list of fitted parameters per model (WHT: c; PL: c, beta;
#'   AR1: c, rho; BPL: c, beta, gamma, fb (after the beta > gamma correction)).
#'   \item \code{model_stats}: data frame with log-likelihood M, AIC, MERIT and
#'   convergence code per model; \code{selected_model}; \code{background_selected}.
#'   \item \code{gamma}: gamma-statistic S/B of the selected model (fit range).
#'   \item \code{thresholds}: data frame with confidence levels (0.5, 0.9, 0.95, 0.99
#'   and \code{conf_level}) and the gamma-test and F-test thresholds.
#'   \item \code{gamma_threshold}, \code{f_threshold}: thresholds at \code{conf_level}.
#'   \item \code{f_peaks}, \code{gamma_peaks}: indices (within the fit range) of F-test
#'   and gamma-test peaks above threshold; \code{gamma_bands}: data frame of gamma
#'   bands (f_start, f_end) wider than NW/(2T).
#'   \item \code{detected}: data frame of periodic signals passing both tests
#'   (frequency, period, ftest, gamma, psd, background).
#'   \item \code{prewhiten}: NULL or list(ar_coef, psd_recoloured).
#'   \item \code{settings}: list of the settings used (nw, k, nfft, dt, n, conf_level).
#' }
#' @export
#' @importFrom stats fft
#' @importFrom stats nlminb
#' @importFrom stats optimize
#' @importFrom stats qf
#' @importFrom stats pgamma
amtm_spectrum <- function(data = NULL,
    nw = 3,
    k = NULL,
    nfft = NULL,
    models = c("WHT", "PL", "AR1", "BPL"),
    selection = "MERIT",
    bkg_model = NULL,
    freq_range = NULL,
    conf_level = 0.99,
    prewhiten = FALSE,
    pw_order = 2,
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

  # internal helpers

  # @title aMTM time-series parameters
  # @description dt, N, Rayleigh frequency and Nyquist frequency exactly as
  # get_tseries_params() in amtm_signaldetect.
  # @keywords internal
  # @noRd
  amtm_tseries_params <- function(tdata){
    dt <- tdata[2] - tdata[1]
    n <- length(tdata)
    list(dt = dt, n = n, df = 1 / (n * dt), fny = 1 / (2 * dt))
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

  # @title One-sided spectrum (Prieto MTSpec.rspec)
  # @keywords internal
  # @noRd
  amtm_rspec <- function(spec, nf, nfft){
    s <- 2 * spec[seq_len(nf)]
    s[1] <- 0.5 * s[1]
    if (nfft %% 2 == 0){
      s[nf] <- 0.5 * s[nf]
    }
    s
  }

  # @title Half degrees of freedom (get_spectra_dofA)
  # @keywords internal
  # @noRd
  amtm_dof_alpha <- function(wt, nf){
    w <- wt[seq_len(nf), , drop = FALSE]
    rowSums(w^2)^2 / rowSums(w^4)
  }

  # @title Harmonic F-test (Prieto utils.ftest)
  # @keywords internal
  # @noRd
  amtm_ftest <- function(vn, yk){
    k <- ncol(vn)
    vn0 <- colSums(vn)
    vn0sq <- sum(vn0^2)
    mu <- as.vector(yk %*% vn0) / vn0sq
    fup <- (k - 1) * Mod(mu)^2 * vn0sq
    resid <- yk - outer(mu, vn0)
    fdw <- rowSums(Mod(resid)^2)
    fup / fdw
  }

  # @title Indices of the fitting frequency range (freqtrim_amtm_specs)
  # @description Python slice [first f >= f_low : first f >= f_high) (upper end
  # excluded).
  # @keywords internal
  # @noRd
  amtm_trim_index <- function(freq, prm, nw, freq_range = NULL){
    if (is.null(freq_range)){
      lo_val <- 2 * nw * prm$df
      up_val <- prm$fny - 2 * nw * prm$df
    } else {
      if (length(freq_range) != 2 || freq_range[1] >= freq_range[2]){
        stop("freq_range must be c(f_low, f_high) with f_low < f_high")
      }
      lo_val <- freq_range[1]
      up_val <- freq_range[2]
    }
    lo <- which(freq >= lo_val)
    up <- which(freq >= up_val)
    if (length(lo) == 0 || length(up) == 0){
      stop("freq_range outside the frequency axis (max = ", max(freq), ")")
    }
    lo <- lo[1]
    up <- up[1]
    if (up <= lo){
      stop("empty fitting frequency range")
    }
    seq(lo, up - 1)
  }

  # @title aMTM log-likelihood M (pl_loglikeM / bpl_loglikeM_ctoo)
  # @keywords internal
  # @noRd
  amtm_loglik <- function(b, s, a){
    m1 <- sum(a * s / b + log(gamma(a) * s))
    m2 <- sum(-a * log(a * s / b))
    2 * (m1 + m2)
  }

  # @title Background model lines
  # @keywords internal
  # @noRd
  amtm_pl_line <- function(f, p){
    p[1] * f^(-p[2])
  }

  # @keywords internal
  # @noRd
  amtm_bpl_line <- function(f, p){
    (p[1] * f^(-p[2])) / (1 + ((f / p[4])^(p[3] - p[2])))
  }

  # @keywords internal
  # @noRd
  amtm_ar1_line <- function(f, p, fny){
    p[1] / (1 - 2 * p[2] * cos(pi * f / fny) + p[2]^2)
  }

  # @title White-noise background (analytical ML solution)
  # @keywords internal
  # @noRd
  amtm_fit_wht <- function(f, s, a){
    c0 <- sum(a * s) / sum(a)
    list(par = c(c = c0), background = rep(c0, length(f)), convergence = 0L)
  }

  # @title Robust objective wrapper: returns a large finite value if invalid
  # @keywords internal
  # @noRd
  amtm_safe_M <- function(b, s, a){
    if (any(!is.finite(b)) || any(b <= 0)){
      return(1e300)
    }
    m <- amtm_loglik(b, s, a)
    if (!is.finite(m)){
      return(1e300)
    }
    m
  }

  # @title Power-law background (get_pl_guess + optimisation of [c, beta])
  # @description Starting values and bounds as in Python (c in [min S, max S],
  # beta in [0, 10]); optimised with nlminb on (log c, beta) instead of SLSQP.
  # @keywords internal
  # @noRd
  amtm_fit_pl <- function(f, s, a){
    n <- length(f)
    if (n < 11){
      bet0 <- log10(s[1] / s[n]) / log10(f[n] / f[1])
    } else {
      bet0 <- mean(log10(s[1:5] / s[(n - 5):(n - 1)])) /
          mean(log10(f[(n - 5):(n - 1)] / f[1:5]))
    }
    c0 <- sum(a * s * (f^bet0)) / sum(a)
    lower <- c(log(min(s)), 0)
    upper <- c(log(max(s)), 10)
    start <- pmin(pmax(c(log(c0), bet0), lower), upper)
    obj <- function(p){
      amtm_safe_M(exp(p[1]) * f^(-p[2]), s, a)
    }
    op <- amtm_nlminb(start, obj, lower, upper)
    par <- c(c = exp(op$par[1]), beta = op$par[2])
    list(par = par, background = amtm_pl_line(f, par), convergence = op$convergence,
        start = c(c0, bet0))
  }

  # @title Bending power-law background (get_bpl_guess + optimisation + correction)
  # @description Starting values and bounds as in Python (c in [min S, max S],
  # beta in [-5, 10], gamma in [0, 15], fb in [f_1, f_n]); optimised with nlminb on
  # (log c, beta, gamma, log fb); if beta > gamma the parameters are swapped
  # (correct_bpl_param).
  # @keywords internal
  # @noRd
  amtm_fit_bpl <- function(f, s, a){
    n <- length(f)
    dmid <- as.integer(n / 2)
    fb0 <- f[dmid + 1]
    if (n < 11){
      bet0 <- log(s[1] / s[dmid + 1]) / log(fb0 / f[1])
      gam0 <- log(s[dmid + 1] / s[n]) / log(f[n] / fb0)
    } else {
      bet0 <- mean(log(s[1:5] / s[(dmid - 4):dmid])) / mean(log(fb0 / f[1:5]))
      gam0 <- mean(log(s[(dmid + 1):(dmid + 5)] / s[(n - 5):(n - 1)])) /
          mean(log(f[(n - 5):(n - 1)] / fb0))
    }
    c0 <- sum(a * s * (f^bet0) * (1 + (f / fb0)^(gam0 - bet0))) / sum(a)
    lower <- c(log(min(s)), -5, 0, log(f[1]))
    upper <- c(log(max(s)), 10, 15, log(f[n]))
    start <- pmin(pmax(c(log(c0), bet0, gam0, log(fb0)), lower), upper)
    obj <- function(p){
      amtm_safe_M(amtm_bpl_line(f, c(exp(p[1]), p[2], p[3], exp(p[4]))), s, a)
    }
    op <- amtm_nlminb(start, obj, lower, upper)
    par <- c(exp(op$par[1]), op$par[2], op$par[3], exp(op$par[4]))
    if (par[2] > par[3]){
      par <- c(par[1] * par[4]^(par[3] - par[2]), par[3], par[2], par[4])
    }
    names(par) <- c("c", "beta", "gamma", "fb")
    list(par = par, background = amtm_bpl_line(f, par), convergence = op$convergence,
        start = c(c0, bet0, gam0, fb0))
  }

  # @title AR(1) background (SPD_MTM spd_mtm_ar1_*): c profiled, rho in [0, 1]
  # @keywords internal
  # @noRd
  amtm_fit_ar1 <- function(f, s, a, fny, rho0 = NULL){
    prof <- function(rho){
      den <- 1 - 2 * rho * cos(pi * f / fny) + rho^2
      cc <- sum(a * s * den) / sum(a)
      c(cc, rho)
    }
    obj <- function(rho){
      amtm_safe_M(amtm_ar1_line(f, prof(rho), fny), s, a)
    }
    op <- stats::optimize(obj, interval = c(0, 1), tol = 1e-12)
    best <- op$minimum
    if (!is.null(rho0) && is.finite(rho0) && rho0 > 0 && rho0 < 1){
      # SPD_MTM starts from the Yule-Walker lag-1 estimate; keep it if Brent's
      # search ends in a worse (local) minimum
      if (obj(rho0) < op$objective){
        best <- rho0
      }
    }
    par <- prof(best)
    names(par) <- c("c", "rho")
    list(par = par, background = amtm_ar1_line(f, par, fny), convergence = 0L)
  }

  # @title Bounded nlminb with a polishing restart
  # @keywords internal
  # @noRd
  amtm_nlminb <- function(start, obj, lower, upper){
    ctrl <- list(eval.max = 5000, iter.max = 3000)
    op <- stats::nlminb(start, obj, lower = lower, upper = upper, control = ctrl)
    op2 <- stats::nlminb(op$par, obj, lower = lower, upper = upper, control = ctrl)
    # the restart only polishes; when it cannot improve, nlminb reports a benign
    # "singular convergence", so the first run (and its code) is kept
    if (op2$objective < op$objective){
      op <- op2
    }
    op
  }

  # @title numpy.arange(start, stop, step) for doubles
  # @keywords internal
  # @noRd
  amtm_np_arange <- function(start, stop, step){
    len <- ceiling((stop - start) / step)
    if (len <= 0){
      return(numeric(0))
    }
    if (len == 1){
      return(start)
    }
    delta <- (start + step) - start
    out <- start + (0:(len - 1)) * delta
    out[2] <- start + step
    out
  }

  # @title Gamma-test thresholds (get_gamtest_confs)
  # @description CDF of gamma = S/B as a histogram-weighted mixture (bin width 0.1,
  # numpy.histogram density, left bin edges used as alpha) of regularised lower
  # incomplete gamma functions evaluated on z = 0, 0.001, ..., 14.999. The threshold
  # is the first z with CDF >= level (NA if never reached).
  # @keywords internal
  # @noRd
  amtm_gamma_thresholds <- function(alpha, levels){
    edges <- amtm_np_arange(min(alpha), max(alpha) + 0.1, 0.1)
    nb <- length(edges) - 1
    if (nb < 1){
      edges <- c(min(alpha), min(alpha) + 0.1)
      nb <- 1
    }
    counts <- numeric(nb)
    for (i in seq_len(nb)){
      if (i < nb){
        counts[i] <- sum(alpha >= edges[i] & alpha < edges[i + 1])
      } else {
        counts[i] <- sum(alpha >= edges[i] & alpha <= edges[i + 1])
      }
    }
    dalpha <- diff(edges)
    palpha <- counts / dalpha / sum(counts)
    z <- amtm_np_arange(0, 15, 0.001)
    cdf <- numeric(length(z))
    for (i in seq_len(nb)){
      cdf <- cdf + stats::pgamma(z * edges[i], shape = edges[i]) * palpha[i] * dalpha[i]
    }
    thr <- numeric(length(levels))
    for (j in seq_along(levels)){
      w <- which(cdf >= levels[j])
      if (length(w) == 0){
        thr[j] <- NA_real_
      } else {
        thr[j] <- z[w[1]]
      }
    }
    thr
  }

  # @title scipy.signal.find_peaks(x, height = h)
  # @description Local maxima (plateaus: middle sample, rounded down) excluding the
  # end points, with x[peak] >= height.
  # @keywords internal
  # @noRd
  amtm_find_peaks <- function(x, height){
    n <- length(x)
    peaks <- integer(0)
    i <- 2
    imax <- n - 1
    while (i <= imax){
      if (x[i - 1] < x[i]){
        ia <- i + 1
        while (ia <= imax && x[ia] == x[i]){
          ia <- ia + 1
        }
        if (x[ia] < x[i]){
          left <- i
          right <- ia - 1
          peaks <- c(peaks, (left + right) %/% 2)
          i <- ia
        }
      }
      i <- i + 1
    }
    peaks[x[peaks] >= height]
  }

  # @title Dual gamma/F-test detection (get_gftest_confpeaks)
  # @description Reproduces the Python logic including np.split behaviour: runs of
  # gamma > threshold are split at the zeros and the first element of every piece is
  # dropped (so a run starting at the first frequency loses its first sample). Bands
  # narrower than NW*df/2 are discarded; F-test peaks inside a band are detections.
  # @keywords internal
  # @noRd
  amtm_detect <- function(freq, gam, ftest, f_crit, g_crit, nw, df){
    n <- length(freq)
    g_peaks <- amtm_find_peaks(gam, g_crit)
    f_peaks <- amtm_find_peaks(ftest, f_crit)
    dummy <- as.numeric(gam > g_crit)
    zeros <- which(dummy == 0)
    starts <- c(1, zeros)
    ends <- c(zeros - 1, n)
    b_start <- numeric(0)
    b_end <- numeric(0)
    for (i in seq_along(starts)){
      len <- ends[i] - starts[i] + 1
      if (len > 1){
        fr <- freq[(starts[i] + 1):ends[i]]
        if ((max(fr) - min(fr)) >= (nw * df / 2)){
          b_start <- c(b_start, min(fr))
          b_end <- c(b_end, max(fr))
        }
      }
    }
    dual <- integer(0)
    for (i in seq_along(b_start)){
      for (j in seq_along(f_peaks)){
        fj <- freq[f_peaks[j]]
        if (fj >= b_start[i] && fj <= b_end[i]){
          dual <- c(dual, f_peaks[j])
        }
      }
    }
    list(f_peaks = f_peaks, g_peaks = g_peaks,
        bands = data.frame(f_start = b_start, f_end = b_end), dual_index = dual)
  }

  # @title |H_AR(f)|^2 of the prewhitening filter (correct_spec_for_pw)
  # @description |1 - sum_j phi_j exp(-i 2 pi f j dt)|^2; Python additionally
  # multiplies by dt (identical for dt = 1); omitted here so that the recoloured PSD
  # keeps the units of the aMTM PSD.
  # @keywords internal
  # @noRd
  amtm_ar_h2 <- function(phi, freq, dt){
    h <- rep(1 + 0i, length(freq))
    for (j in seq_along(phi)){
      h <- h - phi[j] * exp(-1i * 2 * pi * freq * j * dt)
    }
    Mod(h)^2
  }

  # Slepian (DPSS) tapers, MATLAB dpss() conventions (unit energy)

  acy_dpss_cache <- new.env(parent = emptyenv())

  # Main code
  data <- acy_check_two_column(data)
  tdata <- as.numeric(data[, 1])
  xdata <- as.numeric(data[, 2])
  if (any(!is.finite(xdata)) || any(!is.finite(tdata))){
    stop("data must not contain NA/NaN/Inf values")
  }
  if (is.null(k)){
    k <- as.integer(round(2 * nw - 1))
  }
  if (k < 2){
    stop("at least 2 tapers are required (k >= 2)")
  }
  models <- toupper(models)
  for (m in models){
    if (!(m %in% c("WHT", "PL", "AR1", "BPL"))){
      stop("unknown background model: ", m)
    }
  }
  selection <- toupper(selection)
  if (!(selection %in% c("MERIT", "AIC"))){
    stop("selection must be 'MERIT' or 'AIC'")
  }
  if (conf_level <= 0 || conf_level >= 1){
    stop("conf_level must be between 0 and 1")
  }

  pw <- NULL
  if (prewhiten == TRUE){
    pw_res <- amtm_prewhiten(data = cbind(tdata, xdata), order = pw_order)
    tdata <- pw_res$data[, 1]
    xdata <- pw_res$data[, 2]
    pw <- list(ar_coef = pw_res$ar_coef)
  }

  # 1) aMTM spectral data products (get_amtm_specs)
  prm <- amtm_tseries_params(tdata)
  n <- prm$n
  if (is.null(nfft)){
    nfft <- n
  }
  x0 <- xdata - mean(xdata)
  if (verbose == TRUE){
    cat("Nyquist frequency =", prm$fny, ", dt =", prm$dt, ", T =",
        tdata[n] - tdata[1], ", N =", n, "\n")
  }
  mt <- amtm_mtspec(x = x0, nw = nw, k = k, dt = prm$dt, nfft = nfft)
  freq <- mt$freq
  psd <- amtm_rspec(mt$spec, mt$nf, mt$nfft)
  alpha <- amtm_dof_alpha(mt$wt, mt$nf)
  ftest <- amtm_ftest(mt$vn, mt$yk[seq_len(mt$nf), , drop = FALSE])

  # 2) background fits (get_background_psdfit)
  idx <- amtm_trim_index(freq, prm, nw, freq_range)
  f_in <- freq[idx]
  s_in <- psd[idx]
  a_in <- alpha[idx]
  n_in <- length(idx)
  if (n_in < 12){
    stop("too few frequencies in the fitting range (", n_in, ")")
  }
  bkg <- matrix(NA_real_, n_in, length(models))
  colnames(bkg) <- models
  params <- list()
  stats_df <- data.frame(model = models, n_par = NA_integer_, loglik_M = NA_real_,
      AIC = NA_real_, MERIT = NA_real_, convergence = NA_integer_,
      stringsAsFactors = FALSE)
  nfree_tab <- c(WHT = 0, PL = 1, AR1 = 1, BPL = 3)
  rho0 <- sum(x0[-n] * x0[-1]) / sum(x0[-n]^2)
  for (i in seq_along(models)){
    m <- models[i]
    if (m == "WHT"){
      ft <- amtm_fit_wht(f_in, s_in, a_in)
    } else if (m == "PL"){
      ft <- amtm_fit_pl(f_in, s_in, a_in)
    } else if (m == "AR1"){
      ft <- amtm_fit_ar1(f_in, s_in, a_in, fny = prm$fny, rho0 = rho0)
    } else {
      ft <- amtm_fit_bpl(f_in, s_in, a_in)
    }
    bkg[, i] <- ft$background
    params[[m]] <- ft$par
    nfree <- nfree_tab[[m]]
    stats_df$n_par[i] <- length(ft$par)
    stats_df$loglik_M[i] <- amtm_loglik(ft$background, s_in, a_in)
    stats_df$AIC[i] <- stats_df$loglik_M[i] + 2 * nfree
    stats_df$MERIT[i] <- sum(a_in * ((s_in - ft$background) / ft$background)^2) /
        (n_in - nfree)
    stats_df$convergence[i] <- ft$convergence
    if (verbose == TRUE){
      cat(sprintf("  %-4s params: %s\n", m, paste(signif(ft$par, 8), collapse = ", ")))
    }
  }
  score <- stats_df[[selection]]
  score[!is.finite(score)] <- Inf
  selected <- models[which.min(score)]
  if (!is.null(bkg_model)){
    bkg_model <- toupper(bkg_model)
    if (!(bkg_model %in% models)){
      stop("bkg_model must be one of the fitted models")
    }
    selected <- bkg_model
  }
  b_sel <- bkg[, selected]
  if (verbose == TRUE){
    cat("Selected background model (", selection, "): ", selected, "\n", sep = "")
  }

  # 3) confidence thresholds (get_gamtest_confs / get_ftest_confs)
  gam <- s_in / b_sel
  levels <- unique(c(0.5, 0.9, 0.95, 0.99, conf_level))
  g_thr <- amtm_gamma_thresholds(a_in, levels)
  f_thr <- stats::qf(levels, 2, 2 * (k - 1))
  thresholds <- data.frame(conf_level = levels, gamma = g_thr, ftest = f_thr)
  g_crit <- g_thr[match(conf_level, levels)]
  f_crit <- f_thr[match(conf_level, levels)]

  # 4) dual-test detection (get_gftest_confpeaks)
  det <- amtm_detect(f_in, gam, ftest[idx], f_crit, g_crit, nw, prm$df)
  dind <- det$dual_index
  detected <- data.frame(frequency = f_in[dind], period = 1 / f_in[dind],
      ftest = ftest[idx][dind], gamma = gam[dind],
      psd = s_in[dind], background = b_sel[dind])

  if (!is.null(pw)){
    pw$psd_recoloured <- psd / amtm_ar_h2(pw$ar_coef, freq, prm$dt)
  }

  res <- list(freq = freq, psd = psd, alpha = alpha, ftest = ftest,
      fit_index = idx, fit_freq = f_in, background = bkg, params = params,
      model_stats = stats_df, selected_model = selected,
      background_selected = b_sel, gamma = gam, thresholds = thresholds,
      gamma_threshold = g_crit, f_threshold = f_crit,
      f_peaks = det$f_peaks, gamma_peaks = det$g_peaks,
      gamma_bands = det$bands, detected = detected, prewhiten = pw,
      settings = list(nw = nw, k = k, nfft = mt$nfft, dt = prm$dt, n = n,
          conf_level = conf_level, selection = selection))
  if (genplot == TRUE){
    plot_amtm_spectrum(res)
  }
  res
}
