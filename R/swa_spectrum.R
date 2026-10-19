#' @title Power spectrum with SWA background, chi-square and FDR confidence levels
#'
#' @description Calculates a power spectrum, locates the spectral background
#' with the Smoothed Window Average (SWA) method (Kemp, 2016; Weedon, 2020,
#' 2022) and derives chi-square confidence levels and False Discovery Rate
#' (FDR) confidence levels (Benjamini and Hochberg, 1995) following the
#' procedure of Miller et al. (2001), with the correlation of neighbouring
#' spectral estimates accounted for as in Hopkins et al. (2002).
#'
#' Two spectral estimators are available:
#' \describe{
#' \item{\code{method = "mtm"}}{Multi-taper spectrum of the (evenly spaced,
#' zero-mean) series, SWA background from \code{\link{swa_background}}.
#' Degrees of freedom = 2K (K = number of tapers). This is the SWA option
#' of the Acycle spectral analysis window.
#' Translated from Acycle v3.1 spectralswafdr.m (Mingsong Li).}
#' \item{\code{method = "lomb"}}{The program of Weedon (2020, 2022) used by
#' the Acycle SWA window: lag-1 check, linear detrending, 5\% cosine taper,
#' Lomb-Scargle transform (routine PERIOD of Press et al., 1992; allows
#' unevenly spaced data or data with gaps), three passes of a Hanning
#' spectral window (8 degrees of freedom), SWA background, chi-square and
#' FDR confidence levels.
#' Translated from Acycle v3.1 specswafdr.m (Graham Weedon, Mingsong Li).}
#' }
#'
#' @param data Matrix or data frame with two columns: depth/time and proxy
#' value. Rows are sorted on the first column and non-finite rows removed.
#' For \code{method = "mtm"} the data must be evenly spaced; as in Acycle
#' the series is used as is (remove the mean/trend beforehand).
#' @param method \code{"mtm"} or \code{"lomb"} (see Description).
#' \code{Default="mtm"}
#' @param nw Time-bandwidth product of the multi-taper spectrum
#' (\code{method = "mtm"} only). \code{Default=2}
#' @param pad_times Zero-padding factor; nfft = max(n, round(pad_times * n))
#' (\code{method = "mtm"} only). \code{Default=1}
#' @param background Background fit used for \code{method = "lomb"}:
#' \code{"specswa"} uses \code{\link{swa_background}} (current specswa.m
#' algorithm, also used by \code{method = "mtm"}); \code{"legacy"}
#' reproduces the older SWA fit coded inside specswafdr.m, i.e. exactly
#' what the Acycle SWA window reports. \code{Default="specswa"}
#' @param genplot Plot the spectrum, background and confidence levels.
#' \code{Default=FALSE}
#' @param verbose Print progress information. \code{Default=FALSE}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li,
#' translated from the FORTRAN programs of Graham P. Weedon.
#'
#' @references
#' Benjamini, Y., Hochberg, Y., 1995. Controlling the false discovery rate:
#' a practical and powerful approach to multiple testing. Journal of the
#' Royal Statistical Society B 57, 289-300.
#'
#' Hopkins, A.M., Miller, C.J., Connolly, A.J., Genovese, C., Nichol, R.C.,
#' Wasserman, L., 2002. A new source detection algorithm using the
#' false-discovery rate. The Astronomical Journal 123, 1086-1094.
#'
#' Kemp, D.B., 2016. Optimizing significance testing of astronomical forcing
#' in cyclostratigraphy. Paleoceanography 31, 1516-1531.
#'
#' Miller, C.J., Genovese, C., Nichol, R.C., Wasserman, L., Connolly, A.,
#' Reichart, D., Hopkins, A., Schneider, J., Moore, A., 2001. Controlling
#' the false-discovery rate in astrophysical data analysis. The Astronomical
#' Journal 122, 3492-3505.
#'
#' Press, W.H., Teukolsky, S.A., Vetterling, W.T., Flannery, B.P., 1992.
#' Numerical Recipes in FORTRAN, 2nd ed. Cambridge University Press.
#'
#' Weedon, G.P., 2020. Confirmed detection of Palaeogene and Jurassic
#' orbitally-forced sedimentary cycles in the depth domain using False
#' Discovery Rates and Bayesian probability spectra. Boletin Geologico y
#' Minero 131, 207-230.
#'
#' Weedon, G.P., 2022. Problems with the current practice of spectral
#' analysis in cyclostratigraphy: avoiding false detection of regular
#' cyclicity. Earth-Science Reviews 235, 104261.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 600
#' t <- (0:(n - 1)) * 0.5
#' noise <- as.numeric(stats::filter(stats::rnorm(n), 0.7,
#'   method = "recursive"))
#' x <- noise + 1.5 * sin(2 * pi * t / 20) + sin(2 * pi * t / 7.3)
#' dat <- cbind(t, x - mean(x))
#' res_mtm <- swa_spectrum(data = dat, method = "mtm", nw = 2,
#'   genplot = TRUE)
#' res_lomb <- swa_spectrum(data = dat, method = "lomb", genplot = TRUE)
#' res_lomb$fdr
#' }
#'
#' @return A list with elements
#' \describe{
#' \item{freq}{Frequencies.}
#' \item{power}{Spectral power used for testing (MTM PSD, or the
#' Hanning-smoothed Lomb-Scargle spectrum normalised by the variance).}
#' \item{periodogram}{(\code{"lomb"} only) raw normalised Lomb-Scargle
#' periodogram before Hanning smoothing.}
#' \item{swa}{SWA background.}
#' \item{chi2_factor}{Multiplication factors of the background for the
#' 90, 95, 99, 99.9 and 99.99\% chi-square confidence levels.}
#' \item{chi2_cl}{Matrix with the chi-square confidence levels (columns
#' chi90, chi95, chi99, chi999, chi9999).}
#' \item{fdr_cl}{Matrix with the FDR confidence levels (columns fdr10,
#' fdr5, fdr1, fdr01, fdr001 for 10, 5, 1, 0.1 and 0.01\% FDR). As in
#' Acycle a column is NA when the level is not reached or does not differ
#' from the next less strict level; the 10\% column is always NA in Acycle
#' (kept for compatibility).}
#' \item{fdr}{Data frame per FDR level: \code{n_rejected} (number of
#' frequencies declared significant, -1 = none), \code{factor}
#' (multiplication factor of the background), \code{chi2_cl_percent}
#' (equivalent chi-square confidence level) and \code{drawn} (whether the
#' level is present in \code{fdr_cl}).}
#' \item{alpha}{Probability (p-value on the Acycle grid) of every spectral
#' estimate.}
#' \item{variance_ratio}{power / swa.}
#' \item{dof}{Degrees of freedom.}
#' \item{cn}{Hopkins et al. (2002) correction constant.}
#' \item{bandwidth, rayleigh}{Resolution bandwidth and Rayleigh frequency.}
#' \item{n_windows, diagnostics}{Selected SWA window width and fit
#' diagnostics.}
#' \item{rho1}{(\code{"lomb"} only) lag-1 autocorrelation of the input.}
#' \item{irregular}{(\code{"lomb"} only) TRUE if the spacing varies by more
#' than 5\%.}
#' }
#'
#' @export
#' @importFrom stats median var qchisq
#' @importFrom graphics plot lines legend par
swa_spectrum <- function(data = NULL, method = "mtm", nw = 2, pad_times = 1,
                         background = "specswa", genplot = FALSE,
                         verbose = FALSE){
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

  # @title Lomb-Scargle spectrum of Weedon's SWA program (specswafdr.m)
  # @description Spacing check, lag-1 autocorrelation, linear detrending,
  # 5 percent cosine taper, Lomb-Scargle transform (Press et al., 1992) and
  # three Hanning passes, exactly as specswafdr.m.
  # @keywords internal
  # @noRd
  swa_lomb_power <- function(time = NULL, x = NULL, verbose = FALSE){
    ndata <- length(x)
    andata <- ndata
    sp <- swa_spacing(time)
    tdif <- sp$tdif
    if (sp$nrev > 0){
      warning("Equal or reversed time/position steps detected")
    }
    if (sp$nint >= 1){
      nirr <- 1
      if (sp$tintmax / tdif > 0.3){
        warning(paste0("Largest time step > 30% of total interval; analyse ",
                       "the data in shorter segments to ensure stationarity"))
      }
      if (verbose == TRUE){
        cat("The data are not uniformly spaced: Lomb-Scargle on irregular data\n")
      }
    } else {
      nirr <- 0
    }
    varx <- stats::var(x)
    rho1 <- swa_rho(x, varx)
    x1 <- swa_detrend(time, x)
    x1 <- swa_ctaper(x1, time, tdif)

    xfint <- 1.0 / tdif
    tint <- (max(time) - min(time)) / (ndata - 1)
    bw <- 1.0 / (0.225586 * 1.055 * andata * tint)
    if (nirr == 0){
      xnyq <- 1.0 / (2.0 * tint)
      nout <- xnyq / xfint
    } else {
      xnyq <- 1.0 / (2.0 * sp$tintmin)
      nout <- tdif / (2.0 * sp$tintmin)
      bw <- 1.0 / (0.225586 * 1.055 * andata * sp$tintmin)
      xndof <- (2.0 * andata) / (1.055 * 0.225586 * 2.0 * nout)
      ndof <- xndof + 0.5
      if (ndof <= 7){
        nout <- trunc(ndata / 2)
        xnyq <- xfint * nout
        bw <- 1.0 / (0.225586 * 1.055 * andata * tint)
      }
    }
    if (verbose == TRUE){
      cat(sprintf("Time range = %.7f, Rayleigh = %.7f, Nyquist = %.7f, nout = %g\n",
                  tdif, xfint, xnyq, nout))
    }
    avex <- mean(x1)
    varx1 <- stats::var(x1)
    per <- swa_period(x1, time, avex, floor(nout), xfint, varx1)
    pow1 <- swa_hanning3(per$pow)
    list(freq = per$freq, power = pow1, periodogram = per$pow, rho1 = rho1,
         bandwidth = bw, rayleigh = xfint, nyquist = xnyq,
         irregular = (nirr == 1))
  }

  # @title Sample spacing check of specswafdr.m/specbayes.m
  # @keywords internal
  # @noRd
  swa_spacing <- function(time = NULL){
    ndata <- length(time)
    tdif <- time[ndata] - time[1]
    tint1 <- time[2] - time[1]
    tintmin <- tint1
    tintmax <- tint1
    tintva <- tint1 * 1.05
    tintvb <- tint1 * 0.95
    nint <- 0
    nrev <- 0
    if (ndata > 2){
      d <- time[3:ndata] - time[2:(ndata - 1)]
      nrev <- sum(d <= 0)
      hi <- d > tintva
      lo <- d < tintvb
      nint <- sum(hi) + sum(lo)
      if (any(hi)){
        tintmax <- max(c(tintmax, d[hi]))
      }
      if (any(lo)){
        tintmin <- min(c(tintmin, d[lo]))
      }
    }
    list(tdif = tdif, tintmin = tintmin, tintmax = tintmax, nint = nint,
         nrev = nrev)
  }

  # @title Lag-1 autocorrelation (rho.m, AUTOCR of Davis 1973)
  # @keywords internal
  # @noRd
  swa_rho <- function(x = NULL, varx = NULL){
    n <- length(x)
    l <- n - 1
    sx <- sum(x[1:l])
    sy <- sum(x[2:n])
    sxy <- sum(x[1:l] * x[2:n])
    rho1 <- (l * sxy - sx * sy) / (l * (l - 1.0)) / varx
    if (rho1 < 0.0 || rho1 > 1.0){
      stop(sprintf(paste0("Lag-1 autocorrelation RHO = %f is outside [0, 1]; ",
                          "the SWA/Bayes programs of Weedon require red noise"),
                   rho1))
    }
    rho1
  }

  # @title Least-squares linear detrending (lfitls.m)
  # @keywords internal
  # @noRd
  swa_detrend <- function(time = NULL, x = NULL){
    an <- length(x)
    sx <- sum(time)
    sy <- sum(x)
    tt <- time - sx / an
    st2 <- sum(tt * tt)
    bb <- sum(tt * x) / st2
    aa <- (sy - (sx * bb)) / an
    x - ((bb * time) + aa)
  }

  # @title 5 percent split cosine taper (ctaper.m, TAPER of Bloomfield 1976)
  # @keywords internal
  # @noRd
  swa_ctaper <- function(x = NULL, time = NULL, tdif = NULL){
    ndata <- length(x)
    xpi <- 3.1415926535897932
    bit <- tdif / (10.0 * ndata)
    if (time[ndata] - time[1] > 0.0){
      xprop <- (time[ndata] - time + bit) / tdif
    } else {
      xprop <- (time - time[ndata] + bit) / tdif
    }
    sel <- xprop < 0.05 | xprop > 0.95
    x[sel] <- x[sel] * (0.5 - 0.5 * cos(xpi * (xprop[sel] / 0.05)))
    x
  }

  # @title Lomb-Scargle transform (period.m, PERIOD of Press et al. 1992)
  # @description Uses the same trigonometric recurrences as the original.
  # Returns power normalised by the variance at frequencies i * xfint.
  # @keywords internal
  # @noRd
  swa_period <- function(x = NULL, time = NULL, avex = NULL, nout = NULL,
                         xfint = NULL, varx = NULL){
    ndata <- length(x)
    twopid <- 6.2831853071795865
    tave <- 0.5 * (time[1] + time[ndata])
    pnow <- xfint
    arg <- twopid * ((time - tave) * pnow)
    wpr <- -2.0 * sin(0.5 * arg)^2
    wpi <- sin(arg)
    wr <- cos(arg)
    wi <- wpi
    yy <- x - avex
    freq <- numeric(nout)
    pow <- numeric(nout)
    for (i in seq_len(nout)){
      freq[i] <- pnow
      sumsh <- sum(wi * wr)
      sumc <- sum((wr - wi) * (wr + wi))
      wtau <- 0.5 * atan2(2.0 * sumsh, sumc)
      swtau <- sin(wtau)
      cwtau <- cos(wtau)
      ss <- wi * cwtau - wr * swtau
      cc <- wr * cwtau + wi * swtau
      sums <- sum(ss^2)
      sumc <- sum(cc^2)
      sumsy <- sum(yy * ss)
      sumcy <- sum(yy * cc)
      wtemp <- wr
      wr <- (wr * wpr - wi * wpi) + wr
      wi <- (wi * wpr + wtemp * wpi) + wi
      pow[i] <- 0.5 * (sumcy^2 / sumc + sumsy^2 / sums) / varx
      pnow <- pnow + xfint
    }
    list(freq = freq, pow = pow)
  }

  # @title Three passes of a Hanning (1/4, 1/2, 1/4) spectral window
  # @keywords internal
  # @noRd
  swa_hanning3 <- function(p = NULL){
    nout <- length(p)
    mid <- 2:(nout - 1)
    for (k in 1:3){
      ps <- numeric(nout)
      ps[1] <- 0.5 * p[1] + 0.5 * p[2]
      ps[nout] <- 0.5 * p[nout - 1] + 0.5 * p[nout]
      ps[mid] <- 0.25 * p[mid - 1] + 0.5 * p[mid] + 0.25 * p[mid + 1]
      p <- ps
    }
    p
  }

  # @title Acycle probability grid for the FDR test
  # @description alpha = 1e-9 ... 1e-4 in steps of 1e-9, then steps of 1e-7
  # (1.1e6 levels, accumulated sequentially as MATLAB cumsum) and the
  # corresponding chi-square variance ratios qchisq(1 - alpha, dof) / dof.
  # @keywords internal
  # @noRd
  swa_fdr_grid <- function(dof = NULL){
    key <- format(dof, digits = 17)
    if (exists(key, envir = swa_fdr_cache, inherits = FALSE)){
      return(get(key, envir = swa_fdr_cache))
    }
    nlevs <- 1100000
    n1 <- 100001
    alpha1 <- numeric(n1)
    acc <- 0
    for (i in seq_len(n1)){
      acc <- acc + 0.000000001
      alpha1[i] <- acc
    }
    n2 <- nlevs - n1
    alpha2 <- numeric(n2)
    acc <- 0
    for (i in seq_len(n2)){
      acc <- acc + 0.0000001
      alpha2[i] <- acc
    }
    alpha <- c(alpha1, alpha2 + alpha1[n1])
    factor <- stats::qchisq(1 - alpha, dof) / dof
    out <- list(alpha = alpha, factor = factor)
    assign(key, out, envir = swa_fdr_cache)
    out
  }

  # @title FDR confidence levels (Miller et al. 2001) as in spectralswafdr.m
  # @keywords internal
  # @noRd
  swa_fdr <- function(power = NULL, swa = NULL, dof = NULL, cn = NULL,
                      anout = NULL){
    nout <- length(power)
    grid <- swa_fdr_grid(dof)
    alpha <- grid$alpha
    factor <- grid$factor
    nlevs <- length(alpha)
    factob <- power / swa
    # j = last index with factor > factob (factor decreases with alpha)
    if (is.unsorted(-factor) == FALSE){
      j <- findInterval(-factob, -factor, left.open = TRUE)
    } else {
      j <- numeric(nout)
      for (i in seq_len(nout)){
        w <- which(factor > factob[i])
        if (length(w) > 0){
          j[i] <- max(w)
        }
      }
    }
    idx <- rep(0, nout)
    idx[factor[1] < factob] <- 1
    idx[j > 0 & j < nlevs] <- j[j > 0 & j < nlevs] + 1
    idx[j == nlevs] <- nlevs
    if (any(idx == 0)){
      stop("ALPHOB not determined (variance ratio not finite)")
    }
    alphob <- alpha[idx]
    ord <- order(alphob)
    alphsort <- alphob[ord]
    idxsort <- idx[ord]

    levels <- c(0.1, 0.05, 0.01, 0.001, 0.0001)
    maxneg <- rep(-1, 5)
    rank <- seq_len(nout)
    for (k in 1:5){
      jalpha <- rank * levels[k] / (cn * anout)
      neg <- which(alphsort > 0 & (alphsort - jalpha) < 0)
      if (length(neg) > 0){
        maxneg[k] <- max(neg)
      }
    }
    fac <- rep(NA_real_, 5)
    pct <- rep(NA_real_, 5)
    for (k in 1:5){
      if (maxneg[k] > 0){
        fac[k] <- factor[idxsort[maxneg[k]]]
        pct[k] <- 100.0 * (1.0 - alphsort[maxneg[k]])
      }
    }
    # drawing rules of spectralswafdr.m / specswafdr.m
    drawn <- c(maxneg[1] > 0 && maxneg[1] != maxneg[2] && maxneg[2] == 0,
               maxneg[2] > 0,
               maxneg[3] > 0 && maxneg[3] != maxneg[2],
               maxneg[4] > 0 && maxneg[4] != maxneg[3],
               maxneg[5] > 0 && maxneg[5] != maxneg[4])
    fdr_cl <- matrix(NA_real_, nout, 5)
    for (k in 1:5){
      if (drawn[k] == TRUE && nout >= 2){
        fdr_cl[, k] <- swa * fac[k]
      }
    }
    colnames(fdr_cl) <- c("fdr10", "fdr5", "fdr1", "fdr01", "fdr001")
    tab <- data.frame(fdr_percent = c(10, 5, 1, 0.1, 0.01),
                      n_rejected = maxneg, factor = fac,
                      chi2_cl_percent = pct, drawn = drawn)
    list(fdr_cl = fdr_cl, table = tab, alpha = alphob, variance_ratio = factob)
  }

  # @title Chi-square quantile by bisection (getchi2.m, REDFIT)
  # @description Upper-tail chi-square value for probability alpha using the
  # Numerical Recipes incomplete gamma function (tolerance 1e-6 relative on
  # alpha, as in Acycle).
  # @keywords internal
  # @noRd
  swa_getchi2 <- function(dof = NULL, alpha = NULL){
    tol <- 1.0e-6
    itmax <- 10000
    lm <- 0.0
    rm <- 1000.0
    iter <- 0
    if (alpha > 0.5){
      eps <- (1.0 - alpha) * tol
    } else {
      eps <- alpha * tol
    }
    chi2 <- NA_real_
    while (TRUE){
      iter <- iter + 1
      if (iter > itmax){
        warning("getchi2: ITER > ITMAX")
        break
      }
      chi2 <- 0.5 * (lm + rm)
      ac <- 1.0 - swa_gammp(0.5 * dof, 0.5 * chi2)
      if (abs(ac - alpha) < eps){
        break
      }
      if (ac > alpha){
        lm <- chi2
      } else {
        rm <- chi2
      }
    }
    chi2
  }

  # @title Regularised lower incomplete gamma P(a, x) (Numerical Recipes)
  # @keywords internal
  # @noRd
  swa_gammp <- function(a = NULL, x = NULL){
    if (x < 0 || a <= 0){
      stop("Bad arguments in GAMMP")
    }
    if (x < a + 1){
      return(swa_gser(a, x))
    }
    1.0 - swa_gcf(a, x)
  }

  # @title Series representation of P(a, x) (gamser.m)
  # @keywords internal
  # @noRd
  swa_gser <- function(a = NULL, x = NULL){
    itmax <- 100
    eps <- 3.0e-7
    if (x <= 0){
      return(0.0)
    }
    gl <- swa_gammln(a)
    ap <- a
    s <- 1.0 / a
    del <- s
    for (n in seq_len(itmax)){
      ap <- ap + 1.0
      del <- del * x / ap
      s <- s + del
      if (abs(del) < abs(s) * eps){
        break
      }
    }
    s * exp(-x + a * log(x) - gl)
  }

  # @title Continued fraction for Q(a, x) (gammcf.m)
  # @keywords internal
  # @noRd
  swa_gcf <- function(a = NULL, x = NULL){
    itmax <- 100
    eps <- 3.0e-7
    fpmin <- 1.0e-30
    gln <- swa_gammln(a)
    b <- x + 1.0 - a
    cc <- 1.0 / fpmin
    d <- 1.0 / b
    h <- d
    for (i in seq_len(itmax)){
      an <- -i * (i - a)
      b <- b + 2.0
      d <- an * d + b
      if (abs(d) < fpmin){
        d <- fpmin
      }
      cc <- b + an / cc
      if (abs(cc) < fpmin){
        cc <- fpmin
      }
      d <- 1.0 / d
      del <- d * cc
      h <- h * del
      if (abs(del - 1.0) < eps){
        break
      }
    }
    exp(-x + a * log(x) - gln) * h
  }

  # @title log Gamma, Lanczos approximation (gammln.m)
  # @keywords internal
  # @noRd
  swa_gammln <- function(xx = NULL){
    cof <- c(76.18009172947146, -86.50532032941677, 24.01409824083091,
             -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5)
    stp <- 2.5066282746310005
    x <- xx
    y <- x
    tmp <- x + 5.5
    tmp <- (x + 0.5) * log(tmp) - tmp
    ser <- 1.000000000190015
    for (j in 1:6){
      y <- y + 1.0
      ser <- ser + cof[j] / y
    }
    tmp + log(stp * ser / x)
  }

  # @title Legacy SWA fit coded inside specswafdr.m
  # @description Reproduces, including its quirks, the SWA background fit of
  # the Acycle SWA window (specswafdr.m): the width search stops at the first
  # zero-reversal fit, window centres advance winwidth - 1 ordinates, the
  # quadratic tail uses the Davis (1973) normal equations with the shifted
  # frequency index of the original, and at least one Hanning smoothing pass
  # is always applied.
  # @keywords internal
  # @noRd
  swa_background_legacy <- function(freq = NULL, log_power = NULL,
                                     n_data = NULL, verbose = FALSE){
    xlogp <- as.numeric(log_power)
    freq <- as.numeric(freq)
    nout <- length(xlogp)
    pad <- nout + 200
    fitsp <- numeric(pad)
    fitsp1 <- numeric(pad)
    fitsp2 <- numeric(pad)
    winwidth <- 9
    rmsemin <- 1.0e+20
    revmin <- n_data
    winfinal <- NA
    midw <- numeric(0)
    ncycle <- 0
    finish <- 0
    while (revmin > 0 && winwidth < 99){
      winwidth <- winwidth + 2
      wincentre <- floor(winwidth / 2) + 1
      ncycle <- floor(nout / winwidth)
      if (ncycle < 3){
        stop(paste0("Legacy SWA fit failed: fewer than three windows of ",
                    "width ", winwidth, " (the MATLAB code errors here)"))
      }
      midw <- numeric(ncycle)
      k <- 0
      for (icycle in seq_len(ncycle)){
        sump <- 0.0
        for (iwin in seq_len(winwidth)){
          sump <- sump + xlogp[k + iwin]
        }
        midw[icycle] <- sump / winwidth
        k <- k + winwidth
      }
      diff1 <- (midw[2] - midw[1]) / winwidth
      diff2 <- (midw[3] - midw[2]) / winwidth
      delta <- diff2 - diff1
      diff0 <- diff1 - delta
      value1 <- midw[1]
      value2 <- midw[1]
      for (i in (wincentre - 1):1){
        value1 <- value1 - diff1
        fitsp1[i] <- value1
        value2 <- value2 - diff0
        fitsp2[i] <- value2
      }
      ii <- wincentre
      fitsp1[ii] <- midw[1]
      fitsp2[ii] <- midw[1]
      for (k in seq_len(ncycle - 1)){
        d0 <- (midw[k + 1] - midw[k]) / winwidth
        value <- midw[k]
        for (j in seq_len(winwidth - 1)){
          ii <- ii + 1
          value <- value + d0
          fitsp1[ii] <- value
          fitsp2[ii] <- value
        }
        fitsp1[ii] <- midw[k + 1]
        fitsp2[ii] <- midw[k + 1]
      }
      diff1 <- (midw[ncycle] - midw[ncycle - 1]) / winwidth
      diff2 <- (midw[ncycle - 1] - midw[ncycle - 2]) / winwidth
      delta <- diff2 - diff1
      diff0 <- diff1 - delta
      value1 <- midw[ncycle]
      value2 <- midw[ncycle]
      finish <- ii
      for (i in finish:nout){
        value1 <- value1 + diff1
        fitsp1[i] <- value1
        value2 <- value2 + diff0
        fitsp2[i] <- value2
      }
      rmsef1 <- sqrt(sum((xlogp - fitsp1[1:nout])^2) / (nout - 1))
      rmsef2 <- sqrt(sum((xlogp - fitsp2[1:nout])^2) / (nout - 1))
      lim <- ncycle * winwidth
      rev1 <- sum(fitsp1[1:(lim - 1)] <= fitsp1[2:lim])
      rev2 <- sum(fitsp2[1:(lim - 1)] <= fitsp2[2:lim])
      if (revmin > rev1){
        revmin <- rev1
        rmsemin <- rmsef1
        winfinal <- winwidth
        fitsp[1:nout] <- fitsp1[1:nout]
      }
      if (revmin > rev2){
        revmin <- rev2
        rmsemin <- rmsef2
        winfinal <- winwidth
        fitsp[1:nout] <- fitsp2[1:nout]
      }
      if (revmin == 0){
        if (rmsemin > rmsef1){
          rmsemin <- rmsef1
          winfinal <- winwidth
          fitsp[1:nout] <- fitsp1[1:nout]
        }
        if (rmsemin > rmsef2){
          rmsemin <- rmsef2
          winfinal <- winwidth
          fitsp[1:nout] <- fitsp2[1:nout]
        }
      }
    }
    rmsel1 <- rmsemin
    tail_model <- "extrapolated"
    # white noise end (uses the last tested window, as the original)
    fitsp3 <- numeric(pad)
    if (finish > 1){
      fitsp3[1:(finish - 1)] <- fitsp[1:(finish - 1)]
    }
    fitsp3[finish:nout] <- midw[ncycle]
    rmse <- sqrt(sum((xlogp - fitsp3[1:nout])^2) / (nout - 1))
    if (rmse < rmsel1){
      fitsp[1:nout] <- fitsp3[1:nout]
      rmsel1 <- rmse
      tail_model <- "constant"
    }
    # quadratic end fit (Davis 1973), original index handling
    tail_idx <- (finish - floor(winwidth / 2)):nout
    ntail <- length(tail_idx)
    xfr <- c(freq[tail_idx], 0, 0)
    yv <- c(xlogp[tail_idx], 0)
    nest <- ntail + 1
    bmat <- matrix(0, 3, 3)
    cvec <- numeric(3)
    for (i in seq_len(nest)){
      xp <- c(1, xfr[i + 1], xfr[i + 1] * xfr[i + 1])
      for (j in 1:3){
        for (k in 1:3){
          bmat[j, k] <- bmat[j, k] + xp[j] * xp[k]
        }
        cvec[j] <- cvec[j] + xp[j] * yv[i]
      }
    }
    cvec <- swa_solve(bmat, cvec, 1.0e-06)
    ii <- finish - 20
    if (ii + 1 < 1){
      stop("Legacy SWA quadratic tail: index below 1 (the MATLAB code errors here)")
    }
    for (i in seq_len(nest - 2)){
      xxp <- 1.0
      yyp <- 0.0
      for (j in 1:3){
        yyp <- yyp + xxp * cvec[j]
        xxp <- xxp * xfr[i + 1]
      }
      ii <- ii + 1
      if (ii > length(fitsp3)){
        fitsp3 <- c(fitsp3, numeric(ii - length(fitsp3)))
      }
      fitsp3[ii] <- yyp
    }
    fitsp3[nout] <- (fitsp3[nout] - fitsp3[nout - 1]) + fitsp3[nout]
    if (finish > 1){
      fitsp3[1:(finish - 1)] <- fitsp[1:(finish - 1)]
    }
    rmse <- sqrt(sum((xlogp - fitsp3[1:nout])^2) / (nout - 1))
    if (rmse < rmsel1){
      fitsp[1:nout] <- fitsp3[1:nout]
      rmsel1 <- rmse
      tail_model <- "quadratic"
    }
    # Hanning smoothing, legacy stopping rule
    f <- fitsp[1:nout]
    mid <- 2:(nout - 1)
    nsmth <- 1
    while (TRUE){
      ps <- numeric(nout)
      ps[1] <- 0.5 * f[1] + 0.5 * f[2]
      ps[nout] <- 0.5 * f[nout - 1] + 0.5 * f[nout]
      ps[mid] <- 0.25 * f[mid - 1] + 0.5 * f[mid] + 0.25 * f[mid + 1]
      f <- ps
      nsmth <- nsmth + 1
      rmset <- sqrt(sum((xlogp - f)^2) / nout)
      if (!(rmset < rmsel1 * 1.01 && nsmth < 500)){
        break
      }
    }
    if (verbose == TRUE){
      cat(sprintf("Legacy SWA: winfinal = %d, after smoothing RMSE = %.6f\n",
                  winfinal, rmset))
    }
    diagnostics <- list(selected_window = winfinal,
                        initial_nondecreasing_count = revmin,
                        initial_rmse = rmsemin, tail_model = tail_model,
                        n_smoothing = nsmth - 1, final_rmse = rmset,
                        n_data = n_data, n_frequency = nout)
    list(swa = 10^f, n_windows = winfinal, diagnostics = diagnostics)
  }

  # @title Gauss-Jordan solver of the legacy quadratic fit (solve.m)
  # @keywords internal
  # @noRd
  swa_solve <- function(a = NULL, b = NULL, zero = 1.0e-06){
    n <- length(b)
    for (i in seq_len(n)){
      div <- a[i, i]
      if (abs(div) > zero){
        a[i, ] <- a[i, ] / div
        b[i] <- b[i] / div
        for (j in seq_len(n)){
          if (i != j){
            ratio <- a[j, i]
            a[j, ] <- a[j, ] - ratio * a[i, ]
            b[j] <- b[j] - ratio * b[i]
          }
        }
      } else {
        warning("Matrix problem during polynomial fit to spectrum")
        return(b)
      }
    }
    b
  }

  # ---------------------------------------------------------------------------
  # Slepian (DPSS) tapers, MATLAB dpss() conventions (unit energy)
  # ---------------------------------------------------------------------------

  acy_dpss_cache <- new.env(parent = emptyenv())

  swa_fdr_cache <- new.env(parent = emptyenv())

  # ---- main code ----
  data <- acy_check_two_column(data)
  if (!(method %in% c("mtm", "lomb"))){
    stop("method must be 'mtm' or 'lomb'")
  }
  if (!(background %in% c("specswa", "legacy"))){
    stop("background must be 'specswa' or 'legacy'")
  }
  datax <- data[, 1]
  datay <- data[, 2]
  ndata <- length(datax)

  if (method == "mtm"){
    if (background == "legacy"){
      stop("background = 'legacy' is only available for method = 'lomb'")
    }
    if (!is.numeric(nw) || length(nw) != 1 || !is.finite(nw) || nw <= 0){
      stop("nw must be a finite positive number")
    }
    if (nw == 1){
      n_taper <- round(2 * nw)
    } else {
      n_taper <- round(2 * nw) - 1
    }
    if (n_taper < 1){
      stop("nw does not provide a usable Slepian taper")
    }
    dt <- stats::median(diff(datax))
    nzeropad <- max(ndata, round(pad_times * ndata))
    dof <- 2 * n_taper
    sp <- acy_pmtm(datay, nw = nw, nfft = nzeropad, drop_last = (nw != 1))
    pxx <- sp$pxx
    fn <- 1 / (2 * dt)
    freq <- sp$f / pi * fn
    xlogp <- log10(pxx)
    bg <- swa_background(freq, xlogp, ndata, verbose = verbose)
    swa <- bg$swa
    probs <- c(0.90, 0.95, 0.99, 0.999, 0.9999)
    chi2_factor <- stats::qchisq(probs, dof) / dof
    tdif <- max(datax) - min(datax)
    xfint <- 1 / tdif
    bw <- xfint * nw
    ncorr <- bw / xfint
    nout <- length(freq)
    power <- pxx
    extra <- list()
  } else {
    lomb <- swa_lomb_power(datax, datay, verbose = verbose)
    freq <- lomb$freq
    power <- lomb$power
    nout <- length(freq)
    if (lomb$rho1 == 0 || lomb$rho1 == 1){
      stop(paste0("Lag-1 autocorrelation = ", lomb$rho1, ": SWA is not ",
                  "appropriate for series with white noise or power law ",
                  "backgrounds"))
    }
    xlogp <- log10(power)
    if (background == "specswa"){
      bg <- swa_background(freq, xlogp, ndata, verbose = verbose)
    } else {
      bg <- swa_background_legacy(freq, xlogp, ndata, verbose = verbose)
    }
    swa <- bg$swa
    dof <- 8
    chi2_factor <- numeric(5)
    alphas <- c(0.1, 0.05, 0.01, 0.001, 0.0001)
    for (i in 1:5){
      chi2_factor[i] <- swa_getchi2(8, alphas[i]) / 8
    }
    bw <- lomb$bandwidth
    xfint <- lomb$rayleigh
    ncorr <- bw / xfint
    extra <- list(periodogram = lomb$periodogram, rho1 = lomb$rho1,
                  irregular = lomb$irregular)
  }
  names(chi2_factor) <- c("chi90", "chi95", "chi99", "chi999", "chi9999")
  chi2_cl <- outer(swa, chi2_factor)
  colnames(chi2_cl) <- names(chi2_factor)

  cn <- 0
  i <- 1
  while (i <= ncorr){
    cn <- cn + 1 / i
    i <- i + 1
  }
  fdr <- swa_fdr(power, swa, dof, cn, nout)

  out <- list(freq = freq, power = power)
  if (method == "lomb"){
    out$periodogram <- extra$periodogram
  }
  out$swa <- swa
  out$chi2_factor <- chi2_factor
  out$chi2_cl <- chi2_cl
  out$fdr_cl <- fdr$fdr_cl
  out$fdr <- fdr$table
  out$alpha <- fdr$alpha
  out$variance_ratio <- fdr$variance_ratio
  out$dof <- dof
  out$cn <- cn
  out$bandwidth <- bw
  out$rayleigh <- xfint
  out$n_windows <- bg$n_windows
  out$diagnostics <- bg$diagnostics
  if (method == "lomb"){
    out$rho1 <- extra$rho1
    out$irregular <- extra$irregular
  }
  if (verbose == TRUE){
    cat(sprintf("Max variance ratio = %.5f\n", max(fdr$variance_ratio)))
    print(fdr$table)
  }

  if (genplot == TRUE){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    pos <- freq > 0
    yl <- range(c(power[pos], swa[pos], chi2_cl[pos, ]), na.rm = TRUE)
    graphics::plot(freq[pos], power[pos], type = "l", log = "y", ylim = yl,
                   xlab = "Frequency", ylab = "Power",
                   main = paste0("SWA background (", method, ")"))
    graphics::lines(freq[pos], swa[pos], lwd = 2)
    graphics::lines(freq[pos], chi2_cl[pos, 2], col = "red", lwd = 1.5)
    graphics::lines(freq[pos], chi2_cl[pos, 3], col = "blue")
    graphics::lines(freq[pos], chi2_cl[pos, 4], col = "magenta")
    leg <- c("Power", "SWA", "Chi2 95%", "Chi2 99%", "Chi2 99.9%")
    cols <- c("black", "black", "red", "blue", "magenta")
    lty <- c(1, 1, 1, 1, 1)
    lwd <- c(1, 2, 1.5, 1, 1)
    fdr_cols <- c("grey40", "red", "blue", "darkgreen", "black")
    fdr_names <- c("10% FDR", "5% FDR", "1% FDR", "0.1% FDR", "0.01% FDR")
    for (k in 1:5){
      if (!is.na(fdr$fdr_cl[1, k])){
        graphics::lines(freq[pos], fdr$fdr_cl[pos, k], col = fdr_cols[k],
                        lty = 2, lwd = 1.5)
        leg <- c(leg, fdr_names[k])
        cols <- c(cols, fdr_cols[k])
        lty <- c(lty, 2)
        lwd <- c(lwd, 1.5)
      }
    }
    graphics::legend("topright", legend = leg, col = cols, lty = lty,
                     lwd = lwd, bty = "n", cex = 0.8)
  }
  out
}
