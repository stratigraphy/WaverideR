#' @title Bayesian probability spectrum
#'
#' @description Calculates the Bayesian probability spectrum of Gregory
#' (2005) as used by Weedon et al. (2019) and Weedon (2020): the data are
#' pre-whitened with the lag-1 autocorrelation coefficient, linearly
#' detrended, tapered with a 5\% split cosine taper and standardised; the
#' variance-normalised Lomb-Scargle periodogram (routine PERIOD of Press et
#' al., 1992) is then converted into the posterior probability of a single
#' sinusoid at each frequency,
#' \eqn{p(f) \propto (1 - 2 P(f) / (N \overline{d^2}))^{(2 - N)/2}},
#' normalised so that the probabilities sum to one.
#' Translated from Acycle v3.1 specbayes.m (Mingsong Li, after the FORTRAN
#' program of Graham P. Weedon).
#'
#' @param data Record to analyse (matrix or data frame with two columns: depth or time in column 1, proxy value
#' in column 2; further columns are ignored). Rows are sorted on the first column and rows with non-finite values
#' are removed; uneven spacing is allowed (Lomb-Scargle periodogram). Required, no default (\code{Default=NULL}
#' gives an error).
#' @param genplot Plot the periodogram and the Bayesian probability against frequency (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param verbose Print the maximum probability and its frequency and period (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li,
#' translated from the FORTRAN program of Graham P. Weedon.
#'
#' @references
#' Gregory, P.C., 2005. Bayesian Logical Data Analysis for the Physical
#' Sciences. Cambridge University Press.
#'
#' Press, W.H., Teukolsky, S.A., Vetterling, W.T., Flannery, B.P., 1992.
#' Numerical Recipes in FORTRAN, 2nd ed. Cambridge University Press.
#'
#' Weedon, G.P., Page, K.N., Jenkyns, H.C., 2019. Cyclostratigraphy,
#' stratigraphic gaps and the duration of the Hettangian (Jurassic):
#' insights from the Blue Lias Formation of Southern Britain. Geological
#' Magazine 156, 1469-1509.
#'
#' Weedon, G.P., 2020. Confirmed detection of Palaeogene and Jurassic
#' orbitally-forced sedimentary cycles in the depth domain using False
#' Discovery Rates and Bayesian probability spectra. Boletin Geologico y
#' Minero 131, 207-230.
#'
#' @examples
#' \donttest{
#' set.seed(2)
#' n <- 400
#' t <- 1:n
#' noise <- as.numeric(stats::filter(stats::rnorm(n), 0.5,
#'   method = "recursive"))
#' x <- noise + 1.5 * sin(2 * pi * t / 23)
#' res <- bayes_spectrum(data = cbind(t, x), genplot = TRUE)
#' 1 / res$freq_max
#' }
#'
#' @return A list with elements \code{freq} (frequencies), \code{period}
#' (1/freq), \code{power} (variance-normalised Lomb-Scargle periodogram of
#' the pre-whitened, detrended, tapered and standardised data),
#' \code{bayes_prob} (normalised Bayesian probability), \code{max_prob},
#' \code{freq_max} (frequency of maximum probability) and \code{rho1}
#' (lag-1 autocorrelation used for pre-whitening).
#'
#' @export
#' @importFrom graphics plot lines legend par
bayes_spectrum <- function(data = NULL, genplot = FALSE, verbose = FALSE){
  # Internal helper functions (local to this function)

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

  # @title Mean and corrected two-pass variance (meanvar in specbayes.m)
  # @keywords internal
  # @noRd
  bayes_spec_meanvar <- function(x = NULL){
    an <- length(x)
    avex <- sum(x) / an
    d <- x - avex
    ep <- sum(d)
    varx <- (sum(d * d) - (ep * ep) / an) / (an - 1.0)
    list(avex = avex, varx = varx)
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

  # Main code
  data <- acy_check_two_column(data)
  time <- data[, 1]
  x1 <- data[, 2]
  ndata <- length(x1)
  andata <- ndata
  if (ndata < 4){
    stop("data must contain at least 4 rows")
  }
  sp <- swa_spacing(time)
  tdif <- sp$tdif

  # pre-whitening with the lag-1 autocorrelation
  mv <- bayes_spec_meanvar(x1)
  rho1 <- swa_rho(x1, mv$varx)
  if (rho1 > 0.0){
    x1[1:(ndata - 1)] <- x1[1:(ndata - 1)] - rho1 * x1[2:ndata]
  }
  # linear detrending and 5 % cosine taper
  x1 <- swa_detrend(time, x1)
  x1 <- swa_ctaper(x1, time, tdif)

  xfint <- 1.0 / tdif
  if (sp$tintmax == sp$tintmin){
    tint <- time[1] - time[2]
  } else {
    tint <- (time[1] - time[ndata]) / (andata - 1.0)
  }
  if (tint < 0.0){
    tint <- -1 * tint
  }
  xnyq <- 1.0 / (2.0 * tint)
  nout <- floor(xnyq / xfint)

  # standardise
  mv <- bayes_spec_meanvar(x1)
  x1 <- (x1 - mv$avex) / sqrt(mv$varx)
  mv <- bayes_spec_meanvar(x1)

  per <- swa_period(x1, time, mv$avex, nout, xfint, mv$varx)
  pow1 <- per$pow
  freq <- per$freq

  aved2 <- sum(x1 * x1) / andata
  bayr <- (2 - andata) / 2
  base <- 1 - ((2 * pow1) / (andata * aved2))
  if (any(base <= 0)){
    stop("Negative probability: periodogram exceeds N * mean(d^2) / 2")
  }
  bayesprob <- base^bayr
  sumbayes <- sum(bayesprob)
  if (any(!is.finite(bayesprob)) || !is.finite(sumbayes)){
    # overflow in the original: evaluate in log space (same normalised result)
    lb <- bayr * log(base)
    bayesprob <- exp(lb - max(lb))
    sumbayes <- sum(bayesprob)
  }
  bayesprob <- bayesprob / sumbayes
  maxib <- which.max(bayesprob)
  if (verbose == TRUE){
    cat(sprintf("Max. probability = %.3f\nFrequency = %.7f\nPeriod = %.7f\n",
        bayesprob[maxib], freq[maxib], 1 / freq[maxib]))
  }

  if (genplot == TRUE){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    graphics::plot(freq, pow1, type = "l", log = "y",
        ylim = c(min(bayesprob[bayesprob > 0]), max(pow1) * 1.5),
        xlab = "Frequency (cycles/unit)", ylab = "Power",
        main = "Bayesian probability spectrum")
    graphics::lines(freq, bayesprob, col = "red")
    graphics::legend("bottomleft", legend = c("Power", "Bayesian probability"),
        col = c("black", "red"), lty = 1, bty = "n")
  }
  list(freq = freq, period = 1 / freq, power = pow1, bayes_prob = bayesprob,
      max_prob = bayesprob[maxib], freq_max = freq[maxib], rho1 = rho1)
}
