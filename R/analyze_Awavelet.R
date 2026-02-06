#' @title Conduct the adaptive continuous wavelet transform on a time series or signal
#'
#' @description
#' Compute the continuous wavelet transform (CWT) using a complex Morlet
#' wavelet with scale dependent wavelet width. The number of oscillatory
#' cycles in the Morlet wavelet is allowed to vary smoothly with scale,
#' enabling adaptive time frequency resolution similar in spirit to
#' superlet based approaches.
#'
#' @param data Input data, should be a matrix or data frame in which
#' the first column is depth or time and the second column is the proxy record.
#'
#' @param dj Spacing between successive scales. Scales increase by powers
#' of two as \code{2^(j * dj)}. Smaller values of \code{dj} increase
#' frequency resolution at the expense of computational time.
#' Default is \code{1/100}.
#'
#' @param lowerPeriod Lowest period to be analysed. Defines the smallest
#' scale of the transform.
#'
#' @param upperPeriod Highest period to be analysed. Defines the largest
#' scale of the transform.
#'
#' @param verbose Logical. If TRUE, print interpolation diagnostics.
#'
#' @param omega_min Minimum number of oscillatory cycles of the Morlet
#' wavelet.
#'
#' @param omega_max Maximum number of oscillatory cycles of the Morlet
#' wavelet.
#'
#' @param scaling Character string defining how the number of wavelet
#' cycles varies with scale. One of
#' \code{"log2"}, \code{"linear"}, \code{"sqrt"}, \code{"quadratic"},
#' or \code{"power"}.
#'
#' @param alpha Numeric. Exponent used when \code{scaling = "power"}.
#' Values greater than one emphasize high frequency sharpening, whereas
#' values smaller than one emphasize low frequency sharpening.
#'
#' @param n_simulations Number of Monte Carlo simulations. Currently
#' included for interface compatibility.
#'
#' @param run_multicore Logical. Currently included for interface
#' compatibility.
#'
#' @return
#' A list with class \code{"analyze.Awavelet"} containing:
#'
#' \itemize{
#'   \item Wave: complex wavelet coefficients
#'   \item Phase: instantaneous phase
#'   \item Ampl: wavelet amplitude
#'   \item Power: wavelet power spectrum
#'   \item dt: sampling interval
#'   \item dj: scale spacing
#'   \item Power.avg: average power per period
#'   \item Period: physical periods
#'   \item Scale: wavelet scales
#'   \item coi.1: cone of influence x coordinates
#'   \item coi.2: cone of influence y coordinates
#'   \item nc: number of samples
#'   \item nr: number of scales
#'   \item axis.1: x axis values
#'   \item axis.2: log2 scaled periods
#'   \item omega_min: minimum wavelet cycles
#'   \item omega_max: maximum wavelet cycles
#'   \item omega: scale dependent wavelet cycles
#'   \item scaling: scaling method used
#'   \item alpha: power law exponent
#'   \item x: interpolated x values
#'   \item y: interpolated signal values
#' }
#'
#' @author
#' Adapted from the WaveletComp and biwavelet R packages, which are based on
#' the MATLAB wavelet code by Torrence and Compo, with extensions for
#' scale dependent wavelet width.
#'
#' @references
#' Torrence, C., and G. P. Compo (1998). A practical guide to wavelet analysis.
#' Bulletin of the American Meteorological Society, 79, 61–78.
#'
#' Morlet, J., G. Arens, E. Fourgeau, and D. Giard (1982).
#' Wave propagation and sampling theory. Geophysics, 47, 203–236.
#'
#' @examples
#' \donttest{
#'#Example 1. Using the Total Solar Irradiance data set of Steinhilber et al., (2012)
#'TSI_wt <-
#'  analyze_Awavelet(
#'    data = TSI,
#'    dj = 1/200,
#'    lowerPeriod = 16,
#'    upperPeriod = 8192,
#'    verbose = FALSE,
#' omega_min = 6,
#' omega_max = 12,
#' scaling = "log2",
#' alpha = 1
#'  )
#'
#'
#'#Example 2. Using the magnetic susceptibility data set of Pas et al., (2018)
#'mag_wt <-
#'analyze_Awavelet(
#'data = mag,
#'dj = 1/100,
#'lowerPeriod = 0.1,
#'upperPeriod = 254,
#'verbose = FALSE,
#' omega_min = 6,
#' omega_max = 12,
#' scaling = "log2",
#' alpha = 1
#')
#'
#'#Example 3. Using the greyscale data set of Zeeden et al., (2013)
#'grey_wt <-
#'  analyze_Awavelet(
#'    data = grey,
#'    dj = 1/200,
#'    lowerPeriod = 0.02,
#'    upperPeriod = 256,
#'    verbose = FALSE,
#' omega_min = 6,
#' omega_max = 12,
#' scaling = "log2",
#' alpha = 1,
#'  )
#'
#'}
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom stats acf

analyze_Awavelet <-
  function(data = NULL,
           dj = 1 / 100,
           lowerPeriod = 2,
           upperPeriod = 1024,
           verbose = FALSE,
           omega_min = 6,
           omega_max = 12,
           scaling = "log2",
           alpha = 1,
           n_simulations = 10,
           run_multicore = FALSE) {

    scaling <- match.arg(scaling)

    ## -------------------------------------------------------------------------
    ## 1. Data preparation
    ## -------------------------------------------------------------------------

    dat <- as.matrix(data)
    dat <- na.omit(dat)
    dat <- dat[order(dat[, 1], na.last = NA, decreasing = FALSE), ]

    dx <- diff(dat[, 1])
    dt <- median(dx)

    xout <- seq(dat[1, 1], dat[nrow(dat), 1], by = dt)
    interp <- approx(dat[, 1], dat[, 2], xout, method = "linear")
    dat <- as.data.frame(interp)

    if (verbose) {
      cat("dataset interpolated to:", round(dt, 10), "\n")
    }

    x_axis <- dat[, 1]
    x <- dat[, 2]

    series.length <- length(x)
    pot2 <- trunc(log2(series.length) + 0.5)
    pad.length <- 2^(pot2 + 1) - series.length

    ## -------------------------------------------------------------------------
    ## 2. Scale and period definition
    ## -------------------------------------------------------------------------

    omega_ref <- mean(c(omega_min, omega_max))
    fourier.factor.ref <- (2 * pi) / omega_ref

    min.scale <- lowerPeriod / fourier.factor.ref
    max.scale <- upperPeriod / fourier.factor.ref

    J <- as.integer(log2(max.scale / min.scale) / dj)

    scales <- min.scale * 2^((0:J) * dj)
    scales.length <- length(scales)

    periods <- fourier.factor.ref * scales

    ## -------------------------------------------------------------------------
    ## 3. Scale-dependent omega definition
    ## -------------------------------------------------------------------------
    # chose from  "c("log2", "linear", "sqrt", "quadratic", "power")"

    if (scaling == "linear") {
      u <- (periods - min(periods)) /
        (max(periods) - min(periods))
    }

    if(scaling == "log2") {
      u <- (log2(periods) - min(log2(periods))) /
        (max(log2(periods)) - min(log2(periods)))
    }

    if (scaling == "sqrt") {
      u <- sqrt(u)
    }

    if (scaling == "quadratic") {
      u <- u^2
    }

    if (scaling == "power") {
      u <- u^alpha
    }

    omega_vec <- omega_min + u * (omega_max - omega_min)
    fourier.factor.vec <- (2 * pi) / omega_vec

    ## -------------------------------------------------------------------------
    ## 4. FFT preparation
    ## -------------------------------------------------------------------------

    N <- series.length + pad.length
    omega.k <- (1:floor(N / 2)) * (2 * pi) / (N * dt)
    omega.k <- c(0, omega.k, -omega.k[floor((N - 1) / 2):1])

    x_standard <- (x - mean(x)) / sd(x)
    xpad <- c(x_standard, rep(0, pad.length))
    fft.xpad <- fft(xpad)

    wave <- matrix(0, nrow = scales.length, ncol = N)
    wave <- wave + 1i * wave

    ## -------------------------------------------------------------------------
    ## 5. Continuous wavelet transform
    ## -------------------------------------------------------------------------

    for (i in seq_len(scales.length)) {

      my.scale <- scales[i]
      omega0 <- omega_vec[i]

      norm.factor <- pi^(1 / 4) * sqrt(2 * my.scale / dt)
      expnt <- -((my.scale * omega.k - omega0)^2 / 2) * (omega.k > 0)
      daughter <- norm.factor * exp(expnt) * (omega.k > 0)

      wave[i, ] <- fft(fft.xpad * daughter, inverse = TRUE) / N
    }

    Wave <- wave[, 1:series.length]

    ## -------------------------------------------------------------------------
    ## 6. Diagnostics
    ## -------------------------------------------------------------------------

    Power <- Mod(Wave)^2 / matrix(rep(scales, series.length),
                                  nrow = scales.length)

    Phase <- Arg(Wave)

    Ampl <- Mod(Wave) / matrix(rep(sqrt(scales), series.length),
                               nrow = scales.length)

    Power.avg <- rowMeans(Power)

    axis.1 <- x_axis
    axis.2 <- log2(periods)

    coi <- fourier.factor.vec * sqrt(2) * dt *
      c(1e-5,
        1:((series.length + 1) / 2 - 1),
        rev(1:(series.length / 2 - 1)),
        1e-5)

    coi.x <- c(axis.1[1] - dt * 0.5,
               axis.1,
               axis.1[series.length] + dt * 0.5)

    logyint <- axis.2[2] - axis.2[1]
    yl <- c(log2(periods[scales.length]) + 0.5 * logyint,
            log2(periods[1]) - 0.5 * logyint)
    yr <- rev(yl)
    coi.y <- c(yl, log2(coi), yr)

    ## -------------------------------------------------------------------------
    ## 7. Output
    ## -------------------------------------------------------------------------

    output <- list(
      Wave = Wave,
      Phase = Phase,
      Ampl = Ampl,
      Power = Power,
      dt = dt,
      dj = dj,
      Power.avg = Power.avg,
      Period = periods,
      Scale = scales,
      coi.1 = coi.x,
      coi.2 = coi.y,
      nc = series.length,
      nr = scales.length,
      axis.1 = axis.1,
      axis.2 = axis.2,
      omega_min = omega_min,
      omega_max = omega_max,
      omega = omega_vec,
      scaling = scaling,
      alpha = alpha,
      x = dat[, 1],
      y = dat[, 2]
    )

    class(output) <- "analyze.Awavelet"
    return(invisible(output))
  }
