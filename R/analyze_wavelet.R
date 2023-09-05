#' @title Conduct the continuous wavelet transform on a time series/signal
#'
#' @description
#' Compute the continuous wavelet transform (CWT) using a Morlet wavelet
#'
#' @param data Input data, should be a matrix or data frame in which
#' the first column is depth or time and the second column is proxy record.
#' @param dj Spacing between successive scales. The CWT analyses analyses the signal using successive periods
#' which increase by the power of 2 (e.g.2^0=1,2^1=2,2^2=4,2^3=8,2^4=16). To have more resolution
#' in-between these steps the dj parameter exists, the dj parameter specifies how many extra steps/spacing in-between
#' the power of 2 scaled CWT is added. The amount of steps is 1/x with a higher x indicating a smaller spacing.
#' Increasing the increases the computational time of the CWT \code{Default=1/200}.
#' @param lowerPeriod  Lowest period to be analyzed \code{Default=2}.
#' The CWT analyses the signal starting from the lowerPeriod to the upperPeriod so the proper selection these
#' parameters allows to analyze the signal for a specific range of cycles.
#' scaling is done using power 2 so for the best plotting results select a value to the power or 2.
#' @param upperPeriod Upper period to be analyzed \code{Default=1024}.
#' The CWT analyses the signal starting from the lowerPeriod to the upperPeriod so the proper selection these
#' parameters allows to analyze the signal for a specific range of cycles.
#'  scaling is done using power 2 so for the best plotting results select a value to the power or 2.
#' @param verbose Print text \code{Default=FALSE}.
#' @param omega_nr Number of cycles contained within the Morlet wavelet
#'
#' @return
#' The output is a list (wavelet object) which contain 18 objects which are the result of the continuous wavelet transform (CWT).
#'Object 1: Wave - Wave values of the wavelet
#'Object 2: Phase - Phase of the wavelet
#'Object 3: Ampl - Amplitude values of the wavelet
#'Object 4: Power - Power values of the wavelet
#'Object 5: dt - Step size
#'Object 6: dj - Scale size
#'Object 7: Power.avg  - Average power values
#'Object 8: Period - Period values
#'Object 9: Scale - Scale value
#'Object 10: coi.1 - Cone of influence values 1
#'Object 11: coi.2 - Cone of influence values 2
#'Object 12: nc - Number of columns
#'Object 13: nr - Number of rows
#'Object 14: axis.1 - axis values 1
#'Object 15: axis.2 - axis values 2
#'Object 16: omega_nr - Number of cycles in the wavelet
#'Object 17: x - x values of the data set
#'Object 18: y - y values of the data set
#'
#' @author
#' Code based on on the \link[WaveletComp]{analyze.wavelet} function of the 'WaveletComp' R package
#' and \link[biwavelet]{wt} function of the 'biwavelet' R package which are based on the
#' wavelet MATLAB code written by Christopher Torrence and Gibert P. Compo.
#'
#' @references
#'Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
#'Wavelet Analysis. R package version 1.1.
#'\url{https://CRAN.R-project.org/package=WaveletComp}
#'
#'Gouhier TC, Grinsted A, Simko V (2021). R package biwavelet: Conduct Univariate and Bivariate Wavelet Analyses. (Version 0.20.21),
#'\url{https://github.com/tgouhier/biwavelet}
#'
#'Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#'Bulletin of the American Meteorological Society 79:61-78.
#'\url{https://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf}
#'
#'Morlet, Jean, Georges Arens, Eliane Fourgeau, and Dominique Glard.
#'"Wave propagation and sampling theory—Part I: Complex signal and scattering in multilayered media.
#'" Geophysics 47, no. 2 (1982): 203-221.
#' \url{https://pubs.geoscienceworld.org/geophysics/article/47/2/203/68601/Wave-propagation-and-sampling-theory-Part-I}
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#'  Geophysics 1982 47 (2): 222–236. \url{https://pubs.geoscienceworld.org/geophysics/article/47/2/222/68604/Wave-propagation-and-sampling-theory-Part-II}
#'
#' @examples
#' \donttest{
#'#Example 1. Using the Total Solar Irradiance data set of Steinhilver et al., (2012)
#'TSI_wt <-
#'  analyze_wavelet(
#'    data = TSI,
#'    dj = 1/200,
#'    lowerPeriod = 16,
#'    upperPeriod = 8192,
#'    verbose = FALSE,
#'    omega_nr = 6
#'  )
#'
#'
#'#Example 2. Using the magnetic susceptibility data set of Pas et al., (2018)
#'mag_wt <-
#'analyze_wavelet(
#'data = mag,
#'dj = 1/100,
#'lowerPeriod = 0.1,
#'upperPeriod = 254,
#'verbose = FALSE,
#'omega_nr = 10
#')
#'
#'#Example 3. Using the greyscale data set of Zeeden et al., (2013)
#'grey_wt <-
#'  analyze_wavelet(
#'    data = grey,
#'    dj = 1/200,
#'    lowerPeriod = 0.02,
#'    upperPeriod = 256,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'}
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom WaveletComp analyze.wavelet
#' @importFrom biwavelet wt



analyze_wavelet <-
  function(data = NULL,
           dj = 1 / 20,
           lowerPeriod = 2,
           upperPeriod = 1024,
           verbose = FALSE,
           omega_nr = 6) {
    dat <- as.matrix(data)
    dat <- na.omit(dat)
    dat <- dat[order(dat[, 1], na.last = NA, decreasing = F),]
    npts <- length(dat[, 1])
    start <- dat[1, 1]
    end <- dat[length(dat[, 1]), 1]
    x1 <- dat[1:(npts - 1), 1]
    x2 <- dat[2:(npts), 1]
    dx = x2 - x1
    dt = median(dx)
    xout <- seq(start, end, by = dt)
    npts <- length(xout)
    interp <- approx(dat[, 1], dat[, 2], xout, method = "linear",
                     n = npts)
    dat <- as.data.frame(interp)

    if (verbose == TRUE) {
      cat("dataset interpolated to: ", round(dt, 10))
    }

    x_axis <- dat[, 1]
    x <- dat[, 2]

    # Original length and length of zero padding:
    series.length = length(x)
    pot2 = trunc(log2(series.length) + 0.5)
    pad.length = 2 ^ (pot2 + 1) - series.length

    # Define central angular frequency omega0 and fourier factor:
    omega0 = omega_nr
    #   fourier.factor   = (4*pi)/(omega0 + sqrt(2+omega0^2))
    fourier.factor = (2 * pi) / omega0

    # Compute scales and periods:
    min.scale = lowerPeriod / fourier.factor             # Convert lowerPeriod to minimum scale
    max.scale = upperPeriod / fourier.factor             # Convert upperPeriod to maximum scale
    J = as.integer(log2(max.scale / min.scale) / dj)    # Index of maximum scale -1

    scales = min.scale * 2 ^ ((0:J) * dj)        # sequence of scales
    scales.length = length(scales)           # J + 1
    periods = fourier.factor * scales          # sequence of periods

    # Computation of the angular frequencies
    N = series.length + pad.length
    omega.k = 1:floor(N / 2)
    omega.k = omega.k * (2 * pi) / (N * dt)                    # k <= N/2
    omega.k = c(0, omega.k, -omega.k[floor((N - 1) / 2):1])

    ###############################################################################
    ## Define the Morlet wavelet transform function
    ###############################################################################

    x_standard = (x - mean(x)) / sd(x)
    xpad = c(x_standard, rep(0, pad.length))

    # Compute Fast Fourier Transform of xpad
    fft.xpad = fft(xpad)

    # Compute wavelet transform of x
    # Prepare a complex matrix which accomodates the wavelet transform
    wave = matrix(0, nrow = scales.length, ncol = N)
    wave = wave + 1i * wave

    # Computation for each scale...
    # ... simultaneously for all time instances
    for (ind.scale in (1:scales.length)) {
      my.scale = scales[ind.scale]

      norm.factor = pi ^ (1 / 4) * sqrt(2 * my.scale / dt)
      expnt       = -((my.scale * omega.k - omega0) ^ 2 / 2) * (omega.k > 0)
      daughter    = norm.factor * exp(expnt)
      daughter    = daughter * (omega.k > 0)

      wave[ind.scale,] = fft(fft.xpad * daughter, inverse = TRUE) / N
    }

    # Cut out the wavelet transform
    Wave = wave[, 1:series.length]

    # Compute wavelet power
    Power = Mod(Wave) ^ 2 / matrix(rep(scales, series.length), nrow = scales.length)

    # Phase
    Phase = Arg(Wave)

    # Amplitude
    Ampl  = Mod(Wave) / matrix(rep(sqrt(scales), series.length), nrow =
                                 scales.length)



    Power.avg = rowMeans(as.matrix(Power))

    axis.1 = x_axis
    axis.2 = log2(periods)
    fourier.factor = (2 * pi) / omega0
    coi = fourier.factor * sqrt(2) * dt * c(1E-5, 1:((series.length + 1) /
                                                       2 - 1), rev((1:(
                                                         series.length / 2 - 1
                                                       ))), 1E-5)
    coi.x = c(axis.1[c(1, 1)] - dt * 0.5, axis.1, axis.1[c(series.length, series.length)] + dt *
                0.5)
    logyint = axis.2[2] - axis.2[1]
    yl = c(log2(periods[scales.length]) + 0.5 * logyint, log2(periods[1]) - 0.5 *
             logyint)
    yr = rev(yl)
    coi.y = c(yl, log2(coi), yr)
    x <- dat[, 1]
    y <- dat[, 2]

    output <-
      list(
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
        omega_nr = omega_nr,
        x = x,
        y = y
      )


    class(output) = "analyze.wavelet"

    return(invisible(output))

  }
