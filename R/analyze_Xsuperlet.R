#' @title Conduct the superlet transform on a time series or signal
#'
#' @description
#' Compute the cross superlet transform using a complex Morlet wavelet based on
#' the approach of Moca et al. (2021) and  the "analyze.coherency"  function of the WaveletComp
#' package. The superlet transform increases time frequency resolution by
#' combining multiple wavelet orders at each frequency.
#'
#' @param data_1 First input data set as a matrix or data frame. The first column must
#' represent depth or time, and the second column the signal or proxy record.
#'
#' @param data_2 Second  input data as a matrix or data frame. The first column must
#' represent depth or time, and the second column the signal or proxy record.
#'
#' @param upperPeriod Maximum period to be analysed. Controls the lowest
#' analysed frequency.
#'
#' @param lowerPeriod Minimum period to be analysed. Controls the highest
#' analysed frequency.
#'
#' @param Nf Number of frequencies used to construct the superlet spectrum.
#' Frequencies are logarithmically spaced between the limits defined by
#' \code{upperPeriod} and \code{lowerPeriod}.
#'
#' @param c1 Base number of cycles of the Morlet wavelet. Acts as the
#' fundamental wavelet width.
#'
#' @param o numeric vector of length two defining the minimum and
#' maximum superlet order. If NULL, all frequencies are analysed using
#' order one.
#'
#' @param mult Logical. If TRUE, the number of cycles increases
#' multiplicatively with superlet order. If FALSE, cycles increase
#' additively.
#'
#' @param verbose Logical. If TRUE, print progress and interpolation
#' information.
#'
#' @return
#'A list with class \code{"analyze.Xsuperlet"} containing
#'Wave: complex wavelet coefficients
#'Power: time frequency power spectrum
#'dt: sampling interval after interpolation
#'Phase: instantaneous phase of the signal
#'dj: number of frequencies
#'Power.avg: average spectral power
#'Period: physical periods corresponding to frequencies
#'nc: number of columns in the input data
#'nr: number of frequency levels
#'axis.1: x axis values (time or depth)
#'axis.2: y axis values (log2 scaled periods)
#'c1: base number of wavelet cycles
#'o: numeric vector of length two defining the minimum and
#' maximum superlet order. If NULL, all frequencies are analysed using
#' order one.
#'x1: interpolated x values first data set
#'y1: interpolated signal values first data set
#'x2: interpolated x values second data set
#'y2: interpolated signal values second data set
#'
#' @author
#' The "analyze_Xsuperlet" that generates the input for the plotting function
#' is based on the matlab code in Moca et al. (2021) and the the "analyze.coherency" function of the  'WaveletComp' R package
#'
#' @references
#' Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., and Mureșan, R. C. (2021).
#' Time frequency super resolution with superlets.
#' Nature Communications, 12, 337.
#' \url{https://doi.org/10.1038/s41467-020-20539-9}
#' @examples
#'#Example 1. A cross superlet of two etp solutions with noise overprint
#'etp_1 <- etp(
#'  tmin = 0,
#'  tmax = 1500,
#'  dt = 1,
#'  eWt = 1.5,
#'  oWt = 0.75,
#'  pWt = 1,
#'  esinw = TRUE,
#'  standardize = TRUE,
#'  genplot = FALSE,
#'  verbose = FALSE
#')
#'
#'etp_2 <- etp(
#'  tmin = 0,
#'  tmax = 1500,
#'  dt = 1,
#'  eWt = 1,
#'  oWt = 0.5,
#'  pWt = 1.5,
#'  esinw = TRUE,
#'  standardize = TRUE,
#'  genplot = FALSE,
#'  verbose = FALSE
#')
#'
#'etp_1[, 2] <- etp_1[, 2] + colorednoise::colored_noise(
#'  nrow(etp_1),
#'  sd = sd(etp_1[, 2]) / 1.5,
#'  mean = mean(etp_1[, 2]),
#'  phi = 0.9
#')
#'etp_2[, 2] <- etp_2[, 2] + colorednoise::colored_noise(
#'  nrow(etp_2),
#'  sd = sd(etp_2[, 2]) / 1.5,
#'  mean = mean(etp_2[, 2]),
#'  phi = 0.9
#')
#'
#'Xetp <- analyze_Xsuperlet(
#'  data_1 = etp_1,
#'  data_2  = etp_2,
#'  upperPeriod = 1024,
#'  lowerPeriod = 2,
#'  Nf = 128,
#'  c1 = 3,
#'  o = c(1, 10),
#'  mult = TRUE,
#'  verbose = FALSE
#')
#'
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom stats fft


analyze_Xsuperlet <- function(data_1 = NULL,
                              data_2  = NULL,
                              upperPeriod = 2,
                              lowerPeriod = 1024,
                              Nf = 128,
                              c1 = 3,
                              o = c(1, 10),
                              mult = TRUE,
                              verbose = FALSE) {
  data_1 <- data_1[order(data_1[, 1]), ]
  data_2 <- data_2[order(data_2[, 1]), ]

  # Determine overlapping interval
  xmin <- max(min(data_1[, 1], na.rm = TRUE), min(data_2[, 1], na.rm = TRUE))

  xmax <- min(max(data_1[, 1], na.rm = TRUE), max(data_2[, 1], na.rm = TRUE))


  dx1 <- median(diff(data_1[, 1]))
  dx2 <- median(diff(data_2[, 1]))

  dx <- max(dx1, dx2)  # conservative choice

  common_grid <- seq(from = xmin, to = xmax, by = dx)

  data_1_common <- data.frame(x = common_grid,
                              proxy = approx(data_1[, 1], data_1[, 2], xout = common_grid, rule = 1)$y)

  data_2_common <- data.frame(x = common_grid,
                              proxy = approx(data_2[, 1], data_2[, 2], xout = common_grid, rule = 1)$y)

  superlet_1 <- analyze_superlet(
    data = data_1_common,
    upperPeriod = upperPeriod,
    lowerPeriod = lowerPeriod,
    Nf = Nf,
    c1 = c1,
    o = o,
    mult = mult,
    verbose = verbose
  )


  superlet_2 <- analyze_superlet(
    data = data_2_common,
    upperPeriod = upperPeriod,
    lowerPeriod = lowerPeriod,
    Nf = Nf,
    c1 = c1,
    o = o,
    mult = mult,
    verbose = verbose
  )

  #  CROSS-ANALYSIS CALCULATIONS
  # Extract shared parameters from the Superlet output
  dt    <- superlet_1$dt
  Scale <- superlet_1$Scale
  nc    <- superlet_1$nc
  nr    <- superlet_1$nr

  # Raw Cross-Wavelet: Multiply complex coefficients and normalize by Scale
  # This mirrors the logic: Wave.xy = (W1 * Conj(W2)) / Scale
  Wave.xy  <- (superlet_1$Wave * Conj(superlet_2$Wave)) / matrix(rep(Scale, nc), nrow = nr)
  Power.xy <- Mod(Wave.xy)
  Phase.xy <- Arg(Wave.xy)

  output <- list(
    Wave = t(Wave.xy),
    Power = t(Power.xy),
    Phase = t(Phase.xy),
    dt = dt,
    dj = 1 / Nf,
    Power.avg = colMeans(t(Power.xy), na.rm = TRUE),
    Scale = Scale,
    Period = superlet_1$Period,
    nc = nc,
    nr = Nf,
    axis.1 = superlet_1$axis.1,
    axis.2 = superlet_1$axis.2,
    c1 = c1,
    o = o,
    x1 = superlet_1$x,
    y1 = superlet_1$y,
    x2 = superlet_2$x,
    y2 = superlet_2$y
  )
  class(output) <- "analyze.Xsuperlet"
  return(output)
}

