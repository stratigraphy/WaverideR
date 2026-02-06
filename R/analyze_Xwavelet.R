#' @title Conduct the cross wavelet transform on a time series or signal
#'
#' @description
#' Compute the cross wavelet transform using a complex Morlet wavelet based on
#' the "analyze.coherency"  function of the WaveletComp
#' package.
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
#' @param verbose Logical. If TRUE, print progress and interpolation
#' information.
#'
#' @return
#'A list with class \code{"analyze.Xwavelet"} containing
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
#'Xetp <- analyze_xwavelet(
#'  data_1 = etp_1,
#'  data_2  = etp_2,
#'  upperPeriod = 1024,
#'  lowerPeriod = 2,
#'  verbose = FALSE
#')
#'
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom stats fft

analyze_xwavelet <- function(
    data_1,
    data_2,
    dj = 1/100,
    lowerPeriod = 2,
    upperPeriod = 1024,
    verbose = FALSE,
    omega_nr = 8
) {

  # --- PREPROCESSING ---
  data_1 <- data_1[order(data_1[,1]), ]
  data_2 <- data_2[order(data_2[,1]), ]

  xmin <- max(min(data_1[,1]), min(data_2[,1]))
  xmax <- min(max(data_1[,1]), max(data_2[,1]))

  dx <- max(median(diff(data_1[,1])), median(diff(data_2[,1])))
  grid <- seq(xmin, xmax, by = dx)

  d1 <- data.frame(
    x = grid,
    y = approx(data_1[,1], data_1[,2], grid, rule = 1)$y
  )
  d2 <- data.frame(
    x = grid,
    y = approx(data_2[,1], data_2[,2], grid, rule = 1)$y
  )

  # --- WAVELETS ---
  wt1 <- analyze_wavelet(
    d1,
    dj=dj,
    lowerPeriod = lowerPeriod,
    upperPeriod = upperPeriod,
    omega_nr = omega_nr,
    verbose = verbose,
    pval = FALSE,
    n_simulations = 10,
    run_multicore = FALSE
  )

  wt2 <- analyze_wavelet(
    d2,
    dj=dj,
    lowerPeriod = lowerPeriod,
    upperPeriod = upperPeriod,
    omega_nr = omega_nr,
    verbose = verbose,
    pval = FALSE,
    n_simulations = 10,
    run_multicore = FALSE
  )

  nr <- wt1$nr
  nc <- wt1$nc
  Scale <- wt1$Scale

  # --- CROSS-WAVELET ---
  Wave.xy <- (wt1$Wave * Conj(wt2$Wave)) /
    matrix(rep(Scale, nc), nrow = nr)

  Power.xy <- Mod(Wave.xy)
  Phase.xy <- Arg(Wave.xy)

  Power.avg = rowMeans(as.matrix(Power.xy))



  out <- list(
    Wave   = Wave.xy,
    Power  = Power.xy,
    Phase  = Phase.xy,
    Power.avg = Power.avg,
    dt = dt,
    dj = dj,
    Power.avg = colMeans(t(Power.xy), na.rm = TRUE),
    Scale = Scale,
    Period = wt1$Period,
    nc = nc,
    nr = Nf,
    axis.1 = wt1$axis.1,
    axis.2 = wt1$axis.2,
    omega_nr = omega_nr,
    x1 = wt1$x,
    y1 = wt1$y,
    x2 = wt2$x,
    y2 = wt2$y
  )

  class(out) <- "analyze.Xwavelet"
  return(out)
}









