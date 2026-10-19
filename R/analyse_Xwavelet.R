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
#' @param omega_nr . omega number: central angular frequency of the Morlet wavelet
#' @param wavelet Mother wavelet, given as a name (PyWavelets style names with the order appended, such as
#' \code{"db4"}, \code{"sym8"}, \code{"cgau3"}, \code{"gaus2"} or \code{"cmor1.5-1"}, also work) or as a function
#' \code{psi_hat(omega)} that returns the Fourier transform of a user defined mother wavelet.
#' Analytic (complex) wavelets, which give phase and amplitude and can be synchrosqueezed: \code{"morlet"}
#' (default; aliases \code{"gabor"}, \code{"amor"}), \code{"morlet0"} (exactly zero mean Morlet, admissible for a
#' small omega_nr; alias \code{"morlex"}), \code{"paul"} (alias \code{"cauchy"}), \code{"morse"} (generalised Morse
#' wavelet of Olhede and Walden, 2002; aliases \code{"gmw"}, \code{"airy"}), \code{"bump"}, \code{"cmor"} (complex
#' Morlet with bandwidth and centre frequency as in Matlab/PyWavelets), \code{"cgau"} (analytic Gaussian derivative;
#' \code{"cmhat"} is the complex Mexican hat), \code{"shannon"} (complex Shannon, box spectrum), \code{"fbsp"}
#' (frequency B-spline), \code{"meyer"}, \code{"loggabor"} (alias \code{"lognormal"}), \code{"gammatone"} (asymmetric in
#' time: sharp onset, slow decay) and \code{"hann"} (raised cosine band).
#' Real wavelets, which give a real transform (no phase or envelope): \code{"dog"} (derivative of Gaussian;
#' \code{"mexh"}, \code{"ricker"}, \code{"marr"} and \code{"sombrero"} are the Mexican hat, order 2), \code{"morl"} (real
#' Morlet), \code{"hhat"} (Hermitian hat), \code{"poisson"}, and the orthogonal wavelets \code{"haar"}, \code{"db"}
#' (Daubechies), \code{"sym"} (Symlets) and \code{"battle"} (Battle-Lemarie spline wavelets), which are computed in the
#' frequency domain from their filters. The period axis of the orthogonal wavelets is approximate because their
#' spectra are broad and not symmetric. A narrow wavelet (Paul, low-order DOG, Shannon, Haar) resolves position
#' well, a wide one (Morse with a large beta, Morlet with a high omega_nr, bump) resolves period well.
#' @param wavelet_order Parameter(s) of the wavelet: \code{paul}: order m (default 4); \code{dog}: derivative m
#' (default 2); \code{morse}: \code{c(gamma, beta)} or beta alone with gamma = 3 (default \code{c(3, 20)});
#' \code{bump}: \code{c(mu, sigma)} (default \code{c(5, 0.6)}); \code{cmor}: \code{c(bandwidth, centre frequency)}
#' (default \code{c(1.5, 1)}); \code{cgau}: derivative m (default 4); \code{shannon}: \code{c(bandwidth, centre
#' frequency)} (default \code{c(2, 1.5)}); \code{fbsp}: \code{c(m, bandwidth, centre frequency)} (default
#' \code{c(2, 2, 1.5)}); \code{loggabor}: \code{c(centre angular frequency, log-width)} (default \code{c(6, 0.25)});
#' \code{meyer} and \code{hhat} have no parameter; \code{morlet0}: omega_0 (default omega_nr); \code{morl}: omega_0 (default 5);
#' \code{poisson}: order n (default 2); \code{gammatone}: \code{c(n, centre angular frequency, bandwidth)} (default
#' \code{c(4, 6, 1.5)}); \code{hann}: \code{c(centre angular frequency, half width)} (default \code{c(6, 3)});
#' \code{db}: number of vanishing moments 1 to 20 (default 4, \code{"haar"} = \code{"db1"}); \code{sym}: 2 to 20
#' (default 8); \code{battle}: spline order 1, 3, 5 or 7 (default 3). For the Morlet wavelet use \code{omega_nr}.
#'
#' @param dj Spacing between successive scales. The CWT analyses analyses the signal using successive periods
#' which increase by the power of 2 (e.g.2^0=1,2^1=2,2^2=4,2^3=8,2^4=16). To have more resolution
#' in-between these steps the dj parameter exists, the dj parameter specifies how many extra steps/spacing in-between
#' the power of 2 scaled CWT is added. The amount of steps is 1/x with a higher x indicating a smaller spacing.
#' Increasing the increases the computational time of the CWT
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
#' The "analyse_Xsuperlet" that generates the input for the plotting function
#' is based on the matlab code in Moca et al. (2021) and the the "analyze.coherency" function of the  'WaveletComp' R package
#'
#' @references
#' Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., & Mureșan, R. C. (2021).
#' Time-frequency super-resolution with superlets.
#' Nature Communications, 12(1), 337.
#' \doi{10.1038/s41467-020-20539-9}
#' @examples
#' \donttest{
#'#Example 1. A cross superlet of two etp solutions with noise overprint
#'etp_1 <- astrochron::etp(
#'  tmin = 0,
#'  tmax = 500,
#'  dt = 2,
#'  eWt = 1.5,
#'  oWt = 0.75,
#'  pWt = 1,
#'  esinw = TRUE,
#'  standardize = TRUE,
#'  genplot = FALSE,
#'  verbose = FALSE
#')
#'
#'etp_2 <- astrochron::etp(
#'  tmin = 0,
#'  tmax = 500,
#'  dt = 2,
#'  eWt = 1,
#'  oWt = 0.5,
#'  pWt = 1.5,
#'  esinw = TRUE,
#'  standardize = TRUE,
#'  genplot = FALSE,
#'  verbose = FALSE
#')
#'
#'etp_1[, 2] <- etp_1[, 2] + wr_coloured_noise(
#'  nrow(etp_1),
#'  sd = sd(etp_1[, 2]) / 1.5,
#'  mean = mean(etp_1[, 2]),
#'  phi = 0.9
#')
#'etp_2[, 2] <- etp_2[, 2] + wr_coloured_noise(
#'  nrow(etp_2),
#'  sd = sd(etp_2[, 2]) / 1.5,
#'  mean = mean(etp_2[, 2]),
#'  phi = 0.9
#')
#'
#'Xetp <- analyse_Xwavelet(
#'  data_1 = etp_1,
#'  data_2  = etp_2,
#'  upperPeriod = 512,
#'      dj = 1/20,
#'  lowerPeriod = 4,
#'  verbose = FALSE,
#'  omega_nr = 8
#')
#'}
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom stats fft

analyse_Xwavelet <- function(
    data_1,
    data_2,
    dj = 1/100,
    lowerPeriod = 2,
    upperPeriod = 1024,
    verbose = FALSE,
    omega_nr = 4,
    wavelet = "morlet",
    wavelet_order = NULL
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
  wt1 <- analyse_wavelet(
    d1,
    dj=dj,
    lowerPeriod = lowerPeriod,
    upperPeriod = upperPeriod,
    omega_nr = omega_nr,
    wavelet = wavelet,
    wavelet_order = wavelet_order,
    verbose = verbose,
    pval = FALSE,
    n_simulations = 10,
    run_multicore = FALSE
  )

  wt2 <- analyse_wavelet(
    d2,
    dj=dj,
    lowerPeriod = lowerPeriod,
    upperPeriod = upperPeriod,
    omega_nr = omega_nr,
    wavelet = wavelet,
    wavelet_order = wavelet_order,
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
    Wave   = t(Wave.xy),
    Power  = t(Power.xy),
    Phase  = t(Phase.xy),
    Power.avg = Power.avg,
    dt = dx,
    dj = dj,
    Scale = Scale,
    Period = wt1$Period,
    nc = nc,
    nr = nr,
    axis.1 = wt1$axis.1,
    axis.2 = wt1$axis.2,
    omega_nr = omega_nr,
    wavelet = wt1$wavelet,
    wavelet_order = wt1$wavelet_order,
    x1 = wt1$x,
    y1 = wt1$y,
    x2 = wt2$x,
    y2 = wt2$y
  )

  class(out) <- "analyze.Xwavelet"
  return(out)
}









