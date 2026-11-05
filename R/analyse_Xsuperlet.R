#' @title Conduct the cross superlet transform on a time series or signal
#'
#' @description
#' Compute the cross superlet transform using a complex Morlet wavelet based on
#' the approach of Moca et al. (2021) and  the "analyze.coherency"  function of the WaveletComp
#' package. The superlet transform increases time frequency resolution by
#' combining multiple wavelet orders at each frequency.
#'
#' @param data_1 First record: a matrix or data frame with depth/time in column 1 and the proxy value in column 2.
#' Both records are sorted, rows with \code{NA} are removed and they are interpolated to a common, regular grid over
#' their overlapping interval (spacing: the coarser of the two median sampling intervals); at least 10 common
#' points are needed. Required (\code{Default=NULL} gives an error).
#'
#' @param data_2 Second record, in the same format and depth/time unit as \code{data_1} (matrix or data frame:
#' depth/time in column 1, proxy value in column 2). The cross phase is the phase of \code{data_1} minus that of
#' \code{data_2}. Required (\code{Default=NULL} gives an error).
#'
#' @param upperPeriod Longest period to analyse, in the same unit as the depth/time column of the records
#' (numeric, > \code{lowerPeriod}); sets the lowest analysed frequency. \code{Default=1024}.
#'
#' @param lowerPeriod Shortest period to analyse, in the same unit as the depth/time column of the records
#' (numeric, > 0; at least twice the common sampling interval); sets the highest analysed frequency.
#' \code{Default=2}.
#'
#' @param Nf Number of frequencies of the superlet spectrum (integer, >= 1). Frequencies are logarithmically
#' spaced between \code{1/upperPeriod} and \code{1/lowerPeriod}. \code{Default=128}.
#'
#' @param c1 Base number of cycles of the Morlet wavelet of order 1 (numeric, > 0, typically 2 to 5): the
#' fundamental wavelet width. \code{Default=3}.
#'
#' @param o Minimum and maximum superlet order (numeric vector of length 2, values >= 1), from the lowest to the
#' highest frequency; a single number gives a fixed order and \code{NULL} uses order 1 everywhere (see
#' \code{\link{analyse_superlet}}). \code{Default=c(1, 10)}.
#'
#' @param mult How the number of cycles grows with the order i within a superlet (logical, TRUE or FALSE):
#' \code{TRUE} multiplicative (\eqn{c_i = i c_1}); \code{FALSE} additive (\eqn{c_i = c_1 + i - 1}, Moca et al.,
#' 2021). \code{Default=TRUE}.
#'
#' @param order_scaling How the superlet order increases from \code{o[1]} (lowest frequency) to \code{o[2]}
#' (highest frequency) (character): \code{"none"} or \code{"log2"} (linear in log2(frequency)), \code{"linear"}
#' (linear in frequency), \code{"sqrt"}, \code{"quadratic"} or \code{"power"} (see
#' \code{\link{analyse_superlet}}). \code{Default="log2"}.
#'
#' @param order_alpha Exponent of the \code{"power"} order scaling (numeric, > 0); only used when
#' \code{order_scaling="power"}. Values greater than 1 raise the order mainly at the highest frequencies
#' (sharpening short periods), values smaller than 1 already at low frequencies. \code{Default=1}.
#'
#' @param verbose Print the sampling interval to which each record was interpolated and show progress bars
#' (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @param detrend Removal of a trend from both records before the transform, passed to
#' \code{\link{analyse_superlet}} (character): \code{"linear"}, \code{"mean"} or \code{"none"}.
#' \code{Default="linear"}.
#' @param standardize Divide both (detrended) records by their standard deviation before the transform, passed to
#' \code{\link{analyse_superlet}} (logical, TRUE or FALSE). \code{Default=TRUE}.
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
#' The "analyse_Xsuperlet" that generates the input for the plotting function
#' is based on the matlab code in Moca et al. (2021) and the the "analyze.coherency" function of the  'WaveletComp' R package
#'
#' @references
#' Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., & Mureșan, R. C. (2021).
#' Time-frequency super-resolution with superlets.
#' Nature Communications, 12(1), 337.
#' \doi{10.1038/s41467-020-20539-9}
#'@examples
#' \donttest{
#'#Example 1. A cross superlet of two etp solutions with noise overprint
#'etp_1 <- astrochron::etp(
#'  tmin = 0,
#'  tmax = 300,
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
#'etp_2 <- astrochron::etp(
#'  tmin = 0,
#'  tmax = 300,
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
#'Xetp <- analyse_Xsuperlet(
#'  data_1 = etp_1,
#'  data_2  = etp_2,
#'  upperPeriod = 256,
#'  lowerPeriod = 2,
#'  Nf = 32,
#'  c1 = 3,
#'  o = c(1, 10),
#'  mult = TRUE,
#'  order_scaling = "log2",
#'  order_alpha = 1,
#'  verbose = FALSE
#')
#'}
#' @param wavelet Mother wavelet (character string or function), given as a name (PyWavelets style names with the
#' order appended,
#' such as \code{"db4"}, \code{"sym8"}, \code{"cgau3"}, \code{"gaus2"} or \code{"cmor1.5-1"}, also work) or as a
#' function
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
#' \code{"mexh"}, \code{"ricker"}, \code{"marr"} and \code{"sombrero"} are the Mexican hat, order 2), \code{"morl"}
#' (real
#' Morlet), \code{"hhat"} (Hermitian hat), \code{"poisson"}, and the orthogonal wavelets \code{"haar"}, \code{"db"}
#' (Daubechies), \code{"sym"} (Symlets) and \code{"battle"} (Battle-Lemarie spline wavelets), which are computed in the
#' frequency domain from their filters. The period axis of the orthogonal wavelets is approximate because their
#' spectra are broad and not symmetric. A narrow wavelet (Paul, low-order DOG, Shannon, Haar) resolves position
#' well, a wide one (Morse with a large beta, Morlet with a high omega_nr, bump) resolves period well.
#' \code{Default="morlet"}.
#' @param wavelet_order Parameter(s) of the wavelet (see \code{\link{analyse_superlet}}): \code{paul}: order m
#' (default 4); \code{dog}: derivative m
#' (default 2); \code{morse}: \code{c(gamma, beta)} or beta alone with gamma = 3 (default \code{c(3, 20)});
#' \code{bump}: \code{c(mu, sigma)} (default \code{c(5, 0.6)}); \code{cmor}: \code{c(bandwidth, centre frequency)}
#' (default \code{c(1.5, 1)}); \code{cgau}: derivative m (default 4); \code{shannon}: \code{c(bandwidth, centre
#' frequency)} (default \code{c(2, 1.5)}); \code{fbsp}: \code{c(m, bandwidth, centre frequency)} (default
#' \code{c(2, 2, 1.5)}); \code{loggabor}: \code{c(centre angular frequency, log-width)} (default \code{c(6, 0.25)});
#' \code{meyer} and \code{hhat} have no parameter; \code{morlet0}: omega_0 (default omega_nr); \code{morl}: omega_0
#' (default 5);
#' \code{poisson}: order n (default 2); \code{gammatone}: \code{c(n, centre angular frequency, bandwidth)} (default
#' \code{c(4, 6, 1.5)}); \code{hann}: \code{c(centre angular frequency, half width)} (default \code{c(6, 3)});
#' \code{db}: number of vanishing moments 1 to 20 (default 4, \code{"haar"} = \code{"db1"}); \code{sym}: 2 to 20
#' (default 8); \code{battle}: spline order 1, 3, 5 or 7 (default 3).
#' Passed to \code{\link{analyse_superlet}}; for the Morlet wavelet the number of cycles is set by \code{c1}, \code{o}
#' and \code{mult}. \code{Default=NULL}.
#'
#' @param standardise British spelling of \code{standardize} (logical); if given (not \code{NULL}) it overrides
#' \code{standardize}. \code{Default=NULL}.
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom stats fft
#' @importFrom astrochron etp

analyse_Xsuperlet <- function(data_1 = NULL,
    data_2  = NULL,
    upperPeriod = 1024,
    lowerPeriod = 2,
    Nf = 128,
    c1 = 3,
    o = c(1, 10),
    mult = TRUE,
    order_scaling = "log2",
    order_alpha = 1,
    verbose = FALSE,
    detrend = "linear",
    standardize = TRUE,
    standardise = NULL,
    wavelet = "morlet",
    wavelet_order = NULL) {
  if (!is.null(standardise)) {
    standardize <- standardise
  }  # British spelling synonym
  data_1 <- na.omit(data_1[order(data_1[, 1]), ])
  data_2 <- na.omit(data_2[order(data_2[, 1]), ])

  # common regular grid over the overlapping interval (the coarser spacing of the two records)
  xmin <- max(min(data_1[, 1]), min(data_2[, 1]))
  xmax <- min(max(data_1[, 1]), max(data_2[, 1]))
  dx <- max(median(diff(data_1[, 1])), median(diff(data_2[, 1])))
  common_grid <- seq(from = xmin, to = xmax, by = dx)
  y1 <- approx(data_1[, 1], data_1[, 2], xout = common_grid, rule = 1)$y
  y2 <- approx(data_2[, 1], data_2[, 2], xout = common_grid, rule = 1)$y
  keep <- is.finite(y1) & is.finite(y2)
  if (sum(keep) < 10) {
    stop("the two records have (almost) no overlap")
  }
  data_1_common <- data.frame(x = common_grid[keep], proxy = y1[keep])
  data_2_common <- data.frame(x = common_grid[keep], proxy = y2[keep])

  run <- function(d) {
    analyse_superlet(data = d, upperPeriod = upperPeriod, lowerPeriod = lowerPeriod, Nf = Nf, c1 = c1, o = o,
        mult = mult, order_scaling = order_scaling, order_alpha = order_alpha, verbose = verbose,
        detrend = detrend, standardize = standardize,
        wavelet = wavelet, wavelet_order = wavelet_order)
  }
  superlet_1 <- run(data_1_common)
  superlet_2 <- run(data_2_common)

  # Cross spectrum: amplitudes are geometric means over the superlet orders, phases circular means, so
  # the cross phase (phase of record 1 minus phase of record 2) is well defined. The superlet amplitudes are not
  # scale dependent, so there is no division by the period.
  Wave.xy  <- superlet_1$Wave * Conj(superlet_2$Wave)
  Power.xy <- Mod(Wave.xy)
  Phase.xy <- Arg(Wave.xy)

  output <- list(
      Wave = t(Wave.xy),
      Power = t(Power.xy),
      Phase = t(Phase.xy),
      dt = superlet_1$dt,
      dj = 1 / Nf,
      Power.avg = rowMeans(Power.xy, na.rm = TRUE),
      Scale = superlet_1$Scale,
      Period = superlet_1$Period,
      nc = superlet_1$nc,
      nr = Nf,
      axis.1 = superlet_1$axis.1,
      axis.2 = superlet_1$axis.2,
      c1 = c1,
      o = o,
      wavelet = superlet_1$wavelet,
      wavelet_order = wavelet_order,
      coi = superlet_1$coi,
      coi.1 = superlet_1$coi.1,
      coi.2 = superlet_1$coi.2,
      x1 = superlet_1$x,
      y1 = superlet_1$y,
      x2 = superlet_2$x,
      y2 = superlet_2$y
  )
  class(output) <- "analyze.Xsuperlet"
  return(output)
}
