#' @title Cross-wavelet coherence analysis
#'
#' @description
#' Computes the cross-wavelet transform and wavelet coherence between two
#' time/depth series using a superlet-based wavelet approach. The function
#' internally interpolates both input series onto a common grid, computes
#' individual wavelet transforms, and derives cross-wavelet power, phase,
#' and coherence with configurable smoothing.
#'
#' @param data_1 First record: a two-column matrix or data frame with depth/time in column 1 and the proxy value
#' in column 2 (no \code{NA} values). Both records are sorted and linearly interpolated to a common, regular grid
#' over their overlapping interval, with the coarser of the two median sampling intervals. Required, no default.
#' @param data_2 Second record, in the same format and depth/time unit as \code{data_1} (two-column matrix or data
#' frame: depth/time, proxy). The phase is the phase of \code{data_1} minus that of \code{data_2}. Required, no
#' default.
#' @param dj Spacing between successive scales in octaves (numeric, > 0), passed to \code{\link{analyse_wavelet}};
#' \code{dj=1/x} gives x periods per doubling of the period, smaller values give a finer period axis.
#' \code{Default=1/100}.
#' @param lowerPeriod Shortest period to analyse, in the same unit as the depth/time column of the records
#' (numeric, > 0; at least twice the common sampling interval). \code{Default=2}.
#' @param upperPeriod Longest period to analyse, in the same unit as the depth/time column of the records
#' (numeric, > \code{lowerPeriod}). \code{Default=1024}.
#' @param verbose Print the sampling interval to which the records were interpolated, passed to
#' \code{\link{analyse_wavelet}} (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param omega_nr Number of cycles (central angular frequency omega_0) of the Morlet wavelet, passed to
#' \code{\link{analyse_wavelet}} (numeric, > 0); controls the trade-off between period and depth/time
#' resolution. \code{Default=8}.
#' @param wavelet Mother wavelet (character string or function), passed to \code{\link{analyse_wavelet}}, given as a
#' name
#' (PyWavelets style names with the order appended, such as
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
#' \code{"mexh"}, \code{"ricker"}, \code{"marr"} and \code{"sombrero"} are the Mexican hat, order 2), \code{"morl"}
#' (real
#' Morlet), \code{"hhat"} (Hermitian hat), \code{"poisson"}, and the orthogonal wavelets \code{"haar"}, \code{"db"}
#' (Daubechies), \code{"sym"} (Symlets) and \code{"battle"} (Battle-Lemarie spline wavelets), which are computed in the
#' frequency domain from their filters. The period axis of the orthogonal wavelets is approximate because their
#' spectra are broad and not symmetric. A narrow wavelet (Paul, low-order DOG, Shannon, Haar) resolves position
#' well, a wide one (Morse with a large beta, Morlet with a high omega_nr, bump) resolves period well.
#' \code{Default="morlet"}.
#' @param wavelet_order Parameter(s) of the wavelet (numeric scalar or vector; \code{NULL} uses the defaults listed
#' here): \code{paul}: order m (default 4); \code{dog}: derivative m
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
#' (default 8); \code{battle}: spline order 1, 3, 5 or 7 (default 3). For the Morlet wavelet use \code{omega_nr}.
#' Passed to \code{\link{analyse_wavelet}}. \code{Default=NULL}.
#' @param n_simulations Number of simulations (integer). Not used: reserved for a future significance test.
#' \code{Default=10}.
#' @param window.type.t Type of the smoothing window applied in the depth/time direction before the coherence is
#' computed, given as a number or a character code: \code{0} or \code{"none"} (no smoothing), \code{1} or
#' \code{"bar"} (Bartlett), \code{2} or \code{"tri"} (triangular), \code{3} or \code{"box"} (boxcar), \code{4} or
#' \code{"han"} (Hanning), \code{5} or \code{"ham"} (Hamming), \code{6} or \code{"bla"} (Blackman).
#' \code{Default=4} (Hanning).
#' @param window.type.s Type of the smoothing window applied in the period (scale) direction; same options as
#' \code{window.type.t}. \code{Default=4} (Hanning).
#' @param abs_window_t Length of the smoothing window in the depth/time direction, in the same unit as the
#' depth/time column of the records (numeric, > 0); converted to an odd number of samples and capped at about a
#' quarter of the record length (with a warning). Longer windows give smoother but less local coherence.
#' \code{Default=500}.
#' @param abs_window_s Length of the smoothing window in the period (scale) direction, in octaves (numeric, > 0;
#' converted to an odd number of scale steps of size \code{dj}). \code{Default=0.5}.
#'
#' @return
#' A list of class containing:
#'   Wave Smoothed cross-wavelet transform (complex)
#'   Coherence Wavelet coherence matrix
#'   Phase Phase difference between the two signals
#'   Coh.avg Average coherence over time
#'   Period (m)vector corresponding to scales
#'   Scale Wavelet scales
#'   dt Time step of interpolated grid
#'   dj Scale resolution
#'   axis.1 Time/depth axis
#'   axis.2 Period axis
#'   nc, nr Dimensions of the wavelet matrices
#'   x1, y1 Interpolated first dataset
#'   x2, y2 Interpolated second dataset
#'
#' @author
#' plotting code based on the "analyze.coherency" function
#'  of the 'WaveletComp' R package
#'
#' @references
#' Torrence, C., and G. P. Compo (1998).
#' A practical guide to wavelet analysis.
#' \emph{Bulletin of the American Meteorological Society}, 79(1), 61–78.
#'
#' @examples
#' \donttest{
#' # Generate two synthetic signals
#' t <- seq(0, 1000, by = 1)
#' x1 <- sin(2 * pi * t / 100) + rnorm(length(t), 0, 0.5)
#' x2 <- sin(2 * pi * t / 100 + pi/4) + rnorm(length(t), 0, 0.5)
#'
#' data_1 <- cbind(t, x1)
#' data_2 <- cbind(t, x2)
#'
#' coh <- analyse_wavelet_coherence(
#'   data_1 = data_1,
#'   data_2 = data_2,
#'   lowerPeriod = 2,
#'   upperPeriod = 256
#' )
#'
#' }
#'
#' @export
#' @importFrom stats approx

analyse_wavelet_coherence <- function(
    data_1,
    data_2,
    dj = 1/100,
    lowerPeriod = 2,
    upperPeriod = 1024,
    verbose = FALSE,
    omega_nr = 8,
    wavelet = "morlet",
    wavelet_order = NULL,
    n_simulations = 10,
    window.type.t = 4,
    window.type.s = 4,
    abs_window_t = 500,
    abs_window_s = 0.5
) {

  # PREPROCESSING
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

  wt1 <- analyse_wavelet(
      d1,
      dj=dj,
      lowerPeriod = lowerPeriod,
      upperPeriod = upperPeriod,
      omega_nr = omega_nr,
      wavelet = wavelet,
      wavelet_order = wavelet_order,
      verbose = verbose,
      run_multicore = FALSE,
      pval = FALSE)

  wt2 <- analyse_wavelet(
      d2,
      dj=dj,
      lowerPeriod = lowerPeriod,
      upperPeriod = upperPeriod,
      omega_nr = omega_nr,
      wavelet = wavelet,
      wavelet_order = wavelet_order,
      verbose = verbose,
      run_multicore = FALSE,
      pval = FALSE)

  dt <- wt1$dt
  dj <- wt1$dj
  nr <- wt1$nr
  nc <- wt1$nc
  Scale <- wt1$Scale

  # CROSS-WAVELET
  Wave.xy  = (wt1$Wave * Conj(wt2$Wave)) / matrix(rep(Scale, nc), nrow = nr)
  Power.xy = Mod(Wave.xy)
  Phase.xy = Arg(Wave.xy)

  # WINDOWING & SMOOTHING FUNCTIONS (FULL SET)

  window.func <- function(type = 0, n = 0) {

    n <- max(1, floor(n))

    # NONE
    if (is.element(type, c(0, "none"))) {
      window <- 1
    }

    # BARTLETT
    if (is.element(type, c(1, "bar"))) {
      n <- max(n, 3)
      window <- 1 - abs((0:(n - 1)) - (n - 1)/2)/((n - 1)/2)
    }

    # TRIANGULAR
    if (is.element(type, c(2, "tri"))) {
      n <- max(n, 2)
      window <- 1 - abs((0:(n - 1)) - (n - 1)/2)/(n/2)
    }

    # BOXCAR
    if (is.element(type, c(3, "box"))) {
      window <- rep(1, n)
    }

    # HANNING
    if (is.element(type, c(4, "han"))) {
      n <- max(n, 3)
      window <- 0.5 - 0.5 * cos(2 * pi * (0:(n - 1))/(n - 1))
    }

    # HAMMING
    if (is.element(type, c(5, "ham"))) {
      n <- max(n, 2)
      window <- 0.53836 -
          (1 - 0.53836) * cos(2 * pi * (0:(n - 1))/(n - 1))
    }

    # BLACKMAN
    if (is.element(type, c(6, "bla"))) {
      n <- max(n, 2)
      window <- 7938/18608 -
          (9240/18608) * cos(2 * pi * (0:(n - 1))/(n - 1)) +
          (1430/18608) * cos(4 * pi * (0:(n - 1))/(n - 1))
    }

    window/sum(window)
    return(window)

  }

  smooth2D <- function(mat, window2D) {

    mat.pad <- matrix(
        0,
        nrow(mat) + nrow(window2D) - 1,
        ncol(mat) + ncol(window2D) - 1
    )

    window2D.pad <- mat.pad

    mat.pad[1:nrow(mat), 1:ncol(mat)] <- mat
    window2D.pad[1:nrow(window2D), 1:ncol(window2D)] <- window2D

    smooth.mat <- fft(
        fft(mat.pad) * fft(window2D.pad),
        inverse = TRUE
    ) / length(mat.pad)

    smooth.mat[
        floor(nrow(window2D)/2) + (1:nrow(mat)),
        floor(ncol(window2D)/2) + (1:ncol(mat))
    ]
  }

  # WINDOW SIZES
  window.size.t <- 2 * floor(abs_window_t / (2 * dt)) + 1
  window.size.s <- 2 * floor(abs_window_s / (2 * dj)) + 1

  if (window.size.t > nc) {
    window.size.t <- 2 * floor(nc/4) + 1
    warning("Time window too large, capped at 25% of series length.")
  }

  window2D <- window.func(window.type.s, window.size.s) %*%
      t(window.func(window.type.t, window.size.t))

  # SMOOTHED SPECTRA
  sPower.x <- Re(smooth2D(wt1$Power, window2D))
  sPower.y <- Re(smooth2D(wt2$Power, window2D))
  sWave.xy <- smooth2D(Wave.xy, window2D)

  # COHERENCE
  Coherence <- Mod(sWave.xy)^2 / (sPower.x * sPower.y)
  Coherence[sPower.x * sPower.y == 0] <- 0

  output <- list(
      Wave = t(sWave.xy),
      Coherence = t(Coherence),
      Phase = t(Arg(Wave.xy)),
      dt = dt,
      dj = dj,
      Coh.avg = colMeans(t(Coherence), na.rm = TRUE),
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
  class(output) <- "analyze.waveletcoherence"
  return(output)
}
