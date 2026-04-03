#' @title Cross-wavelet coherence analysis
#'
#' @description
#' Computes the cross-wavelet transform and wavelet coherence between two
#' time/depth series using a superlet-based wavelet approach. The function
#' internally interpolates both input series onto a common grid, computes
#' individual wavelet transforms, and derives cross-wavelet power, phase,
#' and coherence with configurable smoothing.
#'
#' @param data_1 First input dataset as a two-column matrix or data.frame.
#' First column must contain time/depth, second column the signal.
#' @param data_2 Second input dataset as a two-column matrix or data.frame.
#' Must have the same structure as \code{data_1}.
#' @param dj Spacing between discrete scales. Smaller values increase resolution.
#' Default is \code{1/100}.
#' @param lowerPeriod Lower bound of the period range to analyze.
#' @param upperPeriod Upper bound of the period range to analyze.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' @param omega_nr Number of oscillations in the Morlet/superlet wavelet.
#' Controls time-frequency resolution trade-off.
#' @param n_simulations Number of simulations (reserved for future significance testing).
#' Currently not used internally.
#' @param window.type.t Type of smoothingw indow applied in the time direction.
#' Options: "none", "bar" (Bartlett), "tri" (triangular), "box" (boxcar),
#' "han" (Hanning, default), "ham" (Hamming) or "bla" (Blackman)
#' @param window.type.s Type of smoothing window applied in the scale direction.
#' Same options as window.type.t.
#' @param abs_window_t Absolute smoothing window size in the time direction
#' (same units as input data).
#' @param abs_window_s Absolute smoothing window size in the scale direction
#' (in units of dj).
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
#' coh <- analyze_wavelet_coherence(
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


analyze_wavelet_coherence <- function(
    data_1,
    data_2,
    dj = 1/100,
    lowerPeriod = 2,
    upperPeriod = 1024,
    verbose = FALSE,
    omega_nr = 8,
    n_simulations = 10,
    window.type.t = 4,
    window.type.s = 4,
    abs_window_t = 500,
    abs_window_s = 0.5
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

  wt1 <- analyze_wavelet(
    d1,
    dj=dj,
    lowerPeriod = lowerPeriod,
    upperPeriod = upperPeriod,
    omega_nr = omega_nr,
    verbose = verbose,
    run_multicore = FALSE,
    pval = FALSE)

  wt2 <- analyze_wavelet(
    d2,
    dj=dj,
    lowerPeriod = lowerPeriod,
    upperPeriod = upperPeriod,
    omega_nr = omega_nr,
    verbose = verbose,
    run_multicore = FALSE,
    pval = FALSE)

  dt <- wt1$dt
  dj <- wt1$dj
  nr <- wt1$nr
  nc <- wt1$nc
  Scale <- wt1$Scale

  # --- CROSS-WAVELET ---
  Wave.xy  = (wt1$Wave * Conj(wt2$Wave)) / matrix(rep(Scale, nc), nrow = nr)
  Power.xy = Mod(Wave.xy)
  Phase.xy = Arg(Wave.xy)


  # ============================================================
  # --- WINDOWING & SMOOTHING FUNCTIONS (FULL SET) ---
  # ============================================================

  window.func <- function(type = 0, n = 0) {

    n <- max(1, floor(n))

    # --- 0. NONE ---
    if (is.element(type, c(0, "none"))) {
      window <- 1
    }

    # --- 1. BARTLETT ---
    if (is.element(type, c(1, "bar"))) {
      n <- max(n, 3)
      window <- 1 - abs((0:(n - 1)) - (n - 1)/2)/((n - 1)/2)
    }

    # --- 2. TRIANGULAR ---
    if (is.element(type, c(2, "tri"))) {
      n <- max(n, 2)
      window <- 1 - abs((0:(n - 1)) - (n - 1)/2)/(n/2)
    }

    # --- 3. BOXCAR ---
    if (is.element(type, c(3, "box"))) {
      window <- rep(1, n)
    }

    # --- 4. HANNING ---
    if (is.element(type, c(4, "han"))) {
      n <- max(n, 3)
      window <- 0.5 - 0.5 * cos(2 * pi * (0:(n - 1))/(n - 1))
    }

    # --- 5. HAMMING ---
    if (is.element(type, c(5, "ham"))) {
      n <- max(n, 2)
      window <- 0.53836 -
        (1 - 0.53836) * cos(2 * pi * (0:(n - 1))/(n - 1))
    }

    # --- 6. BLACKMAN ---
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

  # --- WINDOW SIZES ---
  window.size.t <- 2 * floor(abs_window_t / (2 * dt)) + 1
  window.size.s <- 2 * floor(abs_window_s / (2 * dj)) + 1

  if (window.size.t > nc) {
    window.size.t <- 2 * floor(nc/4) + 1
    warning("Time window too large, capped at 25% of series length.")
  }

  window2D <- window.func(window.type.s, window.size.s) %*%
    t(window.func(window.type.t, window.size.t))

  # --- SMOOTHED SPECTRA ---
  sPower.x <- Re(smooth2D(wt1$Power, window2D))
  sPower.y <- Re(smooth2D(wt2$Power, window2D))
  sWave.xy <- smooth2D(Wave.xy, window2D)

  # --- COHERENCE ---
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
    omega_nr=omega_nr,
    x1 = wt1$x,
    y1 = wt1$y,
    x2 = wt2$x,
    y2 = wt2$y
  )
  class(output) <- "analyze.waveletcoherence"
  return(output)
}
