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

  ?analyze_wavelet

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
  Wave.xy  = (wt_1_wt$Wave * Conj(wt_2_wt$Wave)) / matrix(rep(Scale, nc), nrow = nr)
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

  out <- list(
    Coherence = Coherence,
    Phase.xy = Arg(Wave.xy),

  )

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
