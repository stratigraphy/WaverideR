#' @title Generate a log2 spaced EHA spectra
#'
#' @description
#' Compute a log2 spaced Evolutive Harmonic Analysis (EHA) & Evolutive Power
#' Spectral Analysis
#'This is a wrapper function for the "eha" function of the 'astrochron' R package
#' @param data Input series (matrix or data frame with two columns: depth/time in column 1, proxy value in column
#' 2), preferably evenly sampled (the padded series is linearly interpolated with \code{astrochron::linterp}).
#' Required (\code{Default=NULL}).
#' @param win Length of the moving EHA window, in the depth/time unit of \code{data} (numeric, > 0, shorter than
#' the record). The series is also padded by half a window at both ends (see \code{padding}). Required
#' (\code{Default=NULL}).
#' @param demean Remove the mean of each window (logical). Currently not used: the windows are always demeaned
#' (\code{demean=TRUE} is passed to \code{\link[astrochron]{eha}}). \code{Default=NULL}.
#' @param detrend Remove a linear trend from each window (logical). Currently not used: the windows are always
#' detrended (\code{detrend=TRUE} is passed to \code{\link[astrochron]{eha}}). \code{Default=NULL}.
#' @param tbw Time-bandwidth product of the multitaper (MTM) spectra, passed to \code{\link[astrochron]{eha}}
#' (numeric, <= 10, e.g. 2, 3 or 4); larger values give smoother spectra with a lower frequency resolution.
#' Required (\code{Default=NULL}).
#' @param upperPeriod Longest period to analyse, in the depth/time unit of \code{data} (numeric, > \code{lowerPeriod},
#' at most about \code{win}); sets \code{fmin = 1 / upperPeriod} of the EHA and the top of the log2 period grid.
#' Required (\code{Default=NULL}).
#' @param lowerPeriod Shortest period to analyse, in the depth/time unit of \code{data} (numeric, > 0, at least twice
#' the sampling interval); sets \code{fmax = 1 / lowerPeriod} and the bottom of the log2 period grid. Required
#' (\code{Default=NULL}).
#' @param pad Number of points to which each window is zero-padded (integer or \code{NULL}); odd values are increased
#' by one, values below the number of points per window are raised to it, and \code{5 * pad} is passed to
#' \code{\link[astrochron]{eha}} to obtain a fine low-frequency spacing. \code{NULL} uses
#' \code{2 * 2^ceiling(log2(points per window))}. \code{Default=NULL}.
#' @param padding Values used to extend the series by half a window at both ends, so that the EHA covers the whole
#' record (character): \code{"noise"} (AR(1) red noise with the mean, standard deviation and lag-1 autocorrelation
#' of the end segments; note that this calls \code{set.seed()}), \code{"mean"} (also \code{"average"} or
#' \code{"avg"}; the mean of the series) or \code{"zero"} (zeros). Any other value is treated as \code{"mean"}.
#' \code{Default="noise"}.
#'
#' @return
#' The output is a list (analyze.eha object) which contain
#' 12 objects which are the result of running the eha.
#'Object 1: depth - depth/time of axis
#'Object 2: log2_period - log2 scales period
#'Object 3: pwr - Power values of the EHA run
#'Object 4: amp - Amplitude values of the EHA run
#'Object 5: prob - Probability values of the EHA run
#'Object 6: f_test - F test values of the EHA run Scale size
#'Object 7: Power.avg  - Average power values
#'Object 8: amp.avg -  Average amplitude values
#'Object 9: prob.avg - Average probability value
#'Object 10: f_test.avg - Average f test values
#'Object 11: x - x values of the data set
#'Object 12: y - y values of the data set
#'Object 13: window size
#'
#' @author
#' Code based on on the "eha" function of the 'astrochron' R package
#'
#' @references
#'S.R. Meyers, 2012, Seeing Red in Cyclic Stratigraphy: Spectral Noise
#' Estimation for #'Astrochronology: Paleoceanography, 27, PA3228,
#' <doi:10.1029/2012PA002307>
#'
#'S.R. Meyers, 2019 Cyclostratigraphy and the problem of astrochronologic
#'testing, Earth-Science Reviews,
#'Volume 190, <doi.org/10.1016/j.earscirev.2018.11.015.>
#'
#' @examples
#' \donttest{
#'#Example 1. Total Solar Irradiance
#'# data set of Steinhilber et al., (2012)
#'
#'TSI_eha_log2 <- eha_log2(data = TSI,
#'win =  8192,
#'tbw = 4,
#'demean = TRUE,
#'detrend = TRUE,
#'upperPeriod = 8192,
#'lowerPeriod = 50,
#'pad = NULL,
#'padding = "noise")
#'
#'#Example 2. Magnetic susceptibility data set of Pas et al., (2018)
#'mag_eha_log2 <- eha_log2(data = mag,
#'win =  50,
#'tbw = 4,
#'demean = TRUE,
#'detrend = TRUE,
#'upperPeriod = 50,
#'lowerPeriod = 2,
#'pad = NULL,
#'padding = "noise")
#'
#'#Example 3. Greyscale data set of Zeeden et al., (2013)
#'grey_eha_log2 <- eha_log2(data = grey,
#'win =  20,
#'tbw = 4,
#'demean = TRUE,
#'detrend = TRUE,
#'upperPeriod = 20,
#'lowerPeriod = 2,
#'pad = NULL,
#'padding = "noise")
#'
#'}
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats acf
#' @importFrom astrochron eha


eha_log2 <- function(data = NULL,
    win =  NULL,
    tbw = NULL,
    demean = NULL,
    detrend = NULL,
    upperPeriod = NULL,
    lowerPeriod = NULL,
    pad = NULL,
    padding = "noise") {

  half_win <- (win / 2)
  dz <- mean(diff(data[, 1]), na.rm = TRUE)
  z <- data[, 1]
  x <- data[, 2]


  # Lower padding (before start)
  z_pad_low <- seq(from = z[1] - half_win,
      to   = z[1] - dz,
      by   = dz)

  # Upper padding (after end)
  z_pad_high <- seq(from = z[length(z)] + dz,
      to   = z[length(z)] + half_win,
      by   = dz)

  if (padding == "mean" | padding == "average" | padding == "avg") {
    x_mean <- mean(x, na.rm = TRUE)
    x_pad_low  <- rep(x_mean, length(z_pad_low))
    x_pad_high <- rep(x_mean, length(z_pad_high))
  } else if (padding == "noise") {
    timesteps_data <- length(z_pad_low)
    set.seed(timesteps_data)
    phi_x_up <- wr_autocorrelation(data[1:timesteps_data, 2])
    if (phi_x_up >= 1) {
      phi_x_up <- 0.99
    }
    x_pad_low <- wr_coloured_noise(
        timesteps = timesteps_data,
        mean = mean(data[1:timesteps_data, 2]),
        sd = sd(data[1:timesteps_data, 2]),
        phi = phi_x_up
    )


    wr_autocorrelation(data[(nrow(data) - timesteps_data):nrow(data), 2])

    phi_x_down <- wr_autocorrelation(data[(nrow(data) - timesteps_data):nrow(data), 2])
    if (phi_x_down >= 1) {
      phi_x_down <- 0.99
    }
    x_pad_high <- wr_coloured_noise(
        timesteps = timesteps_data,
        mean = mean(data[(nrow(data) - timesteps_data):nrow(data), 2]),
        sd = sd(data[(nrow(data) - timesteps_data):nrow(data), 2]),
        phi = phi_x_down
    )


  } else if (padding == "zero") {
    x_pad_low  <- rep(0, length(z_pad_low))
    x_pad_high <- rep(0, length(z_pad_high))
  } else{
    x_mean <- mean(x, na.rm = TRUE)
    x_pad_low  <- rep(x_mean, length(z_pad_low))
    x_pad_high <- rep(x_mean, length(z_pad_high))
  }

  pad_low  <- data.frame(depth = z_pad_low, value = x_pad_low)
  names(pad_low) <- names(data)
  pad_high <- data.frame(depth = z_pad_high, value = x_pad_high)
  names(pad_high) <- names(data)

  data_padded <- rbind(pad_low, data, pad_high)
  data_padded <- linterp(data_padded, verbose = FALSE, genplot = FALSE)

  winpts = as.integer(floor(win / dz) + 1)
  defaultPad = 2 * 2^ceiling(log2(winpts))
  if (length(pad) == 0) {
    pad <- defaultPad
  }
  if ((pad) %% 2 != 0) {
    pad = pad + 1
  }
  newpts = as.integer(pad)
  if (pad < winpts) {
    pad <- winpts
  }


  eha_res <- astrochron::eha(
      data_padded,
      tbw = tbw,
      pad = pad*5, #factor 5 is needed to ensure low frequency spacing
      fmin = 1 / upperPeriod,
      fmax = 1 / lowerPeriod,
      win = win,
      demean = T,
      detrend = T,
      siglevel = 0.90,
      sigID = F,
      ydir = 1,
      output = 1,
      pl = 1,
      palette = 6,
      centerZero = T,
      ncolors = 100,
      genplot = 0,
      verbose = FALSE
  )

  depth_x <- as.numeric(sub("^X", "", names(eha_res[[1]])[grepl("^X[0-9]", names(eha_res[[1]]))]))

  power <- eha_res[[1]]
  amplitude <- eha_res[[2]]
  prob <- eha_res[[3]]
  f_test <- eha_res[[4]]

  amp_mat <- as.matrix(amplitude[, -1])
  pwr_mat <- as.matrix(power[, -1])
  prob_mat <- as.matrix(prob[, -1])
  f_mat <- as.matrix(f_test[, -1])

  freq_y <- eha_res[[2]]$freq

  log2_period_y <- log2(1 / freq_y)

  storage.mode(amp_mat) <- "numeric"
  eha_img <- list(depth = depth_x,
      log2_period = log2_period_y,
      amp = amp_mat)

  n_depth <- length(depth_x)

  depth_grid <- seq(
      from = min(depth_x),
      to   = max(depth_x),
      length.out = n_depth
  )

  n_per <- 4 * length(log2_period_y)

  log2p_grid <- seq(
      from = min(log2(lowerPeriod)),
      to   = max(log2(upperPeriod)),
      length.out = n_per
  )


  #interpolate amplitude
  # interpolate along period
  amp_tmp <- matrix(NA_real_,
      nrow = length(log2p_grid),
      ncol = length(depth_x))

  for (i in seq_along(depth_x)) {
    amp_tmp[, i] <- approx(
        x = log2_period_y,
        y = amp_mat[, i],
        xout = log2p_grid,
        rule = 2
    )$y
  }


  # interpolate along depth
  amp_reg <- matrix(NA_real_,
      nrow = length(log2p_grid),
      ncol = length(depth_grid))

  for (j in seq_along(log2p_grid)) {
    amp_reg[j, ] <- approx(
        x = depth_x,
        y = amp_tmp[j, ],
        xout = depth_grid,
        rule = 2
    )$y
  }



  #interpolate power
  # interpolate along period
  pwr_tmp <- matrix(NA_real_,
      nrow = length(log2p_grid),
      ncol = length(depth_x))

  for (i in seq_along(depth_x)) {
    pwr_tmp[, i] <- approx(
        x = log2_period_y,
        y = pwr_mat[, i],
        xout = log2p_grid,
        rule = 2
    )$y
  }


  # interpolate along depth
  pwr_reg <- matrix(NA_real_,
      nrow = length(log2p_grid),
      ncol = length(depth_grid))

  for (j in seq_along(log2p_grid)) {
    pwr_reg[j, ] <- approx(
        x = depth_x,
        y = pwr_tmp[j, ],
        xout = depth_grid,
        rule = 2
    )$y
  }






  #interpolate probabillity
  # interpolate along period
  prob_tmp <- matrix(NA_real_,
      nrow = length(log2p_grid),
      ncol = length(depth_x))

  for (i in seq_along(depth_x)) {
    prob_tmp[, i] <- approx(
        x = log2_period_y,
        y = prob_mat[, i],
        xout = log2p_grid,
        rule = 2
    )$y
  }


  # interpolate along depth
  prob_reg <- matrix(NA_real_,
      nrow = length(log2p_grid),
      ncol = length(depth_grid))

  for (j in seq_along(log2p_grid)) {
    prob_reg[j, ] <- approx(
        x = depth_x,
        y = prob_tmp[j, ],
        xout = depth_grid,
        rule = 2
    )$y
  }






  #interpolate probabillity
  # interpolate along period
  f_tmp <- matrix(NA_real_,
      nrow = length(log2p_grid),
      ncol = length(depth_x))

  for (i in seq_along(depth_x)) {
    f_tmp[, i] <- approx(
        x = log2_period_y,
        y = f_mat[, i],
        xout = log2p_grid,
        rule = 2
    )$y
  }


  # interpolate along depth
  f_reg <- matrix(NA_real_,
      nrow = length(log2p_grid),
      ncol = length(depth_grid))

  for (j in seq_along(log2p_grid)) {
    f_reg[j, ] <- approx(
        x = depth_x,
        y = f_tmp[j, ],
        xout = depth_grid,
        rule = 2
    )$y
  }
  Power.avg <- rowMeans(pwr_reg)
  amp.avg <- rowMeans(amp_reg)
  prob.avg <- rowMeans(prob_reg)
  f_test.avg <- rowMeans(f_reg)

  output <- list(
      depth = depth_grid,
      log2_period = log2p_grid,
      pwr = pwr_reg,
      amp = amp_reg,
      prob = prob_reg,
      f_test = f_reg,
      Power.avg = Power.avg,
      amp.avg = amp.avg,
      prob.avg = prob.avg,
      f_test.avg = f_test.avg,
      x= data[,1],
      y= data[,2],
      win=win
  )

  class(output) <- "analyze.eha"
  return(invisible(output))
}
