#' @title Generate a log2 spaced EHA spectra
#'
#' @description
#' Compute a log2 spaced Evolutive Harmonic Analysis (EHA) & Evolutive Power
#' Spectral Analysis
#'This is a wrapper function for the "eha" function of the 'astrochron' R package
#' @param data Input data, should be a matrix or data frame in which
#' the first column is depth or time and the second column is proxy record.
#' @param data Window size for EHA, in units of space or time.
#' @param demean	Remove mean from data series? (T or F)
#' @param detrend	Remove linear trend from data series? (T or F)
#' @param tbw  MTM time-bandwidth product (<=10)
#' @param upperPeriod Upper period to be analyzed.
#' The CWT analyses the signal starting from the lowerPeriod to the upperPeriod
#' so the proper selection these
#' parameters allows to analyze the signal for a specific range of cycles.
#'  scaling is done using power 2 so for the best plotting results select a
#'  value to the power or 2.
#' @param lowerPeriod  Lowest period to be analyzed.
#' The CWT analyses the signal starting from the lowerPeriod to the upperPeriod
#' so the proper selection these
#' parameters allows to analyze the signal for a specific range of cycles.
#' scaling is done using power 2 so for the best plotting results select a value
#' to the power or 2.
#' @param Pad with zeros to how many points? Must not factor into a prime number
#' >23. Maximum number of points is 200,000.
#' @param padding pad the edges of the data set with half a
#' window length with the following, the  "Mean", "noise" or "zero"
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
#'Object 2: y - y values of the data set
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
#'#Example 1. A plot of a wavelet spectra using the Total Solar Irradiance
#'# data set of Steinhilber et al., (2012)
#'
#'TSI_eha_log2 <- function(data = TSI,
#'win =  8192,
#'tbw = 4,
#'demean = TRUE,
#'detrend = TRUE,
#'upperPeriod = 8192,
#'lowerPeriod = c,
#'pad = NULL,
#'padding = "noise")
#'
#'#Example 2. A plot of a wavelet spectra using the magnetic susceptibility
#'#data set of Pas et al., (2018)
#'mag_eha_log2 <- function(data = mag,
#'win =  254,
#'tbw = 4,
#'demean = TRUE,
#'detrend = TRUE,
#'upperPeriod = 254,
#'lowerPeriod = 0.1,
#'pad = NULL,
#'padding = "noise")
#'
#'#Example 3. A plot of a wavelet spectra using the greyscale
#'# data set of Zeeden et al., (2013)
#'grey_eha_log2 <- function(data = grey,
#'win =  256,
#'tbw = 4,
#'demean = TRUE,
#'detrend = TRUE,
#'upperPeriod = 256,
#'lowerPeriod = 0.02,
#'pad = NULL,
#'padding = "noise")
#'
#'}
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel stopCluster
#' @importFrom stats acf


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
    phi_x_up <- colorednoise::autocorrelation(data[1:timesteps_data, 2])
    if (phi_x_up >= 1) {
      phi_x_up <- 0.99
    }
    x_pad_low <- colorednoise::colored_noise(
      timesteps = timesteps_data,
      mean = mean(data[1:timesteps_data, 2]),
      sd = sd(data[1:timesteps_data, 2]),
      phi = phi_x_up
    )
    phi_x_down <- colorednoise::autocorrelation(data[(nrow(data) - timesteps_data):nrow(data), 2])
    if (phi_x_down >= 1) {
      phi_x_down <- 0.99
    }
    x_pad_high <- colorednoise::colored_noise(
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
  if ((pad) %% 2 != 0)
    pad = pad + 1
  newpts = as.integer(pad)
  if (pad < winpts) {
    pad <- winpts
  }


  eha_res <- eha(
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
    y= data[,2]
  )

  class(output) <- "analyze.eha"
  return(invisible(output))
}

