#' @title Models average spectral power based curves based on a red-noise signal
#' generated using the characteristics of an input signal.
#'
#' @description The \code{\link{model_red_noise_wt}} function is used to generate
#' average spectral power curves based on and input signal and set wavelet settings.
#'
#'@param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#'@param n_simulations Number of red noise simulations.
#'@param run_multicore run simulation using multiple cores \code{Default=FALSE}
#'the simulation is run at x-2 cores to allow the 2 remaining processes to run background processes.
#'@param verbose Print text \code{Default=FALSE}.
#'
#' @author
#' Code based on the "analyze.wavelet" function of the 'WaveletComp' R package
#' and "wt" function of the 'biwavelet' R package which are based on the
#' wavelet 'MATLAB' code written by Christopher Torrence and Gibert P. Compo (1998).
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
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#'  Geophysics 1982 47 (2): 222–236.
#'
#'@examples
#'\donttest{
#'#'#generate average spectral power curves based on red noise curves which are
#'# based on the magnetic susceptibility record of the Sullivan core of Pas et al., (2018)
#'
#'mag_wt <- analyze_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'#increase n_simulations to better define the red noise spectral power curve
#'mag_wt_red_noise <- model_red_noise_wt(wavelet=mag_wt,
#'n_simulations=10, # increase number for better constrained results
#'run_multicore=FALSE,
#'verbose=FALSE)
#'}
#'
#'
#'@return
#'Returns a matrix in which each column represents the average spectral
#'power resulting from a red-noise run.
#'
#' @export
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom foreach foreach
#' @importFrom colorednoise colored_noise
#' @importFrom colorednoise autocorrelation
#' @importFrom stats runif
#' @importFrom stats fft
#' @importFrom stats sd
#' @importFrom foreach %dopar%
#' @importFrom parallel stopCluster
#' @importFrom truncnorm  rtruncnorm


model_red_noise_wt <- function(wavelet = NULL,
                               n_simulations = NULL,
                               run_multicore = FALSE,
                               verbose = FALSE) {
  data <- cbind(wavelet$x, wavelet$y)
  phi_data <-  autocorrelation(data[, 2])
  lowerPeriod <- min(wavelet$Period)
  upperPeriod <- max(wavelet$Period)
  omega_nr <- wavelet$omega_nr
  dj <- wavelet$dj

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
  data <- as.data.frame(interp)

  if (verbose == TRUE) {
    cat("dataset interpolated to: ", round(dt, 10))
  }

  sd_data <- sd(data[, 2])
  timesteps_data <- nrow(data)
  mean_data <- mean(data[, 2])


  if (run_multicore == TRUE) {
    numCores <- detectCores()
    cl <- makeCluster(numCores - 2)
    registerDoSNOW(cl)
  }
  else{
    numCores <- 1
    cl <- makeCluster(numCores)
    registerDoSNOW(cl)
  }

  if (verbose==TRUE){
    pb <- txtProgressBar(max = n_simulations, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)}else{opts=NULL}




  # Original length and length of zero padding:
  series.length = nrow(data)
  pot2 = trunc(log2(series.length) + 0.5)
  pad.length = 2 ^ (pot2 + 1) - series.length

  # Define central angular frequency omega0 and fourier factor:
  fourier.factor = (2 * pi) / omega_nr

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
  omega0 <- omega_nr



  red_noise <-
    foreach (
      j = 1:n_simulations,
      .combine = 'cbind',
      .options.parallel   = opts
    ) %dopar% {
      phi_data_2 <- truncnorm::rtruncnorm(n=1,a=1/100000,b=1-(1/100000),mean=phi_data,sd=0.5) #old option
      x <-
        (
          colorednoise::colored_noise(
            timesteps = timesteps_data,
            mean = mean_data,
            sd = sd_data,
            phi_data_2
          )
        )


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
      Power.avg = rowMeans(as.matrix(Mod(Wave) ^ 2 / matrix(rep(
        scales, series.length
      ), nrow = scales.length)))

      noise_period <- Power.avg

    }

  stopCluster(cl)
  return(red_noise)
}
