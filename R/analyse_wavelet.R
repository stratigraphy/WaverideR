#' @title Conduct the continuous wavelet transform on a time series/signal
#'
#' @description
#' Compute the continuous wavelet transform (CWT) using a Morlet wavelet
#'
#' @param data Input record: a matrix or data frame with depth/time in column 1 and the proxy value in column 2.
#' Rows with \code{NA} are removed, the rows are sorted by depth/time and the record is linearly interpolated to
#' an even spacing equal to the median sampling interval. Required (\code{Default=NULL} gives an error).
#' @param dj Spacing between successive scales in octaves (numeric, > 0). The periods analysed increase by powers
#' of 2 (1, 2, 4, 8, ...), and \code{dj} sets how many steps are made per doubling of the period: \code{dj=1/x}
#' gives x periods per octave. Smaller values give a finer period axis but increase the computation time.
#' \code{Default=1/100}.
#' @param lowerPeriod Shortest period to analyse, in the same unit as the depth/time column of \code{data}
#' (numeric, > 0; at least twice the sampling interval). Together with \code{upperPeriod} it sets the range of
#' cycles analysed; because the period axis is log2-scaled, a power of 2 gives the neatest plots. \code{Default=2}.
#' @param upperPeriod Longest period to analyse, in the same unit as the depth/time column of \code{data}
#' (numeric, > \code{lowerPeriod}; normally not more than the length of the record); a power of 2 gives the
#' neatest plots. \code{Default=1024}.
#' @param verbose Print the sampling interval to which the record was interpolated and, when \code{pval=TRUE},
#' show a progress bar of the simulations (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param omega_nr Number of cycles of the Morlet wavelet, i.e. its central angular frequency omega_0 (numeric,
#' > 0, typically 5 to 15). A larger value gives a longer wavelet with more oscillations: a better period
#' resolution and a worse depth/time resolution. Also used as the default parameter of \code{"morlet0"} and to
#' match the bandwidth of some other wavelets. \code{Default=8}.
#' @param wavelet Mother wavelet (character string or function), given as a name (PyWavelets style names with the
#' order appended, such as
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
#' @param wavelet_order Parameter(s) of the wavelet (numeric scalar or vector; \code{NULL} uses the defaults
#' listed here): \code{paul}: order m (default 4); \code{dog}: derivative m
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
#' The period axis is always the Fourier period, i.e. the period of the sinusoid at which the wavelet has its maximum
#' response, so the results of different wavelets are comparable. \code{Default=NULL}.
#' @param pval Calculate p-values of the wavelet power by Monte Carlo simulation (logical, TRUE or FALSE).
#' \code{Default=FALSE}. Each simulation analyses a surrogate record made of the average of a shuffled copy of
#' the data (white noise) and an AR(1) red-noise series with the mean and standard deviation of the data and a
#' lag-1 autocorrelation drawn around the value estimated in a moving window of half the record length. The
#' shuffled part keeps high-amplitude short-period cycles from becoming significant; the red-noise part gives
#' a spectrum similar to that of the data without being biased by possible astronomical peaks in it. The
#' p-value is the fraction of simulations whose power is at least the observed power (returned as
#' \code{Power.pval} and \code{Power.avg.pval}).
#' @param n_simulations Number of Monte Carlo simulations used for the p-values (integer, > 0); only used when
#' \code{pval=TRUE}. More simulations give more precise p-values but take longer. \code{Default=10}.
#' @param run_multicore Run the Monte Carlo simulations in parallel on several cores (logical, TRUE or FALSE);
#' only used when \code{pval=TRUE}. \code{Default=FALSE}.
#' @param n_cores Number of cores used when \code{run_multicore=TRUE} (integer, > 0); \code{NULL} picks an
#' automatic, safe number: physical cores minus one, at most 4, and limited by the free memory. Set
#' \code{options(waverider.n_cores = ...)} to change that default. A value > 1 also switches on parallel
#' computing. Results do not depend on the number of cores. \code{Default=NULL}.
#' @return
#' The output is a list (wavelet object) which contain 20 objects which are the result of the continuous wavelet transform (CWT).
#'Object 1: Wave - Wave values of the wavelet
#'Object 2: Phase - Phase of the wavelet
#'Object 3: Ampl - Amplitude values of the wavelet
#'Object 4: Power - Power values of the wavelet
#'Object 5: dt - Step size
#'Object 6: dj - Scale size
#'Object 7: Power.avg  - Average power values
#'Object 8: Period - Period values
#'Object 9: Scale - Scale value
#'Object 10: coi.1 - Cone of influence values 1
#'Object 11: coi.2 - Cone of influence values 2
#'Object 12: nc - Number of columns
#'Object 13: nr - Number of rows
#'Object 14: axis.1 - axis values 1
#'Object 15: axis.2 - axis values 2
#'Object 16: omega_nr - Number of cycles in the wavelet
#'Object 17: x - x values of the data set
#'Object 18: y - y values of the data set
#'Object 19: average p value of the spectral power
#'Object 20: p value of spectral power
#'
#' @author
#' Code based on on the "WaveletComp" function of the 'WaveletComp' R package
#' and "wt" function of the 'biwavelet' R package which are based on the
#' wavelet MATLAB code written by Christopher Torrence and Gibert P. Compo.
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
#' @examples
#' \donttest{
#'#Example 1. Using the Total Solar Irradiance data set of Steinhilver et al., (2012)
#'TSI_wt <-
#'  analyse_wavelet(
#'    data = TSI,
#'    dj = 1/200,
#'    lowerPeriod = 16,
#'    upperPeriod = 8192,
#'    verbose = FALSE,
#'    omega_nr = 6,
#'    pval=FALSE,
#'    n_simulations=10,
#'    run_multicore = FALSE
#'  )
#'
#'
#'#Example 2. Using the magnetic susceptibility data set of Pas et al., (2018)
#'mag_wt <-
#'analyse_wavelet(
#'data = mag,
#'dj = 1/100,
#'lowerPeriod = 0.1,
#'upperPeriod = 254,
#'verbose = FALSE,
#'omega_nr = 10,
#'pval=FALSE,
#'n_simulations=10,
#'run_multicore = FALSE
#')
#'
#'#Example 3. Using the greyscale data set of Zeeden et al., (2013)
#'grey_wt <-
#'  analyse_wavelet(
#'    data = grey,
#'    dj = 1/200,
#'    lowerPeriod = 0.02,
#'    upperPeriod = 256,
#'    verbose = FALSE,
#'    omega_nr = 8,
#'    pval=FALSE,
#'    n_simulations=10,
#'    run_multicore = FALSE
#'  )
#'
#'}
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats acf

analyse_wavelet <-
    function(data = NULL,
        dj = 1 / 100,
        lowerPeriod = 2,
        upperPeriod = 1024,
        verbose = FALSE,
        omega_nr = 8,
        wavelet = "morlet",
        wavelet_order = NULL,
        pval=FALSE,
        n_simulations=10,
        run_multicore = FALSE, n_cores = NULL) {

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
      dat <- as.data.frame(interp)

      if (verbose == TRUE) {
        cat("dataset interpolated to: ", round(dt, 10))
      }

      x_axis <- dat[, 1]
      x <- dat[, 2]
      x_data <- x
      # Original length and length of zero padding:
      series.length = length(x)
      pot2 = trunc(log2(series.length) + 0.5)
      pad.length = 2 ^ (pot2 + 1) - series.length

      # Define central angular frequency omega0 and fourier factor:
      wdef <- .wr_wavelet(wavelet, wavelet_order, omega_nr)
      omega0 = if (is.numeric(wdef$order)) { wdef$order[1] } else { wdef$peak_omega }
      fourier.factor = wdef$fourier_factor

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

      ## Define the Morlet wavelet transform function

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

        daughter = wdef$daughter(my.scale, omega.k, dt)

        wave[ind.scale,] = fft(fft.xpad * daughter, inverse = TRUE) / N
      }

      # Cut out the wavelet transform
      Wave = wave[, 1:series.length]

      # Compute wavelet power
      Power = Mod(Wave) ^ 2 / matrix(rep(scales, series.length), nrow = scales.length)

      # Phase
      Phase = Arg(Wave)

      # Amplitude
      Ampl  = Mod(Wave) / matrix(rep(sqrt(scales), series.length), nrow =
          scales.length)

      Power.avg = rowMeans(as.matrix(Power))

      axis.1 = x_axis
      axis.2 = log2(periods)
      fourier.factor = wdef$fourier_factor
      coi = fourier.factor * .wr_coi_factor(wdef) * dt * c(1E-5, 1:((series.length + 1) /
          2 - 1), rev((1:(
          series.length / 2 - 1
      ))), 1E-5)
      coi.x = c(axis.1[c(1, 1)] - dt * 0.5, axis.1, axis.1[c(series.length, series.length)] + dt *
          0.5)
      logyint = axis.2[2] - axis.2[1]
      yl = c(log2(periods[scales.length]) + 0.5 * logyint, log2(periods[1]) - 0.5 *
          logyint)
      yr = rev(yl)
      coi.y = c(yl, log2(coi), yr)
      x <- dat[, 1]
      y <- dat[, 2]

      Power.pval <- NULL
      Power.avg.pval <- NULL

      if(pval==TRUE){

        Power.pval <- NULL
        Power.avg.pval <- NULL

        # data_test <-
        # wr_coloured_noise(
        #   timesteps = timesteps_data,
        #   mean = 500,
        #   sd = 50,
        #   phi=0.1
        # )

        #calc sd_phi

        # Define window size (half the length of the data)
        window_size <- floor(length(data[, 2]) / 2)
        n <- length(data[, 2])

        # Initialize result vector
        acf_vals <- numeric(n - window_size + 1)

        # Loop over the data with a sliding window
        for (i in 1:(n - window_size + 1)) {
          window_data <- data[i:(i + window_size - 1), 2]
          acf_result <- acf(window_data, lag.max = 1, plot = FALSE, type = "correlation")
          acf_vals[i] <- acf_result$acf[2]  # lag 1 is at position 2
        }

        # result <- zoo::rollapply(data[,2], width = length(data[,2])/2, FUN=acf,
        #                          lag.max = 1,type = "correlation",plot = FALSE)
        #
        #
        # acf_result <-   matrix(data=0,nrow= length(result),ncol=1)

        # for (i in 1:length(result)){
        #   acf_res <- unlist(result[i])
        #
        #   if( is.null(acf_res) != TRUE){
        #     acf_result[i,1]  <- as.numeric(acf_res[2])
        #   }else{acf_result[i,1]  <- NA}
        #
        # }
        acf_result<- na.omit(acf_vals)
        acf_result <- acf_result[acf_result<1]

        phi_data <- mean(acf_result)
        sd_data_2 <-  sd(acf_result)
        sd_data <- sd(data[,2])
        timesteps_data <- length(data[, 2])
        mean_data <- mean(data[, 2])

        fits <- NULL

        # rough memory needed by one simulation (complex wavelet matrix scales x N and a few copies of it),
        # used to limit the number of cores to what fits in the free memory
        est_scales <- floor(log2(upperPeriod / lowerPeriod) / dj) + 1
        est_N <- 2 ^ (trunc(log2(nrow(data)) + 0.5) + 1)
        est_mem <- 6 * est_scales * est_N * 16

        fits <- wr_lapply(seq_len(n_simulations), function(i) {
          #i <- 1
          set.seed(i)

          x_1 <- sample(x_data, length(x_data))

          if(phi_data+sd_data_2 > 1){
            max_phi <- 1-(1/10^9)
            min_phi <- phi_data-(1-phi_data)
          }else{
            max_phi <- phi_data+sd_data_2
            min_phi <- phi_data-sd_data_2

          }
          #?colorednoise::colored_noise
          phi_data_2 <- wr_rtruncnorm(n=1,a=1/100000,b=1-(1/10^9),mean=phi_data,sd=sd_data_2) #old option
          x_2 <-
              (
                  wr_coloured_noise(
                      timesteps = timesteps_data,
                      mean = mean_data,
                      sd = sd_data,
                      phi=phi_data_2
                  )
              )

          # n <- length(x_data)
          # z <- fft(x_data)
          # if (n%%2 == 0) {
          #   ph <- 2 * pi * runif(n/2 - 1)
          #   ph <- c(0, ph, 0, -rev(ph))
          # }
          # if (n%%2 != 0) {
          #   ph <- 2 * pi * runif((n - 1)/2)
          #   ph <- c(0, ph, -rev(ph))
          #   }
          # ph <- complex(imaginary = ph)
          # z <- z * exp(ph)
          # x_3 <- Re(fft(z, inverse = TRUE)/n)
          #
          # data_2 <- (x_2+x_1+x_3)/3
          #
          #

          data_2 <- (x_2+x_1)/2

          dat <- as.matrix(cbind(x_axis,data_2))
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
          dat <- as.data.frame(interp)

          x_axis <- dat[, 1]
          x <- dat[, 2]

          # Original length and length of zero padding:
          series.length = length(x)
          pot2 = trunc(log2(series.length) + 0.5)
          pad.length = 2 ^ (pot2 + 1) - series.length

          # Define central angular frequency omega0 and fourier factor:
          wdef <- .wr_wavelet(wavelet, wavelet_order, omega_nr)
          omega0 = if (is.numeric(wdef$order)) { wdef$order[1] } else { wdef$peak_omega }
          fourier.factor = wdef$fourier_factor

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

          ## Define the Morlet wavelet transform function

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

            daughter = wdef$daughter(my.scale, omega.k, dt)

            wave[ind.scale,] = fft(fft.xpad * daughter, inverse = TRUE) / N
          }

          # Cut out the wavelet transform
          Wave = wave[, 1:series.length]

          # Compute wavelet power
          Power.sim = Mod(Wave) ^ 2 / matrix(rep(scales, series.length), nrow = scales.length)

          Power.avg.sim = rowMeans(as.matrix(Power.sim))

          Power.pval = matrix(0, nrow = scales.length, ncol = nrow(dat))
          Power.avg.pval = rep(0, scales.length)

          Power.avg.sim = rowMeans(Power.sim)

          Power.pval[Power.sim >= Power] = Power.pval[Power.sim >=
              Power] + 1

          Power.avg.pval[Power.avg.sim >= Power.avg] = Power.avg.pval[Power.avg.sim >=
              Power.avg] + 1

          #plot(Power.avg)
          #lines(Power.avg.sim)
          #lines(Power.avg.pval)
          #Power.pval <- Power.sim
          #Power.avg.pval <- Power.avg.sim
          return(cbind(Power.pval,Power.avg.pval))

        }, n_cores = n_cores, run_multicore = run_multicore, combine = "sum", progress = verbose,
            mem_per_worker = est_mem)

        Power.pval = (fits[,1:(ncol(fits)-1)])/n_simulations
        Power.avg.pval <- fits[,ncol(fits)]/n_simulations
        #Power.avg.pval <- (rowMeans(Power.pval)+Power.avg.pval)/2

      }

      output <-
          list(
              Wave = Wave,
              Phase = Phase,
              Ampl = Ampl,
              Power = Power,
              dt = dt,
              dj = dj,
              Power.avg = Power.avg,
              Period = periods,
              Scale = scales,
              coi.1 = coi.x,
              coi.2 = coi.y,
              nc = series.length,
              nr = scales.length,
              axis.1 = axis.1,
              axis.2 = axis.2,
              omega_nr = omega_nr,
              wavelet = wdef$name,
              wavelet_order = wdef$order,
              fourier_factor = wdef$fourier_factor,
              psi0 = wdef$psi0,
              x = x,
              y = y,
              Power.avg.pval = Power.avg.pval,
              Power.pval = Power.pval
          )

      class(output) = "analyze.wavelet"

      return(invisible(output))

    }
