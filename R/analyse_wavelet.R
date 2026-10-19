#' @title Conduct the continuous wavelet transform on a time series/signal
#'
#' @description
#' Compute the continuous wavelet transform (CWT) using a Morlet wavelet
#'
#' @param data Input data, should be a matrix or data frame in which
#' the first column is depth or time and the second column is proxy record.
#' @param dj Spacing between successive scales. The CWT analyses analyses the signal using successive periods
#' which increase by the power of 2 (e.g.2^0=1,2^1=2,2^2=4,2^3=8,2^4=16). To have more resolution
#' in-between these steps the dj parameter exists, the dj parameter specifies how many extra steps/spacing in-between
#' the power of 2 scaled CWT is added. The amount of steps is 1/x with a higher x indicating a smaller spacing.
#' Increasing the increases the computational time of the CWT \code{Default=1/200}.
#' @param lowerPeriod  Lowest period to be analysed \code{Default=2}.
#' The CWT analyses the signal starting from the lowerPeriod to the upperPeriod so the proper selection these
#' parameters allows to analyse the signal for a specific range of cycles.
#' scaling is done using power 2 so for the best plotting results select a value to the power or 2.
#' @param upperPeriod Upper period to be analysed \code{Default=1024}.
#' The CWT analyses the signal starting from the lowerPeriod to the upperPeriod so the proper selection these
#' parameters allows to analyse the signal for a specific range of cycles.
#'  scaling is done using power 2 so for the best plotting results select a value to the power or 2.
#' @param verbose Print text \code{Default=FALSE}.
#' @param omega_nr The omega number: central angular frequency of the Morlet wavelet. A larger value gives a longer wavelet with more oscillations, a better frequency resolution and a worse resolution in depth/time
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
#' The period axis is always the Fourier period, i.e. the period of the sinusoid at which the wavelet has its maximum
#' response, so the results of different wavelets are comparable.
#' @param pval calculate the P-value  \code{Default=FALSE}. The p-value is based on
#' Monte Carlo modelling runs on surrogate data generated based on autocorrelated noise (red noise) the calculated using a windowed
#' (the window is half the size of the data set) temporal autocorrelation
#' and on shuffling the data set resulting in a random data sets which has similar spectral characteristics
#' to the original data set.The shuffling of the data set creates white noise which ensures that high amplitude high frequency/short
#' period cycles do not result in statistical significant peaks. The part of the data generated using the  autocorrelated noise (red noise)
#' based on the windowed  (the window is half the size of the data set) temporal autocorrelation represent a spectral signature similar to
#' to that of the original data. The original data might include spectral peaks which are the result of astronomical
#' forcing. The result is that the spectral power profile is biased towards rejecting the 0-hypothesis (e.g. no astronomical forcing).
#' By combining the shuffling of the data set with autocorrelated noise a surrogate data set is created which rejects
#' high amplitude high frequency/short period cycles and a reduced biased towards towards rejecting the 0-hypothesis if the data was
#' solely the result of autocorrelated noise
#' @param n_simulations Number of simulation to be ran to generate the p-value
#' @param run_multicore Run p-value calculation with one core or multiple cores
#'
#'
#'
#' @param n_cores Number of cores used when `run_multicore = TRUE` (`Default = NULL`: an
#'   automatic, safe number: physical cores minus one, at most 4, and limited by the free memory).
#'   Set `options(waverider.n_cores = ...)` to change the default. A value > 1 also switches on
#'   parallel computing. Results do not depend on the number of cores.
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
    omega0 = if (is.numeric(wdef$order)) wdef$order[1] else wdef$peak_omega
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
        }else{max_phi <- phi_data+sd_data_2
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
        omega0 = if (is.numeric(wdef$order)) wdef$order[1] else wdef$peak_omega
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

