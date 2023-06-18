#' @title Extract amplitude from a signal
#'
#' @description Extracts the amplitude from a signal using the continuous wavelet transform using a Morlet wavelet.
#' The extraction of the amplitude is useful for cyclostratigraphic studies because the amplitude of
#' an astronomical cycle is modulated by higher order astronomical cycles.
#'
#'@param signal Input signal from which the amplitude is extracted any signal in which the first column is
#'depth/time and the second column is the proxy record from which the amplitude is extracted
#'@param pts The pts parameter specifies how many points to the left/right up/down the peak detect algorithm goes in detecting
#'a peak. The peak detecting algorithm works by comparing the values left/right up/down of it, if the values are both higher or lower
#'then the value a peak. To deal with error produced by this algorithm the pts parameter can be changed which can
#'aid in peak detection. Usually increasing the pts parameter means more peak certainty, however it also means that minor peaks might not be
#'picked up by the algorithm \code{Default=3}
#'@param genplot If set to TRUE a plot with extracted amplitude will be displayed \code{Default=FALSE}.
#'@param remean Prior to analysis the mean is subtracted from the data set to re-mean set \code{Default=TRUE}.
#'@param ver_results To verify the amplitude extraction is representative of the amplitude
#'extracted using the \code{\link{extract_amplitude}} function the results can be compared to the amplitude extracted
#'using the \code{\link{Hilbert_transform}} if the mean difference is more then 5% one might need to reconsider
#'whether the input contains a reliable enough signal with high a enough amplitude modulation to actually extract an amplitude from. \code{Default=FALSE}.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#' @author
#' Code based on the \link[WaveletComp]{reconstruct} function of the 'WaveletComp' R package
#' which is based on the wavelet 'MATLAB' code written by Christopher Torrence and Gibert P. Compo.
#' The assignment of the standard deviation of the uncertainty of the wavelet
#' is based on the work of Gabor (1946) and Russell et al., (2016)
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
#' <\doi{doi:10.1190/1.1441328}>
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#'  Geophysics 1982 47 (2): 222–236. <\doi{doi:10.1190/1.1441329}>
#'
#'@examples
#'\donttest{
#'#Extract amplitude of the 405 kyr eccentricity cycle from the the magnetic
#'# susceptibility data set of De pas et al., (2018)
#'#Perform the CWT on the magnetic susceptibility data set of Pas et al., (2018)
#'
#'mag_wt <- analyze_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'#Track the 405 kyr eccentricity cycle in a wavelet spectra
#'
#'#mag_track <- track_period_wavelet(astro_cycle = 405,
#'#                                   wavelet=mag_wt,
#'#                                   n.levels = 100,
#'#                                   periodlab = "Period (metres)",
#'#                                   x_lab = "depth (metres)")
#'
#'#Instead of tracking, the tracked solution data set mag_track_solution
#'#is used
#'mag_track <- mag_track_solution
#'
#' mag_track_complete <- completed_series(
#'   wavelet = mag_wt,
#'   tracked_curve = mag_track,
#'   period_up = 1.2,
#'   period_down = 0.8,
#'   extrapolate = TRUE,
#'   genplot = FALSE
#' )
#'
#'#Smooth the completed tracking of the 405 kyr eccentricity cycle in the wavelet spectra
#'
#'mag_track_complete <- loess_auto(time_series = mag_track_complete,
#' genplot = FALSE, print_span = FALSE)
#'
#'
#'
#'
#'mag_405_ecc  <- extract_signal(
#'tracked_cycle_curve = mag_track_complete,
#'wavelet = mag_wt,
#'period_up = 1.2,
#'period_down = 0.8,
#'add_mean = TRUE,
#'tracked_cycle_period = 405,
#'extract_cycle = 405,
#'tune = FALSE,
#'plot_residual = FALSE
#')
#'
#'#extract the amplitude  of the 405 kyr eccentricity cycle
#'mag_ampl <- extract_amplitude(
#'signal = mag_405_ecc,
#'pts=3,
#'genplot = FALSE,
#'ver_results = FALSE,
#'keep_editable=FALSE)
#'
#'}
#'@return
#'Returns a matrix with 2 columns.
#'The first column is depth/time.
#'The second column is the extracted amplitude
#'
#' @export
#' @importFrom WaveletComp reconstruct

extract_amplitude <- function(signal = NULL,
                              pts=3,
                              genplot = FALSE,
                              remean = TRUE,
                              ver_results = FALSE,
                              keep_editable = FALSE) {
  mean_signal <- mean(signal[, 2])
  signal[, 2] <- signal[, 2] - mean_signal
  LP <- (signal[2, 1] - signal[1, 1])
  UP <- (signal[nrow(signal), 1] - signal[1, 1])
  hilb <- Hilbert_transform(signal)

  wavelet <- analyze_wavelet(
    data = signal,
    dj = 1 / 200,
    lowerPeriod = LP,
    upperPeriod = UP,
    verbose = FALSE,
    omega_nr = 0.1
  )

  my.w <- wavelet
  my.data <- cbind(wavelet$x, wavelet$y)
  filtered_cycle <- my.data[, 1]
  filtered_cycle <- as.data.frame(filtered_cycle)
  filtered_cycle$value <- NA
  Wave = my.w$Ampl
  Power = my.w$Power
  nc = my.w$nc
  nr = my.w$nr
  dt = my.w$dt
  dj = my.w$dj
  Scale = my.w$Scale
  Period = my.w$Period
  rec.waves = matrix(0, nrow = nr, ncol = nc)
  for (s.ind in seq_len(nr)) {
    rec.waves[s.ind, ] = (Wave[s.ind, ] / sqrt(Scale[s.ind])) *
      dj * sqrt(dt) / (pi ^ (-1 / 4))
  }

  #re-scaling of the data against the standard deviation of the maximum peaks of the input data
  peaks_max <- max_detect(signal,pts=pts)
  value <- colSums(rec.waves, na.rm = T)
  value <- as.numeric(value)
  filtered_cycle[, 2] <- value
  filtered_cycle[, 2] <-
    (filtered_cycle[, 2]) * (sd(peaks_max[, 2]) / sd(value))

  if (remean == TRUE) {
    filtered_cycle[, 2] <- filtered_cycle[, 2] + mean_signal
    my.data[, 2] <- my.data[, 2] + mean_signal
  }

  ampl <- filtered_cycle

  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    par(mfrow = c(2, 1))
    par(mar = c(4, 4, 1, 1))
    plot(my.data, type = "l")
    lines(ampl)
    plot(ampl, type = "l")

  }

  if (ver_results == TRUE) {
    hilb <- Hilbert_transform(my.data)
    if (mean(ampl[, 2] / hilb[, 2]) > 1.05 |
        mean(ampl[, 2] / hilb[, 2]) < 0.95) {
      cat(
        "Mean difference between the amplitude extracted via the Hilbert Transform and the amplitude extracted
      using wavelets is larger then 5% therefore results might not be representative"
      )
      if ((mean(ampl[, 2] / hilb[, 2]) < 1.05) &
          (mean(ampl[, 2] / hilb[, 2]) > 0.95)) {
        cat(
          "Mean difference between the amplitude extracted via the Hilbert Transform and the amplitude extracted
      using wavelets is smaller then 5% hence results are representative"
        )
      }
    }
  }

  return(ampl)
}
