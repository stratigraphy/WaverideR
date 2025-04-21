#' @title Perform a Hilbert transform on a signal
#'
#' @description Extract the amplitude modulation using the Hilbert transform.
#'
#'@param data Input is a time series with the first column being depth or time and the second column being a proxy.
#'@param demean Remove the mean from the time series.
#'@param nr_pad nr of points added tot the top and bottom of the data set
#'to mitigate the edging effect of the Hilbert transform.
#'
#'
#'@author
#'Based on the the inst.pulse function of the 'DecomposeR' R package.
#'@references
#'Wouters, S., Crucifix, M., Sinnesael, M., Da Silva, A.C., Zeeden, C., Zivanovic, M., Boulvain, F.,
#'Devleeschouwer, X., 2022, "A decomposition approach to cyclostratigraphic signal processing".
#'Earth-Science Reviews 225 (103894). <doi:10.1016/j.earscirev.2021.103894>
#'
#'Huang, Norden E., Zhaohua Wu, Steven R. Long, Kenneth C. Arnold, Xianyao Chen, and Karin Blank. 2009.
#'"On Instantaneous Frequency". Advances in Adaptive Data Analysis 01 (02): 177–229. <doi:10.1142/S1793536909000096>
#'
#'@examples
#'#Example in which the Hilbert transform (eg. amplitude modulation) of the ~210yr
#'#de Vries cycle is extracted from the Total Solar Irradiance data set of
#'#Steinhilber et al., (2012)
#'
#'#Perform the CWT
#'TSI_wt <-
#'analyze_wavelet(
#'data = TSI,
#'dj = 1/200,
#'lowerPeriod = 16,
#'upperPeriod = 8192,
#'    verbose = FALSE,
#'    omega_nr = 6
#'  )
#'
#'#Extract the 210 yr de Vries cycle from the wavelet spectra
#'de_Vries_cycle <- extract_signal_stable(wavelet=TSI_wt,
#'cycle=210,
#'period_up =1.25,
#'period_down = 0.75,
#'add_mean=TRUE,
#'plot_residual=FALSE)
#'
#'#Perform the Hilbert transform on the amplitude record of the 210 yr de Vries
#'# cycle which was extracted from the wavelet spectra
#'
#'de_Vries_cycle_hilbert <- Hilbert_transform(data=de_Vries_cycle,demean=TRUE)
#'
#'@return
#'Returns a matrix with 2 columns.
#'The first column is depth/time.
#'The second column is the Hilbert transform of the signal.
#'
#'@export
#' @importFrom DecomposeR inst.pulse


Hilbert_transform <- function(data = NULL,
                              demean = TRUE,
                              nr_pad = 100) {

  mean_dat <- mean(data[, 2])
  data[, 2] <- data[, 2] - mean_dat

  step <- (data[2, 1] - data[1, 1])
  top <- seq(from = data[1, 1] - (nr_pad * step),
             by = step ,
             length.out = nr_pad)
  top_dat <- rep(data[1, 2], times = nr_pad)
  top <- cbind(top, top_dat)
  colnames(top) <- colnames(data)

  bottom <- seq(from = data[nrow(data), 1] + step ,
                by = step ,
                length.out = nr_pad)
  bottom_dat <- rep(data[nrow(data), 2], times = nr_pad)
  bottom <- cbind(bottom, bottom_dat)
  colnames(bottom) <- colnames(data)

  my.data <- rbind(top, data, bottom)
  my.hilbert <- as.data.frame(my.data)


  hilb_result <-
    inst.pulse(
      imf = my.data[, 2],
      dt = my.data[, 1],
      method = "HT",
      plot = FALSE
    )

  my.hilbert[, 2] <- hilb_result$a

  if (demean != TRUE) {
    my.hilbert[, 2] <- my.hilbert[, 2]+mean_dat
  }

  dat_2 <- my.hilbert[c((nr_pad + 1):(nrow(data) + nr_pad)),]

  return(dat_2)
}
