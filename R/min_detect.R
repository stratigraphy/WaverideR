#' @title Detect and filter out all minima in a signal
#'
#' @description The \code{\link{min_detect}} function is used to detect and
#' filter out local minima in a sinusoidal signal
#'
#'
#' @param data Signal in which the local extremes are detected: a two-column matrix or data frame with depth or
#' time in column 1 (in increasing order) and the proxy value in column 2. Required (\code{Default=NULL}).
#' @param pts Number of points to each side that a point is compared with (integer, >= 2): a point is a local minimum
#' when it is above (maximum) or below (minimum) the value \code{pts} points further on and the value
#' \code{pts - 1} points back. Of consecutive extremes of the same type only the most extreme one is kept. Larger
#' values give more robust extremes but may miss minor ones.
#' \code{Default=3}.
#'
#'@examples
#'#Example in which the ~210yr de Vries cycle is extracted from the Total Solar
#'#Irradiance data set of Steinhilber et al., (2012)
#'#after which all minima are extracted
#'
#'TSI_wt <-
#'analyse_wavelet(
#'data = TSI,
#'dj = 1/200,
#'lowerPeriod = 16,
#'upperPeriod = 8192,
#'    verbose = FALSE,
#'    omega_nr = 6
#'  )
#'
#'de_Vries_cycle <- extract_signal_stable(wavelet=TSI_wt,
#'cycle=210,
#'period_up =1.25,
#'period_down = 0.75,
#'add_mean=TRUE,
#'plot_residual=FALSE)
#'
#'
#'min_de_Vries_cycle <- min_detect(de_Vries_cycle)
#'
#'@return
#'#Returns a matrix with 2 columns
#'first column is depth/time
#'the second column are local minima values
#'
#' @export

min_detect <- function(data = NULL,pts=3) {
  astro_mindetect <- as.data.frame(data)
  astro_mindetect$min <- 0
  for (i in pts:(nrow(data) - pts)) {
    if ((data[i, 2] - data[(i + pts), 2] < 0) &
        (data[i, 2] - data[(i - (pts-1)), 2] < 0))
    {
      astro_mindetect[i, 3] <- 1
    }
  }

  astro_mindetect_error_corr <- astro_mindetect
  astro_mindetect_error_corr <-
      astro_mindetect_error_corr[astro_mindetect_error_corr$min == 1 ,]

  astro_maxdetect <- as.data.frame(data)
  astro_maxdetect$max <- 0
  for (i in pts:(nrow(data) - pts)) {
    if ((data[i, 2] - data[(i + pts), 2] > 0) &
        (data[i, 2] - data[(i - (pts-1)), 2]  > 0))
    {
      astro_maxdetect[i, 3] <- 1
    }
  }

  astro_maxdetect_error_corr <- astro_maxdetect
  astro_maxdetect_error_corr <-
      astro_maxdetect_error_corr[astro_maxdetect_error_corr$max == 1 ,]

  max <- astro_maxdetect_error_corr
  colnames(max) <- c("A", "B", "C")
  min <- astro_mindetect_error_corr
  colnames(min) <- c("A", "B", "C")

  min[, 3] <- -1
  peaks <- rbind(max, min)

  peaks <- peaks[order(peaks[, 1]), ]
  i <- 1
  res_rownr <- nrow(peaks)

  while (i < res_rownr) {
    if ((i < res_rownr) & (peaks[i, 3] == peaks[(i + 1), 3])) {
      if ((i < res_rownr) &(peaks[i, 3]  == 1 & peaks[(i + 1), 3] == 1) &
          (peaks[i, 2] > peaks[(i + 1), 2])) {
        peaks[(i + 1), ] <- NA
        peaks <- na.omit(peaks)
        res_rownr <- res_rownr - 1
      }
      if ((i < res_rownr) &(peaks[i, 3]  == 1 & peaks[(i + 1), 3] == 1) &
          (peaks[i, 2] < peaks[(i + 1), 2])) {
        peaks[i, ] <- NA
        peaks <- na.omit(peaks)
        res_rownr <- res_rownr - 1
      }
      if ((i < res_rownr) &(peaks[i, 3] == -1 & peaks[(i + 1), 3] == -1) &
          (peaks[i, 2] < peaks[(i + 1), 2])) {
        peaks[(i + 1), ] <- NA
        peaks <- na.omit(peaks)
        res_rownr <- res_rownr - 1
      }
      if ((i < res_rownr) &(peaks[i, 3] == -1 & peaks[(i + 1), 3] == -1) &
          (peaks[i, 2] > peaks[(i + 1), 2])) {
        peaks[i, ] <- NA
        peaks <- na.omit(peaks)
        res_rownr <- res_rownr - 1

      }
    }
    if ((peaks[i, 3] != peaks[(i + 1), 3]) |
        is.na(peaks[i, 3] != peaks[(i + 1), 3])) {
      i <- i + 1
    }

  }

  peaks_min <- peaks[peaks[, 3] < 0,]
  return(peaks_min)

}
