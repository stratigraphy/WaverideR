#' @title Extract a signal/cycle from a wavelet spectra using a set period and boundaries
#'
#' @description Extracts a cycle from the wavelet object created using the \code{\link{analyze_wavelet}}
#' function using a fixed period and fixed period boundaries defined as factors of the original cycle
#'
#'@param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#'@param cycle Period of the cycle which needs to be extracted.
#'@param period_up Upper period as a factor of the to be extracted cycle \code{Default=1.2}.
#'@param period_down Lower period as a factor of the to be extracted cycle \code{Default=0.8}.
#'@param add_mean Add mean to the extracted cycle \code{Default=TRUE}.
#'@param plot_residual plot the residual signal after extraction of set cycle \code{Default=FALSE}.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#' @author
#' Code based on the \link[WaveletComp]{reconstruct} function of the 'WaveletComp' R package
#' which is based on the wavelet 'MATLAB' code written by Christopher Torrence and Gibert P. Compo (1998).
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
#'@examples
#'#Example in which the ~210yr de Vries cycle is extracted from the Total Solar
#'#Irradiance data set of Steinhilber et al., (2012)
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
#'plot_residual=FALSE,
#'keep_editable=FALSE)
#'
#'@return
#'#'Returns a matrix with 2 columns.
#'The first column is time/depth.
#'The second column is the extracted signal/cycle.
#' @export
#' @importFrom Hmisc approxExtrap
#' @importFrom stats na.omit
#' @importFrom graphics par
#' @importFrom graphics hist
#' @importFrom graphics lines
#' @importFrom WaveletComp reconstruct

extract_signal_stable <- function(wavelet = NULL,
                                  cycle = NULL,
                                  period_up = 1.2,
                                  period_down = 0.8,
                                  add_mean = TRUE,
                                  plot_residual = FALSE,
                                  keep_editable = FALSE) {
  my.w <- wavelet
  my.data <- cbind(wavelet$x, wavelet$y)
  extracted_cycle <- my.data[, 1]
  extracted_cycle <- as.data.frame(extracted_cycle)
  extracted_cycle$value <- NA


  Wave = my.w$Wave
  Power = my.w$Power

  nc = my.w$nc
  nr = my.w$nr
  dt = my.w$dt
  dj = my.w$dj

  Scale = my.w$Scale
  Period = my.w$Period
  loess.span = my.w$loess.span
  rec.waves = matrix(0, nrow = nr, ncol = nc)


  for (s.ind in seq_len(nr)) {
    rec.waves[s.ind, ] = (Re(Wave[s.ind, ]) / sqrt(Scale[s.ind])) *
      dj * sqrt(dt) / (pi ^ (-1 / 4))
  }

  cycle_high <- cycle * period_up
  cycle_low <- cycle * period_down


  row_nr_high <- Closest(Period[], cycle_high, which = TRUE)
  row_nr_low <- Closest(Period[], cycle_low, which = TRUE)

  row_nr_high <- row_nr_high[1]
  row_nr_low <- row_nr_low[1]

  value <- rec.waves[c(row_nr_low:row_nr_high),]
  value  <- colSums(value, na.rm = T)
  value <- as.numeric(value)
  extracted_cycle[, 2] <- value

  rec_value  <- colSums(rec.waves, na.rm = T)

  extracted_cycle[, 2] <-
    (extracted_cycle[, 2]) * sd(my.data[, 2]) / sd(rec_value)


  if (add_mean == TRUE) {
    extracted_cycle[, 2] <- extracted_cycle[, 2] + mean(my.data[, 2])
  }

  if (plot_residual == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    residual <- extracted_cycle[, 2] - my.data[, 2]
    layout.matrix <- matrix(c(1, 2), nrow = 2, ncol = 1)
    graphics::layout(mat = layout.matrix,
                     heights = c(1, 1),
                     # Heights of the two rows
                     widths = c(1, 1))
    par(mar = c(4, 4, 1, 1))
    plot(x = extracted_cycle[, 1], y = residual, xlab = "depth m")
    hist(residual)
  }

  colnames(extracted_cycle) <- c("x", "value")

  return(extracted_cycle)

}
