#' @title Extract signal from a wavelet spectrum using a upper and lower period boundary
#'
#' @description Extract a signal from the wavelet using a upper and lower period boundary
#'
#' @param wavelet wavelet object created using the \code{\link{analyze_wavelet}} function.
#' @param period_max Maximum period (upper boundary) to be used to extract a cycle.
#' @param period_min Minimum period (lower boundary) to be used to extract a cycle.
#' @param add_mean Add mean to the extracted cycle \code{Default=TRUE}.
#' @param plot_residual Plot the signal from which the extracted cycle is subtracted \code{Default=FALSE}.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#' @author
#' Code based on  ased on the \link[WaveletComp]{reconstruct} function of the WaveletComp R package
#' which is based on the wavelet MATLAB code written by Christopher Torrence and Gibert P. Compo.
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
#'@return Signal extracted from the wavelet spectra.
#'Output is a matrix with the first column being depth/time
#'and the second column is the cycle extracted from the proxy record.
#'
#'@examples
#'#Example in which the ~210yr de Vries cycle is extracted from the Total Solar
#'# Irradiance data set of Steinhilber et al., (2012)
#'
#'TSI_wt <-
#'analyze_wavelet(
#'data = TSI,
#'dj = 1/200,
#'lowerPeriod = 16,
#'upperPeriod = 8192,
#'    verbose = TRUE,
#'    omega_nr = 6
#'  )
#'
#'de_Vries_cycle <- extract_signal_stable_V2(wavelet=TSI_wt,
#'period_max = 240,
#'period_min = 180,
#'add_mean=TRUE,
#'plot_residual=FALSE,
#'keep_editable=FALSE)
#'
#'
#' @export
#' @importFrom graphics par
#' @importFrom graphics hist
#' @importFrom graphics lines
#' @importFrom WaveletComp reconstruct


extract_signal_stable_V2 <- function(wavelet = NULL,
                                     period_max = NULL,
                                     period_min = NULL,
                                     add_mean = TRUE,
                                     plot_residual = FALSE,
                                     keep_editable = FALSE) {
  my.w <- wavelet
  my.data <- cbind(wavelet$x, wavelet$y)
  filtered_cycle <- my.data[, 1]
  filtered_cycle <- as.data.frame(filtered_cycle)
  filtered_cycle$value <- NA


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

  row_nr_high <- Closest(Period[], period_max, which = TRUE)
  row_nr_low <- Closest(Period[], period_min, which = TRUE)

  row_nr_high <- row_nr_high[1]
  row_nr_low <- row_nr_low[1]

  value <- rec.waves[c(row_nr_low:row_nr_high),]
  value  <- colSums(value, na.rm = T)
  value <- as.numeric(value)
  filtered_cycle[, 2] <- value

  rec_value  <- colSums(rec.waves, na.rm = T)

  filtered_cycle[, 2] <-
    (filtered_cycle[, 2]) * sd(my.data[, 2]) / sd(rec_value)


  if (add_mean == TRUE) {
    filtered_cycle[, 2] <- filtered_cycle[, 2] + mean(my.data[, 2])
  }

  if (plot_residual == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    residual <- filtered_cycle[, 2] - my.data[, 2]
    layout.matrix <- matrix(c(1, 2), nrow = 2, ncol = 1)
    graphics::layout(mat = layout.matrix,
                     heights = c(1, 1),
                     # Heights of the two rows
                     widths = c(1, 1))
    par(mar = c(4, 4, 1, 1))
    plot(x = filtered_cycle[, 1], y = residual, xlab = "depth m")
    hist(residual)
  }

  filtered_cycle
}
