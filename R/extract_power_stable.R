#' @title Extract power from a wavelet spectra by using a constant period/duration
#'
#' @description Extract spectral power from the wavelet using a constant period/duration and
#' boundaries as selection criteria. The extraction of spectral is useful for cyclostratigraphic studies because the spectral power of an
#' astronomical cycle is modulated by higher order astronomical cycles.
#'The spectral power record from an astronomical cycle can thus be used as a proxy for
#' amplitude modulating cycles. The spectral power is extracted from a wavelet spectra
#' which was created using the \code{\link{analyze_wavelet}}
#'function for a given, \code{cycle}, \code{period_up} and \code{period_down}
#'
#'
#'@param scalogram Wavelet or superlet object created using the \code{\link{analyze_wavelet}} or
#'or \code{\link{track_period}} function
#'@param cycle Period of cycle for which the power will be extracted from the record.
#'@param period_up Species the upper period of the to be extracted power \code{Default=1.2}.
#'@param period_down specifies the lower period of the to be extracted power \code{Default=0.8}.
#'
#' @author
#' Code based on the reconstruct function of the 'WaveletComp' R package
#' which is based on the wavelet 'MATLAB' code written by Christopher Torrence and Gibert P. Compo (1998).
#' The functionality of this function is is inspired by the
#' \link[astrochron]{integratePower} function of the 'astrochron' R package
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
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#'@examples
#'#Extract the spectral power of the 210 yr de Vries cycle from the Total Solar
#'#Irradiance data set of Steinhilber et al., (2012).
#'
#'TSI_wt <-
#'  analyze_wavelet(
#'    data = TSI,
#'    dj = 1/200,
#'    lowerPeriod = 16,
#'    upperPeriod = 8192,
#'    verbose = FALSE,
#'    omega_nr = 6
#'  )
#'TSI_wt_pwr_de_Vries_cycle <-  extract_power_stable(
#'  scalogram = TSI_wt,
#'  cycle = 210,
#'  period_up = 1.2,
#'  period_down = 0.8
#')
#'
#'
#'@return
#'Returns a matrix with 3 columns.
#'The first column is depth/time.
#'The second column is extracted power.
#'The third column is extracted power/total power.
#'
#' @export
#' @importFrom DescTools Closest

extract_power_stable <- function(scalogram = NULL,
                                 cycle = NULL,
                                 period_up = 1.2,
                                 period_down = 0.8,...) {
  # Alias: wavelet -> scalogram
  if (is.null(scalogram) && !is.null(dots$wavelet)) {
    scalogram <- dots$wavelet
  }

  my.w <- scalogram
  my.data <- cbind(scalogram$x, scalogram$y)

  Period <- as.data.frame(my.w$Period)
  extract_power_high <- cycle * period_up
  extract_power_low <- cycle * period_down
  low_rownr <- Closest(Period[, 1], extract_power_low, which = TRUE)
  high_rownr <-
    Closest(Period[, 1], extract_power_high, which = TRUE)

  if(inherits(scalogram,"analyze.superlet") == TRUE){
    Power <- t(my.w$Power)}
  if(inherits(scalogram,"analyze.wavelet") == TRUE ){
    Power <- my.w$Power}

  Power_sel <- Power[(low_rownr:high_rownr),]
  Power_sel <- as.data.frame(colSums(Power_sel))
  filtered_power <- cbind(my.w$axis.1 , Power_sel)
  filtered_power <- as.data.frame(filtered_power)
  total_pwr <- colSums(Power, na.rm = TRUE)
  filtered_power$total <- total_pwr
  filtered_power$pwr_div_total <-
    filtered_power[, 2] / filtered_power[, 3]

  colnames(filtered_power) <- c("x_axis","sum_pwr","total_pwr","pwr_div_total")

  return(filtered_power)

}
