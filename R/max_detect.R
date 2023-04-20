#' @title Detect and filter out all maxima in a signal
#'
#' @description The \code{\link{max_detect}} function is used
#' to detect and filter out local maxima in a sinusoidal signal.
#'
#' @param data Matrix or data frame with the first column being depth or time
#'  and the second column being a proxy
#'
#'
#' @examples
#'#Example in which the ~210yr de Vries cycle is extracted from the Total Solar
#'#Irradiance data set of Steinhilver et al., (2012)\cr
#'#after which all maxima are extracted
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
#'the second column are local maxima values
#'
#' @export


 max_detect <- function(data=NULL){
  astro_maxdetect <- data
  astro_maxdetect$max <- 0
  for(i in 3:(nrow(data)-2)){
    if ((data[i,2]- data[(i+1),2] > 0) & (data[i,2] - data[(i-1),2]  > 0))
    {astro_maxdetect[i,3] <- 1}
  }

  astro_maxdetect_error_corr <- astro_maxdetect

  for(i in 1:(nrow(astro_maxdetect)-1)){
    if ((astro_maxdetect[i,3]==1) & (astro_maxdetect[i+1,3] == 1))
    {astro_maxdetect_error_corr[i+1,3] <- 0}
  }

  for(i in 1:(nrow(astro_maxdetect)-2)){
    if (((astro_maxdetect[i,3]==1) & (astro_maxdetect[i+2,3] == 1)))
    {astro_maxdetect_error_corr[i+2,3] <- 0}
  }

  for(i in 1:(nrow(astro_maxdetect)-3)){
    if ((astro_maxdetect[i,3]==1) & (astro_maxdetect[i+3,3] == 1))
    {astro_maxdetect_error_corr[i+3,3] <- 0}
  }
  for(i in 1:(nrow(astro_maxdetect)-4)){
    if ((astro_maxdetect[i,3]==1) & (astro_maxdetect[i+4,3] == 1))
    {astro_maxdetect_error_corr[i+4,3] <- 0}
  }
  for(i in 1:(nrow(astro_maxdetect)-5)){
    if ((astro_maxdetect[i,3]==1) & (astro_maxdetect[i+5,3] == 1))
    {astro_maxdetect_error_corr[i+5,3] <- 0}
  }
  for(i in 1:(nrow(astro_maxdetect)-6)){
    if ((astro_maxdetect[i,3]==1) & (astro_maxdetect[i+6,3] == 1))
    {astro_maxdetect_error_corr[i+6,3] <- 0}
  }

  for(i in 1:(nrow(astro_maxdetect)-7)){
    if ((astro_maxdetect[i,3]==1) & (astro_maxdetect[i+7,3] == 1))
    {astro_maxdetect_error_corr[i+7,3] <- 0}
  }
  for(i in 1:(nrow(astro_maxdetect)-8)){
    if ((astro_maxdetect[i,3]==1) & (astro_maxdetect[i+8,3] == 1))
    {astro_maxdetect_error_corr[i+8,3] <- 0}
  }
  astro_maxdetect_error_corr <- astro_maxdetect_error_corr[astro_maxdetect_error_corr$max == 1 , ]
  astro_maxdetect_error_corr <- astro_maxdetect_error_corr[astro_maxdetect_error_corr[,2] > mean(data[,2]), ]
  astro_maxdetect_error_corr
}

