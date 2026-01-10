#' @title exponentially ramped  sedimentation rate model to convert time to stratigraphy
#'
#' @description Apply a exponentially increasing increasing sedimentation rate model
#'  to convert time to stratigraphy.
#'
#'@param data Input data set  should consist of a matrix with 2 columns with the first being time
#'and the second being the proxy
#'@param sr_start Initial sedimentation rate (in cm/kyr).
#'@param sr_end Final sedimentation rate (in cm/kyr).
#'
#' @author
#'Based on the \link[astrochron]{sedRamp}
#'function of the 'astrochron' R package.
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#'
#'@examples
#'
#'data_1 <- etp(
#'tmin = 0,
#'tmax = 4000,
#'dt = 1,
#'eWt = 1.5,
#'oWt = 0.75,
#'pWt = 1,
#'esinw = T,
#'solution = NULL,
#'standardize = T,
#'genplot = T,
#'verbose = T
#')
#'data_ramped_ls <- expSedRamp(data_1, sr_start = 1, sr_end = 5)
#' @return
#'Returns a list which contains 10 elements
#'element 1: time series in the depth (m) domain
#'element 2: sedimentation rate curve
#' @export


expSedRamp <- function(data, sr_start, sr_end) {
  #time and proxy
  t <- data[, 1] #kyr
  x <- data[, 2]
  # convert SR from cm/kyr to m/kyr
  sr_start <- sr_start / 100
  sr_end <- sr_end / 100
  # time step (robust to non-constant dt)
  dt <- mean(diff(t), na.rm = TRUE) # kyr
  # ornmalized time
  Ttot <- max(t) - min(t)
  tau <- (t - min(t)) / Ttot
  # exponential SR(t) in m/kyr, log2-linear increase
  k <- log2(sr_end / sr_start)
  sr_t <- sr_start * 2^(k * tau)
  # integrate to depth (meters)
  depth <- cumsum(sr_t * dt)
  # interpolate SR onto depth domain
  sr_depth <- approx(depth, sr_t, xout = depth)$y
  # output: depth–proxy and depth–SR
  out <- list(depth_proxy = cbind(depth, x),
              sr_depth = cbind(depth, sr_depth * 100))
  return(out)
}
