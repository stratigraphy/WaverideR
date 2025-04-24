#' @title Windowed timeOpt sedimentation rate estimation
#'
#' @description The \code{\link{win_timeOpt}} function for conducts a widowed
#' timeOpt sedimentation rate estimation
#'This function is based on the \code{\link[astrochron:eTimeOpt]{eTimeOpt}} but allows for
#' multithreaded analysis speeding up the
#'process of conducting a Windowed timeOpt sedimentation rate estimation
#'@param data Input data set  should consist of a matrix with 2 columns with the
#'first column being depth and the second column being a proxy \code{Default=NULL}
#'@param window_size size of the moving window in metres \code{Default=15}
#'@param  sedmin Minimum sedimentation rate for investigation (cm/ka). \code{Default=0.1}
#'@param  sedmax Maximum sedimentation rate for investigation (cm/ka). \code{Default=1}
#'@param  numsed Number of sedimentation rates to investigate
#' in optimization grid. \code{Default=100}
#'@param  limit Limit evaluated sedimentation rates to region in which full
#'target signal can be recovered? .\code{Default=FALSE}
#'@param  fit Test for (1) precession amplitude modulation or (2) short
#'eccentricity amplitude modulation? \code{Default=2}
#'@param  fitModPwr Include the modulation periods
#'in the spectral fit? \code{Default=TRUE}
#'@param  flow 	Low frequency cut-off for
#'Taner bandpass (half power point in cycles/ka) \code{Default=TRUE}
#'@param  fhigh High frequency cut-off for
#'Taner bandpass (half power point; in cycles/ka) \code{Default=NULL}
#'@param  roll Taner filter roll-off rate, in dB/octave. \code{Default=c(10^6)}
#'@param  targetE A vector of eccentricity periods to evaluate (in ka).
#'These must be in order of decreasing period, with a first value of 405 ka.
#' \code{Default= "c(405.7, 130.7, 123.8, 98.9, 94.9)"}
#'@param  targetP A vector of precession periods to evaluate (in ka).
#' These must be in order of decreasing period. \code{Default=c(20.9, 19.9, 17.1, 17.2)}
#'@param  detrend Remove linear trend from data series? \code{Default=TRUE}
#'@param  normalize normalize the r2 curves of individual timeOpt runs \code{Default=TRUE}
#'@param  linLog Use linear or logarithmic scaling for sedimentation
#'rate grid spacing? (0=linear, 1=log; default value is 1) \code{Default=1}
#'@param  run_multicore Run function using multiple cores \code{Default=FALSE}
#'@param  verbose print text \code{Default=FALSE}

#' @author
#'Based on the \code{\link[astrochron:eTimeOpt]{eTimeOpt}}
#'function of the 'astrochron' R package.
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#'
#'@examples
#'\donttest{
#'#Conduct a windowed timeOpt on the magnetic susceptibility record
#'#of the Sullivan core of Pas et al., (2018).
#'mag_win_timeOpt <-win_timeOpt(
#'data = mag,
#'window_size = 15,
#'sedmin = 0.1,
#'sedmax = 1,
#'numsed = 100,
#'limit = FALSE,
#'fit = 2,
#'fitModPwr = TRUE,
#'flow = NULL,
#'fhigh = NULL,
#'roll = 10 ^ 6,
#'targetE = c(405.7, 130.7, 123.8, 98.9, 94.9),
#'targetP = c(20.9, 19.9, 17.1, 17.2),
#'detrend = TRUE,
#'normalize =TRUE,
#'linLog = 1,
#'run_multicore =FALSE,
#'verbose=FALSE)
#'}
#'
#' @return
#'Returns a list which contains 10 elements
#'element 1: r_2_envelope matrix
#'element 2: r_2_power matrix
#'element 3: r_2_opt matrix
#'element 4: r_2_envelope_avg
#'element 5: r_2_opt_avg
#'element 6: depth
#'element 7: y_axis
#'element 8: linLog value
#'
#' @export
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel stopCluster
#' @importFrom astrochron timeOpt
#' @importFrom DescTools Closest

win_timeOpt <- function(data = NULL,
                        window_size = 10,
                        sedmin = 0.5,
                        sedmax = 2,
                        numsed = 100,
                        limit = FALSE,
                        fit = 2,
                        fitModPwr = TRUE,
                        flow = NULL,
                        fhigh = NULL,
                        roll = 10 ^ 6,
                        targetE = c(405.7, 130.7, 123.8, 98.9, 94.9),
                        targetP = c(20.9, 19.9, 17.1, 17.2),
                        detrend = TRUE,
                        normalize = TRUE,
                        linLog = 1,
                        run_multicore = FALSE,
                        verbose=FALSE ) {

  if (run_multicore == TRUE) {
    numCores <- detectCores()
    cl <- parallel::makeCluster(numCores - 2)
    registerDoSNOW(cl)
  } else{
    numCores <- 1
    cl <- parallel::makeCluster(numCores)
    registerDoSNOW(cl)
  }

  output = 1
  title = NULL
  genplot = FALSE
  check = FALSE
  verbose_2 = FALSE
  proxy <- data
  n_simulations <- nrow(proxy)

  if (verbose==TRUE){
    pb <- txtProgressBar(max = n_simulations, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)}else{opts=NULL}

  op <- 1
  win_timeOpt_res <-
    foreach (
      op = 1:n_simulations,
      .options.snow = opts,
      .combine = 'cbind',
      .packages = c("astrochron", "DescTools"),
      .errorhandling = c("pass")
    ) %dopar% {

      d <- proxy
      dt <- d[2, 1] - d[1, 1]

      left <- d[op, 1] - (window_size / 2)
      right <- d[op, 1] + (window_size / 2)

      row_nr_1 <- DescTools::Closest(d[, 1], left, which = TRUE)
      row_nr_1 <- row_nr_1[1]

      row_nr_2 <- DescTools::Closest(d[, 1], right, which = TRUE)
      row_nr_2 <- row_nr_2[1]

      d_subsel <- d[row_nr_1:row_nr_2, ]
      timeOpt_res <- astrochron::timeOpt(
        dat = d_subsel,
        sedmin = sedmin,
        sedmax = sedmax,
        numsed = numsed,
        linLog = 1,
        limit = limit,
        fit = fit,
        fitModPwr = fitModPwr,
        flow = flow,
        fhigh = fhigh,
        roll = roll,
        targetE = targetE,
        targetP = targetP,
        detrend = detrend,
        output = output,
        title = title,
        genplot = genplot,
        check = check,
        verbose = verbose_2,
      )


      if (normalize == TRUE) {
        timeOpt_res[, 2] <- timeOpt_res[, 2] / max(timeOpt_res[, 2])
        timeOpt_res[, 3] <- timeOpt_res[, 3] / max(timeOpt_res[, 3])
        timeOpt_res[, 4] <- timeOpt_res[, 4] / max(timeOpt_res[, 4])
        timeOpt_res[!is.finite(as.matrix(timeOpt_res)), ] <- 0
      }
      timeOpt_res <- c(timeOpt_res[, 2], timeOpt_res[, 3], timeOpt_res[, 4])


    }



  stopCluster(cl)


  if (linLog == 0) {
    sedinc = (sedmax - sedmin) / (numsed - 1)
    sedrate = sedmin + 0:(numsed - 1) * sedinc
  }
  if (linLog == 1) {
    sedinc = (log10(sedmax) - log10(sedmin)) / (numsed -
                                                  1)
    sedrate = log10(sedmin) + 0:(numsed - 1) * sedinc
    sedrate = 10 ^ sedrate
  }

  r_2_envelope <- win_timeOpt_res[1:numsed, ]
  r_2_power <- win_timeOpt_res[numsed + 1:numsed * 2, ]
  r_2_opt <- win_timeOpt_res[(numsed * 2 + 1):nrow(win_timeOpt_res), ]

  r_2_envelope_avg = rowMeans(r_2_envelope)
  r_2_power_avg = rowMeans(r_2_power)
  r_2_opt_avg = rowMeans(r_2_opt)


  output <-
    list(
      r_2_envelope = r_2_envelope,
      r_2_power = r_2_power,
      r_2_opt = r_2_opt,
      r_2_envelope_avg = r_2_envelope_avg,
      r_2_power_avg = r_2_power_avg,
      r_2_opt_avg = r_2_opt_avg,
      sedrate = sedrate,
      x = proxy[, 1],
      y = proxy[, 2],
      linLog = linLog
    )


  class(output) = "win_timeOpt.result"
  return(invisible(output))
}

