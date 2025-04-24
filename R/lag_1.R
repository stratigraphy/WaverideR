#' @title lag-1 autocorrelation coefficient
#'
#' @description The \code{\link{lag_1}} function calculates the lag-1 autocorrelation coefficient using a windowed analysis
#' monte carlo analysis
#'
#'@param data Input data set  should consist of a matrix with 2 columns with first column being depth and the second column being a proxy
#'@param n_sim number of simulations to be ran
#'@param run_multicore Run function using multiple cores \code{Default="FALSE"}
#'@param win_max maximum window size
#'@param win_min minimum window size
#'@param verbose print text
#' @author
#'Michiel Arts
#'
#'@examples
#'\donttest{
#'#The example uses the magnetic susceptibility data set of Pas et al., (2018).
#'# perform the CWT
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
#'#                                   periodlab = "Period (meters)",
#'#                                   x_lab = "depth (meters)")
#'
#'#Instead of tracking, the tracked solution data set mag_track_solution is used
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
#'# smooth the tracking of the 405 kyr eccentricity cycle
#' mag_track_complete <- loess_auto(time_series = mag_track_complete,
#' genplot = FALSE, print_span = FALSE)
#'#convert period in meters to sedrate depth vs time
#'mag_track_time<- curve2tune(data=mag,
#'                            tracked_cycle_curve=mag_track_complete,
#'                            tracked_cycle_period=405,
#'                            genplot = FALSE,
#'                            keep_editable=FALSE)
#'
#'mag_lag_1 <- lag_1(data = mag_track_time,n_sim = 10,
#'run_multicore = FALSE,
#'win_max = 505,
#'win_min = 150,
#'verbose=FALSE)
#'
#'}
#' @return
#'Returns a matrix which contains 3 columns
#'column 1: depth/time matrix
#'column 2: mean autocorrelation coefficient
#'column 3: sd autocorrelation coefficient
#'
#' @export
#' @importFrom Matrix rowMeans
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom foreach foreach
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom foreach %dopar%
#' @importFrom parallel stopCluster
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats approx
#' @importFrom DescTools Closest
#' @importFrom matrixStats rowSds
#' @importFrom stats acf


lag_1 <- function(data = NULL,
                  n_sim = 10,
                  run_multicore = FALSE,
                  win_max = NULL,
                  win_min = NULL,
                  verbose = FALSE) {

  if (run_multicore == TRUE) {
    numCores <- detectCores()
    cl <- parallel::makeCluster(numCores - 2)
    registerDoSNOW(cl)
  } else{
    numCores <- 1
    cl <- parallel::makeCluster(numCores)
    registerDoSNOW(cl)
  }


  if (verbose == TRUE) {
    pb <- txtProgressBar(max = n_sim, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else{
    opts = NULL
  }


  dat <- data.frame(data)
  dt2 <- dat[2, 1] - dat[1, 1]
  dat <- dat[order(dat[, 1], na.last = NA, decreasing = F), ]
  npts <- length(dat[, 1])

  start <- dat[1, 1]
  end <- dat[length(dat[, 1]), 1]

  x1 <- dat[1:(npts - 1), 1]
  x2 <- dat[2:(npts), 1]
  dx = x2 - x1
  dt = mean(dx)
  sdt = sd(dx)

  xout <- seq(start, end, by = dt)
  npts <- length(xout)
  interp <- approx(dat[, 1], dat[, 2], xout, method = "linear", n = npts)
  d <- as.data.frame(interp)



  mat_sim <- matrix(data = NA,
                    nrow = nrow(d),
                    ncol = n_sim)



  xout_vals <- d[, 1]

  i <- 1 # needed to assign 1 to ijk to avoid note
  npts <- length(xout_vals)
  new_sampling_rate <- NULL


  fit <-
    foreach (i = 1:n_sim,
             .combine = 'cbind',
             .options.parallel = opts) %dopar% {
               #i <- 1
               if ((dt - (2 * sdt)) > 0) {
                 new_sampling_rate <-
                   truncnorm::rtruncnorm(
                     n = 1,
                     a = dt - 2 * sdt,
                     b = dt + (2 * sdt),
                     mean = dt,
                     sd = sdt
                   )
               } else {
                 new_sampling_rate <-
                   truncnorm::rtruncnorm(
                     n = 1,
                     a = dt / 2,
                     b = dt +
                       (2 * sdt),
                     mean = dt,
                     sd = sdt
                   )
               }
               win_size <- stats::runif(n = 1,
                                        min = min(c(win_min, win_max)),
                                        max = max(c(win_min, win_max)))
               dt_new <- new_sampling_rate
               new_sampling_rate


               xout_new <- seq(from=start, to=end, by = dt_new)
               npts_new <- length(xout_new)
               interp_new <-
                 approx(dat[, 1], dat[, 2], xout_new, method = "linear", n = npts)
               d_new <- as.data.frame(interp_new)
               d_new[, 3] <- NA

               for (k in 1:nrow(d_new)) {
                 row_nr_1 <-
                   DescTools::Closest(d_new[, 1], d_new[k, 1] - (win_size / 2), which = TRUE)

                 row_nr_2 <-
                   DescTools::Closest(d_new[, 1], d_new[k, 1] + (win_size / 2), which = TRUE)

                 data_sel <- d_new[row_nr_1[1]:row_nr_2[1], ]
                 corr <- acf(data_sel[, 2], plot = F)

                 a <- as.numeric(unlist(corr[1])[1])

                 d_new[k, 3] <- a
               }




               yleft_comp <- d_new[1, 3]
               yright_com <- d_new[nrow(d_new), 3]

               app <-
                 approx(
                   d_new[, 1],
                   d_new[, 3],
                   xout_vals,
                   method = "linear",
                   n = npts,
                   yleft = yleft_comp,
                   yright = yright_com
                 )


               app_res <- cbind(app$y)

               return(app_res)
             }

  stopCluster(cl)

  mat_sim_mean <- rowMeans(fit)
  mat_sim_sd  <- rowSds(fit)
  results <- cbind(xout_vals, mat_sim_mean, mat_sim_sd)


  return(results)
}
