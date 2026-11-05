#' @title lag-1 autocorrelation coefficient
#'
#' @description The \code{\link{lag_1}} function calculates the lag-1 autocorrelation coefficient using a windowed analysis
#' monte carlo analysis
#'
#' @param data Two-column matrix or data frame: depth or time in column 1 and the proxy value in column 2; it is
#' sorted and linearly interpolated to its mean sampling interval. Required (\code{Default=NULL}).
#' @param n_sim Number of Monte Carlo simulations (integer, > 0). In every simulation the record is resampled at a
#' random sampling interval (truncated normal around the mean interval) and a random window size between
#' \code{win_min} and \code{win_max} is used; the output is the mean and standard deviation over the simulations.
#' \code{Default=10}.
#' @param run_multicore Run the simulations in parallel on several cores (logical, TRUE or FALSE). Results do not
#' depend on it. \code{Default=FALSE}.
#' @param n_cores Number of cores used when \code{run_multicore = TRUE} (integer or \code{NULL}). \code{NULL} gives an
#' automatic, safe number (physical cores minus one, at most 4, limited by the free memory); change this default with
#' \code{options(waverider.n_cores = ...)}. A value > 1 also switches on parallel computing. Results do not depend
#' on the number of cores. \code{Default=NULL}.
#' @param win_max Largest width of the moving window in which the lag-1 autocorrelation is computed, in the unit of
#' the first column of \code{data} (numeric, > 0). Required (\code{Default=NULL}).
#' @param win_min Smallest width of the moving window, in the unit of the first column of \code{data} (numeric, > 0;
#' the window width of each simulation is drawn uniformly between \code{win_min} and \code{win_max}). Required
#' (\code{Default=NULL}).
#' @param verbose Show a progress bar of the simulations (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @author
#'Michiel Arts
#'
#'@examples
#'\donttest{
#'#The example uses the magnetic susceptibility data set of Pas et al., (2018).
#'# perform the CWT
#'mag_wt <- analyse_wavelet(data = mag,
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
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats approx
#' @importFrom stats acf


lag_1 <- function(data = NULL,
    n_sim = 10,
    run_multicore = FALSE,
    n_cores = NULL,
    win_max = NULL,
    win_min = NULL,
    verbose = FALSE) {





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
      wr_lapply(seq_len(n_sim), function(i) {
        #i <- 1
        if ((dt - (2 * sdt)) > 0) {
          new_sampling_rate <-
              wr_rtruncnorm(
                  n = 1,
                  a = dt - 2 * sdt,
                  b = dt + (2 * sdt),
                  mean = dt,
                  sd = sdt
              )
        } else {
          new_sampling_rate <-
              wr_rtruncnorm(
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
              wr_closest(d_new[, 1], d_new[k, 1] - (win_size / 2), which = TRUE)

          row_nr_2 <-
              wr_closest(d_new[, 1], d_new[k, 1] + (win_size / 2), which = TRUE)

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
      }, n_cores = n_cores, run_multicore = run_multicore, combine = "cbind", progress = verbose)


  mat_sim_mean <- rowMeans(fit)
  mat_sim_sd  <- wr_row_sds(fit)
  results <- cbind(xout_vals, mat_sim_mean, mat_sim_sd)


  return(results)
}
