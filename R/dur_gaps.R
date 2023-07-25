#' @title calculate the duration of stratigraphic gaps using astronomical cycles
#'
#' @description calculate the duration of stratigraphic gaps using the duration
#' of stable astronomical cycles
#'
#'
#'@param proxies list of proxies which were used to create a astrochronological
#'age model and which are used to calculate the duration of the gap
#'@param retracked_period_1 A matrix of 3 columns in which the first column
#' is depth/height.The second column is the period of the tracked cycle.
#' The thirds column is uncertainty given as 1 standard deviation for the
#' period of the tracked cycle. The gap to be modeled should be located
#' in between retracked_period_1 and retracked_period_2
#'@param retracked_period_2 A matrix of 3 columns in which the first column
#' is depth/height.The second column is the period of the tracked cycle.
#' The thirds column is uncertainty given as 1 standard deviation for the
#' period of the tracked cycle. The gap to be modeled should be located
#' in between retracked_period_1 and retracked_period_2
#'@param min_max list of "min" or "max" indicating whether time should be
#'estimated between minima or maxima for each proxy
#'@param n_simulations number of gap duration to calculate
#'@param tracked_cycle_period period in time of the tracked cycle
#'@param tracked_cycle_period_unc uncertainty in the period of the tracked cycle
#'@param tracked_cycle_period_unc_dist distribution of the uncertainty of the
#'tracked cycle value need to be either "u" for uniform distribution or
#'"n" for normal distribution  \code{Default="u"}
#' @param pts the pts parameter specifies how many points to the left/right up/down the peak detect algorithm goes in detecting
#'a peak. The peak detecting algorithm works by comparing the values left/right up/down of it, if the values are both higher or lower
#'then the value a peak. To deal with error produced by this algorithm the pts parameter can be changed which can
#'aid in peak detection. Usually increasing the pts parameter means more peak certainty, however it also means that minor peaks might not be
#'picked up by the algorithm \code{Default=5}#'
#' @param dj Spacing between successive scales. The CWT analyses analyses the signal using successive periods
#' which increase by the power of 2 (e.g.2^0=1,2^1=2,2^2=4,2^3=8,2^4=16). To have more resolution
#' in-between these steps the dj parameter exists, the dj parameter specifies how many extra steps/spacing in-between
#' the power of 2 scaled CWT is added. The amount of steps is 1/x with a higher x indicating a smaller spacing.
#' Increasing the increases the computational time of the CWT \code{Default=1/200}.
#' @param lowerPeriod  Lowest period to be analyzed \code{Default=2}.
#' The CWT analyses the signal starting from the lowerPeriod to the upperPeriod so the proper selection these
#' parameters allows to analyze the signal for a specific range of cycles.
#' scaling is done using power 2 so for the best plotting results select a value to the power or 2.
#' @param upperPeriod Upper period to be analyzed \code{Default=1024}.
#' The CWT analyses the signal starting from the lowerPeriod to the upperPeriod so the proper selection these
#' parameters allows to analyze the signal for a specific range of cycles.
#'  scaling is done using power 2 so for the best plotting results select a value to the power or 2.
#' @param period_max Maximum period (upper boundary) to be used to extract a cycle.
#' @param period_min Minimum period (lower boundary) to be used to extract a cycle.
#'@param dur_between_pts duration in time between the minima or maxima
#'@param run_multicore Run function using multiple cores \code{Default="FALSE"}
#'
#'@return
#'a vector with all the calculated gap durations
#'
#' @export
#' @importFrom stats quantile
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom tcltk setTkProgressBar
#' @importFrom tcltk setTkProgressBar
#' @importFrom stats runif
#' @importFrom stats rnorm



dur_gaps <- function(proxies = NULL,
                     retracked_period_1 = NULL,
                     retracked_period_2 = NULL,
                     min_max =  NULL,
                     n_simulations = 10,
                     tracked_cycle_period = NULL,
                     tracked_cycle_period_unc =NULL,
                     tracked_cycle_period_unc_dist = "u",
                     pts = 5,
                     dj = 1 / 200,
                     lowerPeriod = 1,
                     upperPeriod = 3200,
                     period_max = NULL,
                     period_min = NULL,
                     dur_between_pts = NULL,
                     run_multicore = FALSE) {
  if (run_multicore == TRUE) {
    j <- 1
    numCores <- detectCores()
    cl <- makeCluster(numCores - 2)
    registerDoSNOW(cl)

    pb <- txtProgressBar(max = n_simulations, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)


    dur_gaps <-
      foreach::foreach (
        j = 1:(n_simulations),
        .options.snow = opts,
        .errorhandling = "remove",
        .packages = c("WaverideR", "stats"),
        .combine = "cbind"
      ) %dopar% {
        new_curve_1 <- matrix(
          data = NA,
          nrow = nrow(retracked_period_1),
          ncol = 2
        )
        new_curve_1[, 1] <- retracked_period_1[, 1]

        new_curve_2 <-
          matrix(
            data = NA,
            nrow = nrow(retracked_period_2),
            ncol = 2
          )
        new_curve_2[, 1] <-  retracked_period_2[, 1]

        proxy_nr <- floor(runif(1, min = 1, max = length(proxies)))


        val_1 <-
          rnorm(1, mean = retracked_period_1[1, 2], sd = retracked_period_1[1, 3])
        pnorm_val_1 <-pnorm(val_1, mean = retracked_period_1[1, 2], sd = retracked_period_1[1, 3])
        for (k in 1:nrow(new_curve_1)) {
          new_curve_1[k, 2] <-
            qnorm(pnorm_val_1, mean = retracked_period_1[k, 2], sd = retracked_period_1[k, 3])
        }

        val_2 <-
          rnorm(1, mean = retracked_period_2[1, 2], sd = retracked_period_2[1, 3])
        pnorm_val_2 <-
          pnorm(val_2, mean = retracked_period_2[1, 2], sd = retracked_period_2[1, 3])
        for (k in 1:nrow(new_curve_2)) {
          new_curve_2[k, 2] <-
            qnorm(pnorm_val_2, mean = retracked_period_2[k, 2], sd = retracked_period_2[k, 3])
        }


        proxy_1 <- proxies[[proxy_nr]]
        proxy_1 <- proxy_1[proxy_1[, 1] >= new_curve_1[1, 1], ]
        proxy_1 <-
          proxy_1[proxy_1[, 1] <= new_curve_1[nrow(new_curve_1), 1], ]

        proxy_2 <- proxies[[proxy_nr]]
        proxy_2 <- proxy_2[proxy_2[, 1] >= new_curve_2[1, 1], ]
        proxy_2 <-
          proxy_2[proxy_2[, 1] <= new_curve_2[nrow(new_curve_2), 1], ]


        if (tracked_cycle_period_unc_dist == "u"){
          tracked_cycle_period_new <- runif(1, min = tracked_cycle_period-tracked_cycle_period_unc, max = tracked_cycle_period+tracked_cycle_period_unc)
        }

        if (tracked_cycle_period_unc_dist == "n"){
          tracked_cycle_period_new <- rnorm(1, mean = tracked_cycle_period, sd = tracked_cycle_period_unc)
        }

        tuned_1 <- WaverideR::curve2tune(
          data = proxy_1,
          tracked_cycle_curve = new_curve_1,
          tracked_cycle_period = tracked_cycle_period_new
        )

        tuned_2 <- WaverideR::curve2tune(
          data = proxy_2,
          tracked_cycle_curve = new_curve_2,
          tracked_cycle_period = tracked_cycle_period_new
        )


        tuned_1_wt <-
          analyze_wavelet(
            tuned_1,
            dj = dj,
            lowerPeriod = lowerPeriod,
            upperPeriod = upperPeriod
          )
        tuned_1_wt_cycle <- extract_signal_stable_V2(
          wavelet = tuned_1_wt,
          period_max = period_max,
          period_min = period_min,
          add_mean = FALSE,
          plot_residual = FALSE,
          keep_editable = FALSE
        )

        tuned_1_wt <- NULL

        tuned_2_wt <-
          analyze_wavelet(
            tuned_2,
            dj = dj,
            lowerPeriod = lowerPeriod,
            upperPeriod = upperPeriod
          )
        tuned_2_wt_cycle <- extract_signal_stable_V2(
          wavelet = tuned_2_wt,
          period_max = period_max,
          period_min = period_min,
          add_mean = FALSE,
          plot_residual = FALSE,
          keep_editable = FALSE
        )
        tuned_2_wt <- NULL

        peak_opt <- min_max[[proxy_nr]]

        if (peak_opt == "min") {
          pts_tuned_1_wt_cycle <- min_detect(data = tuned_1_wt_cycle, pts = pts)
          pts_tuned_2_wt_cycle <-
            min_detect(data = tuned_2_wt_cycle, pts = pts)
        }

        if (peak_opt == "max") {
          pts_tuned_1_wt_cycle <- max_detect(data = tuned_1_wt_cycle, pts = pts)
          pts_tuned_2_wt_cycle <-
            max_detect(data = tuned_2_wt_cycle, pts = pts)
        }


        dur_gap <-
          (dur_between_pts) - ((pts_tuned_2_wt_cycle[1, 1] + max(tuned_1_wt_cycle[, 1])) -
                                 max(pts_tuned_1_wt_cycle[, 1]))
      }
  } else {
    dur_gaps <- NA
    for (j in 1:n_simulations) {
      new_curve_1 <- matrix(
        data = NA,
        nrow = nrow(retracked_period_1),
        ncol = 2
      )
      new_curve_1[, 1] <- retracked_period_1[, 1]

      new_curve_2 <-
        matrix(
          data = NA,
          nrow = nrow(retracked_period_2),
          ncol = 2
        )
      new_curve_2[, 1] <-  retracked_period_2[, 1]

      proxy_nr <- floor(runif(1, min = 1, max = length(proxies)))


      val_1 <-
        rnorm(1, mean = retracked_period_1[1, 2], sd = retracked_period_1[1, 3])
      pnorm_val_1 <-
        pnorm(val_1, mean = retracked_period_1[1, 2], sd = retracked_period_1[1, 3])
      for (k in 1:nrow(new_curve_1)) {
        new_curve_1[k, 2] <-
          qnorm(pnorm_val_1, mean = retracked_period_1[k, 2], sd = retracked_period_1[k, 3])
      }

      val_2 <-
        rnorm(1, mean = retracked_period_2[1, 2], sd = retracked_period_2[1, 3])
      pnorm_val_2 <-
        pnorm(val_2, mean = retracked_period_2[1, 2], sd = retracked_period_2[1, 3])
      for (k in 1:nrow(new_curve_2)) {
        new_curve_2[k, 2] <-
          qnorm(pnorm_val_2, mean = retracked_period_2[k, 2], sd = retracked_period_2[k, 3])
      }


      proxy_1 <- proxies[[proxy_nr]]
      proxy_1 <- proxy_1[proxy_1[, 1] >= new_curve_1[1, 1], ]
      proxy_1 <-
        proxy_1[proxy_1[, 1] <= new_curve_1[nrow(new_curve_1), 1], ]

      proxy_2 <- proxies[[proxy_nr]]
      proxy_2 <- proxy_2[proxy_2[, 1] >= new_curve_2[1, 1], ]
      proxy_2 <-
        proxy_2[proxy_2[, 1] <= new_curve_2[nrow(new_curve_2), 1], ]


      if (tracked_cycle_period_unc_dist == "u"){
        tracked_cycle_period_new <- runif(1, min = tracked_cycle_period-tracked_cycle_period_unc, max = tracked_cycle_period+tracked_cycle_period_unc)
      }

      if (tracked_cycle_period_unc_dist == "n"){
        tracked_cycle_period_new <- rnorm(1, mean = tracked_cycle_period, sd = tracked_cycle_period_unc)
      }

      tuned_1 <- WaverideR::curve2tune(
        data = proxy_1,
        tracked_cycle_curve = new_curve_1,
        tracked_cycle_period = tracked_cycle_period_new
      )

      tuned_2 <- WaverideR::curve2tune(
        data = proxy_2,
        tracked_cycle_curve = new_curve_2,
        tracked_cycle_period = tracked_cycle_period_new
      )


      tuned_1_wt <-
        analyze_wavelet(
          tuned_1,
          dj = dj,
          lowerPeriod = lowerPeriod,
          upperPeriod = upperPeriod
        )
      tuned_1_wt_cycle <- extract_signal_stable_V2(
        wavelet = tuned_1_wt,
        period_max = period_max,
        period_min = period_min,
        add_mean = FALSE,
        plot_residual = FALSE,
        keep_editable = FALSE
      )

      tuned_1_wt <- NULL

      tuned_2_wt <-
        analyze_wavelet(
          tuned_2,
          dj = dj,
          lowerPeriod = lowerPeriod,
          upperPeriod = upperPeriod
        )
      tuned_2_wt_cycle <- extract_signal_stable_V2(
        wavelet = tuned_2_wt,
        period_max = period_max,
        period_min = period_min,
        add_mean = FALSE,
        plot_residual = FALSE,
        keep_editable = FALSE
      )
      tuned_2_wt <- NULL

      peak_opt <- min_max[[proxy_nr]]

      if (peak_opt == "min") {
        pts_tuned_1_wt_cycle <- min_detect(data = tuned_1_wt_cycle, pts = pts)
        pts_tuned_2_wt_cycle <-
          min_detect(data = tuned_2_wt_cycle, pts = pts)
      }

      if (peak_opt == "max") {
        pts_tuned_1_wt_cycle <- max_detect(data = tuned_1_wt_cycle, pts = pts)
        pts_tuned_2_wt_cycle <-
          max_detect(data = tuned_2_wt_cycle, pts = pts)
      }


      dur_gap <-
        (dur_between_pts) - ((pts_tuned_2_wt_cycle[1, 1] + max(tuned_1_wt_cycle[, 1])) -
                               max(pts_tuned_1_wt_cycle[, 1]))


      dur_gaps <- c(dur_gaps, dur_gap)
    }

    dur_gaps <- dur_gaps[-c(1)]

  }

  return(dur_gaps)
}
