#' @title calculate the duration of stratigraphic gaps using astronomical cycles
#'
#' @description calculate the duration of stratigraphic gaps using the duration
#' of stable astronomical cycles
#'
#'
#' @param proxies Proxy records that were used to build the astrochronological age model and from which the gap
#' duration is estimated (list of two-column matrices or data frames: depth/height in column 1, proxy value in
#' column 2). In each simulation one proxy of the list is drawn at random. Required (\code{Default=NULL}).
#' @param retracked_period_1 Tracked cycle curve of the interval below/before the gap (three-column matrix or data
#' frame: depth/height in column 1, period of the tracked cycle in column 2, its 1 standard deviation uncertainty
#' in column 3, e.g. the output of \code{\link{retrack_wt_MC}}). In each simulation a curve is drawn from this
#' normal uncertainty envelope and used to tune the proxy with \code{\link{curve2tune}}. The gap must lie between
#' \code{retracked_period_1} and \code{retracked_period_2}. Required (\code{Default=NULL}).
#' @param retracked_period_2 Tracked cycle curve of the interval above/after the gap, in the same format as
#' \code{retracked_period_1}. Required (\code{Default=NULL}).
#' @param min_max For each proxy in \code{proxies}, whether the gap duration is measured between minima
#' (\code{"min"}, peaks found with \code{\link{min_detect}}) or maxima (\code{"max"}, \code{\link{max_detect}}) of
#' the extracted cycle (list or character vector of the same length as \code{proxies}). Required
#' (\code{Default=NULL}).
#' @param n_simulations Number of Monte Carlo simulations, i.e. number of gap durations calculated (integer, > 0).
#' \code{Default=10}.
#' @param tracked_cycle_period Duration in time of the tracked cycle (numeric, e.g. in kyr; e.g. 405 for the long
#' eccentricity cycle). The tuned records and the gap duration are in this time unit. Required
#' (\code{Default=NULL}).
#' @param tracked_cycle_period_unc Uncertainty of \code{tracked_cycle_period}, in the same unit (numeric, >= 0):
#' half-width of the range for \code{"u"}, standard deviation for \code{"n"} (see
#' \code{tracked_cycle_period_unc_dist}). Required (\code{Default=NULL}).
#' @param tracked_cycle_period_unc_dist Distribution from which the tracked cycle period is drawn in each
#' simulation (character): \code{"u"} uniform between \code{tracked_cycle_period -/+ tracked_cycle_period_unc}
#' or \code{"n"} normal with mean \code{tracked_cycle_period} and standard deviation
#' \code{tracked_cycle_period_unc}. \code{Default="u"}.
#' @param pts Number of points on each side that a point must exceed (maximum) or be below (minimum) to be
#' detected as a peak in the extracted cycle; passed to \code{\link{max_detect}}/\code{\link{min_detect}}
#' (integer, >= 1). Larger values give more certain peaks but may miss minor ones. \code{Default=5}.
#' @param dj Spacing between successive scales of the wavelet transform of the tuned records, in octaves
#' (numeric, > 0), passed to \code{\link{analyse_wavelet}}. Smaller values give a finer period resolution but a
#' longer computation time. \code{Default=1/200}.
#' @param lowerPeriod Shortest period of the wavelet transform of the tuned records, in the time unit of
#' \code{tracked_cycle_period} (numeric, > 0), passed to \code{\link{analyse_wavelet}}. \code{Default=1}.
#' @param upperPeriod Longest period of the wavelet transform of the tuned records, in the time unit of
#' \code{tracked_cycle_period} (numeric, > \code{lowerPeriod}), passed to \code{\link{analyse_wavelet}}.
#' \code{Default=3200}.
#' @param period_max Upper period boundary (time unit) of the cycle extracted from the tuned records with
#' \code{\link{extract_signal_stable_V2}} (numeric, > \code{period_min}). Required (\code{Default=NULL}).
#' @param period_min Lower period boundary (time unit) of the cycle extracted from the tuned records with
#' \code{\link{extract_signal_stable_V2}} (numeric, > 0). Required (\code{Default=NULL}).
#' @param missing_cycle_dur Duration of one missing cycle in the gap, in the time unit (numeric). If it equals
#' \code{tracked_cycle_period}, the value drawn for the tracked cycle is also used for the missing cycle.
#' Required (\code{Default=NULL}).
#' @param n_cycles_missing Number of cycles of duration \code{missing_cycle_dur} that are missing in the gap
#' (numeric, >= 0). \code{Default=1}.
#' @param missing_cycle_unc Uncertainty of \code{missing_cycle_dur} in the time unit (numeric, >= 0): half-width
#' for \code{"u"}, standard deviation for \code{"n"}. Only used when \code{missing_cycle_dur} differs from
#' \code{tracked_cycle_period}. Required in that case (\code{Default=NULL}).
#' @param missing_cycle_unc_dist Distribution from which the missing cycle duration is drawn (character):
#' \code{"u"} uniform or \code{"n"} normal (see \code{missing_cycle_unc}). \code{Default="u"}.
#' @param seed_nr Seed of the random number generator (integer), set with \code{set.seed()} so that the Monte
#' Carlo results can be reproduced. \code{Default=1337}.
#' @param run_multicore Run the simulations in parallel on several cores (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param n_cores Number of cores used when \code{run_multicore=TRUE} (integer or \code{NULL}). \code{NULL} picks
#' an automatic, safe number: physical cores minus one, at most 4, and limited by the free memory; set
#' \code{options(waverider.n_cores = ...)} to change that default. Results do not depend on the number of cores.
#' \code{Default=NULL}.
#'@return
#'a vector with all the calculated gap durations
#'
#' @export
#' @importFrom stats quantile
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats runif
#' @importFrom stats rnorm



dur_gaps <- function(proxies = NULL,
    retracked_period_1 = NULL,
    retracked_period_2 = NULL,
    min_max =  NULL,
    n_simulations = 10,
    tracked_cycle_period = NULL,
    tracked_cycle_period_unc = NULL,
    tracked_cycle_period_unc_dist = "u",
    pts = 5,
    dj = 1 / 200,
    lowerPeriod = 1,
    upperPeriod = 3200,
    period_max = NULL,
    period_min = NULL,
    missing_cycle_dur = NULL,
    n_cycles_missing = 1,
    missing_cycle_unc = NULL,
    missing_cycle_unc_dist = "u",
    seed_nr=1337,
    run_multicore = FALSE, n_cores = NULL) {


  if (run_multicore == TRUE) {
    set.seed(seed_nr)

    dur_gaps <-
        wr_lapply(seq_len(n_simulations), function(j) {
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


          val_1 <- abs(
              rnorm(1, mean = retracked_period_1[1, 2], sd = retracked_period_1[1, 3]))

          pnorm_val_1 <-  pnorm(
              val_1,
              mean =retracked_period_1[1, 2],
              sd = retracked_period_1[1, 3],
              lower.tail = FALSE
          )


          for (k in 1:nrow(new_curve_1)) {
            new_curve_1[k, 2] <-
                1/(qnorm(pnorm_val_1, mean = retracked_period_1[k, 2], sd = retracked_period_1[k, 3]))
          }

          val_2 <- abs(
              rnorm(1, mean = retracked_period_2[1, 2], sd = retracked_period_2[1, 3]))

          pnorm_val_2 <-  pnorm(
              val_2,
              mean = retracked_period_2[1, 2],
              sd = retracked_period_2[1, 3],
              lower.tail = FALSE
          )

          for (k in 1:nrow(new_curve_2)) {
            new_curve_2[k, 2] <-
                1/(qnorm(pnorm_val_2, mean = retracked_period_2[k, 2], sd = retracked_period_2[k, 3]))
          }


          proxy_1 <- proxies[[proxy_nr]]
          proxy_1 <- proxy_1[proxy_1[, 1] >= new_curve_1[1, 1], ]
          proxy_1 <-
              proxy_1[proxy_1[, 1] <= new_curve_1[nrow(new_curve_1), 1], ]

          proxy_2 <- proxies[[proxy_nr]]
          proxy_2 <- proxy_2[proxy_2[, 1] >= new_curve_2[1, 1], ]
          proxy_2 <-
              proxy_2[proxy_2[, 1] <= new_curve_2[nrow(new_curve_2), 1], ]


          if (tracked_cycle_period_unc_dist == "u") {
            tracked_cycle_period_new <-
                runif(
                    1,
                    min = tracked_cycle_period - tracked_cycle_period_unc,
                    max = tracked_cycle_period + tracked_cycle_period_unc
                )
          }

          if (tracked_cycle_period_unc_dist == "n") {
            tracked_cycle_period_new <-
                rnorm(1, mean = tracked_cycle_period, sd = tracked_cycle_period_unc)
          }

          if (tracked_cycle_period == missing_cycle_dur) {
            missing_cycle_dur_new <- tracked_cycle_period_new
          } else{
            if (missing_cycle_unc_dist == "u") {
              missing_cycle_dur_new <-
                  runif(1,
                      min = missing_cycle_dur - missing_cycle_unc,
                      max = tracked_cycle_period + missing_cycle_unc)
            }
            if (missing_cycle_unc_dist == "n") {
              missing_cycle_dur_new <-
                  rnorm(1, mean = missing_cycle_dur, sd = missing_cycle_unc)
            }
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
              analyse_wavelet(
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
              analyse_wavelet(
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
            pts_tuned_1_wt_cycle <-
                min_detect(data = tuned_1_wt_cycle, pts = pts)
            pts_tuned_2_wt_cycle <-
                min_detect(data = tuned_2_wt_cycle, pts = pts)
          }

          if (peak_opt == "max") {
            pts_tuned_1_wt_cycle <-
                max_detect(data = tuned_1_wt_cycle, pts = pts)
            pts_tuned_2_wt_cycle <-
                max_detect(data = tuned_2_wt_cycle, pts = pts)
          }



          dur_gap <-
              (missing_cycle_dur_new * n_cycles_missing) - ((pts_tuned_2_wt_cycle[1, 1] + max(tuned_1_wt_cycle[, 1])) -
                  max(pts_tuned_1_wt_cycle[, 1]))
        }, n_cores = n_cores, run_multicore = run_multicore, combine = "cbind", progress = TRUE)
  } else {
    dur_gaps <- NA
    set.seed(seed_nr)
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


      val_1 <- abs(
          rnorm(1, mean = 1/retracked_period_1[1, 2], sd = 1/retracked_period_1[1, 3]))


      pnorm_val_1 <-  pnorm(
          1 / val_1,
          mean = retracked_period_1[1, 2],
          sd = retracked_period_1[1, 3],
          lower.tail = FALSE
      )

      for (k in 1:nrow(new_curve_1)) {
        new_curve_1[k, 2] <-
            1/(qnorm(pnorm_val_1, mean = retracked_period_1[k, 2], sd = retracked_period_1[k, 3]))
      }

      val_2 <- abs(
          rnorm(1, mean = 1/retracked_period_2[1, 2], sd = 1/retracked_period_2[1, 3]))

      pnorm_val_2 <-  pnorm(
          1 / val_1,
          mean = retracked_period_2[1, 2],
          sd = retracked_period_2[1, 3],
          lower.tail = FALSE
      )

      for (k in 1:nrow(new_curve_2)) {
        new_curve_2[k, 2] <-
            1/(qnorm(pnorm_val_2, mean = retracked_period_2[k, 2], sd = retracked_period_2[k, 3]))
      }


      proxy_1 <- proxies[[proxy_nr]]
      proxy_1 <- proxy_1[proxy_1[, 1] >= new_curve_1[1, 1], ]
      proxy_1 <-
          proxy_1[proxy_1[, 1] <= new_curve_1[nrow(new_curve_1), 1], ]

      proxy_2 <- proxies[[proxy_nr]]
      proxy_2 <- proxy_2[proxy_2[, 1] >= new_curve_2[1, 1], ]
      proxy_2 <-
          proxy_2[proxy_2[, 1] <= new_curve_2[nrow(new_curve_2), 1], ]


      if (tracked_cycle_period_unc_dist == "u") {
        tracked_cycle_period_new <-
            runif(
                1,
                min = tracked_cycle_period - tracked_cycle_period_unc,
                max = tracked_cycle_period + tracked_cycle_period_unc
            )
      }

      if (tracked_cycle_period_unc_dist == "n") {
        tracked_cycle_period_new <-
            rnorm(1, mean = tracked_cycle_period, sd = tracked_cycle_period_unc)
      }

      if (tracked_cycle_period == missing_cycle_dur) {
        missing_cycle_dur_new <- tracked_cycle_period_new
      } else{
        if (missing_cycle_unc_dist == "u") {
          missing_cycle_dur_new <-
              runif(1,
                  min = missing_cycle_dur - missing_cycle_unc,
                  max = tracked_cycle_period + missing_cycle_unc)
        }
        if (missing_cycle_unc_dist == "n") {
          missing_cycle_dur_new <-
              rnorm(1, mean = missing_cycle_dur, sd = missing_cycle_unc)
        }
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
          analyse_wavelet(
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
          analyse_wavelet(
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
        pts_tuned_1_wt_cycle <-
            min_detect(data = tuned_1_wt_cycle, pts = pts)
        pts_tuned_2_wt_cycle <-
            min_detect(data = tuned_2_wt_cycle, pts = pts)
      }

      if (peak_opt == "max") {
        pts_tuned_1_wt_cycle <-
            max_detect(data = tuned_1_wt_cycle, pts = pts)
        pts_tuned_2_wt_cycle <-
            max_detect(data = tuned_2_wt_cycle, pts = pts)
      }



      dur_gap <-
          (missing_cycle_dur_new * n_cycles_missing) - ((pts_tuned_2_wt_cycle[1, 1] + max(tuned_1_wt_cycle[, 1])) -
              max(pts_tuned_1_wt_cycle[, 1]))


      dur_gaps <- c(dur_gaps, dur_gap)
    }

    dur_gaps <- dur_gaps[-c(1)]

  }

  return(dur_gaps)
}
