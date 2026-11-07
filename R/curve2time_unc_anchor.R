#' @title Convert the re-tracked curve results to a
#'  depth time space with uncertainty
#'
#' @description Converts the re-tracked curve results from
#' \code{\link{retrack_wt_MC}} function to a depth time space using an anchor date
#' while also taking into account the uncertainty of the tracked astronomical cycle
#'
#'
#'@param age_constraint Age constraints (anchor points, e.g. radiometric dates) for the modelling run: a data frame
#'with 7 columns and one row per constraint. Column 1: ID name (character); column 2: age (numeric, usually in kyr,
#'in the unit of the output ages); column 3: age uncertainty in the same unit (the standard deviation for a normal
#'distribution, the full width of the range for a uniform one); column 4: distribution of the age, \code{"n"}
#'(normal) or \code{"u"} (uniform); column 5: depth of the constraint (in the depth unit of
#'\code{tracked_cycle_curve}); column 6: depth uncertainty (full width of the range for a uniform distribution);
#'column 7: distribution of the depth, \code{"u"} (uniform), \code{"n"} (normal) or \code{"t"} (trapezoidal, with
#'column 5 then holding the string \code{"min mode1 mode2 max"}). Required (\code{Default=NULL} gives an error).
#' @param tracked_cycle_curve Re-tracked cycle with its uncertainty (matrix or data frame, three columns): depth in
#' column 1 (metres), the mean tracked frequency in cycles per metre in column 2 and its standard deviation (1 sigma)
#' in column 3, e.g. the output of \code{\link{retrack_wt_MC}}. Rows with \code{NA} are dropped and the curve is
#' interpolated to its median depth spacing. Required (\code{Default=NULL} gives an error).
#'@param tracked_cycle_period Duration of the tracked cycle in kyr (single numeric value, > 0, e.g. 405 for the long
#'eccentricity cycle), used to turn the tracked period into a sedimentation rate with \code{\link{curve2time}}.
#'Required (\code{Default=NULL} gives an error).
#'@param tracked_cycle_period_unc Uncertainty of \code{tracked_cycle_period} in kyr (single numeric value, >= 0):
#'the standard deviation when \code{tracked_cycle_period_unc_dist="n"}, the half-width of the range when
#'\code{"u"}. Required (\code{Default=NULL} gives an error).
#'@param tracked_cycle_period_unc_dist Distribution from which the duration of the tracked cycle is drawn in every
#'simulation (character): \code{"n"} normal or \code{"u"} uniform. \code{Default="n"}.
#'@param n_simulations Number of anchored age models to simulate (integer, > 0). \code{Default=20}.
#'@param gap_constraints Hiatuses (gaps) to insert in the age models: a data frame with 4 columns and one row per
#'gap. Column 1: duration of the gap (kyr); column 2: its uncertainty (kyr; the standard deviation for a normal
#'distribution, the half-width of the range for a uniform one; negative draws are set to 0); column 3: depth of the
#'gap (depth unit of \code{tracked_cycle_curve}), below which all ages are shifted by the drawn duration; column 4:
#'distribution, \code{"n"} (normal) or \code{"u"} (uniform). \code{NULL} inserts no gaps. \code{Default=NULL}.
#'@param max_runs Number of rejected simulations (in which the anchored astronomical ages do not match the age
#'constraints) after which the age constraint that disagrees most often is dropped (integer, > 0).
#'\code{Default=1000}.
#'@param keep_nr Minimum number of age constraints that are kept when constraints are dropped after
#'\code{max_runs} rejections (integer, >= 1). \code{Default=2}.
#'@param run_multicore Run the simulations in parallel on several cores (logical, TRUE or FALSE). Note that only the
#'parallel version also checks the tuned proxy data for the cycles in \code{cycles_check}. \code{Default=FALSE}.
#' @param n_cores Number of cores used when \code{run_multicore=TRUE} (integer); \code{NULL} uses an automatic, safe
#' number (physical cores minus one, at most 4, and limited by the free memory), which can be changed with
#' \code{options(waverider.n_cores = ...)}. Ignored when \code{run_multicore=FALSE}. Results do not depend on the
#' number of cores. \code{Default=NULL}.
#'@param verbose Show a progress bar (logical, TRUE or FALSE). \code{Default=FALSE}.
#'@param genplot Plot the age constraints (with their 2-sigma range) and the mean simulated age model with its
#'2-sigma envelope against depth (logical, TRUE or FALSE). \code{Default=FALSE}.
#'@param keep_all_time_curves Also return all generated age curves, including the ones rejected during the
#'modelling run, as the last element of the output list (logical, TRUE or FALSE). \code{Default=FALSE}.
#'@param proxy_data Proxy records used to check that the astronomical cycles are preserved after tuning (list of
#'matrices or data frames, each with depth in column 1, in the unit of \code{tracked_cycle_curve}, and the proxy in
#'column 2); in every run one record is picked at random, tuned with the simulated age model and analysed with
#'\code{\link{analyse_wavelet}}. Only used when \code{run_multicore=TRUE}. \code{Default=NULL}.
#'@param cycles_check Periods in kyr of the astronomical cycles (numeric vector) that must appear as peaks in the
#'average wavelet spectrum of the tuned \code{proxy_data} for a simulation to be accepted. Only used when
#'\code{run_multicore=TRUE}. \code{Default=NULL}.
#'@param uncer_cycles_check Tolerance of each period in \code{cycles_check}, in kyr (numeric vector of the same
#'length): a peak counts if it lies within \code{cycles_check +/- uncer_cycles_check}. Only used when
#'\code{run_multicore=TRUE}. \code{Default=NULL}.
#' @param dj Spacing of the periods of the wavelet transform of the tuned proxy data, in octaves (numeric, > 0;
#' e.g. \code{1/200} gives 200 periods per doubling); smaller values give a finer period axis but take longer.
#' Passed to \code{\link{analyse_wavelet}}; only used when \code{run_multicore=TRUE}. \code{Default=1/200}.
#' @param lowerPeriod Shortest period of the wavelet transform of the tuned proxy data, in kyr (numeric, > 0; at
#' least twice the sample spacing in time). Passed to \code{\link{analyse_wavelet}}; only used when
#' \code{run_multicore=TRUE}. \code{Default=1}.
#' @param upperPeriod Longest period of the wavelet transform of the tuned proxy data, in kyr (numeric, larger than
#' \code{lowerPeriod}). Passed to \code{\link{analyse_wavelet}}; only used when \code{run_multicore=TRUE}.
#' \code{Default=4600}.
#' @param omega_nr Central angular frequency (omega_0) of the Morlet wavelet used on the tuned proxy data (numeric,
#' > 0; larger values give a better period resolution and a poorer depth resolution). Passed to
#' \code{\link{analyse_wavelet}}; only used when \code{run_multicore=TRUE}. \code{Default=6}.
#'@param seed_nr Seed of the random number generator (integer), so that the Monte Carlo simulations can be
#'reproduced. \code{Default=1337}.
#'@param dir Direction of time with depth (logical): \code{TRUE} time increases with depth (younger at the top),
#'\code{FALSE} time decreases with depth (the age model is reversed). \code{Default=TRUE}.
#'
#' @author
#'Part of the code is based on the \link[astrochron]{sedrate2time}
#'function of the 'astrochron' R package
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#'@examples
#'\dontrun{
#'
#'
#'Bisciaro_al <- Bisciaro_XRF[, c(1, 61)]
#'Bisciaro_al <-
#'  astrochron::sortNave(Bisciaro_al, verbose = FALSE, genplot = FALSE)
#'Bisciaro_al <-
#'  astrochron::linterp(Bisciaro_al,
#'                      dt = 0.01,
#'                      verbose = FALSE,
#'                      genplot = FALSE)
#'Bisciaro_al <- Bisciaro_al[Bisciaro_al[, 1] > 2, ]
#'
#'Bisciaro_al_wt <-
#'  analyse_wavelet(
#'    data = Bisciaro_al,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'# Bisciaro_al_wt_track <-
#'#   track_period_wavelet(
#'#     astro_cycle = 110,
#'#     wavelet = Bisciaro_al_wt,
#'#     n.levels = 100,
#'#     periodlab = "Period (metres)",
#'#     x_lab = "depth (metres)"
#'#   )
#'#
#'# Bisciaro_al_wt_track <- completed_series(
#'#   wavelet = Bisciaro_al_wt,
#'#   tracked_curve = Bisciaro_al_wt_track,
#'#   period_up = 1.2,
#'#   period_down = 0.8,
#'#   extrapolate = TRUE,
#'#   genplot = FALSE,
#'#   keep_editable = FALSE
#'# )
#'#
#'# Bisciaro_al_wt_track <-
#'#   loess_auto(
#'#     time_series = Bisciaro_al_wt_track,
#'#     genplot = FALSE,
#'#     print_span = FALSE,
#'#     keep_editable = FALSE
#'#   )
#'
#'
#'
#'Bisciaro_ca <- Bisciaro_XRF[, c(1, 55)]
#'Bisciaro_ca <-
#'  astrochron::sortNave(Bisciaro_ca, verbose = FALSE, genplot = FALSE)
#'Bisciaro_ca <-
#'  astrochron::linterp(Bisciaro_ca,
#'                      dt = 0.01,
#'                      verbose = FALSE,
#'                      genplot = FALSE)
#'Bisciaro_ca <- Bisciaro_ca[Bisciaro_ca[, 1] > 2, ]
#'
#'Bisciaro_ca_wt <-
#'  analyse_wavelet(
#'    data = Bisciaro_ca,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'
#'#
#'# Bisciaro_ca_wt_track <-
#'#   track_period_wavelet(
#'#     astro_cycle = 110,
#'#     wavelet = Bisciaro_ca_wt,
#'#     n.levels = 100,
#'#     periodlab = "Period (metres)",
#'#     x_lab = "depth (metres)"
#'#   )
#'#
#'# Bisciaro_ca_wt_track <- completed_series(
#'#   wavelet = Bisciaro_ca_wt,
#'#   tracked_curve = Bisciaro_ca_wt_track,
#'#   period_up = 1.2,
#'#   period_down = 0.8,
#'#   extrapolate = TRUE,
#'#   genplot = FALSE,
#'#   keep_editable = FALSE
#'# )
#'#
#'# Bisciaro_ca_wt_track <-
#'#   loess_auto(
#'#     time_series = Bisciaro_ca_wt_track,
#'#     genplot = FALSE,
#'#     print_span = FALSE,
#'#     keep_editable = FALSE
#'#   )
#'
#'
#'Bisciaro_sial <- Bisciaro_XRF[, c(1, 64)]
#'Bisciaro_sial <-
#'  astrochron::sortNave(Bisciaro_sial, verbose = FALSE, genplot = FALSE)
#'Bisciaro_sial <-
#'  astrochron::linterp(Bisciaro_sial,
#'                      dt = 0.01,
#'                      verbose = FALSE,
#'                      genplot = FALSE)
#'Bisciaro_sial <- Bisciaro_sial[Bisciaro_sial[, 1] > 2, ]
#'
#'Bisciaro_sial_wt <-
#'  analyse_wavelet(
#'    data = Bisciaro_sial,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'#Bisciaro_sial_wt_track <-
#'#   track_period_wavelet(
#'#     astro_cycle = 110,
#'#     wavelet = Bisciaro_sial_wt,
#'#     n.levels = 100,
#'#     periodlab = "Period (metres)",
#'#     x_lab = "depth (metres)"
#'#   )
#'#
#'#
#'# Bisciaro_sial_wt_track <- completed_series(
#'#   wavelet = Bisciaro_sial_wt,
#'#   tracked_curve = Bisciaro_sial_wt_track,
#'#   period_up = 1.2,
#'#   period_down = 0.8,
#'#   extrapolate = TRUE,
#'#   genplot = FALSE,
#'#   keep_editable = FALSE
#'# )
#'#
#'# Bisciaro_sial_wt_track <-
#'#   loess_auto(
#'#     time_series = Bisciaro_sial_wt_track,
#'#     genplot = FALSE,
#'#     print_span = FALSE,
#'#     keep_editable = FALSE
#'#   )
#'
#'
#'Bisciaro_Mn <- Bisciaro_XRF[, c(1, 46)]
#'Bisciaro_Mn <-
#'  astrochron::sortNave(Bisciaro_Mn, verbose = FALSE, genplot = FALSE)
#'Bisciaro_Mn <-
#'  astrochron::linterp(Bisciaro_Mn,
#'                      dt = 0.01,
#'                      verbose = FALSE,
#'                      genplot = FALSE)
#'Bisciaro_Mn <- Bisciaro_Mn[Bisciaro_Mn[, 1] > 2, ]
#'
#'Bisciaro_Mn_wt <-
#'  analyse_wavelet(
#'    data = Bisciaro_Mn,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'# Bisciaro_Mn_wt_track <-
#'#   track_period_wavelet(
#'#     astro_cycle = 110,
#'#     wavelet = Bisciaro_Mn_wt,
#'#     n.levels = 100,
#'#     periodlab = "Period (metres)",
#'#     x_lab = "depth (metres)"
#'#   )
#'#
#'#
#'# Bisciaro_Mn_wt_track <- completed_series(
#'#   wavelet = Bisciaro_Mn_wt,
#'#   tracked_curve = Bisciaro_Mn_wt_track,
#'#   period_up = 1.2,
#'#   period_down = 0.8,
#'#   extrapolate = TRUE,
#'#   genplot = FALSE,
#'#   keep_editable = FALSE
#'# )
#'# Bisciaro_Mn_wt_track <-
#'#   loess_auto(
#'#     time_series = Bisciaro_Mn_wt_track,
#'#     genplot = FALSE,
#'#     print_span = FALSE,
#'#     keep_editable = FALSE
#'#   )
#'
#'Bisciaro_Mg <- Bisciaro_XRF[, c(1, 71)]
#'Bisciaro_Mg <-
#'  astrochron::sortNave(Bisciaro_Mg, verbose = FALSE, genplot = FALSE)
#'Bisciaro_Mg <-
#'  astrochron::linterp(Bisciaro_Mg,
#'                      dt = 0.01,
#'                      verbose = FALSE,
#'                      genplot = FALSE)
#'Bisciaro_Mg <- Bisciaro_Mg[Bisciaro_Mg[, 1] > 2, ]
#'
#'Bisciaro_Mg_wt <-
#'  analyse_wavelet(
#'    data = Bisciaro_Mg,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'# Bisciaro_Mg_wt_track <-
#'#   track_period_wavelet(
#'#     astro_cycle = 110,
#'#     wavelet = Bisciaro_Mg_wt,
#'#     n.levels = 100,
#'#     periodlab = "Period (metres)",
#'#     x_lab = "depth (metres)"
#'#   )
#'#
#'#
#'# Bisciaro_Mg_wt_track <- completed_series(
#'#   wavelet = Bisciaro_Mg_wt,
#'#   tracked_curve = Bisciaro_Mg_wt_track,
#'#   period_up = 1.2,
#'#   period_down = 0.8,
#'#   extrapolate = TRUE,
#'#   genplot = FALSE,
#'#   keep_editable = FALSE
#'# )
#'#
#'# Bisciaro_Mg_wt_track <-
#'#   loess_auto(
#'#     time_series = Bisciaro_Mg_wt_track,
#'#     genplot = FALSE,
#'#     print_span = FALSE,
#'#     keep_editable = FALSE
#'#   )
#'
#'
#'
#'
#'wt_list_bisc <- list(Bisciaro_al_wt,
#'                     Bisciaro_ca_wt,
#'                     Bisciaro_sial_wt,
#'                     Bisciaro_Mn_wt,
#'                     Bisciaro_Mg_wt)
#'
#'
#'data_track_bisc <- cbind(
#'  Bisciaro_al_wt_track[, 2],
#'  Bisciaro_ca_wt_track[, 2],
#' Bisciaro_sial_wt_track[, 2],
#'  Bisciaro_Mn_wt_track[, 2],
#'  Bisciaro_Mg_wt_track[, 2]
#')
#'
#'x_axis_bisc <- Bisciaro_al_wt_track[, 1]
#'
#'bisc_retrack <- retrack_wt_MC(
#'wt_list = wt_list_bisc,
#'  data_track = data_track_bisc,
#'  x_axis = x_axis_bisc,
#'  nr_simulations = 500,
#'  seed_nr = 1337,
#'  verbose = TRUE,
#'  genplot = FALSE,
#'  keep_editable = FALSE,
#'  create_GIF = FALSE,
#'  plot_GIF = FALSE,
#'  width_plt =  600,
#'  height_plt = 450,
#'  period_up  =  1.5,
#'  period_down = 0.5,
#'  plot.COI = TRUE,
#'  n.levels = 100,
#'  palette_name = "rainbow",
#'  color_brewer = "grDevices",
#'periodlab = "Period (metres)",
#'x_lab = "depth (metres)",
#'add_avg = FALSE,
#'time_dir = TRUE,
#'file_name = "TEST",
#'run_multicore = TRUE,
#'output = 5,
#'  n_imgs = 50,
#'  plot_horizontal = TRUE,
#'  empty_folder = FALSE
#')
#'
#'proxy_list_bisc <- list(Bisciaro_al,
#'                     Bisciaro_ca,
#'                     Bisciaro_sial,
#'                     Bisciaro_Mn,
#'                     Bisciaro_Mg)
#'
#'
#'
#'
#'id <- c("CCT18_322", "CCT18_315", "CCT18_311")
#'ages <- c(20158, 20575, 20857)
#'ageSds <- c(28, 40, 34)
#'ages_unc_dist <- c("n", "n", "n")
#'position <- c(13.3, 7.25, 3.2)
#'anchor_thick <- c(0.2, 0.1, 0.1)
#'anchor_thick_unc_dist <- c("u", "u", "u")
#'
#'ash_Bisc <-
#'  as.data.frame(
#'    cbind(
#'      id,
#'      ages,
#'      ageSds,
#'      ages_unc_dist,
#'      position,
#'      anchor_thick,
#'      anchor_thick_unc_dist
#'    )
#'  )
#'
#'gap_dur = c(10, 20)
#'gap_unc = c(3, 10)
#'gap_depth = c(2.5, 9)
#'gap_unc_dist = c("n", "n")
#'
#'
#'gap_constraints_Bisc <-
#'  as.data.frame(cbind(gap_dur, gap_unc, gap_depth, gap_unc_dist))
#'
#'cycles_checks <- c(110,40,22)
#'uncer_cycles_checks <- c(20,5,7)
#'
#'curve2time_unc_anchor_res <-
#'  curve2time_unc_anchor(
#'  age_constraint =  ash_Bisc,
#'   tracked_cycle_curve =  bisc_retrack,
#'   tracked_cycle_period = 110,
#'    tracked_cycle_period_unc = ((135 - 110) + (110 - 95)) / 2,
#'   tracked_cycle_period_unc_dist = "n",
#'    n_simulations = 20,
#'    gap_constraints = gap_constraints_Bisc,
#'    proxy_data = proxy_list_bisc,
#'    cycles_check = NULL,
#'    uncer_cycles_check = NULL,
#'    cycles_check = cycles_checks,
#'    uncer_cycles_check = uncer_cycles_checks,
#'   max_runs = 1000,
#'    run_multicore = FALSE,
#'    verbose = FALSE,
#'    genplot = FALSE,
#'    keep_nr = 2,
#'    keep_all_time_curves = FALSE,
#'    dj = 1/200,
#'    lowerPeriod =1,
#'    upperPeriod =2500,
#'    omega_nr = 6,
#'    seed_nr=1337,
#'    dir=TRUE
#'  )
#'
#'}
#' @return
#'The output is a list of 3 or 4 elements
#'if keep_all_time_curves is set to TRUE
#'then the list consist of the x-axis, all the fitted curves in a matrix format,
#'the astrochronologically fitted age of the anchor, all the generated depth time curves
#'if keep_all_time_curves is set to TRUE then the list consists of the x-axis,
#' all the fitted curves in a matrix format and the astrochronologically fitted age of the anchor
#'If \code{genplot=TRUE} then 3 plots stacked on top of each other will be plotted.
#'Plot 1: the original data set.
#'Plot 2: the depth time plot.
#'Plot 3: the data set in the time domain.
#'#'
#'
#' @export
#' @importFrom astrochron sedrate2time
#' @importFrom stats quantile
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats pnorm
#' @importFrom astrochron linterp
#' @importFrom stats quantile
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics text
#' @importFrom graphics box
#' @importFrom graphics polygon
#' @importFrom graphics title
#' @importFrom grDevices rgb
#' @importFrom astrochron mtm
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats runif


curve2time_unc_anchor <- function (age_constraint = NULL, tracked_cycle_curve = NULL,
    tracked_cycle_period = NULL, tracked_cycle_period_unc = NULL,
    tracked_cycle_period_unc_dist = "n", n_simulations = 20,
    gap_constraints = NULL, proxy_data = NULL, cycles_check = NULL,
    uncer_cycles_check = NULL, max_runs = 1000, run_multicore = FALSE,
    n_cores = NULL,
    verbose = FALSE, genplot = FALSE, keep_nr = 2, keep_all_time_curves = FALSE,
    dj = 1/200, lowerPeriod = 1, upperPeriod = 4600, omega_nr = 6,
    seed_nr = 1337, dir = TRUE)
{
  dat <- as.matrix(tracked_cycle_curve[, ])
  dat <- na.omit(dat)
  dat <- dat[order(dat[, 1], na.last = NA, decreasing = F), ]
  npts <- length(dat[, 1])
  start <- dat[1, 1]
  end <- dat[length(dat[, 1]), 1]
  x1 <- dat[1:(npts - 1), 1]
  x2 <- dat[2:(npts), 1]
  dx = x2 - x1
  dt = median(dx)
  xout <- seq(start, end, by = dt)
  npts <- length(xout)
  interp <- approx(dat[, 1], dat[, 2], xout, method = "linear", n = npts)
  interp_2 <- approx(dat[, 1], dat[, 3], xout, method = "linear", n = npts)
  tracked_cycle_curve_2 <- cbind(interp[[1]], interp[[2]], interp_2[[2]])
  out <- matrix(
      data = NA,
      nrow = nrow(tracked_cycle_curve_2),
      ncol = n_simulations
  )
  multi_tracked <- tracked_cycle_curve_2
  age_curve <- tracked_cycle_curve_2[, c(1, 2)]
  x_axis <- tracked_cycle_curve_2[, c(1)]
  res_matrix <- matrix(data = NA,
      nrow = nrow(age_curve),
      ncol = n_simulations)
  if (verbose == TRUE) {
    pb <- txtProgressBar(max = n_simulations, style = 3)
  }
  set.seed(seed_nr)
  if (nrow(age_constraint) > 0) {
    if (run_multicore == FALSE) {
      res_list <- list()
      for (i in 1:n_simulations) {
        check_astro_age <- matrix(
            data = NA,
            ncol = 1,
            nrow = nrow(age_constraint)
        )
        anchor_depth <- matrix(
            data = NA,
            ncol = 1,
            nrow = nrow(age_constraint)
        )
        anchor_age <- matrix(
            data = NA,
            ncol = 1,
            nrow = nrow(age_constraint)
        )
        if (keep_all_time_curves == TRUE) {
          time_curve_comb <- matrix(data = NA,
              ncol = 0,
              nrow = length(x_axis))
        }
        for (klm in 1:nrow(age_constraint)) {
          if (age_constraint[klm, 4] == "u") {
            # the age of the constraint (column 3 is the full width of the range), not the duration of the cycle
            check_astro_age[klm] <- runif(
                1,
                min = as.numeric(age_constraint[klm, 2]) - as.numeric(age_constraint[klm, 3]) / 2,
                max = as.numeric(age_constraint[klm, 2]) + as.numeric(age_constraint[klm, 3]) / 2
            )
          }
          if (age_constraint[klm, 4] == "n") {
            check_astro_age[klm] <- rnorm(1,
                mean = as.numeric(age_constraint[klm, 2]),
                sd = as.numeric(age_constraint[klm, 3]))
          }
        }
        anchor_astr <- 1
        anchor_radio <- 0
        sel_rws <- seq(from = 1,
            to = nrow(age_constraint),
            by = 1)
        runs <- 0
        anchor_diff <- matrix(data = NA,
            ncol = 0,
            nrow = length(sel_rws))
        while (anchor_astr > anchor_radio) {
          age_constraint_2 <- age_constraint[c(sel_rws), ]
          check_astro_age_2 <- check_astro_age[c(sel_rws), ]
          anchor_depth_2 <- anchor_depth[c(sel_rws), ]
          anchor_age_2 <- anchor_age[c(sel_rws), ]
          validator <- 1
          while (validator == 1) {
            new_curve <- multi_tracked[, c(1, 2)]
            # standard-normal draw (as rnorm(1, mean, sd)): no Inf from qnorm(1) when the uncertainty is 0
            pnorm_val <- pnorm(rnorm(1))
            for (j in 1:nrow(new_curve)) {
              new_curve[j, 2] <- 1 / (qnorm(pnorm_val, mean = multi_tracked[j, 2], sd = multi_tracked[j, 3]))
            }
            if (tracked_cycle_period_unc_dist == "u") {
              tracked_cycle_period_new <- runif(
                  1,
                  min = tracked_cycle_period -
                  tracked_cycle_period_unc,
                  max = tracked_cycle_period +
                  tracked_cycle_period_unc
              )
            }
            if (tracked_cycle_period_unc_dist == "n") {
              tracked_cycle_period_new <- rnorm(1, mean = tracked_cycle_period, sd = tracked_cycle_period_unc)
            }
            time_curve <- WaverideR::curve2time(
                tracked_cycle_curve = new_curve,
                tracked_cycle_period = tracked_cycle_period_new,
                genplot = FALSE,
                keep_editable = FALSE
            )
            dif_mat <- time_curve[2:(nrow(time_curve)), 2] - time_curve[1:(nrow(time_curve) -
                1), 2]
            dif_mat_min <- min(dif_mat)
            if (dif_mat_min > 0) {
              validator <- 0
            }
          }
          if (dir == FALSE) {
            time_curve[, 2] <- max(time_curve[, 2]) -
                time_curve[, 2]
          }
          gaps_dur <- 0
          if (is.null(gap_constraints) == FALSE) {
            gaps_dur <- matrix(
                data = NA,
                nrow = nrow(gap_constraints[, ]),
                ncol = 1
            )
            for (qx in 1:nrow(gap_constraints)) {
              if (gap_constraints[qx, 4] == "u") {
                if (as.numeric(gap_constraints[qx, 1]) -
                    as.numeric(gap_constraints[qx, 2]) <
                    0) {
                  a <- 0
                } else {
                  a <- as.numeric(gap_constraints[qx, 1]) - as.numeric(gap_constraints[qx, 2])
                }
                gap_dur_new <- runif(
                    1,
                    min = a,
                    max = as.numeric(gap_constraints[qx, 1]) + as.numeric(gap_constraints[qx, 2])
                )
              }
              if (gap_constraints[qx, 4] == "n") {
                gap_dur_new <- rnorm(
                    1,
                    mean = as.numeric(gap_constraints[qx, 1]),
                    sd = as.numeric(gap_constraints[qx, 2])
                )
                if (gap_dur_new < 0) {
                  gap_dur_new <- 0
                }
              }
              gaps_dur[qx, 1] <- gap_dur_new
              row_nr <- wr_closest(time_curve[, 1],
                  as.numeric(gap_constraints[qx, 3]),
                  which = TRUE)
              out_3 <- time_curve[, 2]
              out_4 <- out_3[row_nr:length(out_3)]
              if (dir == FALSE) {
                out_4 <- out_4 - gap_dur_new
                time_curve[, 2] <- c(out_3[1:(row_nr -
                    1)], (out_3[row_nr:length(out_3)] -
                    gap_dur_new))
              } else {
                out_4 <- (out_4 + gap_dur_new)
                time_curve[, 2] <- c(out_3[1:(row_nr -
                    1)], (out_3[row_nr:length(out_3)] +
                    gap_dur_new))
              }
            }
          }
          if (keep_all_time_curves == TRUE) {
            time_curve_comb <- cbind(time_curve_comb, time_curve[, 2])
          }
          anchor_depth_2 <- matrix(data = NA,
              ncol = 1,
              nrow = length(sel_rws))
          anchor_age_2 <- matrix(data = NA,
              ncol = 1,
              nrow = length(sel_rws))
          anchor_astro_age <- matrix(data = NA,
              ncol = 1,
              nrow = length(sel_rws))
          for (klm in 1:nrow(age_constraint_2)) {
            if (age_constraint_2[klm, 7] == "u") {
              anchor_depth_new <- runif(
                  1,
                  min = as.numeric(age_constraint_2[klm, 5]) - as.numeric(age_constraint_2[klm, 6]) /
                  2,
                  max = as.numeric(age_constraint_2[klm, 5]) + as.numeric(age_constraint_2[klm, 6]) /
                  2
              )
            }
            if (age_constraint_2[klm, 7] == "n") {
              anchor_depth_new <- rnorm(
                  1,
                  mean = as.numeric(age_constraint_2[klm, 4]),
                  sd = as.numeric(age_constraint_2[klm, 4])
              )
            }
            if (age_constraint_2[klm, 7] == "t") {
              trap_par <- as.numeric(unlist(strsplit(age_constraint_2[klm, 5], " +")))
              anchor_depth_new <- wr_rtrapezoid(
                  1,
                  min = trap_par[1],
                  mode1 = trap_par[2],
                  mode2 = trap_par[3],
                  max = trap_par[3]
              )
            }
            if (age_constraint_2[klm, 4] == "u") {
              anchor_age_new <- runif(
                  1,
                  min = as.numeric(age_constraint_2[klm, 2]) - as.numeric(age_constraint_2[klm, 3]) /
                  2,
                  max = as.numeric(age_constraint_2[klm, 2]) + as.numeric(age_constraint_2[klm, 3]) /
                  2
              )
            }
            if (age_constraint_2[klm, 4] == "n") {
              anchor_age_new <- rnorm(
                  1,
                  mean = as.numeric(age_constraint_2[klm, 2]),
                  sd = as.numeric(age_constraint_2[klm, 3])
              )
            }
            anchor_depth_2[klm] <- anchor_depth_new
            anchor_age_2[klm] <- anchor_age_new
            row_nr <- wr_closest(time_curve[, 1], anchor_depth_2[klm], which = TRUE)
            anchor_astro_age[klm] <- time_curve[row_nr, 2]
          }
          ages_radio <- anchor_age_2
          ages_astro <- anchor_astro_age
          ages_astro_anchored <- ages_astro
          ages_sim <- cbind(ages_radio, ages_astro)
          if(nrow(ages_sim)>=2){
            ages_sim <- ages_sim[order(ages_sim[, 1]), ]
          }
          check_astro_age_2 <- check_astro_age_2[order(check_astro_age_2)]
          time_curve_anchored <- time_curve


          if (ages_sim[1, 2] > ages_sim[nrow(ages_sim), 2]) {
            a <- mean(ages_radio) + mean(ages_astro)
            ages_astro_anchored <- a - ages_astro
            time_curve_anchored[, 2] <- a - time_curve[, 2]
          } else {
            a <- mean(ages_radio) - mean(ages_astro)
            ages_astro_anchored <- a + ages_astro
            time_curve_anchored[, 2] <- a + time_curve[, 2]
          }
          dif_radio <- abs(check_astro_age_2 - ages_radio)
          dif_astro <- abs(check_astro_age_2 - ages_astro_anchored)
          rown_nr <- order(dif_astro, decreasing = TRUE)[1]
          dif_astro_2 <- dif_astro
          dif_astro_2[, ] <- 0
          dif_astro_2[rown_nr] <- 1
          anchor_diff <- cbind(anchor_diff, dif_astro_2)
          if ((sum(sign(dif_astro - dif_radio)) * -1) ==
              nrow(age_constraint_2)) {
            anchor_astro_age <- cbind(age_constraint[c(sel_rws), 1], ages_astro_anchored)
            anchor_astr <- -1
          } else {
            runs <- runs + 1
          }
          if (runs == max_runs) {
            anchor_diff_rw_sum <- rowSums(as.matrix(anchor_diff))
            row_nr <- wr_closest(anchor_diff_rw_sum,
                max(anchor_diff_rw_sum),
                which = TRUE)
            if (length(sel_rws) <= keep_nr) {
              anchor_diff <- anchor_diff[c(sel_rws)]
              runs <- 0
            } else {
              sel_rws <- sel_rws[-c(row_nr)]
              anchor_diff <- anchor_diff[c(sel_rws)]
              runs <- 0
            }
          }
        }
        if (keep_all_time_curves == TRUE) {
          result <- list(time_curve_anchored[, 2],
              anchor_astro_age,
              gaps_dur,
              time_curve_comb)
        } else {
          (result <- list(time_curve_anchored[, 2], anchor_astro_age, gaps_dur))
        }
        res_list <- c(res_list, list(result))
        if (verbose == TRUE) {
          setTxtProgressBar(pb, i)
        }
      }
    }
    if (run_multicore == TRUE) {
      fit <- wr_lapply(seq_len(n_simulations), function(i) {


        check_astro_age <- matrix(
            data = NA,
            ncol = 1,
            nrow = nrow(age_constraint)
        )
        anchor_depth <- matrix(
            data = NA,
            ncol = 1,
            nrow = nrow(age_constraint)
        )
        anchor_age <- matrix(
            data = NA,
            ncol = 1,
            nrow = nrow(age_constraint)
        )
        if (keep_all_time_curves == TRUE) {
          time_curve_comb <- matrix(data = NA,
              ncol = 0,
              nrow = length(x_axis))
        }
        for (klm in 1:nrow(age_constraint)) {
          if (age_constraint[klm, 4] == "u") {
            # the age of the constraint (column 3 is the full width of the range), not the duration of the cycle
            check_astro_age[klm] <- runif(
                1,
                min = as.numeric(age_constraint[klm, 2]) - as.numeric(age_constraint[klm, 3]) / 2,
                max = as.numeric(age_constraint[klm, 2]) + as.numeric(age_constraint[klm, 3]) / 2
            )
          }
          if (age_constraint[klm, 4] == "n") {
            check_astro_age[klm] <- rnorm(1,
                mean = as.numeric(age_constraint[klm, 2]),
                sd = as.numeric(age_constraint[klm, 3]))
          }
        }
        anchor_astr <- 1
        anchor_radio <- 0
        sel_rws <- seq(from = 1,
            to = nrow(age_constraint),
            by = 1)
        runs <- 0
        anchor_diff <- matrix(data = NA,
            ncol = 0,
            nrow = length(sel_rws))
        while (anchor_astr > anchor_radio) {
          age_constraint_2 <- age_constraint[c(sel_rws), ]
          check_astro_age_2 <- check_astro_age[c(sel_rws), ]
          anchor_depth_2 <- anchor_depth[c(sel_rws), ]
          anchor_age_2 <- anchor_age[c(sel_rws), ]
          validator <- 1
          while (validator == 1) {
            new_curve <- multi_tracked[, c(1, 2)]
            # standard-normal draw (as rnorm(1, mean, sd)): no Inf from qnorm(1) when the uncertainty is 0
            pnorm_val <- pnorm(rnorm(1))
            for (j in 1:nrow(new_curve)) {
              new_curve[j, 2] <- 1 / (qnorm(pnorm_val, mean = multi_tracked[j, 2], sd = multi_tracked[j, 3]))
            }
            if (tracked_cycle_period_unc_dist == "u") {
              tracked_cycle_period_new <- runif(
                  1,
                  min = tracked_cycle_period - tracked_cycle_period_unc,
                  max = tracked_cycle_period + tracked_cycle_period_unc
              )
            }
            if (tracked_cycle_period_unc_dist == "n") {
              tracked_cycle_period_new <- rnorm(1, mean = tracked_cycle_period, sd = tracked_cycle_period_unc)
            }
            time_curve <- WaverideR::curve2time(
                tracked_cycle_curve = new_curve,
                tracked_cycle_period = tracked_cycle_period_new,
                genplot = FALSE,
                keep_editable = FALSE
            )
            dif_mat <- time_curve[2:(nrow(time_curve)), 2] - time_curve[1:(nrow(time_curve) -
                1), 2]
            dif_mat_min <- min(dif_mat)
            if (dif_mat_min > 0) {
              validator <- 0
            }
          }
          if (dir == FALSE) {
            time_curve[, 2] <- max(time_curve[, 2]) -
                time_curve[, 2]
          }
          gaps_dur <- 0
          if (is.null(gap_constraints) == FALSE) {
            gaps_dur <- matrix(
                data = NA,
                nrow = nrow(gap_constraints[, ]),
                ncol = 1
            )
            for (qx in 1:nrow(gap_constraints)) {
              if (gap_constraints[qx, 4] == "u") {
                if (as.numeric(gap_constraints[qx, 1]) - as.numeric(gap_constraints[qx, 2]) < 0) {
                  a <- 0
                } else {
                  a <- as.numeric(gap_constraints[qx, 1]) - as.numeric(gap_constraints[qx, 2])
                }
                gap_dur_new <- runif(
                    1,
                    min = a,
                    max = as.numeric(gap_constraints[qx, 1]) + as.numeric(gap_constraints[qx, 2])
                )
              }
              if (gap_constraints[qx, 4] == "n") {
                gap_dur_new <- rnorm(
                    1,
                    mean = as.numeric(gap_constraints[qx, 1]),
                    sd = as.numeric(gap_constraints[qx, 2])
                )
                if (gap_dur_new < 0) {
                  gap_dur_new <- 0
                }
              }
              gaps_dur[qx, 1] <- gap_dur_new
              row_nr <- wr_closest(time_curve[, 1],
                  as.numeric(gap_constraints[qx, 3]),
                  which = TRUE)
              out_3 <- time_curve[, 2]
              out_4 <- out_3[row_nr:length(out_3)]
              if (dir == FALSE) {
                out_4 <- out_4 - gap_dur_new
                time_curve[, 2] <- c(out_3[1:(row_nr -
                    1)], (out_3[row_nr:length(out_3)] -
                    gap_dur_new))
              } else {
                out_4 <- (out_4 + gap_dur_new)
                time_curve[, 2] <- c(out_3[1:(row_nr -
                    1)], (out_3[row_nr:length(out_3)] +
                    gap_dur_new))
              }
            }
          }
          if (keep_all_time_curves == TRUE) {
            time_curve_comb <- cbind(time_curve_comb, time_curve[, 2])
          }
          anchor_depth_2 <- matrix(data = NA,
              ncol = 1,
              nrow = length(sel_rws))
          anchor_age_2 <- matrix(data = NA,
              ncol = 1,
              nrow = length(sel_rws))
          anchor_astro_age <- matrix(data = NA,
              ncol = 1,
              nrow = length(sel_rws))
          for (klm in 1:nrow(age_constraint_2)) {
            if (age_constraint_2[klm, 7] == "u") {
              anchor_depth_new <- runif(
                  1,
                  min = as.numeric(age_constraint_2[klm, 5]) - as.numeric(age_constraint_2[klm, 6]) /
                  2,
                  max = as.numeric(age_constraint_2[klm, 5]) + as.numeric(age_constraint_2[klm, 6]) /
                  2
              )
            }
            if (age_constraint_2[klm, 7] == "n") {
              anchor_depth_new <- rnorm(
                  1,
                  mean = as.numeric(age_constraint_2[klm, 4]),
                  sd = as.numeric(age_constraint_2[klm, 4])
              )
            }
            if (age_constraint_2[klm, 7] == "t") {
              trap_par <- as.numeric(unlist(strsplit(age_constraint_2[klm, 5], " +")))
              anchor_depth_new <- wr_rtrapezoid(
                  1,
                  min = trap_par[1],
                  mode1 = trap_par[2],
                  mode2 = trap_par[3],
                  max = trap_par[3]
              )
            }
            if (age_constraint_2[klm, 4] == "u") {
              anchor_age_new <- runif(
                  1,
                  min = as.numeric(age_constraint_2[klm, 2]) - as.numeric(age_constraint_2[klm, 3]) /
                  2,
                  max = as.numeric(age_constraint_2[klm, 2]) + as.numeric(age_constraint_2[klm, 3]) /
                  2
              )
            }
            if (age_constraint_2[klm, 4] == "n") {
              anchor_age_new <- rnorm(
                  1,
                  mean = as.numeric(age_constraint_2[klm, 2]),
                  sd = as.numeric(age_constraint_2[klm, 3])
              )
            }
            anchor_depth_2[klm] <- anchor_depth_new
            anchor_age_2[klm] <- anchor_age_new
            row_nr <- wr_closest(time_curve[, 1], anchor_depth_2[klm], which = TRUE)
            anchor_astro_age[klm] <- time_curve[row_nr, 2]
          }
          ages_radio <- anchor_age_2
          ages_astro <- anchor_astro_age
          ages_astro_anchored <- ages_astro
          ages_sim <- cbind(ages_radio, ages_astro)
          if(nrow(ages_sim)>=2){
            ages_sim <- ages_sim[order(ages_sim[, 1]), ]
          }
          check_astro_age_2 <- check_astro_age_2[order(check_astro_age_2)]
          time_curve_anchored <- time_curve

          if (ages_sim[1, 2] > ages_sim[nrow(ages_sim), 2]) {
            a <- mean(ages_radio) + mean(ages_astro)
            ages_astro_anchored <- a - ages_astro
            time_curve_anchored[, 2] <- a - time_curve[, 2]
          } else {
            a <- mean(ages_radio) - mean(ages_astro)
            ages_astro_anchored <- a + ages_astro
            time_curve_anchored[, 2] <- a + time_curve[, 2]
          }
          dif_radio <- abs(check_astro_age_2 - ages_radio)
          dif_astro <- abs(check_astro_age_2 - ages_astro_anchored)
          rown_nr <- order(dif_astro, decreasing = TRUE)[1]
          dif_astro_2 <- dif_astro
          dif_astro_2[, ] <- 0
          dif_astro_2[rown_nr] <- 1
          anchor_diff <- cbind(anchor_diff, dif_astro_2)
          # every proxy equally likely (round(runif(1, 1, n)) drew the first and last half as often)
          sel_proxy_nr <- sample.int(length(proxy_data), 1)
          my.data <- proxy_data[[sel_proxy_nr]]
          out <- time_curve_anchored
          completed_series <- na.omit(out)
          yleft_comp <- completed_series[1, 2]
          yright_com <- completed_series[nrow(completed_series), 2]
          app <- approx(
              x = out[, 1],
              y = out[, 2],
              xout = my.data[, 1],
              method = "linear",
              yleft = yleft_comp,
              yright = yright_com
          )
          completed_series <- as.data.frame(cbind(app$y, my.data[, 2]))
          dat <- as.matrix(completed_series)
          dat <- na.omit(dat)
          dat <- dat[order(dat[, 1], na.last = NA, decreasing = F), ]
          npts <- length(dat[, 1])
          start <- dat[1, 1]
          end <- dat[length(dat[, 1]), 1]
          x1 <- dat[1:(npts - 1), 1]
          x2 <- dat[2:(npts), 1]
          dx = x2 - x1
          dt = median(dx)
          xout <- seq(start, end, by = dt)
          npts <- length(xout)
          interp <- approx(dat[, 1], dat[, 2], xout, method = "linear", n = npts)
          completed_series <- as.data.frame(interp)
          wt_res <- WaverideR::analyse_wavelet(
              data = completed_series,
              dj = dj,
              lowerPeriod = lowerPeriod,
              upperPeriod = upperPeriod,
              verbose = FALSE,
              omega_nr = omega_nr
          )
          avg_wt <- cbind(wt_res$Period, wt_res$Power.avg)
          avg_wt <- WaverideR::max_detect(data = avg_wt, pts = 5)
          mtm_res_test <- is.na(avg_wt)
          mtm_res_test <- mtm_res_test[1]
          if (mtm_res_test == TRUE) {
            runs <- runs + 1
          } else if ((sum(sign(dif_astro - dif_radio)) *
              -1) == nrow(age_constraint_2) |
              nrow(age_constraint_2) ==
              1) {
            mtm_per <- avg_wt[, 1]
            high_vals <- cycles_check + uncer_cycles_check
            low_vals <- cycles_check - uncer_cycles_check
            check <- matrix(
                data = NA,
                nrow = length(cycles_check),
                ncol = 1
            )
            for (i in 1:length(cycles_check)) {
              check[i, 1] <- any(mtm_per < high_vals[i] &
                  mtm_per > low_vals[i])
            }
            if (sum(check) == length(cycles_check)) {
              anchor_astro_age <- cbind(age_constraint[c(sel_rws), 1], ages_astro_anchored)
              anchor_astr <- -1
            } else {
              runs <- runs + 1
              print(runs)
            }
          } else {
            runs <- runs + 1
            print(runs)
          }
          if (runs == max_runs) {
            anchor_diff_rw_sum <- rowSums(as.matrix(anchor_diff))
            row_nr <- wr_closest(anchor_diff_rw_sum,
                max(anchor_diff_rw_sum),
                which = TRUE)
            if (length(sel_rws) <= keep_nr) {
              anchor_diff <- anchor_diff[c(sel_rws)]
              runs <- 0
            } else {
              sel_rws <- sel_rws[-c(row_nr)]
              anchor_diff <- anchor_diff[c(sel_rws)]
              runs <- 0
            }
          }
        }
        time_curve_anchored <- time_curve_anchored[, 2]
        if (keep_all_time_curves == TRUE) {
          result <- list(time_curve_anchored,
              anchor_astro_age,
              gaps_dur,
              time_curve_comb)
        } else {
          (result <- list(time_curve_anchored, anchor_astro_age, gaps_dur))
        }
      }, n_cores = n_cores, run_multicore = run_multicore, combine = "list", progress = verbose)
      res_list <- fit
    }
  }
  result_matrix <- matrix(data = NA,
      nrow = nrow(age_curve),
      ncol = 0)
  res_radio <- matrix(data = NA, nrow = 0, ncol = 2)
  colnames(res_radio) <- c("name", "age")
  res_time_runs <- matrix(data = NA,
      nrow = length(x_axis),
      ncol = 0)
  if (is.null(gap_constraints) == FALSE & keep_all_time_curves ==
      FALSE) {
    res_gap <- matrix(data = NA,
        nrow = 0,
        ncol = nrow(gap_constraints))
    colnames(res_gap) <- c(gap_constraints[, 3])
    for (i in 1:length(res_list)) {
      time_curve <- res_list[[i]][[1]]
      new_anchor_dates <- res_list[[i]][[2]]
      new_gap_dur <- res_list[[i]][[3]]
      result_matrix <- cbind(result_matrix, time_curve)
      colnames(new_anchor_dates) <- c("name", "age")
      res_radio <- rbind(res_radio, new_anchor_dates)
      new_gap_dur <- t(new_gap_dur)
      colnames(new_gap_dur) <- c(gap_constraints[, 3])
      res_gap <- rbind(as.matrix(res_gap), new_gap_dur)
      res_radio_split <- split(x = res_radio[, 2], f = as.factor(unique(res_radio[, 1])))
      for (i in 1:length(res_radio_split)) {
        res_radio_split[[i]] <- as.numeric(res_radio_split[[i]])
      }
    }
    result <- list(x_axis, result_matrix, res_gap, res_radio_split)
  }
  if (is.null(gap_constraints) == TRUE & keep_all_time_curves ==
      FALSE) {
    for (i in 1:length(res_list)) {
      time_curve <- res_list[[i]][[1]]
      new_anchor_dates <- res_list[[i]][[2]]
      result_matrix <- cbind(result_matrix, time_curve)
      colnames(new_anchor_dates) <- c("name", "age")
      res_radio <- rbind(res_radio, new_anchor_dates)
      res_radio_split <- split(x = res_radio[, 2], f = as.factor(unique(res_radio[, 1])))
      for (i in 1:length(res_radio_split)) {
        res_radio_split[[i]] <- as.numeric(res_radio_split[[i]])
      }
    }
    result <- list(x_axis, result_matrix, res_radio_split)
  }
  if (is.null(gap_constraints) == FALSE & keep_all_time_curves ==
      TRUE) {
    res_gap <- matrix(data = NA,
        nrow = 0,
        ncol = nrow(gap_constraints))
    colnames(res_gap) <- c(gap_constraints[, 3])
    for (i in 1:length(res_list)) {
      time_curve <- res_list[[i]][[1]]
      new_anchor_dates <- res_list[[i]][[2]]
      new_gap_dur <- res_list[[i]][[3]]
      time_curve_all <- res_list[[i]][[4]]
      result_matrix <- cbind(result_matrix, time_curve)
      colnames(new_anchor_dates) <- c("name", "age")
      res_radio <- rbind(res_radio, new_anchor_dates)
      new_gap_dur <- t(new_gap_dur)
      colnames(new_gap_dur) <- c(gap_constraints[, 3])
      res_gap <- rbind(as.matrix(res_gap), new_gap_dur)
      res_time_runs <- cbind(res_time_runs, time_curve_all)
      res_radio_split <- split(x = res_radio[, 2], f = as.factor(unique(res_radio[, 1])))
      for (i in 1:length(res_radio_split)) {
        res_radio_split[[i]] <- as.numeric(res_radio_split[[i]])
      }
    }
    result <- list(x_axis,
        result_matrix,
        res_gap,
        res_radio_split,
        res_time_runs)
  }
  if (is.null(gap_constraints) == TRUE & keep_all_time_curves ==
      TRUE) {
    for (i in 1:length(res_list)) {
      time_curve <- res_list[[i]][[1]]
      new_anchor_dates <- res_list[[i]][[2]]
      time_curve_all <- res_list[[i]][[4]]
      result_matrix <- cbind(result_matrix, time_curve)
      colnames(new_anchor_dates) <- c("name", "age")
      res_radio <- rbind(res_radio, new_anchor_dates)
      res_time_runs <- cbind(res_time_runs, time_curve_all)
      res_radio_split <- split(x = res_radio[, 2], f = as.factor(unique(res_radio[, 1])))
      for (i in 1:length(res_radio_split)) {
        res_radio_split[[i]] <- as.numeric(res_radio_split[[i]])
      }
    }
    result <- list(x_axis, result_matrix, res_radio_split, res_time_runs)
  }
  if (genplot == TRUE) {
    plot(
        age_constraint[, 5],
        age_constraint[, 2],
        col = "black",
        cex = 2,
        type = "p",
        ylim = c(min(result[[2]]), max(result[[2]])),
        xlim = c(max(x_axis), min(x_axis))
    )
    points(
        age_constraint[, 5],
        as.numeric(age_constraint[, 2]) + 2 * as.numeric(age_constraint[, 3]),
        col = "red",
        cex = 1,
        pch = 6
    )
    points(
        age_constraint[, 5],
        as.numeric(age_constraint[, 2]) - 2 * as.numeric(age_constraint[, 3]),
        col = "blue",
        cex = 1,
        pch = 2
    )
    res <- result[[2]]
    sds <- wr_row_sds(result[[2]])
    lines(x_axis, rowMeans(result[[2]]), col = "black")
    lines(x_axis, rowMeans(result[[2]]) + 2 * sds, col = "red")
    lines(x_axis, rowMeans(result[[2]]) - 2 * sds, col = "blue")
  }
  return(result)
}
