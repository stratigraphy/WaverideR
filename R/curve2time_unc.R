#' @title Convert the re-tracked curve results to a
#'  depth time space with uncertainty
#'
#' @description Converts the re-tracked curve results from
#' \code{\link{retrack_wt_MC}} function to a depth time space while also
#' taking into account the uncertainty of the tracked astronomical cycle
#'
#' @param tracked_cycle_curve  Curve of the cycle tracked using the
#' \code{\link{retrack_wt_MC}} function \cr
#' Any input (matrix or data frame) with 3 columns in which column 1 is the
#' x-axis, column 2 is the  mean tracked frequency (in cycles/metres) column 3
#' 1 standard deviation
#'@param tracked_cycle_period Period of the tracked curve in kyr.
#'@param tracked_cycle_period_unc uncertainty in the period of the tracked cycle
#'@param tracked_cycle_period_unc_dist distribution of the uncertainty of the
#'tracked cycle value need to be either "u" for uniform distribution or
#'"n" for normal distribution  \code{Default="n"}
#'@param n_simulations number of time series to be modeled
#'@param output If output = 1 a matrix with the predicted ages given the input for each run
#'is given. If output = 2 a matrix with 6 columns is generated,
#'the first column is depth/height, the other columns are the quantile
#'(0.025,0.373,0.5,0.6827,0.975) age values of the runs
#'if output = 3 a matrix with 4 columns is generated with the first column
#'being depth/height, column 2 is the mean tracked duration (in kyrs) column 3
#'is mean duration + 1 standard deviation and column 4  is mean duration -  1
#'standard deviation
#'
#' @author
#'Based on the \link[astrochron]{sedrate2time}
#'function of the 'astrochron' R package
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#' @examples
#' \donttest{
#'# Re-track the 110kyr eccentricity cycle in the wavelet scalogram
#'# from the XRF record of the Bisciaro data set of Arts (2014) and then
#'# add generate and age model including uncertainty
#'
#'Bisciaro_al <- Bisciaro_XRF[, c(1, 61)]
#'Bisciaro_al <- astrochron::sortNave(Bisciaro_al,verbose=FALSE,genplot=FALSE)
#'Bisciaro_al <- astrochron::linterp(Bisciaro_al, dt = 0.01,verbose=FALSE,genplot=FALSE)
#'Bisciaro_al <- Bisciaro_al[Bisciaro_al[, 1] > 2, ]
#'
#'Bisciaro_al_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_al,
#'    dj = 1 /200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
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
#'Bisciaro_ca <- Bisciaro_XRF[, c(1, 55)]
#'Bisciaro_ca <- astrochron::sortNave(Bisciaro_ca,verbose=FALSE,genplot=FALSE)
#'Bisciaro_ca <- astrochron::linterp(Bisciaro_ca, dt = 0.01,verbose=FALSE,genplot=FALSE)
#'Bisciaro_ca <- Bisciaro_ca[Bisciaro_ca[, 1] > 2, ]
#'
#'Bisciaro_ca_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_ca,
#'    dj = 1 /200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
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
#'#     keep_editable = FALSE)
#'
#'Bisciaro_sial <- Bisciaro_XRF[,c(1,64)]
#'Bisciaro_sial <- astrochron::sortNave(Bisciaro_sial,verbose=FALSE,genplot=FALSE)
#'Bisciaro_sial <- astrochron::linterp(Bisciaro_sial, dt = 0.01,verbose=FALSE,genplot=FALSE)
#'Bisciaro_sial <- Bisciaro_sial[Bisciaro_sial[, 1] > 2, ]
#'
#'Bisciaro_sial_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_sial,
#'    dj = 1 /200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'# Bisciaro_sial_wt_track <-
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
#'Bisciaro_Mn <- Bisciaro_XRF[,c(1,46)]
#'Bisciaro_Mn <- astrochron::sortNave(Bisciaro_Mn,verbose=FALSE,genplot=FALSE)
#'Bisciaro_Mn <- astrochron::linterp(Bisciaro_Mn, dt = 0.01,verbose=FALSE,genplot=FALSE)
#'Bisciaro_Mn <- Bisciaro_Mn[Bisciaro_Mn[, 1] > 2, ]
#'
#'Bisciaro_Mn_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_Mn,
#'    dj = 1 /200 ,
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
#'Bisciaro_Mg <- Bisciaro_XRF[,c(1,71)]
#'Bisciaro_Mg <- astrochron::sortNave(Bisciaro_Mg,verbose=FALSE,genplot=FALSE)
#'Bisciaro_Mg <- astrochron::linterp(Bisciaro_Mg, dt = 0.01,verbose=FALSE,genplot=FALSE)
#'Bisciaro_Mg <- Bisciaro_Mg[Bisciaro_Mg[, 1] > 2, ]
#'
#'Bisciaro_Mg_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_Mg,
#'    dj = 1 /200 ,
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
#'#     keep_editable = FALSE)
#'
#'
#'
#'
#'wt_list_bisc <- list(Bisciaro_al_wt,
#'                Bisciaro_ca_wt,
#'                Bisciaro_sial_wt,
#'                Bisciaro_Mn_wt,
#'                Bisciaro_Mg_wt)
#'
#'#Instead of tracking, the tracked solution data sets Bisciaro_al_wt_track,
#'#Bisciaro_ca_wt_track, Bisciaro_sial_wt_track, Bisciaro_Mn_wt_track,
#'# Bisciaro_Mn_wt_track and Bisciaro_Mg_wt_track are used
#'
#'data_track_bisc <- cbind(Bisciaro_al_wt_track[,2],
#'                      Bisciaro_ca_wt_track[,2],
#'                      Bisciaro_sial_wt_track[,2],
#'                      Bisciaro_Mn_wt_track[,2],
#'                      Bisciaro_Mg_wt_track[,2])
#'
#'x_axis_bisc <- Bisciaro_al_wt_track[,1]
#'
#'
#'bisc_retrack <- retrack_wt_MC(wt_list = wt_list_bisc,
#'              data_track = data_track_bisc,
#'              x_axis = x_axis_bisc,
#'              nr_simulations = 20,
#'              seed_nr = 1337,
#'              verbose = FALSE,
#'              genplot = FALSE,
#'              keep_editable = FALSE,
#'              create_GIF = FALSE,
#'              plot_GIF = FALSE,
#'              width_plt =  600,
#'              height_plt = 450,
#'             period_up  =  1.5,
#'              period_down = 0.5,
#'              plot.COI = TRUE,
#'              n.levels = 100,
#'              palette_name = "rainbow",
#'              color_brewer = "grDevices",
#'              periodlab = "Period (metres)",
#'              x_lab = "depth (metres)",
#'              add_avg = FALSE,
#'              time_dir = TRUE,
#'              file_name = NULL,
#'              run_multicore = FALSE,
#'              output = 5,
#'              n_imgs = 50,
#'              plot_horizontal = TRUE,
#'              empty_folder = FALSE)
#'
#'bisc_retrack_age_incl_unc <- curve2time_unc(tracked_cycle_curve = bisc_retrack,
#'tracked_cycle_period = 110,
#'tracked_cycle_period_unc = ((135-110)+(110-95))/2,
#'tracked_cycle_period_unc_dist = "n",
#'n_simulations = 20,
#'output = 1)
#'
#'
#'}
#'@return
#'If output = 1 a matrix with the predicted ages given the input for each run
#'is given
#'If output = 2 a matrix with 6 columns is generated, the first column is
#'depth/height, the other columns are the quantile
#'(0.02275, 0.373, 0.5, 0.6827, 0.97725) age values of the runs
#'if output = 3 a matrix with 4 columns is generated with the first column
#'being depth/height, column 2 is the mean tracked duration (in kyrs) column 3
#'is mean duration + 1 standard deviation and column 4  is mean duration -  1
#'standard deviation
#' @export
#' @importFrom astrochron sedrate2time
#' @importFrom stats quantile
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom matrixStats rowSds
#' @importFrom matrixStats colMins
#' @importFrom matrixStats rowQuantiles
#' @importFrom stats pnorm

curve2time_unc <- function(tracked_cycle_curve = NULL,
                           tracked_cycle_period = NULL,
                           tracked_cycle_period_unc = NULL,
                           tracked_cycle_period_unc_dist = "n",
                           n_simulations = NULL,
                           output = 1) {
  dat <- as.matrix(tracked_cycle_curve[,])
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
  interp <- approx(dat[, 1], dat[, 2], xout, method = "linear",
                   n = npts)

  interp_2 <- approx(dat[, 1], dat[, 3], xout, method = "linear",
                     n = npts)

  tracked_cycle_curve_2 <-
    cbind(interp[[1]], interp[[2]], interp_2[[2]])
  out <-
    matrix(
      data = NA,
      nrow = nrow(tracked_cycle_curve_2),
      ncol = n_simulations
    )


  for (i in 1:n_simulations) {
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


    new_curve <- tracked_cycle_curve_2[, c(1, 2)]
    val <-
      rnorm(1, mean = tracked_cycle_curve_2[1, 2], sd = tracked_cycle_curve_2[1, 3])
    pnorm_val <- 1 - pnorm(val,
                           mean = tracked_cycle_curve_2[1, 2],
                           sd = tracked_cycle_curve_2[1, 3],
                           lower.tail = FALSE)

    for (j in 1:nrow(new_curve)) {
      new_curve[j, 2] <-
        1 / (qnorm(pnorm_val, mean = tracked_cycle_curve_2[j, 2], sd = tracked_cycle_curve_2[j, 3]))
    }


    tracked_cycle_curve_3 <- new_curve
    tracked_cycle_curve_3[, 2] <-
      (new_curve[, 2] / (tracked_cycle_period_new / 100))



    dat <- as.matrix(tracked_cycle_curve_3[, c(1, 2)])
    dat <- na.omit(dat)
    dat <- dat[order(dat[, 1], na.last = NA, decreasing = F), ]
    interp <-
      approx(dat[, 1],
             dat[, 2],
             tracked_cycle_curve_2[, 1],
             method = "linear",
             n = npts)
    sedrates <- as.data.frame(interp)
    npts <- length(sedrates[, 1])
    sedrates[1] = sedrates[1] * 100
    sedrates[2] = 1 / sedrates[2]
    dx = sedrates[2, 1] - sedrates[1, 1]
    midptx = (sedrates[2:npts, 1] + sedrates[1:(npts - 1), 1]) / 2
    slope = (sedrates[2:npts, 2] - sedrates[1:(npts - 1), 2]) / dx
    yint = sedrates[2:npts, 2] - (slope * sedrates[2:npts, 1])
    midpty = (slope * midptx) + yint
    hsum = cumsum(midpty * dx)
    hsum = append(0, hsum)
    out[, i] <- hsum
  }

  dif_mat <- out[2:(nrow(out)),] - out[1:(nrow(out) - 1),]
  dif_mat_mins <- colMins(dif_mat)
  out_2 <- out[, dif_mat_mins > 0]

  if (output == 1) {
    res <- out_2

  }

  if (output == 2) {
    res <- tracked_cycle_curve_2[, c(1, 2)]
    res <-
      cbind(res, rowQuantiles(
        out_2,
        probs = c(0.02275, 0.373, 0.5, 0.6827, 0.97725),
        na.rm = TRUE
      ))
  }

  if (output == 3) {
    res <- tracked_cycle_curve_2[, c(1)]
    mean <- rowMeans(out_2)
    sd <- rowSds(out_2)

    res <- cbind(res, mean, sd)
  }

  return(res)

}
