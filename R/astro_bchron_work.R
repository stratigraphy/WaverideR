# library(WaverideR)
# #install.packages("DecomposeR")
# library("astrochron")
# #remove.packages("ggplot2")
# #install.packages("ggplot2")
# #install.packages("vctrs")
# library("vctrs")
# library("ggplot2")
# library("gganimate")
# library("tidyr")
# #install.packages("DecomposeR")
# #library("DecomposeR")
# #install.packages("ggrepel")
# #install.packages("rlang")
# library(rlang)
# library("ggrepel")
# library(tidyverse)
# library(plyr)
# library(plotly)
# #library(ggpubr)
# library(WaveletComp)
# library(WaverideR)
# library(parallel)
# library(foreach)
# library(iterators)
# library(doParallel)
# library(astrochron)
# library(ggplot2)
# library(tidyr)
# library(doSNOW)
# library(progress)
# library(tcltk)
# library(matrixStats)
# library(tidyr)
# library(colorednoise)
# library(stats)
# library(WaveletComp)
# library(reshape2)
# library(viridis)
# library(readxl)
#
# # Re-track the 110kyr eccentricity cycle in the wavelet scalogram
# # from the XRF record of the Bisciaro data set of Arts (2014)
#
# Bisciaro_al <- Bisciaro_XRF[, c(1, 61)]
# Bisciaro_al <-
#   astrochron::sortNave(Bisciaro_al, verbose = FALSE, genplot = FALSE)
# Bisciaro_al <-
#   astrochron::linterp(Bisciaro_al,
#                       dt = 0.01,
#                       verbose = FALSE,
#                       genplot = FALSE)
# Bisciaro_al <- Bisciaro_al[Bisciaro_al[, 1] > 2.8,]
#
# Bisciaro_al_wt <-
#   analyze_wavelet(
#     data = Bisciaro_al,
#     dj = 1 / 200 ,
#     lowerPeriod = 0.01,
#     upperPeriod = 50,
#     verbose = FALSE,
#     omega_nr = 8
#   )
#
#
# Bisciaro_ca <- Bisciaro_XRF[, c(1, 55)]
# Bisciaro_ca <-
#   astrochron::sortNave(Bisciaro_ca, verbose = FALSE, genplot = FALSE)
# Bisciaro_ca <-
#   astrochron::linterp(Bisciaro_ca,
#                       dt = 0.01,
#                       verbose = FALSE,
#                       genplot = FALSE)
# Bisciaro_ca <- Bisciaro_ca[Bisciaro_ca[, 1] > 2.8,]
#
# Bisciaro_ca_wt <-
#   analyze_wavelet(
#     data = Bisciaro_ca,
#     dj = 1 / 200 ,
#     lowerPeriod = 0.01,
#     upperPeriod = 50,
#     verbose = FALSE,
#     omega_nr = 8
#   )
#
# Bisciaro_sial <- Bisciaro_XRF[, c(1, 64)]
# Bisciaro_sial <-
#   astrochron::sortNave(Bisciaro_sial, verbose = FALSE, genplot = FALSE)
# Bisciaro_sial <-
#   astrochron::linterp(Bisciaro_sial,
#                       dt = 0.01,
#                       verbose = FALSE,
#                       genplot = FALSE)
# Bisciaro_sial <- Bisciaro_sial[Bisciaro_sial[, 1] > 2.8,]
#
# Bisciaro_sial_wt <-
#   analyze_wavelet(
#     data = Bisciaro_sial,
#     dj = 1 / 200 ,
#     lowerPeriod = 0.01,
#     upperPeriod = 50,
#     verbose = FALSE,
#     omega_nr = 8
#   )
#
# # Bisciaro_sial_wt_track <-
# #   track_period_wavelet(
# #     astro_cycle = 110,
# #     wavelet = Bisciaro_sial_wt,
# #     n.levels = 100,
# #     periodlab = "Period (metres)",
# #     x_lab = "depth (metres)"
# #   )
# #
# #
#
#
# Bisciaro_Mn <- Bisciaro_XRF[, c(1, 46)]
# Bisciaro_Mn <-
#   astrochron::sortNave(Bisciaro_Mn, verbose = FALSE, genplot = FALSE)
# Bisciaro_Mn <-
#   astrochron::linterp(Bisciaro_Mn,
#                       dt = 0.01,
#                       verbose = FALSE,
#                       genplot = FALSE)
# Bisciaro_Mn <- Bisciaro_Mn[Bisciaro_Mn[, 1] > 2.8,]
#
# Bisciaro_Mn_wt <-
#   analyze_wavelet(
#     data = Bisciaro_Mn,
#     dj = 1 / 200 ,
#     lowerPeriod = 0.01,
#     upperPeriod = 50,
#     verbose = FALSE,
#     omega_nr = 8
#   )
#
# # Bisciaro_Mn_wt_track <-
# #   track_period_wavelet(
# #     astro_cycle = 110,
# #     wavelet = Bisciaro_Mn_wt,
# #     n.levels = 100,
# #     periodlab = "Period (metres)",
# #     x_lab = "depth (metres)"
# #   )
# #
# #
#
# Bisciaro_Mg <- Bisciaro_XRF[, c(1, 71)]
# Bisciaro_Mg <-
#   astrochron::sortNave(Bisciaro_Mg, verbose = FALSE, genplot = FALSE)
# Bisciaro_Mg <-
#   astrochron::linterp(Bisciaro_Mg,
#                       dt = 0.01,
#                       verbose = FALSE,
#                       genplot = FALSE)
# Bisciaro_Mg <- Bisciaro_Mg[Bisciaro_Mg[, 1] > 2.8,]
#
# Bisciaro_Mg_wt <-
#   analyze_wavelet(
#     data = Bisciaro_Mg,
#     dj = 1 / 200 ,
#     lowerPeriod = 0.01,
#     upperPeriod = 50,
#     verbose = FALSE,
#     omega_nr = 8
#   )
#
# # Bisciaro_Mg_wt_track <-
# #   track_period_wavelet(
# #     astro_cycle = 110,
# #     wavelet = Bisciaro_Mg_wt,
# #     n.levels = 100,
# #     periodlab = "Period (metres)",
# #     x_lab = "depth (metres)"
# #   )
# #
# #
#
#
#
#
# Bisciaro_ca_wt <-
#   analyze_wavelet(
#     data = Bisciaro_ca,
#     dj = 1 / 200 ,
#     lowerPeriod = 0.01,
#     upperPeriod = 50,
#     verbose = FALSE,
#     omega_nr = 6
#   )
#
# # Bisciaro_ca_wt_track <-
# #   track_period_wavelet(
# #     astro_cycle = 110,
# #     wavelet = Bisciaro_al_wt,
# #     n.levels = 100,
# #     periodlab = "Period (metres)",
# #     x_lab = "depth (metres)"
# #   )
#
# #  Bisciaro_ca_wt_track <- completed_series(
# #    wavelet = Bisciaro_al_wt,
# #    tracked_curve = Bisciaro_ca_wt_track,
# #    period_up = 1.3,
# #    period_down = 0.7,
# #    extrapolate = TRUE,
# #    genplot = TRUE,
# #    keep_editable = FALSE
# #  )
# #
# #
# #
# #
# #  Bisciaro_ca_wt_track <-
# #    loess_auto(
# #      time_series = Bisciaro_ca_wt_track,
# #      genplot = FALSE,
# #      print_span = FALSE,
# #      keep_editable = FALSE)
# #
# #
# # Bisciaro_ca_wt_track_110 <- extract_signal(
# #   tracked_cycle_curve = Bisciaro_ca_wt_track,
# #   wavelet = Bisciaro_al_wt,
# #   period_up = 1.25,
# #   period_down = 0.8,
# #   add_mean = TRUE,
# #   tracked_cycle_period = 110,
# #   extract_cycle = 110,
# #   tune = FALSE,
# #   plot_residual = FALSE
# # )
#
# # library(astrochron)
# # la11_opt <- getLaskar(sol = "la11")
# #
# # plot(la11_opt, type = "l")
# # la11_opt[, 1] <- la11_opt / 1000
# #
# # la11_opt <- la11_opt[la11_opt[, 1] > 19.6, ]
# # la11_opt <- la11_opt[la11_opt[, 1] < 21.2, ]
#
#
#
# # astro_anchor_points <-
# #   astro_anchor(
# #     astro_solution = la11_opt,
# #     proxy_signal = Bisciaro_ca_wt_track_110,
# #     proxy_min_or_max = "max",
# #     clip_astrosolution = FALSE,
# #     astrosolution_min_or_max = "max",
# #     clip_high = NULL,
# #     clip_low = NULL,
# #     extract_astrosolution = FALSE,
# #     astro_period_up = 1.2,
# #     astro_period_down = 0.8,
# #     astro_period_cycle = NULL,
# #     extract_proxy_signal = FALSE,
# #     proxy_period_up = 1.2,
# #     proxy_period_down = 0.8,
# #     proxy_period_cycle = NULL,
# #     pts = 10,
# #     verbose = FALSE,
# #     time_dir = FALSE
# #   )
# #
# # anchor_points_Bisciaro_al <- astro_anchor_points
# # save(anchor_points_Bisciaro_al, file = "D:/Phd/R/packages/WaverideR/data/anchor_points_Bisciaro_al.RData")
# #
# #
# # plot_astro_anchor(
# #   astro_solution = la11_opt,
# #   proxy_signal = Bisciaro_ca_wt_track_110,
# #   anchor_points = astro_anchor_points,
# #   time_dir = FALSE,
# #   keep_editable = FALSE
# # )
# #
# #
# # Bisciaro_ca_time <-
# #   anchor2time(
# #     anchor_points = astro_anchor_points,
# #     data = Bisciaro_al,
# #     genplot = TRUE,
# #     keep_editable = FALSE
# #   )
# #
# #
# #
# # la11_opt[,2]*100000
# #
# # graphics.off()
# # plot(Bisciaro_ca_time,type="l",xaxt='n')
# # axis(side = 1, at= seq(from=19.60,to=22,by=0.05),las = 2)
# # lines(la11_opt[,1],la11_opt[,2]*500000+20000)
# #
# #
# # Bisciaro_ca_time <- linterp(Bisciaro_ca_time)
# # Bisciaro_ca_time[,1] <- Bisciaro_ca_time[,1]*1000
# #
# #
# # graphics.off()
# # ecc_110 <- taner(dat=Bisciaro_ca_time,fhigh=1/90,flow=1/135,xmax=1/50)
# # prec <- taner(dat=Bisciaro_ca_time,fhigh=1/25,flow=1/17,xmax=1/10)
# # prec_hilbert <- hilbert(prec)
# # prec_hilbert_110 <- taner(dat=prec_hilbert,fhigh=1/90,flow=1/135,xmax=1/50)
# #
# # plot(prec_hilbert_110,type="l")
# # lines(ecc_110[,1],ecc_110[,2]/10-800,col="blue")
# #
# # obl <- taner(dat=Bisciaro_ca_time,fhigh=1/45,flow=1/35,xmax=1/10)
# # obl_hilbert <- hilbert(obl)
# # obl_hilbert_filt <- taner(dat=obl_hilbert,fhigh=1/1200,flow=1/800,xmax=1/400)
# #
# #
# # taner(dat=Bisciaro_ca_time,fhigh=1/505,flow=1/305,xmax=1/100)
# #
# #
# #
# #
# #
# # colnames(Bisciaro_XRF)
# #
# # plot(Bisciaro_XRF[,c(1,26)],type="l")
# # abline(v=c(10.2,6.1,3.2),col="red")
# #
# # plot(Bisciaro_XRF[,c(1,55)],type="l")
# # abline(v=c(10.2,6.1,3.2),col="red")
# #
# #
# # astro_anchor_points_linterp <- linterp(astro_anchor_points[,c(1,2)],0.01)
# #
# # rown_nr <-DescTools::Closest(astro_anchor_points_linterp[, 1],
# #                              10.2, which = TRUE)
# # rown_nr <- rown_nr[1]
# # astro_anchor_points_linterp[rown_nr,2]
# # 20.158 - astro_anchor_points_linterp[rown_nr,2]
# #
# #
# # rown_nr <-DescTools::Closest(astro_anchor_points_linterp[, 1],
# #                              6.15, which = TRUE)
# # rown_nr <- rown_nr[1]
# # astro_anchor_points_linterp[rown_nr,2]
# # 20.572 - astro_anchor_points_linterp[rown_nr,2]
# #
# #
# # rown_nr <-DescTools::Closest(astro_anchor_points_linterp[, 1],
# #                              3.6, which = TRUE)
# # rown_nr <- rown_nr[1]
# # astro_anchor_points_linterp[rown_nr,2]
# # 20.857 - astro_anchor_points_linterp[rown_nr,2]
# #
#
#
#
#
# wt_list_bisc <- list(Bisciaro_al_wt,
#                      Bisciaro_ca_wt,
#                      Bisciaro_sial_wt,
#                      Bisciaro_Mn_wt,
#                      Bisciaro_Mg_wt)
#
# Bisciaro_al_wt_track <-
#   Bisciaro_al_wt_track[Bisciaro_al_wt_track[, 1] > 2.8,]
# Bisciaro_ca_wt_track <-
#   Bisciaro_ca_wt_track[Bisciaro_ca_wt_track[, 1] > 2.8,]
# Bisciaro_sial_wt_track <-
#   Bisciaro_sial_wt_track[Bisciaro_sial_wt_track[, 1] > 2.8,]
# Bisciaro_Mn_wt_track <-
#   Bisciaro_Mn_wt_track[Bisciaro_Mn_wt_track[, 1] > 2.8,]
# Bisciaro_Mg_wt_track <-
#   Bisciaro_Mg_wt_track[Bisciaro_Mg_wt_track[, 1] > 2.8,]
#
#
# data_track_bisc <- cbind(
#   Bisciaro_al_wt_track[, 2],
#   Bisciaro_ca_wt_track[, 2],
#   Bisciaro_sial_wt_track[, 2],
#   Bisciaro_Mn_wt_track[, 2],
#   Bisciaro_Mg_wt_track[, 2]
# )
#
# x_axis_bisc <- Bisciaro_al_wt_track[, 1]
#
# getwd()
# setwd("D:/Phd/R/R_results")
#
#
# library(magick)
# retrack_wt_MC_2 <- function(wt_list = NULL,
#                             data_track = NULL,
#                             x_axis = NULL,
#                             nr_simulations = 50,
#                             seed_nr = 1337,
#                             verbose = FALSE,
#                             genplot = FALSE,
#                             keep_editable = FALSE,
#                             create_GIF = FALSE,
#                             plot_GIF = FALSE,
#                             width_plt =  600,
#                             height_plt = 450,
#                             period_up  =  1.5,
#                             period_down = 0.5,
#                             plot.COI = TRUE,
#                             n.levels = 100,
#                             palette_name = "rainbow",
#                             color_brewer = "grDevices",
#                             periodlab = "Period (metres)",
#                             x_lab = "depth (metres)",
#                             add_avg = FALSE,
#                             time_dir = TRUE,
#                             file_name = "TEST_1",
#                             run_multicore = FALSE,
#                             output = 1,
#                             n_imgs = 50,
#                             plot_horizontal = TRUE,
#                             empty_folder = FALSE) {
#   simulations = nr_simulations
#   img_animated <- NULL
#
#   if (run_multicore == TRUE) {
#     numCores <- detectCores()
#     cl <- parallel::makeCluster(numCores - 2)
#     registerDoSNOW(cl)
#   } else{
#     numCores <- 1
#     cl <- makeCluster(numCores)
#     registerDoSNOW(cl)
#   }
#
#
#   if (empty_folder == TRUE & create_GIF == TRUE) {
#     f <-
#       list.files(
#         file_name,
#         include.dirs = F,
#         full.names = T,
#         recursive = T
#       )
#     file.remove(f)
#   }
#
#
#   if (verbose == TRUE) {
#     pb <- txtProgressBar(max = simulations, style = 3)
#     progress <- function(n)
#       setTxtProgressBar(pb, n)
#     opts <- list(progress = progress)
#   } else{
#     opts = NULL
#   }
#
#   if (create_GIF == TRUE) {
#     dir.create(
#       file_name,
#       showWarnings = TRUE,
#       recursive = FALSE,
#       mode = "0777"
#     )
#   }
#
#   n_curves <- ncol(data_track)
#
#   j <- 1 # needed to assign 1 to ijk to avoid note
#   set.seed(seed_nr)
#
#
#   fit <-
#     foreach (
#       j = 1:simulations,
#       .options.snow = opts,
#       .combine = 'cbind'
#     ) %dopar% {
#       sel_curve <- sample(1:n_curves, 1, replace = F)
#       n <- n_curves
#       x <- runif(n, 0, 1)
#       y <- x / sum(x)
#
#
#       vals <-
#         matrix(rep(t(y), nrow(data_track)),
#                ncol = ncol(data_track),
#                byrow = TRUE)
#       fractions <- data_track * vals
#       fractions <- rowSums(fractions)
#       fractions <- cbind(x_axis, fractions)
#
#       wt_sel <- rlist::list.extract(wt_list, sel_curve)
#
#       completed_curve <- WaverideR::completed_series(
#         wavelet =  wt_sel,
#         tracked_curve =  fractions[, c(1, 2)],
#         period_up  = period_up,
#         period_down  = period_down,
#         extrapolate = TRUE,
#         genplot = FALSE
#       )
#
#       completed_curve <-
#         astrochron::linterp(
#           completed_curve,
#           dt = x_axis[2] - x_axis[1],
#           start = x_axis[1],
#           genplot = FALSE,
#           verbose = FALSE
#         )
#       completed_curve <- completed_curve[1:length(x_axis), ]
#       completed_curve <- WaverideR::loess_auto(completed_curve)
#       completed_curve <- completed_curve[, c(1, 2)]
#
#
#       if (create_GIF == TRUE) {
#         png(
#           filename = paste0(file_name, "/", file_name, "_", j, ".jpeg"),
#           type = "cairo",
#           width = width_plt,
#           height = height_plt
#         )
#
#         WaverideR::plot_wavelet(
#           wavelet = wt_sel,
#           plot.COI = plot.COI,
#           n.levels = n.levels,
#           useRaster = TRUE,
#           palette_name = palette_name,
#           color_brewer = color_brewer,
#           periodlab = periodlab,
#           x_lab = x_lab,
#           add_lines = completed_curve,
#           add_avg = add_avg,
#           dev_new = FALSE,
#           time_dir = time_dir,
#           plot_horizontal = plot_horizontal
#         )
#         dev.off()
#       }
#       completed_curve <- completed_curve[, 2]
#
#     }
#
#   stopCluster(cl)
#
#   sims_2 <- fit
#   sims_mean_2 <- rowMeans(sims_2)
#   sims_2 <-  as.matrix(sims_2)
#   sims_sd_2 <- matrixStats::rowSds(sims_2)
#
#   sed_run <-
#     cbind(x_axis,
#           sims_mean_2,
#           sims_mean_2 - sims_sd_2,
#           sims_mean_2 + sims_sd_2)
#
#   if (genplot == TRUE) {
#     if (keep_editable == FALSE) {
#       oldpar <- par(no.readonly = TRUE)
#       on.exit(par(oldpar))
#     }
#     layout.matrix <- matrix(c(1), nrow = 1, ncol = 1)
#     graphics::layout(mat = layout.matrix,
#                      heights = c(1),
#                      # Heights of the two rows
#                      widths = c(1))
#     par(mar = c(4, 4, 1, 1))
#     plot(
#       x = sed_run[, 1],
#       y = sed_run[, 2],
#       type = "l",
#       ylim = c(min(sed_run[, 3]),
#                max(sed_run[, 4])),
#       col = "green",
#       lwd = 2,
#       xlab = x_lab,
#       ylab = periodlab
#     )
#     lines(
#       x = sed_run[, 1],
#       y = sed_run[, 3],
#       col = "red",
#       lwd = 2
#     )
#     lines(
#       x = sed_run[, 1],
#       y = sed_run[, 4],
#       col = "blue",
#       lwd = 2
#     )
#   }
#
#   if (create_GIF == TRUE) {
#     imgs <- list.files(file_name, full.names = TRUE)
#
#     if (n_imgs > nr_simulations) {
#       n_imgs <- nr_simulations
#     }
#
#     imgs <- imgs[1:n_imgs]
#
#     img_list <- lapply(imgs, image_read)
#     img_joined <- image_join(img_list)
#     img_animated <- image_animate(img_joined, fps = 5)
#
#     if (plot_GIF == TRUE) {
#       img_animated
#     }
#
#     image_write(image = img_animated,
#                 path = paste0(file_name, "/", file_name, ".gif"))
#   }
#
#
#   sed_run[, 3] <- sed_run[, 2] - sed_run[, 3]
#   sed_run <- sed_run[, c(1:3)]
#   colnames(sed_run) <- c("depth", "mean_period", "sd")
#
#
#   if (output == 1) {
#     res <- list(sed_run, fit, img_animated)
#
#   }
#
#   if (output == 2) {
#     res <- list(sed_run, fit)
#
#   }
#
#   if (output == 3) {
#     res <- list(sed_run, img_animated)
#
#   }
#
#
#   if (output == 4) {
#     res <- list(fit, img_animated)
#
#   }
#
#
#   if (output == 5) {
#     res <- sed_run
#
#   }
#
#   if (output == 6) {
#     res <- fit
#
#   }
#
#   if (output == 7) {
#     res <- img_animated
#
#   }
#
#
#   return(res)
# }
#
#
#
#
#
# retrack <- retrack_wt_MC_2(
#   wt_list = wt_list_bisc,
#   data_track = data_track_bisc,
#   x_axis = x_axis_bisc,
#   nr_simulations = 1000,
#   seed_nr = 1337,
#   verbose = TRUE,
#   genplot = TRUE,
#   keep_editable = FALSE,
#   create_GIF = TRUE,
#   plot_GIF = TRUE,
#   width_plt =  600,
#   height_plt = 450,
#   period_up  =  1.5,
#   period_down = 0.5,
#   plot.COI = TRUE,
#   n.levels = 100,
#   palette_name = "rainbow",
#   color_brewer = "grDevices",
#   periodlab = "Period (metres)",
#   x_lab = "depth (metres)",
#   add_avg = FALSE,
#   time_dir = TRUE,
#   file_name = "TEST_1",
#   run_multicore = TRUE,
#   output = 1,
#   n_imgs = 50,
#   plot_horizontal = TRUE,
#   empty_folder = TRUE
# )
#
#
# retrack_1 <- retrack[[1]]
#
# retrack_1
#
#
# plot(retrack_1[, 1], retrack_1[, 2], type = "l")
# lines(retrack_1[, 1], retrack_1[, 2] + retrack_1[, 3], col = "red")
# lines(retrack_1[, 1], retrack_1[, 2] - retrack_1[, 3], col = "blue")
#
#
#
# n_simulations <- 22 * 500
# multi_tracked <- retrack_1
# age_curves <-
#   matrix(data = NA,
#          nrow = nrow(multi_tracked),
#          ncol = n_simulations)
# new_curve <- multi_tracked[, c(1, 2)]
#
#
#
# for (i in 1:ncol(age_curves)) {
#   val <-  rnorm(1, mean = multi_tracked[1, 2], sd = multi_tracked[1, 3])
#   pnorm_val <-
#     pnorm(val, mean = multi_tracked[1, 2], sd = multi_tracked[1, 3])
#   for (j in 1:nrow(new_curve)) {
#     new_curve[j, 2] <-
#       qnorm(pnorm_val, mean = multi_tracked[j, 2], sd = multi_tracked[j, 3])
#   }
#
#
#   time_curve <- WaverideR::curve2time(
#     tracked_cycle_curve = new_curve,
#     tracked_cycle_period = 110,
#     genplot = FALSE,
#     keep_editable = FALSE
#   )
#
#   age_curves[, i] <-
#     time_curve[, 2]
# }#<- max(time_curve[,2])-time_curve[,2] }
#
#
#
#
# str(age_curves)
#
#
# id <- c("CCT18_322", "CCT18_315", "CCT18_311")
# ages <- c(20158, 20575, 20857)
# ageSds <- c(28, 40, 34)
# position <- c(10.2, 6.1, 3.2)
# thickness <- c(0.2, 0.1, 0.1)
# calCurves <- rep("normal", 3)
# outlierProbs <- rep(0.01, 3)
#
#
#
#
# ash <-
#   as.data.frame(cbind(id, ages, ageSds, position, thickness, calCurves, outlierProbs))
#
# ash$postime <- as.numeric(ash$ages) - as.numeric(ash$ages[1])
#
# retrack_1_time <- curve2time(
#   tracked_cycle_curve = multi_tracked[, c(1, 2)],
#   tracked_cycle_period = 110,
#   genplot = TRUE,
#   keep_editable = FALSE
# )
#
# head(ash)
# ash$pos_time_rand <- NA
#
#
# ash
#
# for (i in 1:nrow(ash)) {
#   rown_nr <-
#     DescTools::Closest(retrack_1_time[, 1], as.numeric(ash[i, 4]), which = TRUE)
#   rown_nr <- rown_nr[1]
#   ash[i, 9] <- retrack_1_time[rown_nr, 2]
# }
#
# ash_2 <- ash
# pos_thick_rand <-
#   matrix(data = NA, nrow(ash_2), ncol = n_simulations)
# pos_time_rand <-
#   matrix(data = NA, nrow(ash_2), ncol = n_simulations)
# pos_rownr_rand <-
#   matrix(data = NA, nrow(ash_2), ncol = n_simulations)
# pos_rownr_min <-
#   matrix(data = NA, nrow(ash_2), ncol = n_simulations)
# pos_rownr_plus <-
#   matrix(data = NA, nrow(ash_2), ncol = n_simulations)
#
#
# for (i in 1:ncol(pos_time_rand)) {
#   for (k in 1:nrow(ash_2)) {
#     unc <- as.numeric(ash_2[k, 5])
#     position <-  as.numeric(ash_2[k, 4])
#     rand_pos_plus <-
#       DescTools::Closest(multi_tracked[, 1], position + (unc / 2), which = TRUE)
#     rand_pos_min <-
#       DescTools::Closest(multi_tracked[, 1], position - (unc / 2), which = TRUE)
#     rand_pos_plus <- rand_pos_plus[1]
#     rand_pos_min <- rand_pos_min[1]
#     pos_rownr_min[k, i] <- rand_pos_min
#     pos_rownr_plus[k, i] <- rand_pos_plus
#     pos_thick_rand[k, i] <-
#       age_curves[rand_pos_plus, i] - age_curves[rand_pos_min, i]
#     rown_nr <-
#       DescTools::Closest(multi_tracked[, 1], position, which = TRUE)
#     rown_nr <- rown_nr[1]
#     pos_rownr_rand[k, i] <- rown_nr
#     pos_time_rand[k, i] <- age_curves[rown_nr, i]
#   }
# }
#
#
# ash$outlierProbs <- 1 / 10 ^ 4
#
# j <- 1
#
# numCores <- detectCores()
# cl <- makeCluster(numCores - 2)
# registerDoSNOW(cl)
#
# pb <- txtProgressBar(max = n_simulations, style = 3)
# progress <- function(n)
#   setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
#
#
#
# ages <-
#   foreach::foreach (
#     j = 1:(n_simulations),
#     .options.snow = opts,
#     .errorhandling = "remove",
#     .packages = c("Bchron", "matrixStats"),
#     .combine = "cbind"
#   ) %dopar% {
#     #stretch squeeze age-model with factor N
#     #run Bayesian
#     #save the mean ages of the boundary
#     #repeat x times until happy
#
#
#     #
#     # defaultW <- getOption("warn")
#     # options(warn = -1)
#
#
#     bentonites_out = Bchron::Bchronology(
#       ages = as.numeric(ash_2[, 2]),
#       ageSds = as.numeric(ash_2[, 3]),
#       calCurves = ash_2[, 6],
#       positions = as.numeric(pos_time_rand[, j]),
#       positionThicknesses =
#         as.numeric(pos_thick_rand[, j]),
#       ids = ash_2[, 1],
#       predictPositions = age_curves[, j],
#       outlierProbs = as.numeric(ash_2[, 7]),
#       artificialThickness = 1,
#       allowOutside = TRUE,
#       thin = 4,
#       iterations = 10000,
#       burn = 5000,
#       thetaMhSd = 0.5,
#       muMhSd = 0.1,
#       psiMhSd = 0.1,
#       maxExtrap = 10000,
#       ageScaleVal = 10000,
#       positionEps = 1e-04,
#       positionNormalise = TRUE
#     )
#
#     # str(Altajme_bentonites_out)
#     #plot(Altajme_bentonites_out)
#     #plot(Altajme_bentonites_out$thetaPredict[1,],type="l")
#     #
#     # for (i in 1:ncol(Altajme_bentonites_out$thetaPredict[,])){
#     # lines(Altajme_bentonites_out$thetaPredict[i,]) }
#     #
#     # lines(colMeans(Altajme_bentonites_out$thetaPredict[,]),col="red")
#
#     #save(Altajme_bentonites_out, file="Altajme_bentonites_out.RData")
#     #options(warn = defaultW)
#     #save(Altajme_bentonites_out, file="Altajme_bentonites_out.RData")
#     #load("Altajme_bentonites_out.RData")
#     #str(Altajme_bentonites_out)
#
#     # Altajme_bentonites_out_sum=summary(Altajme_bentonites_out) # Default is for quantiles of ages at predictPosition values
#     # summary(Altajme_bentonites_out, type='convergence') # Check model convergence
#     # summary(Altajme_bentonites_out, type='outliers') # Look at outlier probabilities
#     #
#     # plot(Altajme_bentonites_out,main="Altajme_bentonites",xlab='Age (Ma)',ylab='Depth (m)',las=1)
#     # str(Altajme_bentonites_out)
#     # head(Altajme_bentonites_out)
#     #
#     # write.csv(Altajme_bentonites_out_sum,"Altajme_bentonites_out_sum.csv",row.names = F)
#     # Altajme_bentonites_out_sum2=data.frame(Altajme_bentonites_out_sum$Depth,Altajme_bentonites_out_sum$`50%`)
#     # Altajme_bentonites_out_sum2=delPts(Altajme_bentonites_out_sum2)
#
#     # str(Altajme_bentonites_out)
#     #
#     #
#     # DevChronOut_astro <- c()
#     # DevChronOut_astro_temp <- c()
#     # rhs <-paste("DevChronOut_astro_temp<-Altajme_bentonites_out$thetaPredict",
#     #         sep = "")
#     # eval(parse(text = rhs))
#     # DevChronOut_astro = as.matrix(rbind(DevChronOut_astro, DevChronOut_astro_temp))
#     #
#     # Altajme_bentonites_out$thetaPredict
#     #
#     # col_means <- matrixStats::colMedians(DevChronOut_astro)
#     # mean_ages <- cbind(seq(0,max(age_curves[,]), by = 1), col_means)
#     #
#     #
#     #
#     # new_ages <- matrix(
#     #   data = NA,
#     #   ncol = 1,
#     #   nrow = nrow(pred_age_rand)
#     # )
#     #
#     # for (i in 1:nrow(pred_age_rand)) {
#     #   rown_nr_3 <-DescTools::Closest(mean_ages[, 1], pred_age_rand[i, j], which = TRUE)
#     #   if(is.na(rown_nr_3)){
#     #     new_ages[i, 1] <- NA
#     #   }else{
#     #     rown_nr_3 <- rown_nr_3[1]
#     #   new_ages[i, 1] <- mean_ages[rown_nr_3, 2]}
#     # }
#     #
#     # new_ages <- new_age
#
#
#     bentonites_out_mean <- colMeans(bentonites_out$thetaPredict)
#     # thetaPredict <-  Altajme_bentonites_out$thetaPredict
#     # thetaPredictdiv <- thetaPredict
#     #
#     # for (xyz in 1: nrow(thetaPredict)){
#     #  val <-  sqrt((thetaPredict[xyz,] - Altajme_bentonites_out_mean)^2)
#     #  thetaPredictdiv[xyz,] <- val }
#     #
#     #   thetaPredictdiv_sum <- rowSums(thetaPredictdiv)
#     # rown_nr <-DescTools::Closest(thetaPredictdiv_sum, min(thetaPredictdiv_sum), which = TRUE)
#     #
#     # plot(Altajme_bentonites_out_mean)
#     # lines(thetaPredict[rown_nr,])
#
#
#     bentonites_out <- bentonites_out_mean
#   }
#
#
# stopCluster(cl)
#
#
#
#
#
#
#
# ages_2 <- ages
#
# ages_2 <- as.matrix(ages_2)
#
#
# rown_nr <- DescTools::Closest(new_curve[, 1], 10.2, which = TRUE)
# rown_nr <- rown_nr[1]
# rown_nr
#
# hist(ages_2[rown_nr, ])
# mean(ages_2[rown_nr, ])
# sd(ages_2[rown_nr, ])
#
# rown_nr <- DescTools::Closest(new_curve[, 1], 6.1, which = TRUE)
# rown_nr <- rown_nr[1]
# rown_nr
# hist(ages_2[rown_nr, ])
# mean(ages_2[rown_nr, ])
# sd(ages_2[rown_nr, ])
#
#
# hist(ages_2[c(pos_rownr_min[1, 1]:pos_rownr_plus[1, 1]), ])
# mean(ages_2[c(pos_rownr_min[1, 1]:pos_rownr_plus[1, 1]), ])
# sd(ages_2[c(pos_rownr_min[1, 1]:pos_rownr_plus[1, 1]), ])
#
#
# hist(ages_2[c(pos_rownr_min[2, 1]:pos_rownr_plus[2, 1]), ])
# mean(ages_2[c(pos_rownr_min[2, 1]:pos_rownr_plus[2, 1]), ])
# sd(ages_2[c(pos_rownr_min[2, 1]:pos_rownr_plus[2, 1]), ])
#
# hist(ages_2[c(pos_rownr_min[3, 1]:pos_rownr_plus[3, 1]), ])
# mean(ages_2[c(pos_rownr_min[3, 1]:pos_rownr_plus[3, 1]), ])
# sd(ages_2[c(pos_rownr_min[3, 1]:pos_rownr_plus[3, 1]), ])
#
#
#
# ages_sds <- matrixStats::rowSds(ages_2)
# ages_median <- rowMeans(ages_2)
#
# custom_colour_pal <-
#   colorRampPalette(c('#d7191c', '#fdae61', '#ffffbf', '#abd9e9', '#2c7bb6'))
# col_list <- custom_colour_pal(ncol(ages_2))
#
#
# graphics.off()
# plot(ages_median, type = "l")
# for (i in 1:ncol(ages_2)) {
#   lines(ages_2[, i], col = col_list[i])
# }
# lines(ages_median, type = "l", lwd = 2)
#
#
#
#
# #lines(1.97*ages_sds+ages_median,col="red")
# #lines(ages_median-1.97*ages_sds,col="blue")
#
# plot(ages_sds * 1.97,
#      type = "l",
#      log = "y",
#      ylim = c(2, max(ages_sds * 1.97)))
# abline(v = c(40, 330, 740))
#
#
#
#
#
# res <- matrix(data = NA,
#               ncol = nrow(ages),
#               nrow = 0)
#
#
# graphics.off()
# ages_means <- apply(res, 2, median)
# sd <- apply(res, 2, sd)
# sd * 2
# pred_ages
#
#
# for (i in 1:ncol(res)) {
#   hist(res[, i], main = pred_ages[i, 1], xlim = c((ages_means[i] - sd[i] *
#                                                      7.5),
#                                                   (ages_means[i] + sd[i] *
#                                                      7.5)), )
#   abline(v = pred_ages[i, 2],
#          col = "red",
#          lwd = 5)
#   abline(v = ages_means[i], col = "green")
#   abline(v = ages_means[i] + sd[i] * 2, col = "blue")
#   abline(v = ages_means[i] - sd[i] * 2, col = "blue")
#
#
# }
#
#
#
#
#
# library(modifiedBChron)
#
#
# pos_time_rand
# pos_thick_rand
# age_curves
#
#
# j <- 1
#
# n_simulations <- 22 * 10
#
# numCores <- detectCores()
# cl <- makeCluster(numCores - 2)
# registerDoSNOW(cl)
#
# pb <- txtProgressBar(max = n_simulations, style = 3)
# progress <- function(n)
#   setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
#
#
#
# ages <-
#   foreach::foreach (
#     j = 1:(n_simulations),
#     .options.snow = opts,
#     .errorhandling = "remove",
#     .packages = c("Bchron", "matrixStats"),
#     .combine = "cbind"
#   ) %dopar% {
#     ages_1 <- as.numeric(ash_2[, 2])
#     ageSds_1 <- as.numeric(ash_2[, 3])
#     positions_1 <- as.numeric(pos_time_rand[, j])
#     positionThicknesses_1 <- as.numeric(pos_thick_rand[, j])
#     ids_1 <- ash_2[, 1]
#     distTypes_1 <- rep("G", length(ages_1))
#     predictPositions_1 <-   age_curves[, j]
#     #rownrs <- pos_rownr_rand[, j]
#
#
#     bchron_modeled <- modifiedBChron::bchron_model(
#       ages = ages_1,
#       ageSds = ageSds_1,
#       positions = positions_1,
#       positionThicknesses = positionThicknesses_1,
#       ids = ids_1,
#       distTypes = distTypes_1,
#       iterations = 5000,
#       burn = 1000,
#       probability = 0.95,
#       predictPositions = predictPositions_1,
#       truncateUp = 0,
#       extrapUp = 1000,
#       truncateDown = 1e+9,
#       extrapDown = 1000
#     )
#
#
#     # model_res <- bchron_modeled$model
#     # model_res_dif <- model_res
#     # bchron_modeled_50 <- bchron_modeled$HDI$`50%`
#     #
#     # for (i in 1:ncol(model_res_dif)){
#     #   model_res_dif[,i] <- sqrt((model_res_dif[,i]-bchron_modeled_50)^2)}
#     # model_res_dif_sums <- colSums(model_res_dif)
#     # col_nr <-DescTools::Closest(min(model_res_dif_sums),model_res_dif_sums, which = TRUE)
#     # brchron_out <- model_res[,col_nr[1]]
#     brchron_out <- bchron_modeled$HDI$`50%`
#
#
#   }
#
# stopCluster(cl)
#
#
#
# str(ages)
#
# ages_2 <- ages
# ages_2 <- as.matrix(ages_2)
#
#
# hist(ages_2[c(pos_rownr_min[1, 1]:pos_rownr_plus[1, 1]), ])
# mean(ages_2[c(pos_rownr_min[1, 1]:pos_rownr_plus[1, 1]), ])
# sd(ages_2[c(pos_rownr_min[1, 1]:pos_rownr_plus[1, 1]), ])
#
# hist(ages_2[c(pos_rownr_min[2, 1]:pos_rownr_plus[2, 1]), ])
# mean(ages_2[c(pos_rownr_min[2, 1]:pos_rownr_plus[2, 1]), ])
# sd(ages_2[c(pos_rownr_min[2, 1]:pos_rownr_plus[2, 1]), ])
#
# hist(ages_2[c(pos_rownr_min[3, 1]:pos_rownr_plus[3, 1]), ])
# mean(ages_2[c(pos_rownr_min[3, 1]:pos_rownr_plus[3, 1]), ])
# sd(ages_2[c(pos_rownr_min[3, 1]:pos_rownr_plus[3, 1]), ])
#
#
#
# ages_sds <- matrixStats::rowSds(ages_2)
# ages_median <- rowMeans(ages_2)
#
# plot(ages_sds, log = "y", type = "l")
# abline(v = c(40, 330, 740))
# abline(v = c(25, 55, 315, 345, 725, 755))
#
#
# custom_colour_pal <-
#   colorRampPalette(c('#d7191c', '#fdae61', '#ffffbf', '#abd9e9', '#2c7bb6'))
# col_list <- custom_colour_pal(ncol(ages_2))
#
#
# graphics.off()
# plot(ages_median, type = "l")
# for (i in 1:ncol(ages_2)) {
#   lines(ages_2[, i], col = col_list[i])
# }
# lines(ages_median, type = "l", lwd = 2)
# abline(v = c(40, 330, 740))
#
#
# positionThicknesses_1 / (predictPositions_1[2] - predictPositions_1[1])
#
#
# str(bchron_modeled)
#
# modelPlot(bchron_modeled)
#
# plot(bchron_modeled$HDI$predictPositions,
#      bchron_modeled$HDI$`50%`,
#      type = "l")
# lines(bchron_modeled$HDI$predictPositions,
#       bchron_modeled$HDI$predictPositions)
# lines(bchron_modeled$HDI$predictPositions,
#       bchron_modeled$HDI$`25%`)
# lines(bchron_modeled$HDI$predictPositions,
#       bchron_modeled$HDI$`97.5`)
# lines(bchron_modeled$HDI$predictPositions,
#       bchron_modeled$HDI$`2.5`)
#
#
# hist(bchron_modeled$thetas$CCT18_322)
# hist(bchron_modeled$thetas$CCT18_315)
# hist(bchron_modeled$thetas$CCT18_311)
#
#
#
#
# ages_1 <- as.numeric(ash_2[, 2])
# ageSds_1 <- as.numeric(ash_2[, 3])
# positions_1 <- as.numeric(pos_time_rand[, j])
# positionThicknesses_1 <- as.numeric(pos_thick_rand[, j])
# ids_1 <- ash_2[, 1]
# distTypes_1 <- rep("G", length(ages_1))
# predictPositions_1 <-   age_curves[, j]
# rownrs <- pos_rownr_rand[, j]
#
# ages = ages_1
# ageSds = ageSds_1
# positions = positions_1
# positionThicknesses = positionThicknesses_1
# ids = ids_1
# distTypes = distTypes_1
# iterations = 5000
# burn = 1000
# probability = 0.95
# predictPositions = predictPositions_1
# truncateUp = 0
# extrapUp = 1000
# truncateDown = 1e+9
# extrapDown = 1000
#
#
#
#
#
#
# truncatedWalk = function(old, sd, low, high) {
#   if (isTRUE(all.equal(low, high, tolerance = 1e-10)))
#     return(list(new = low, rat = 1))
#   new = .C(
#     "truncatedWalk",
#     as.double(old),
#     as.double(sd),
#     as.double(low),
#     as.double(high),
#     as.double(0),
#     PACKAGE = "WaverideR"
#   )[5][[1]]
#   rat = .C(
#     "truncatedRat",
#     as.double(old),
#     as.double(sd),
#     as.double(low),
#     as.double(high),
#     as.double(new),
#     as.double(0),
#     PACKAGE = "WaverideR"
#   )[6][[1]]
#   return(list(new = new, rat = rat))
# }
# dtweediep1 = Vectorize(function(y, p, mu, phi) {
#   return(
#     .C(
#       "dtweediep1",
#       as.double(y),
#       as.double(p),
#       as.double(mu),
#       as.double(phi),
#       as.double(0),
#       PACKAGE = "WaverideR"
#     )[5][[1]]
#   )
# })
# predictInterp = function(alpha,
#                          lambda,
#                          beta,
#                          predictPositions,
#                          diffPositionj,
#                          currPositionsj,
#                          currPositionsjp1,
#                          thetaj,
#                          thetajp1) {
#   return(
#     .C(
#       "predictInterp",
#       as.double(alpha),
#       as.double(lambda),
#       as.double(beta),
#       as.double(predictPositions),
#       as.integer(length(predictPositions)),
#       as.double(diffPositionj),
#       as.double(currPositionsj),
#       as.double(currPositionsjp1),
#       as.double(thetaj),
#       as.double(thetajp1),
#       as.double(rep(0, length(predictPositions))),
#       PACKAGE = "WaverideR"
#     )[11][[1]]
#   )
# }
# predictExtrapUp = function(alpha,
#                            lambda,
#                            beta,
#                            predictPositions,
#                            currPositions1,
#                            theta1,
#                            maxExtrap,
#                            extractDate) {
#   return(
#     .C(
#       "predictExtrapUp",
#       as.double(alpha),
#       as.double(lambda),
#       as.double(beta),
#       as.double(predictPositions),
#       as.integer(length(predictPositions)),
#       as.double(currPositions1),
#       as.double(theta1),
#       as.integer(maxExtrap),
#       as.double(extractDate),
#       as.double(rep(0, length(predictPositions))),
#       PACKAGE = "WaverideR"
#     )[10][[1]]
#   )
# }
# predictExtrapDown = function(alpha,
#                              lambda,
#                              beta,
#                              predictPositions,
#                              currPositions1,
#                              theta1,
#                              maxExtrap,
#                              extractDate) {
#   return(
#     .C(
#       "predictExtrapDown",
#       as.double(alpha),
#       as.double(lambda),
#       as.double(beta),
#       as.double(predictPositions),
#       as.integer(length(predictPositions)),
#       as.double(currPositions1),
#       as.double(theta1),
#       as.integer(maxExtrap),
#       as.double(extractDate),
#       as.double(rep(0, length(predictPositions))),
#       PACKAGE = "WaverideR"
#     )[10][[1]]
#   )
# }
#
# compoundProb <- function(ages, sigs, distType, x) {
#   interval <- matrix(0, nrow = length(x), ncol = length(ages))
#   for (i in 1:length(ages)) {
#     if (distType[i] == "G") {
#       interval[, i] <- dnorm(x, ages[i], sigs[i]) / length(ages) *
#         mean(diff(x))
#     }
#     else if (distType[i] == "U") {
#       interval[, i] <- dunif(x, ages[i] - sigs[i],
#                              ages[i] + sigs[i]) / length(ages) * mean(diff(x))
#     }
#   }
#   G <- apply(interval, 1, sum)
#   return(G)
# }
#
# o = order(positions)
# if (any(positions[o] != positions)) {
#   #warning("positions not given in order - re-ordering")
#   ages <- ages[o]
#   ageSds <- ageSds[o]
#   positions <- positions[o]
#   distTypes <- distTypes[o]
#   positionThicknesses <- positionThicknesses[o]
#   ids = ids[o]
# }
# nSamples <- length(unique(ids))
# masterPositions <- vector()
# nNames <- unique(ids)
#
# for (i in 1:nSamples) {
#   masterPositions[i] <- positions[ids == nNames[i]][1]
# }
#
# nThicknesses <- vector(length = nSamples)
#
# for (i in 1:nSamples) {
#   nThicknesses[i] <- unique(positionThicknesses[positions ==
#                                                   masterPositions[i]])
# }
# nAges <- rep(NA, length(nSamples))
#
# for (i in 1:nSamples) {
#   nAges[i] <- length(ids[ids == nNames[i]])
# }
# prob <- matrix(0, nrow = 1e+05, ncol = nSamples)
# ageGrid <- seq(min(ages - ageSds * 10), max(ages + ageSds *
#                                               10), length.out = 1e+05)
# for (j in 1:nSamples) {
#   prob[, j] <- compoundProb(ages[ids == nNames[j]], ageSds[ids ==
#                                                              nNames[j]], distType = distTypes[ids == nNames[j]],
#                             x = ageGrid)
# }
# rm(j)
# thetaStore <- matrix(NA, nrow = iterations, ncol = nSamples)
# muStore <- vector(length = iterations)
# psiStore <- muStore
# modelStore <-
#   matrix(NA, ncol = iterations, nrow = length(predictPositions))
# positionStore <- matrix(NA, nrow = iterations, ncol = nSamples)
# predictStore <-
#   matrix(NA, ncol = iterations, nrow = length(predictPositions))
# psiSDStore <- psiStore
# muSDStore <- psiStore
# mhSDStore <- thetaStore
# thetas <- vector(length = nSamples)
# for (i in 1:nSamples) {
#   thetas[i] <- ageGrid[which.max(prob[, i])] + rnorm(1,
#                                                      0, 1e-05)
# }
# currPositions <- masterPositions
# mhSD <- runif(length(unique(ids)))
# psiSD <- runif(1)
# muSD <- runif(1)
# mu <-
#   abs(rnorm(1, mean = mean(diff(thetas)) / mean(diff(currPositions)),
#             sd = muSD))
# psi <- abs(rnorm(1, 1, sd = psiSD))
# p = 1.2
# alpha <- (2 - p) / (p - 1)
#
# #pb <-  progress::progress_bar$new(total = iterations, format = "[:bar] :percent eta: :eta")
#
# for (n in 1:iterations) {
#   #pb$tick()
#   lambda <- (mu ^ (2 - p)) / (psi * (2 - p))
#   beta <- 1 / (psi * (p - 1) * mu ^ (p - 1))
#   currPositions <-
#     runif(nSamples,
#           masterPositions - nThicknesses,
#           masterPositions + nThicknesses)
#   do <- order(currPositions)
#   diffPositions <- diff(currPositions[do])
#   thetas[do] <- sort(thetas, decreasing = T)
#   positionStore[n, ] <- currPositions
#   for (j in 1:(nSamples - 1)) {
#     inRange <- predictPositions >= currPositions[do[j]] &
#       predictPositions <= currPositions[do[j + 1]]
#     predval <- predictInterp(
#       alpha = alpha,
#       lambda = lambda,
#       beta = beta,
#       predictPositions = predictPositions[inRange],
#       diffPositionj = diffPositions[j],
#       currPositionsj = currPositions[do[j]],
#       currPositionsjp1 = currPositions[do[j + 1]],
#       thetaj = thetas[do[j]],
#       thetajp1 = thetas[do[j +
#                              1]]
#     )
#     modelStore[inRange, n] <- predval
#   }
#   if (any(predictPositions >= currPositions[do[nSamples]])) {
#     inRange <- predictPositions >= currPositions[do[nSamples]]
#     predval <- predictExtrapUp(
#       alpha = alpha,
#       lambda = lambda,
#       beta = beta,
#       predictPositions = predictPositions[inRange],
#       theta1 = thetas[do[nSamples]],
#       currPositions1 = currPositions[do[nSamples]],
#       maxExtrap = extrapUp,
#       extractDate = truncateUp
#     )
#     modelStore[inRange, n] <- predval
#   }
#   if (any(predictPositions <= currPositions[do[1]])) {
#     inRange <- predictPositions <= currPositions[do[1]]
#     predval <- predictExtrapDown(
#       alpha = alpha,
#       lambda = lambda,
#       beta = beta,
#       predictPositions = predictPositions[inRange],
#       theta1 = thetas[do[1]],
#       currPositions1 = currPositions[do[1]],
#       maxExtrap = extrapDown,
#       extractDate = truncateDown
#     )
#     modelStore[inRange, n] <- predval
#   }
#   for (i in 1:nSamples) {
#     thetaCurrent <- thetas[do[i]]
#     thetaProposed <- truncatedWalk(
#       old = thetaCurrent,
#       sd = mhSD[do[i]],
#       low = ifelse(i == nSamples,
#                    truncateUp, thetas[do[i + 1]] - 1e-10),
#       high = ifelse(i ==
#                       1, 1e+10, thetas[do[i - 1]] + 1e-10)
#     )
#     pProposed <- log(prob[, do[i]][which.min(abs(ageGrid -
#                                                    thetaProposed$new))])
#     pProposed <- max(pProposed, -1e+06, na.rm = T)
#     pCurrent <- log(prob[, do[i]][which.min(abs(ageGrid -
#                                                   thetaCurrent))])
#     pCurrent <- max(pCurrent, -1e+06, na.rm = T)
#     priorProposed <-
#       ifelse(i == 1, 0, log(
#         dtweediep1(
#           thetas[do[i -
#                       1]] - thetaProposed$new,
#           p = p,
#           mu = mu * (diffPositions[i -
#                                      1]),
#           phi = psi * (diffPositions[i - 1]) ^ (p -
#                                                   1)
#         )
#       )) + ifelse(i == nSamples, 0, log(
#         dtweediep1(
#           thetaProposed$new -
#             thetas[do[i + 1]],
#           p = p,
#           mu = mu * (diffPositions[i]),
#           phi = psi * (diffPositions[i]) ^
#             (p - 1)
#         )
#       ))
#     priorCurrent <-
#       ifelse(i == 1, 0, log(
#         dtweediep1(
#           thetas[do[i -
#                       1]] - thetaCurrent,
#           p = p,
#           mu = mu * (diffPositions[i -
#                                      1]),
#           phi = psi * (diffPositions[i - 1]) ^ (p -
#                                                   1)
#         )
#       )) + ifelse(i == nSamples, 0, log(
#         dtweediep1(
#           thetaCurrent -
#             thetas[do[i + 1]],
#           p = p,
#           mu = mu * (diffPositions[i]),
#           phi = psi * (diffPositions[i]) ^
#             (p - 1)
#         )
#       ))
#     priorCurrent <- max(priorCurrent, -1e+06, na.rm = T)
#     priorProposed <- max(priorProposed, -1e+06, na.rm = T)
#     logRtheta <- priorProposed - priorCurrent + pProposed -
#       pCurrent + log(thetaProposed$rat)
#     logRtheta <- max(logRtheta, -1e+06, na.rm = T)
#     if (runif(1, 0, 1) < min(1, exp(logRtheta))) {
#       thetas[do[i]] <- thetaProposed$new
#     }
#   }
#   thetaStore[n, ] <- thetas
#   muCurrent <- mu
#   muProposed <- truncatedWalk(
#     old = mu,
#     sd = muSD,
#     low = 1e-10,
#     high = 1e+10
#   )
#   priorMu <- sum(log(
#     dtweediep1(
#       abs(diff(thetas[do])),
#       p = p,
#       mu = muProposed$new * diffPositions,
#       phi = psi / diffPositions ^ (p -
#                                      1)
#     )
#   )) - sum(log(
#     dtweediep1(
#       abs(diff(thetas[do])),
#       p = p,
#       mu = mu * diffPositions,
#       phi = psi / diffPositions ^ (p -
#                                      1)
#     )
#   )) + log(muProposed$rat)
#   mu <- ifelse(runif(1, 0, 1) < min(1, exp(priorMu)),
#                muProposed$new, muCurrent)
#   muStore[n] <- mu
#   psiCurrent <- psi
#   psiProposed <- truncatedWalk(
#     old = psi,
#     sd = psiSD,
#     low = 1e-10,
#     high = 1e+10
#   )
#   priorPsi <- sum(log(
#     dtweediep1(
#       abs(diff(thetas[do])),
#       p = p,
#       mu = mu * diffPositions,
#       phi = psiProposed$new / diffPositions ^ (p -
#                                                  1)
#     )
#   )) - sum(log(dtweediep1(
#     abs(diff(thetas[do])),
#     p = p,
#     mu = mu * diffPositions,
#     phi = psi / (diffPositions ^ (p -
#                                     1))
#   ))) + log(psiProposed$rat)
#   psi <- ifelse(runif(1, 0, 1) < min(1, exp(priorPsi)),
#                 psiProposed$new, psiCurrent)
#   psiStore[n] <- psi
#   h = 200
#   if (n %% h == 0) {
#     cd <- 2.4 / sqrt(1)
#     psiK <- psiStore[(n - h + 1):n] - mean(psiStore[(n -
#                                                        h + 1):n])
#     psiRt <- var(psiK)
#     psiSD <- sqrt((cd ^ 2) * psiRt)
#     psiSD <- ifelse(is.na(psiSD), cd * sd(psiStore,
#                                           na.rm = T), psiSD)
#     psiSD <- ifelse(psiSD == 0, cd * sd(psiStore, na.rm = T),
#                     psiSD)
#     muK <- muStore[(n - h + 1):n] - mean(muStore[(n -
#                                                     h + 1):n])
#     muRt <- var(muK)
#     muSD <- sqrt(cd ^ 2 * muRt)
#     muSD <- ifelse(is.na(muSD), cd * sd(muStore, na.rm = T),
#                    muSD)
#     muSD <- ifelse(muSD == 0, cd * sd(muStore, na.rm = T),
#                    muSD)
#     for (q in 1:nSamples) {
#       thetaK <- thetaStore[(n - h + 1):n, q] - mean(thetaStore[(n -
#                                                                   h + 1):n, q])
#       thetaRt <- var(thetaK)
#       mhSD[q] <- ifelse(is.na(sqrt(cd ^ 2 * thetaRt)),
#                         cd * sd(thetaStore[, q], na.rm = T),
#                         sqrt(cd ^ 2 *
#                                thetaRt))
#     }
#   }
#   psiSDStore[n] <- psiSD
#   muSDStore[n] <- muSD
#   mhSDStore[n, ] <- mhSD
# }
#
# # suppressWarnings(isOutlier <- function(probability) {
# #   outlier <- vector(length = ncol(thetaStore))
# #   for (i in 1:ncol(thetaStore)) {
# #     f <- approxfun(x = cumsum(prob[, i]), y = ageGrid)
# #     likelihood <- f(c((1 - probability) / 2, 0.5, (1 +
# #                                                      probability) / 2))
# #     posterior <- as.numeric(quantile(thetaStore[burn:iterations,
# #                                                 i], prob = c((1 - probability) /
# #                                                                2, 0.5, (1 +
# #                                                                           probability) /
# #                                                                2
# #                                                 )))
# #     outlier[i] <- !(any(posterior > likelihood[1]) &
# #                       any(posterior < likelihood[3]))
# #   }
# #   return(outlier)
# # })
# #
# # suppressWarnings(outliers <- isOutlier(probability = probability))
#
# #
# # if (any(outliers)) {
# #   for (i in 1:length(ids[outliers])) {
# #     warning(
# #       paste(
# #         "sample",
# #         ids[outliers][i],
# #         "posterior distribution",
# #         probability * 100,
# #         "% HDI",
# #         "does not overlap the input likelihood. Consider screening for outliers."
# #       )
# #     )
# #   }
# # }
#
#
# HDI <- as.vector(cbind(t(apply(modelStore[, burn:iterations], 1, quantile, c(0.5)))))
#
#
# HDI <- as.data.frame(cbind(t(apply(
#   modelStore[, burn:iterations],
#   1, quantile, c((1 - probability) /
#                    2, 0.5, (1 + probability) / 2)
# )),
# predictPositions))
# thetaStore <- setNames(as.data.frame(thetaStore), nNames)
# inputs <-
#   data.frame(ids, ages, ageSds, positions, positionThicknesses,
#              distTypes)
# prob <- setNames(as.data.frame(prob), nNames)
# output <- list(
#   HDI = HDI,
#   inputs = inputs,
#   model = modelStore,
#   thetas = thetaStore,
#   positionStore = positionStore,
#   psi = psiStore,
#   mu = muStore,
#   ageGrid = ageGrid,
#   likelihoods = prob,
#   nAges = nAges,
#   masterPositions = masterPositions,
#   ids = nNames,
#   psiSDStore = psiSDStore,
#   muSDStore = muSDStore,
#   mhSDStore = mhSDStore,
#   burn = burn,
#   iterations = iterations,
#   outliers = outliers,
#   probability = probability,
#   positionThicknesses = positionThicknesses
# )
# class(output) <- "bchron_model_run"
#
# output
#
# return(output)
#
# #}
#
#
#
#
#
#
#
#
#
