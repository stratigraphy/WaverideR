#' @title Fit linear models to spectral peaks extracted from the wavelet spectra to astronomical cycles multiplied by sedimentation rate x
#'
#' @description The \code{\link{flmw}} function is used calculate the linear correlation
#' for a list of astronomical cycles transformed using a range of sedimentation rates and then compared
#' to spectral peaks of a wavelet spectra
#'
#'@param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function
#'@param sedrate_low Minimum sedimentation rate (cm/kyr)for which the sum of maximum spectral power is calculated for.
#'@param sedrate_high Maximum sedimentation rate (cm/kyr) for which the sum of maximum spectral power is calculated  for.
#'@param spacing Spacing (cm/kyr) between sedimentation rates
#'@param cycles Astronomical cycles (in kyr) for which the combined sum of maximum spectral power is calculated for
#'@param x_lab label for the y-axis \code{Default="depth"}
#'@param y_lab label for the y-axis \code{Default="sedrate"}
#'@param run_random run multiple simulation to calculate percentile against the 0 hypothesis
#'@param rand_simulations nr of simulations to calculate percentile against the 0 hypothesis
#'@param run_multicore run simulation using multiple cores \code{Default=FALSE}
#'the simulation is run at x-2 cores to allow the 2 remaining processes to run background processes
#'@param genplot Generate plot \code{Default="FALSE"}
#'@param plot_res options 1-8 option 1: slope coefficient, option 2: r squared,
#'option 3: nr of components, option 4: difference to the  origin , option 5: slope coefficient percentile
#'option 6: r squared percentile, option 7: nr of components percentile,
#'option 8: difference to the origin percentile \code{Default=2}
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#' @param verbose Print text \code{Default=FALSE}.
#'
#' @author
#'Based on the \link[astrochron]{eAsm} function of the 'astrochron' R package and the 'eCOCO' and 'COCO' function of the 'Acycle' software
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis \doi{<doi:10.1016/j.earscirev.2018.11.015>}
#'
#'Acycle: Time-series analysis software for paleoclimate research and education,
#'Mingsong Li, Linda Hinnov, Lee Kump,
#'Computers & Geosciences,Volume 127,2019,Pages 12-22,ISSN 0098-3004,
#'\doi{<doi:10.1016/j.cageo.2019.02.011>}
#'
#'Tracking variable sedimentation rates and astronomical forcing in Phanerozoic paleoclimate proxy series with evolutionary correlation coefficients and hypothesis testing,
#'Mingsong Li, Lee R. Kump, Linda A. Hinnov, Michael E. Mann,
#'Earth and Planetary Science Letters,Volume 501,
#'T2018,Pages 165-179,ISSN 0012-821X,\doi{<doi:10.1016/j.epsl.2018.08.041>}
#'
#'@examples
#'\donttest{
#'#estimate sedimentation rate for the magnetic susceptibility record
#'# of the Sullivan core of Pas et al., (2018).
#'
#'mag_wt <- analyze_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'sedrates <- flmw(wavelet = mag_wt,
#'     sedrate_low = 0.5,
#'     sedrate_high = 4,
#'     spacing = 0.05,
#'     cycles = c(2376,1600,1180,696,406,110),
#'     x_lab = "depth",
#'     y_lab = "sedrate",
#'     run_random = FALSE,
#'     rand_simulations = 1000,
#'     run_multicore = FALSE,
#'     genplot = FALSE,
#'     plot_res = 2,
#'     keep_editable=FALSE,
#'     verbose=FALSE)
#'}
#'
#' @return
#'Returns a list which contains 10 elements
#'element 1: slope coefficient
#'element 2: r squared
#'element 3: nr of components
#'element 4: difference to the origin
#'element 5: slope coefficient percentile
#'element 6: r squared percentile
#'element 7: nr of components percentile,
#'element 8: difference to the origin percentile
#'element 9: y-axis values of the matrices which is sedimentation rate
#'element 10: x-axis values of the matrices which is depth
#'
#' @export
#' @importFrom matrixStats rowSds
#' @importFrom Matrix rowMeans
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom tcltk setTkProgressBar
#' @importFrom tcltk setTkProgressBar
#' @importFrom foreach foreach
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats lm
#' @importFrom stats ecdf
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics text
#' @importFrom graphics box
#' @importFrom graphics polygon
#' @importFrom grDevices rgb
#' @importFrom foreach %dopar%
#' @importFrom graphics layout
#' @importFrom parallel stopCluster
#' @importFrom truncnorm rtruncnorm
#' @importFrom astrochron asm
#' @importFrom astrochron eAsm

flmw  <- function(wavelet = NULL,
                  sedrate_low = NULL,
                  sedrate_high = NULL,
                  spacing = NULL,
                  cycles = c(NULL),
                  x_lab = "depth",
                  y_lab = "sedrate",
                  run_random = FALSE,
                  rand_simulations = 1000,
                  run_multicore = FALSE,
                  genplot = FALSE,
                  plot_res = 2,
                  keep_editable = FALSE,
                  verbose=FALSE) {
  my.data <- cbind(wavelet$x, wavelet$y)
  my.w <- wavelet

  Power <- matrix(unlist(as.data.frame(my.w$Power)))
  Power <-
    as.data.frame(matrix(
      unlist(as.data.frame(my.w$Power)),
      ncol = my.w$nc,
      nrow = my.w$nr
    ))
  Power <- t(Power)
  Powert <- t(Power)

  testsedrates <-
    seq(from = sedrate_low, to = sedrate_high, by = spacing)
  testsedrates <- as.data.frame(testsedrates)

  if (run_multicore == TRUE) {
    numCores <- detectCores()
    cl <- parallel::makeCluster(numCores - 2)
    registerDoSNOW(cl)
  } else{
    numCores <- 1
    cl <- makeCluster(numCores)
    registerDoSNOW(cl)
  }


  simulations <- ncol(Powert)

  if (verbose==TRUE){
    pb <- txtProgressBar(max = simulations, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)}else{opts=NULL}


  ijk <- 1 # needed to assign 1 to ijk to avoid note

  if (run_random == TRUE) {
    #generate random simulation shifts
    randomize_mat <-
      matrix(
        data = truncnorm::rtruncnorm(
          n = (nrow(Powert) / 2) * rand_simulations,
          a = 1 / 10,
          b = 10,
          mean = 1,
          sd = 1
        ),
        ncol = rand_simulations
      )
  }


  fit <-
    foreach (ijk = 1:simulations, .options.snow = opts) %dopar% {
      fits <- matrix(data = NA,
                     nrow = nrow(testsedrates),
                     ncol = 8)
      slope_coeff_random_mat <-
        matrix(data = NA,
               nrow = nrow(testsedrates),
               ncol = rand_simulations)
      corr_coeff_random_mat  <-
        matrix(data = NA,
               nrow = nrow(testsedrates),
               ncol = rand_simulations)
      nr_compenents_random_mat <-
        matrix(data = NA,
               nrow = nrow(testsedrates),
               ncol = rand_simulations)
      slope_orig_random_mat <-
        matrix(data = NA,
               nrow = nrow(testsedrates),
               ncol = rand_simulations)
      nearest_cycle <-
        as.data.frame(matrix(
          data = NA,
          nrow = length(cycles),
          ncol = 3
        ))

      data <- as.data.frame(cbind(my.w$Period, Powert[, ijk]))

      astro_mindetect <- as.data.frame(data)
      astro_mindetect$min <- 0
      for (i in 3:(nrow(data) - 2)) {
        if ((data[i, 2] - data[(i + 2), 2] < 0) &
            (data[i, 2] - data[(i - 2), 2] < 0))
        {
          astro_mindetect[i, 3] <- 1
        }
      }

      astro_mindetect_error_corr <- astro_mindetect
      astro_mindetect_error_corr <-
        astro_mindetect_error_corr[astro_mindetect_error_corr$min == 1 ,]

      astro_maxdetect <- as.data.frame(data)
      astro_maxdetect$max <- 0
      for (i in 3:(nrow(data) - 2)) {
        if ((data[i, 2] - data[(i + 2), 2] > 0) &
            (data[i, 2] - data[(i - 2), 2]  > 0))
        {
          astro_maxdetect[i, 3] <- 1
        }
      }

      astro_maxdetect_error_corr <- astro_maxdetect
      astro_maxdetect_error_corr <-
        astro_maxdetect_error_corr[astro_maxdetect_error_corr$max == 1 ,]

      max <- astro_maxdetect_error_corr
      colnames(max) <- c("A", "B", "C")
      min <- astro_mindetect_error_corr
      colnames(min) <- c("A", "B", "C")

      min[, 3] <- -1
      peaks <- rbind(max, min)

      peaks <- peaks[order(peaks[, 1]), ]
      i <- 1
      res_rownr <- nrow(peaks)

      while (i < res_rownr) {
        if (peaks[i, 3] == peaks[(i + 1), 3]) {
          if ((peaks[i, 3]  == 1 & peaks[(i + 1), 3] == 1) &
              (peaks[i, 2] > peaks[(i + 1), 2])) {
            peaks[(i + 1), ] <- NA
            peaks <- na.omit(peaks)
            res_rownr <- res_rownr - 1
          }
          if ((peaks[i, 3]  == 1 & peaks[(i + 1), 3] == 1) &
              (peaks[i, 2] < peaks[(i + 1), 2])) {
            peaks[i, ] <- NA
            peaks <- na.omit(peaks)
            res_rownr <- res_rownr - 1
          }
          if ((peaks[i, 3] == -1 & peaks[(i + 1), 3] == -1) &
              (peaks[i, 2] < peaks[(i + 1), 2])) {
            peaks[(i + 1), ] <- NA
            peaks <- na.omit(peaks)
            res_rownr <- res_rownr - 1
          }
          if ((peaks[i, 3] == -1 & peaks[(i + 1), 3] == -1) &
              (peaks[i, 2] > peaks[(i + 1), 2])) {
            peaks[i, ] <- NA
            peaks <- na.omit(peaks)
            res_rownr <- res_rownr - 1

          }
        }
        if ((peaks[i, 3] != peaks[(i + 1), 3]) |
            is.na(peaks[i, 3] != peaks[(i + 1), 3])) {
          i <- i + 1
        }
      }

      peaks <- peaks[peaks[, 3] > 0, c(1, 2)]


      for (ij in 1:nrow(testsedrates)) {
        cycles_in_depth <- ((cycles * testsedrates[ij, ]) / 100)
        nearest_cycle[, 1] <- cycles_in_depth

        for (jj in 1:length(cycles)) {
          row_nr <-
            DescTools::Closest(peaks[, 1], cycles_in_depth[jj], which = TRUE)
          nearest_cycle[jj, 2] <- peaks[row_nr[1], 1]
        }

        nearest_cycle[, 3] <-
          sqrt((nearest_cycle[, 2] - nearest_cycle[, 1]) ^ 2)
        nearest_cycle <-
          nearest_cycle[order(-nearest_cycle[, 2], -nearest_cycle[, 3]),]
        sel_cycles <-
          nearest_cycle[!duplicated(nearest_cycle[, 2], fromLast = TRUE),]

        if (nrow(sel_cycles) <= 1) {
          fits[ij, 1] <- 0
          fits[ij, 2] <- 0
          fits[ij, 3] <- nrow(sel_cycles)
          fits[ij, 4] <- 0
        } else{
          lmodel <- summary(lm(sel_cycles[, 1] ~ sel_cycles[, 2]))
          fits[ij, 1] <- sqrt((1 - lmodel$coefficients[2, 1]) ^ 2)
          fits[ij, 2] <- lmodel$r.squared
          fits[ij, 3] <- nrow(sel_cycles)
          fits[ij, 4] <- 1 - sqrt((lmodel$coefficients[1, 1] ^ 2))

        }

      }

      if (run_random == TRUE) {
        for (rt in 1:rand_simulations) {
          random_peaks <-
            sort(peaks[, 1] * (randomize_mat[1:nrow(peaks), rt]))

          for (ij in 1:nrow(testsedrates)) {
            cycles_in_depth <- ((cycles * testsedrates[ij, ]) / 100)
            nearest_cycle[, 1] <- cycles_in_depth

            for (pl in 1:length(cycles)) {
              row_nr <-
                DescTools::Closest(random_peaks[], cycles_in_depth[pl], which = TRUE)
              nearest_cycle[pl, 2] <- random_peaks[row_nr[1]]
            }

            nearest_cycle[, 3] <-
              sqrt((nearest_cycle[, 2] - nearest_cycle[, 1]) ^ 2)
            nearest_cycle <-
              nearest_cycle[order(-nearest_cycle[, 2],-nearest_cycle[, 3]),]
            sel_cycles <-
              nearest_cycle[!duplicated(nearest_cycle[, 2], fromLast = TRUE),]

            if (nrow(sel_cycles) <= 1) {
              slope_coeff_random_mat[ij, rt] <- 0
              corr_coeff_random_mat[ij, rt] <- 0
              nr_compenents_random_mat[ij, rt] <- nrow(sel_cycles)
              slope_orig_random_mat[ij, rt] <- 0
            } else{
              lmodel_random <- summary(lm(sel_cycles[, 1] ~ sel_cycles[, 2]))
              slope_coeff_random_mat[ij, rt]  <-
                sqrt((1 - lmodel_random$coefficients[2, 1]) ^ 2)
              corr_coeff_random_mat[ij, rt] <-
                lmodel_random$r.squared
              nr_compenents_random_mat[ij, rt] <- nrow(sel_cycles)
              slope_orig_random_mat[ij, rt] <-
                1 - sqrt((lmodel_random$coefficients[1, 1] ^ 2))
            }

          }
        }


        for (yu in 1:nrow(testsedrates)) {
          slope_coeff_random_per <- ecdf(slope_coeff_random_mat[yu, ])
          corr_coeff_random_per <- ecdf(corr_coeff_random_mat[yu, ])
          nr_compenents_random_per <-
            ecdf(nr_compenents_random_mat[yu, ])
          slope_orig_random_per <- ecdf(slope_orig_random_mat[yu, ])

          fits[yu, 5] <- 1 - slope_coeff_random_per(fits[yu, 1])
          fits[yu, 6] <- corr_coeff_random_per(fits[yu, 2])
          fits[yu, 7] <- nr_compenents_random_per(fits[yu, 3])
          fits[yu, 8] <- slope_orig_random_per(fits[yu, 4])

        }
      }


      fits <- fits

    }


  stopCluster(cl)

  fit2 <- fit


  slope_coeff_mat <-
    matrix(
      data = NA,
      ncol = length(my.w$x),
      nrow = nrow(testsedrates)
    )
  corr_coeff_mat <-
    matrix(
      data = NA,
      ncol = length(my.w$x),
      nrow = nrow(testsedrates)
    )
  nr_compenents_mat <-
    matrix(
      data = NA,
      ncol = length(my.w$x),
      nrow = nrow(testsedrates)
    )
  orig_mat <-
    matrix(
      data = NA,
      ncol = length(my.w$x),
      nrow = nrow(testsedrates)
    )

  slope_coeff_mat_per <-
    matrix(
      data = NA,
      ncol = length(my.w$x),
      nrow = nrow(testsedrates)
    )
  corr_coeff_mat_per <-
    matrix(
      data = NA,
      ncol = length(my.w$x),
      nrow = nrow(testsedrates)
    )
  nr_compenents_mat_per <-
    matrix(
      data = NA,
      ncol = length(my.w$x),
      nrow = nrow(testsedrates)
    )
  orig_mat_per <-
    matrix(
      data = NA,
      ncol = length(my.w$x),
      nrow = nrow(testsedrates)
    )


  for (kk in 1:length(fit2)) {
    extract <- as.data.frame(fit2[[kk]])
    slope_coeff_mat[, kk] <- extract[, 1]
    corr_coeff_mat[, kk] <- extract[, 2]
    nr_compenents_mat[, kk] <- extract[, 3]
    orig_mat[, kk] <- extract[, 4]
    slope_coeff_mat_per[, kk]  <- extract[, 5]
    corr_coeff_mat_per[, kk]  <- extract[, 6]
    nr_compenents_mat_per[, kk]  <- extract[, 7]
    orig_mat_per[, kk]  <- extract[, 8]

  }


  depth <- (my.data[, 1])
  y_axis <- unlist(testsedrates)

  results <- list(
    slope_coeff_mat = slope_coeff_mat,
    corr_coeff_mat = corr_coeff_mat,
    nr_compenents_mat = nr_compenents_mat,
    orig_mat = orig_mat,
    slope_coeff_mat_per = slope_coeff_mat_per,
    corr_coeff_mat_per = corr_coeff_mat_per,
    nr_compenents_mat_per = nr_compenents_mat_per,
    orig_mat_per = orig_mat_per,
    depth = my.data[, 1],
    y_axis = as.numeric(unlist(testsedrates))
  )


  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    dev.new(width = 14, height = 7)
    layout.matrix <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
    layout(
      mat = layout.matrix,
      heights = c(1, 1, 1),
      # Heights of the two rows
      widths = c(8, 2, 2)
    )

    par(mar = c(4, 4, 2, 4))
    pmax_avg <- t(as.matrix(results[[plot_res]]))

    n.levels = 100


    power_max_mat.levels = quantile(pmax_avg, probs = seq(
      from = 0,
      to = 1,
      length.out = n.levels + 1
    ))


    image.plt = par()$plt
    color.palette = "rainbow(n.levels, start = 0, end = 0.7)"
    key.cols = rev(eval(parse(text = color.palette)))


    depth <- (my.data[, 1])
    y_axis <- unlist(testsedrates)
    depth <- as.numeric(depth)
    y_axis <- as.numeric(y_axis)
    image(
      x = depth,
      y = y_axis,
      z = (pmax_avg),
      col = key.cols,
      breaks = power_max_mat.levels,
      xlab = x_lab,
      ylab = y_lab,
      useRaster = TRUE
    )

    r_sum <- colMeans(pmax_avg)
    plot(
      y = y_axis,
      x = r_sum,
      type = "l",
      ylim = c(min(y_axis), max(y_axis)),
      yaxs = "i",
      xlab = "normalized power",
      ylab = y_lab
    )

    lwd.axis = 1
    n.ticks = 6
    label.digits = 6
    label.format = "f"
    width = 1.2
    lab.line = 2.5
    lab = NULL

    key.marks = round(seq(
      from = 0,
      to = 1,
      length.out = n.ticks
    ) *
      n.levels)
    key.labels = formatC(as.numeric(power_max_mat.levels),
                         digits = label.digits,
                         format = label.format)[key.marks +
                                                  1]


    image(
      1,
      seq(from = 0, to = n.levels),
      matrix(power_max_mat.levels,
             nrow = 1),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = ""
    )


    axis(
      4,
      lwd = lwd.axis,
      at = key.marks,
      labels = NA,
      tck = 0.02,
      tcl = (par()$usr[2] - par()$usr[1]) *
        width - 0.04
    )

    mtext(
      key.labels,
      side = 4,
      at = key.marks,
      line = 0.5,
      las = 2,
      font = par()$font.axis,
      cex = par()$cex.axis
    )

    box(lwd = lwd.axis)
  }
return(results)

}
