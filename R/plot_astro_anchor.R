#'@title Plot proxy record anchored to an astronomical solution
#'
#' @description Plot the results of the anchoring  the extracted signal to an astronomical solution using
#' which was conducted using the \code{\link{astro_anchor}}
#'
#'@param astro_solution Input is an astronomical solution with with the the proxy record was be anchored to,
#' the input should be a matrix or data frame with the first column being
#' age and the second column should be a insolation/angle/value
#'@param proxy_signal Input is the proxy data set which will which was
#' anchored to an astronomical solution, the input should be a matrix or
#' data frame with the first column being  depth/time and the second column should be a proxy value.
#'@param anchor_points Anchor points generated using the \code{\link{astro_anchor}} function
#'@param time_dir The direction of the proxy record which was assumed during anchoring if time increases with increasing depth/time values
#'(e.g. bore hole data which gets older with increasing depth ) then time_dir should be set to TRUE
#'if time decreases with depth/time values (eg stratospheric logs where 0m is the bottom of the section)
#'then time_dir should be set to FALSE \code{time_dir=TRUE}
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#'@return
#'The output is a set of 2 plots connected by lines
#'The top plot is the proxy record with anchor points on top of it
#'The bottom plot is the astronomical solution
#'The lines connect the anchor points
#'
#'@examples
#'\donttest{
#'# Use the grey_track example tracking points to anchor the grey scale data set
#'# of Zeeden et al., (2013) to the p-0.5t la2004 solution
#'
#'grey_wt <-
#'  analyze_wavelet(
#'    data = grey,
#'    dj = 1/200,
#'    lowerPeriod = 0.02,
#'    upperPeriod = 256,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'#Use the pretracked grey_track curve which traced the precession cycle
#'grey_track <- completed_series(
#'  wavelet = grey_wt,
#'  tracked_curve  = grey_track,
#'  period_up  = 1.25,
#'  period_down  = 0.75,
#'  extrapolate = TRUE,
#' genplot = FALSE
#')
#
#
#'# Extract precession, obliquity and eccentricity to create a synthetic insolation curve
#'
#'grey_prec <- extract_signal(
#'tracked_cycle_curve = grey_track[,c(1,2)],
#'wavelet = grey_wt,
#'period_up = 1.2,
#'period_down = 0.8,
#'add_mean = FALSE,
#'tracked_cycle_period = 22,
#'extract_cycle = 22,
#'tune = FALSE,
#'plot_residual = FALSE
#')
#'
#'grey_obl <- extract_signal(
#'  tracked_cycle_curve = grey_track[,c(1,2)],
#'  wavelet = grey_wt,
#'  period_up = 1.2,
#'  period_down = 0.8,
#'  add_mean = FALSE,
#'  tracked_cycle_period = 22,
#'  extract_cycle = 110,
#'  tune = FALSE,
#'  plot_residual = FALSE
#')
#'
#'grey_ecc <- extract_signal(
#'  tracked_cycle_curve = grey_track[,c(1,2)],
#'  wavelet = grey_wt,
#'  period_up = 1.25,
#'  period_down = 0.75,
#'  add_mean = FALSE,
#'  tracked_cycle_period = 22,
#'  extract_cycle = 40.8,
#'  tune = FALSE,
#'  plot_residual = FALSE
#')
#'
#'insolation_extract <- cbind(grey_ecc[,1],grey_prec[,2]+grey_obl[,2]+grey_ecc[,2]+mean(grey[,2]))
#'insolation_extract <- as.data.frame(insolation_extract)
#'insolation_extract_mins <- min_detect(insolation_extract,pts=3)
#'
#'#use the astrosignal_example to tune to which is an \cr
#'# ETP solution (p-0.5t la2004 solution).
#'
#'astrosignal_example <- na.omit(astrosignal_example)
#'astrosignal_example[,2] <- -1*astrosignal_example[,2]
#'astrosignal <- as.data.frame(astrosignal_example)
#'
#'#anchor the synthetic insolation curve extracted from the
#'# grey scale record to the insolation curve.
#'#use the anchor_points_grey data set to plot the
#'#result of using the astro_anchor function
#'
#'#anchor_points_grey <- astro_anchor(
#'#astro_solution = astrosignal,
#'#proxy_signal = insolation_extract,
#'#proxy_min_or_max = "min",
#'#clip_astrosolution = FALSE,
#'#astrosolution_min_or_max = "min",
#'#clip_high = NULL,
#'#clip_low = NULL,
#'#extract_astrosolution  = FALSE,
#'#astro_period_up  = NULL,
#'#astro_period_down  = NULL,
#'#astro_period_cycle  = NULL,
#'#extract_proxy_signal  = FALSE,
#'#proxy_period_up  = NULL,
#'#proxy_period_down  = NULL,
#'#proxy_period_cycle  = NULL,
#'#pts=3,
#'#verbose=FALSE,
#'#genplot=FALSE # set verbose to TRUE to allow for anchoring using text feedback commands
#'#)
#'
#'
#'plot_astro_anchor(astro_solution = astrosignal,
#'proxy_signal = insolation_extract,
#'anchor_points = anchor_points_grey,
#'time_dir = FALSE,
#'keep_editable = FALSE)
#'
#'}
#'
#' @export
#' @importFrom graphics points
#' @importFrom grDevices dev.new
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom grDevices graphics.off
#' @importFrom graphics segments
#' @importFrom graphics plot.new
#' @importFrom grDevices dev.size
#' @importFrom graphics plot.window


plot_astro_anchor <- function(astro_solution = NULL,
                              proxy_signal = NULL,
                              anchor_points = NULL,
                              time_dir = TRUE,
                              keep_editable = FALSE) {
  anchor_pts <- anchor_points
  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  if (time_dir == TRUE) {
    dev.new(width = 20,
            height = 10,
            noRStudioGD = TRUE)
    par(
      mfrow = c(2, 1),
      mai = c(0.5 / 2.54, 2 / 2.54, 2 / 2.54, 1 / 2.54),
      mgp = c(2, 1, 0)
    )
    plot(
      x = proxy_signal[, 1],
      y = proxy_signal[, 2],
      type = "l",
      xlab = "",
      ylab = "Proxy value",
      xaxt = "n",
      xaxs = "i",
      lwd = 2,
      ylim = c(min(proxy_signal[, 2]), max(proxy_signal[, 2])),
      xlim = (c(
        min(proxy_signal[, 1]), max(proxy_signal[, 1])
      ))
    )
    mtext(text = "Depth (metres)",
          side = 3,
          #side 2 = left
          line = 2)
    box(lwd = 2)
    axis(3)
    segments(
      x0 = anchor_pts[, 1],
      y0 = rep(min(proxy_signal[, 2]) * 0.2, times = (nrow(anchor_pts))),
      x1 = anchor_pts[, 1],
      y1 = anchor_pts[, 3],
      ylim = c(min(proxy_signal[, 2]), max(proxy_signal[, 2])),
      xlim = (c(
        min(proxy_signal[, 1]), max(proxy_signal[, 1])
      )),
      col = "black",
      lwd = 2
    )

    points(
      x = anchor_pts[, 1],
      y = anchor_pts[, 3],
      col = "green",
      pch = 19,
      ylim = c(min(proxy_signal[, 2]), max(proxy_signal[, 2])),
      xlim = (c(
        min(proxy_signal[, 1]), max(proxy_signal[, 1])
      )),
      cex = 2
    )




    par(
      new = FALSE,
      mai = c(2 / 2.54, 2 / 2.54, 0.5 / 2.54 , 1 / 2.54),
      mgp = c(2, 1, 0)
    )


    plot(
      astro_solution,
      type = "l",
      xlab = "Time",
      ylab = "tie point value",
      xaxs = "i",
      yaxs = "i",
      lwd = 2
    )
    box(lwd = 2)

    segments(
      x0 = anchor_pts[, 2],
      y0 = anchor_pts[, 4],
      x1 = anchor_pts[, 2],
      y1 = rep(max(astro_solution[, 2]) * 1.2, times = (nrow(anchor_pts))),
      col = "black",
      lwd = 2
    )
    points(
      x = anchor_pts[, 2],
      y = anchor_pts[, 4],
      col = "red",
      pch = 19,
      cex = 2
    )



    par(
      new = TRUE,
      mfrow = c(1, 1),
      mai = c(2 / 2.54, 2 / 2.54, 2 / 2.54, 1 / 2.54)
    )
    plot.new()
    settings_dev <- dev.size("in")

    plot.window(
      xlim = c(min(astro_solution[, 1]) , max(astro_solution[, 1])),
      ylim = c(2 / 2.54, settings_dev[2] - (3 / 2.54)),
      yaxs = "i",
      xaxs = "i"
    )

    y0_val <- ((settings_dev[2] - ((5) / 2.54)) / 2) + 1.5 / 2.54
    y1_val <-  ((settings_dev[2] - ((5) / 2.54)) / 2) + (2.5) / 2.54


    fact <- (anchor_pts[, 1] - min(proxy_signal[, 1])) /
      (max(proxy_signal[, 1]) - min(proxy_signal[, 1]))
    x1_vals <-
      ((max(astro_solution[, 1]) - min(astro_solution[, 1])) * (fact)) + min(astro_solution[, 1])


    segments(
      x0 = anchor_pts[, 2],
      y0 = rep(y0_val, times = (nrow(anchor_pts)))
      ,
      x1 = x1_vals,
      y1 = rep(y1_val, times = (nrow(anchor_pts))),
      col = "black",
      lwd = 2,
      lend = 1
    )

  } else{
    dev.new(width = 20,
            height = 10,
            noRStudioGD = TRUE)
    par(
      mfrow = c(2, 1),
      mai = c(0.5 / 2.54, 2 / 2.54, 2 / 2.54, 1 / 2.54),
      mgp = c(2, 1, 0)
    )
    plot(
      x = proxy_signal[, 1],
      y = proxy_signal[, 2],
      type = "l",
      xlab = "",
      ylab = "Proxy value",
      xaxt = "n",
      xaxs = "i",
      lwd = 2,
      ylim = c(min(proxy_signal[, 2]), max(proxy_signal[, 2])),
      xlim = rev(c(
        min(proxy_signal[, 1]), max(proxy_signal[, 1])
      ))
    )
    mtext(text = "Depth (metres)",
          side = 3,
          #side 2 = left
          line = 2)
    box(lwd = 2)
    axis(3)


    segments(
      x0 = anchor_pts[, 1],
      y0 = rep(min(proxy_signal[, 2]) * 0.2, times = (nrow(anchor_pts))),
      x1 = anchor_pts[, 1],
      y1 = anchor_pts[, 3],
      ylim = c(min(proxy_signal[, 2]), max(proxy_signal[, 2])),
      xlim = rev(c(
        min(proxy_signal[, 1]), max(proxy_signal[, 1])
      )),
      col = "black",
      lwd = 2
    )

    points(
      x = anchor_pts[, 1],
      y = anchor_pts[, 3],
      col = "green",
      pch = 19,
      ylim = c(min(proxy_signal[, 2]), max(proxy_signal[, 2])),
      xlim = rev(c(
        min(proxy_signal[, 1]), max(proxy_signal[, 1])
      )),
      cex = 2
    )




    par(
      new = FALSE,
      mai = c(2 / 2.54, 2 / 2.54, 0.5 / 2.54 , 1 / 2.54),
      mgp = c(2, 1, 0)
    )


    plot(
      astro_solution,
      type = "l",
      xlab = "Time",
      ylab = "tie point value",
      xaxs = "i",
      yaxs = "i",
      lwd = 2
    )
    box(lwd = 2)

    segments(
      x0 = anchor_pts[, 2],
      y0 = anchor_pts[, 4],
      x1 = anchor_pts[, 2],
      y1 = rep(max(astro_solution[, 2]) * 1.2, times = (nrow(anchor_pts))),
      col = "black",
      lwd = 2
    )
    points(
      x = anchor_pts[, 2],
      y = anchor_pts[, 4],
      col = "red",
      pch = 19,
      cex = 2
    )

    par(
      new = TRUE,
      mfrow = c(1, 1),
      mai = c(2 / 2.54, 2 / 2.54, 2 / 2.54, 1 / 2.54)
    )
    plot.new()
    settings_dev <- dev.size("in")

    plot.window(
      xlim = c(min(astro_solution[, 1]) , max(astro_solution[, 1])),
      ylim = c(2 / 2.54, settings_dev[2] - (3 / 2.54)),
      yaxs = "i",
      xaxs = "i"
    )

    y0_val <- ((settings_dev[2] - ((5) / 2.54)) / 2) + 1.5 / 2.54
    y1_val <-  ((settings_dev[2] - ((5) / 2.54)) / 2) + (2.5) / 2.54


    fact <- (anchor_pts[, 1] - min(proxy_signal[, 1])) /
      (max(proxy_signal[, 1]) - min(proxy_signal[, 1]))
    x1_vals <-
      ((max(astro_solution[, 1]) - min(astro_solution[, 1])) * (1 - fact)) +
      min(astro_solution[, 1])


    segments(
      x0 = anchor_pts[, 2],
      y0 = rep(y0_val, times = (nrow(anchor_pts)))
      ,
      x1 = x1_vals,
      y1 = rep(y1_val, times = (nrow(anchor_pts))),
      col = "black",
      lwd = 2,
      lend = 1
    )
  }
}
