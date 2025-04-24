#' @title plot the windowed timeOpt sedimentation rate estimation
#'
#' @description The \code{\link{plot_win_timeOpt}} function plots a widowed
#'  timeOpt sedimentation rate estimation
#'This function is based on the \code{\link[astrochron:eTimeOpt]{eTimeOpt}}function
#'@param win_timeOpt_result result of the \code{\link{win_timeOpt}}function
#'that needs to be used as input \code{Default=NULL}
#'@param proxy_name the name of the used proxy record \code{Default=NULL}
#'@param  abline_h Add horizontal lines to the plot. Specify the lines as a
#'vector e.g. 2,3,5,6  \code{Default=NULL}
#'@param  abline_v Add vertical lines to the plot. Specify the lines as a
#' vector e.g. 2,3,5,6  \code{Default=NULL}
#'@param  add_lines Add  lines to the wavelet plot input should be matrix
#'with first axis being depth/time the columns after that
#'should be period values  \code{Default=NULL}
#'@param  fig_lts Add a text box \code{Default=NULL}
#'@param  xlab add a label to x-axis \code{Default="depth (m)"}
#'@param  ylab add a label to y-axis  \code{Default="sedrate (cm/kyr)"}
#'@param  sel_parameter select one of the three returns of the
#' \code{\link{win_timeOpt}}function
#'element 1: r_2_envelope matrix
#'element 2: r_2_power matrix
#'element 3: r_2_opt matrix	 \code{Default=3}
#'@param n.levels  Number of color levels \code{Default=100}.
#'
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
#'#plot the windowed timeOpt of the magnetic susceptibility record
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
#'run_multicore = FALSE,
#'verbose=FALSE)
#'
#'plot_win_timeOpt(win_timeOpt_result = mag_win_timeOpt,
#'proxy_name= "mag",
#'abline_h=NULL,
#'abline_v = NULL,
#'add_lines=NULL,
#'fig_lts = NULL,
#'xlab="depth (m)",
#'ylab= "sedrate (cm/kyr)",
#'sel_parameter=3,
#'n.levels=100)
#'}
#'
#' @return
#'The output is a plot of the average spectral power of a windowed timeOpt
#'
#' @export
#' @importFrom stats quantile
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics text
#' @importFrom graphics box
#' @importFrom graphics polygon
#' @importFrom graphics layout
#' @importFrom graphics title
#' @importFrom grDevices rgb
#' @importFrom DescTools Closest
#' @importFrom graphics abline
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom colorRamps blue2green
#' @importFrom colorRamps blue2green2red
#' @importFrom colorRamps blue2red
#' @importFrom colorRamps blue2yellow
#' @importFrom colorRamps cyan2yellow
#' @importFrom colorRamps green2red
#' @importFrom colorRamps magenta2green
#' @importFrom colorRamps matlab.like
#' @importFrom colorRamps matlab.like2
#' @importFrom colorRamps ygobb
#' @importFrom viridis viridis
#' @importFrom viridis magma
#' @importFrom viridis plasma
#' @importFrom viridis inferno
#' @importFrom viridis cividis
#' @importFrom viridis mako
#' @importFrom viridis rocket
#' @importFrom viridis turbo
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom colorRamps blue2green
#' @importFrom colorRamps blue2green2red
#' @importFrom colorRamps blue2red
#' @importFrom colorRamps blue2yellow
#' @importFrom colorRamps cyan2yellow
#' @importFrom colorRamps green2red
#' @importFrom colorRamps magenta2green
#' @importFrom colorRamps matlab.like
#' @importFrom colorRamps matlab.like2
#' @importFrom colorRamps ygobb
#' @importFrom viridis viridis
#' @importFrom viridis magma
#' @importFrom viridis plasma
#' @importFrom viridis inferno
#' @importFrom viridis cividis
#' @importFrom viridis mako
#' @importFrom viridis rocket
#' @importFrom viridis turbo
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors

plot_win_timeOpt <- function(win_timeOpt_result = NULL,
                             proxy_name= NULL,
                             abline_h=NULL,
                             abline_v = NULL,
                             add_lines=NULL,
                             fig_lts = NULL,
                             xlab="depth (m)",
                             ylab= "sedrate (cm/kyr)",
                             sel_parameter=3,
                             n.levels=100){
  plot.new()
  win_timeOpt_res <- win_timeOpt_result
  proxy <- cbind(win_timeOpt_result$x,win_timeOpt_result$y)
  win_timeOpt_means <- win_timeOpt_result[[sel_parameter+3]]
  win_timeOpt_result <- win_timeOpt_result[[sel_parameter]]

  maximum.level = max(win_timeOpt_result)
  power_max_mat.levels = quantile(as.matrix(win_timeOpt_result), probs = seq(from = 0,
                                                                             to = 1, length.out = n.levels + 1))
  color_brewer_Sel <- "grDevices::rainbow(n=n.levels, start = 0, end = 0.7)"
  key.cols <- rev(eval(parse(text = color_brewer_Sel)))
  periodtck = 0.02
  periodtcl = 0.5
  main = NULL
  lwd = 2
  lwd.axis = 1
  legend.params = list(width = 1.2, shrink = 0.9, mar = 5.1,
                       n.ticks = 6, label.digits = 3, label.format = "f", lab = NULL,
                       lab.line = 2.5)
  key.marks = round(seq(from = 0, to = 1, length.out = legend.params$n.ticks) *
                      n.levels)
  key.labels = formatC(as.numeric(power_max_mat.levels), digits = legend.params$label.digits,
                       format = legend.params$label.format)[key.marks + 1]

  sedrate <- win_timeOpt_res[[7]]
  y_axis <- win_timeOpt_res[[7]]
  pmax_avg_sel <- as.matrix(win_timeOpt_result)
  depth <- proxy[,1]
  xlim_vals = c(min(depth), max(depth))

  lowerPeriod <- min(y_axis)

  upperPeriod <- max(y_axis)

  ylim_vals = c(lowerPeriod, upperPeriod)


  layout.matrix <- matrix(c(1, 2, 4, 3), nrow = 2, ncol = 2,
                          byrow = TRUE)
  graphics::layout(mat = layout.matrix, heights = c(0.25,
                                                    1), widths = c(8, 2))
  power_max_mat.levels = quantile(pmax_avg_sel, probs = seq(from = 0,
                                                            to = 1, length.out = n.levels + 1))
  par(mar = c(0, 4, 2, 0))


  plot(y =proxy[,2], x = proxy[,1], type = "l", xaxs = "i",
       xlab = "", ylab = proxy_name, xaxt = "n", xlim = xlim_vals)

  text(x=min(proxy[,1])+2.5, y=max(proxy[,2])*0.9, fig_lts, col = "black", cex = 3)

  par(new = FALSE, mar = c(3, 2, 2, 2), mgp = c(2, 1,
                                                0))
  image(x = seq(from = 0, to = n.levels), y = 1, z = t(matrix(power_max_mat.levels,
                                                              nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
        useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "R^2",
        ylab = "")


  axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
       tck = 0.02, tcl = 1.2)
  mtext(key.labels, side = 1, at = key.marks, line = 0.1,
        las = 2, cex = 0.75)
  box(lwd = lwd.axis)
  par(new = FALSE, mar = c(4, 0, 0, 0.5))


  if(win_timeOpt_res$linLog==1){
    plot(x = win_timeOpt_means, y = sedrate, log = "y",
         type = "l", yaxs = "i", yaxt = "n", xlab = "avg. R^2",
         xaxs = "i", ylim = ylim_vals)}else{
           plot(x = win_timeOpt_means, y = sedrate,
                type = "l", yaxs = "i", yaxt = "n", xlab = "avg. R^2",
                xaxs = "i", ylim = ylim_vals)}


  par(new = FALSE, mar = c(4, 4, 0, 0), xpd = FALSE)


  if(win_timeOpt_res$linLog==1){
    image(
      x = proxy[,1],
      y = log10(sedrate),
      z = t(win_timeOpt_result),
      col = key.cols, breaks = power_max_mat.levels,
      useRaster = TRUE, yaxt = "n", xlab = xlab,
      ylab = ylab)

    period.tick = unique(round(sedrate,1))

    #period.tick[period.tick < log10(sedrate[1])] = NA
    period.tick = na.omit(period.tick)
    period.tick.label = (period.tick)
    axis(2, lwd = lwd.axis, at = log10(period.tick), labels = NA,
         tck = periodtck, tcl = periodtcl)

    mtext(period.tick.label, side = 2, at = log10(period.tick),
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)

    if (is.null(add_lines) != TRUE) {
      for (i  in 2:ncol(add_lines))
        lines(x=add_lines[, 1], y=log10(add_lines[, i]))
    }

  }else{
    image(
      x = proxy[,1],
      y = sedrate,
      z = t(win_timeOpt_result),
      col = key.cols, breaks = power_max_mat.levels,
      useRaster = TRUE, yaxt = "n", xlab = "depth(m)",
      ylab = "")

    period.tick = unique(round(sedrate,1))

    #period.tick[period.tick < log10(sedrate[1])] = NA
    period.tick = na.omit(period.tick)
    period.tick.label = (period.tick)
    axis(2, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)

    mtext(period.tick.label, side = 2, at = period.tick,
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)

    if (is.null(add_lines) != TRUE) {
      for (i  in 2:ncol(add_lines))
        lines(x=add_lines[, 1], y=(add_lines[, i]))
    }



  }

}
















