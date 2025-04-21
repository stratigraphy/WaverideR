#' @title Add a wavelet plot
#'
#' @description Generates a plot of a wavelet scalogram which can be integrated into
#' a larger composite plot
#'
#'@param wavelet wavelet object created using the \code{\link{analyze_wavelet}} function.
#'@param lowerPeriod Lowest period value which will be plotted
#'@param upperPeriod Highest period value which will be plotted
#'@param lower_depth_time lowest depth/time value which will be plotted
#'@param upper_depth_time Highest depth/time value which will be plotted
#'@param plot.COI Option to plot the cone of influence \code{Default=TRUE}.
#'@param n.levels  Number of color levels \code{Default=100}.
#'@param palette_name Name of the color palette which is used for plotting.
#'The color palettes than can be chosen depends on which the R package is specified in
#'the color_brewer parameter. The included R packages from which palettes can be chosen
#'from are; the 'RColorBrewer', 'grDevices', 'ColorRamps' and 'Viridis' R packages.
#'There are many options to choose from so please
#'read the documentation of these packages \code{Default=rainbow}.
#'The R package 'viridis' has the color palette options: “magma”, “plasma”,
#'“inferno”, “viridis”, “mako”, and “rocket”  and “turbo”
#'To see the color palette options of the The R pacakge 'RColorBrewer' run
#'the RColorBrewer::brewer.pal.info() function
#'The R package 'colorRamps' has the color palette options:"blue2green",
#'"blue2green2red", "blue2red",	"blue2yellow", "colorRamps",	"cyan2yellow",
#'"green2red", "magenta2green", "matlab.like", "matlab.like2" and	"ygobb"
#'The R package 'grDevices' has the built in  palette options:"rainbow",
#'"heat.colors", "terrain.colors","topo.colors" and "cm.colors"
#'To see even more color palette options of the The R pacakge 'grDevices' run
#'the grDevices::hcl.pals() function
#'@param color_brewer Name of the R package from which the color palette is chosen from.
#'The included R packages from which palettes can be chosen
#'are; the RColorBrewer, grDevices, ColorRamps and Viridis R packages.
#'There are many options to choose from so please
#'read the documentation of these packages. "\code{Default=grDevices}
#'@param plot_dir The direction of the proxy record which is assumed for tuning if time increases with increasing depth/time values
#'(e.g. bore hole data which gets older with increasing depth ) then plot_dir should be set to TRUE
#'if time decreases with depth/time values (eg stratospheric logs where 0m is the bottom of the section)
#'then plot_dir should be set to FALSE \code{plot_dir=TRUE}
#'@param add_lines Add  lines to the wavelet plot input should be matrix with first axis being depth/time the columns after that
#'should be period values  \code{Default=NULL}
#'@param add_points Add points to the wavelet plot input should be matrix with first axis being depth/time and columns after that
#'should be period values \code{Default=NULL}
#'@param add_abline_h Add horizontal lines to the plot. Specify the lines as a vector e.g. c(2,3,5,6)  \code{Default=NULL}
#'@param add_abline_v Add vertical lines to the plot. Specify the lines as a vector e.g. c(2,3,5,6)  \code{Default=NULL}
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis is depth or y axis power \code{Default=TRUE}
#'@param period_ticks tick mark spacing 1 is all tickmarks and higher value removes tick marks by
#'the fraction of the tick mark spacing value, the opposite is true for value lower than 1 which
#'will add aditional tickmarks
#'@param periodlab lable for the the period column
#'@param main main title
#'@param yaxt turn on of off the yaxis "s" is on "n" is off \code{Default="s"}
#'@param xaxt turn on of off the xaxis "s" is on "n" is off \code{Default="s"}
#'@param depth_time_lab lable for the the depth/time column
#'
#' @author
#' Code based on the "analyze.wavelet" and "wt.image" functions of the 'WaveletComp' R package
#' and  the "wt" function of the 'biwavelet' R package which are based on the
#' wavelet MATLAB code written by Christopher Torrence and Gibert P. Compo (1998).
#' The MTM analysis is from the astrochron R package of Meyers et al., (2012)
#'@references
#'Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
#'Wavelet Analysis. R package version 1.1.
#'\url{https://CRAN.R-project.org/package=WaveletComp}
#'
#'Gouhier TC, Grinsted A, Simko V (2021). R package biwavelet: Conduct Univariate and Bivariate Wavelet Analyses. (Version 0.20.21),
#'\url{https://github.com/tgouhier/biwavelet}
#'
#'Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#'Bulletin of the American Meteorological Society 79:61-78.
#'\url{https://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf}
#'
#'Morlet, Jean, Georges Arens, Eliane Fourgeau, and Dominique Glard.
#'"Wave propagation and sampling theory—Part I: Complex signal and scattering in multilayered media.
#'" Geophysics 47, no. 2 (1982): 203-221.
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#'  Geophysics 1982 47 (2): 222–236.
#'
#' @examples
#' \donttest{
#'#generate a plot for the magnetic susceptibility data set of Pas et al., (2018)
#'
#'plot.new()
#'layout.matrix <- matrix(c(rep(0, 2), 1, 0,0,seq(2, 6, by = 1)),
#'                        nrow = 2,
#'                       ncol = 5 ,
#'                        byrow = TRUE)
#'graphics::layout(mat = layout.matrix,
#'                 heights = c(0.25, 1),
#'                 # Heights of the two rows
#'                 widths = c(rep(c(1, 2, 4,2,2), 2)))
#'
#'par(mar = c(0, 0.5, 1, 0.5))
#'
#'
#'mag_wt <-
#'  analyze_wavelet(
#'    data = mag,
#'    dj = 1 / 100,
#'   lowerPeriod = 0.1,
#'    upperPeriod = 254,
#'    verbose = FALSE,
#'    omega_nr = 10
#'  )
#'
#'  add_wavelet_avg(
#'  wavelet = mag_wt,
#'  plot_horizontal = TRUE,
#'  add_abline_h = NULL,
#'  add_abline_v = NULL,
#'  lowerPeriod = 0.15,
#'  upperPeriod = 80
#')
#'
#'par(mar = c(4, 4, 0, 0.5))
#'
#'
#'plot(
#'  x = c(0, 1),
#'  y = c(max(mag[, 1]), min(mag[, 1])),
#'  col = "white",
#'  xlab = "",
#'  ylab = "Time (Ma)",
#'  xaxt = "n",
#'  xaxs = "i",
#'  yaxs = "i",
#'  ylim = rev(c(max(mag[, 1]), min(mag[, 1])))
#')            # Draw empty plot
#'
#'
#'polygon(
#'  x = c(0, 1, 1, 0),
#'  y = c(max(mag[, 1]), max(mag[, 1]), min(mag[, 1]), min(mag[, 1])),
#'  col = geo_col("Famennian")
#')
#'
#'text(
#'  0.5,
#'  (max(mag[, 1]) - min(mag[, 1])) / 2,
#'  "Fammenian",
#'  cex = 1,
#'  col = "black",
#'  srt = 90
#')
#'par(mar = c(4, 0.5, 0, 0.5))
#'
#'
#'plot(
#'  mag[, 2],
#'  mag[, 1],
#'  type = "l",
#'  ylim = rev(c(max(mag[, 1]), min(mag[, 1]))),
#'  yaxs = "i",
#'  yaxt = "n",
#'  xlab = "Mag. suc.",
#'  ylab = ""
#')
#'
#'
#'add_wavelet(
#'  wavelet = mag_wt,
#'  lowerPeriod = 0.15,
#'  upperPeriod = 80,
#'  lower_depth_time = NULL,
#'  upper_depth_time = NULL,
#'  n.levels = 100,
#'  plot.COI = TRUE,
#'  color_brewer = "grDevices",
#'  palette_name = "rainbow",
#'  plot_dir = FALSE,
#'  add_lines = NULL,
#'  add_points = NULL,
#'  add_abline_h = NULL,
#'  add_abline_v = NULL,
#'  plot_horizontal = TRUE,
#'  period_ticks = 1,
#'  periodlab = "period (m)",
#'  main = NULL,
#'  yaxt = "n",
#'  xaxt = "s",
#'  depth_time_lab = ""
#')
#'
#'
#'lines(log2(mag_track_solution[,2]),mag_track_solution[,1],lwd=4,lty=4)
#'
#'mag_405 <- extract_signal(
#'  tracked_cycle_curve = mag_track_solution,
#'  wavelet = mag_wt,
#'  period_up = 1.2,
#'  period_down = 0.8,
#'  add_mean = TRUE,
#'  tracked_cycle_period = 405,
#'  extract_cycle = 405,
#'  tune = FALSE,
#'  plot_residual = FALSE
#')
#'
#'plot(mag_405[,2],mag_405[,1],type="l",
#'     yaxt="n", yaxs = "i",
#'     xlab="405-kyr ecc")
#'
#'
#'mag_110 <- extract_signal(
#'  tracked_cycle_curve = mag_track_solution,
#'  wavelet = mag_wt,
#'  period_up = 1.25,
#'  period_down = 0.75,
#'  add_mean = TRUE,
#'  tracked_cycle_period = 405,
#'  extract_cycle = 110,
#'  tune = FALSE,
#'  plot_residual = FALSE
#')
#'
#'mag_110_hil <- Hilbert_transform(mag_110,demean=FALSE)
#'
#'plot(mag_110[,2],mag_110[,1],type="l",
#'     yaxt="n", yaxs = "i",
#'     xlab="110-kyr ecc")
#'
#'lines(mag_110_hil[,2],mag_110_hil[,1])
#'}
#'@return
#'returns a plot of a wavelet scalogram
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
#' @importFrom astrochron mtm
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

add_wavelet <- function(wavelet = NULL,
                        lowerPeriod = NULL,
                        upperPeriod = NULL,
                        lower_depth_time = NULL,
                        upper_depth_time = NULL,
                        n.levels = 100,
                        plot.COI = TRUE,
                        color_brewer = "grDevices",
                        palette_name = "rainbow",
                        plot_dir = FALSE,
                        add_lines = NULL,
                        add_points = NULL,
                        add_abline_h = NULL,
                        add_abline_v = NULL,
                        plot_horizontal = TRUE,
                        period_ticks = 1,
                        periodlab = "period (m)",
                        main = NULL,
                        yaxt = "s",
                        xaxt = "s",
                        depth_time_lab = "depth (m)"){
  levels = quantile(wavelet$Power, probs = seq(
    from = 0,
    to = 1,
    length.out = n.levels + 1
  ))


  if (is.null(lower_depth_time) == TRUE) {
    lower_depth_time <- min(wavelet$x)
  }
  if (is.null(upper_depth_time) == TRUE) {
    upper_depth_time <- max(wavelet$x)
  }

  xlim_vals <- c(lower_depth_time, upper_depth_time)


  if (is.null(lowerPeriod) == TRUE) {
    lowerPeriod <- min(wavelet$Period)
  }
  if (is.null(upperPeriod) == TRUE) {
    upperPeriod <- max(wavelet$Period)
  }

  ylim_vals <- c(lowerPeriod, upperPeriod)

  if (plot_dir == FALSE) {
    xlim_vals <- rev(xlim_vals)
  }


  if (color_brewer == "RColorBrewer") {
    key.cols <-
      rev(colorRampPalette(brewer.pal(brewer.pal.info[palette_name, 1], palette_name))(n.levels))

  }

  if (color_brewer == "colorRamps") {
    color_brewer_Sel <-
      paste("colorRamps::", palette_name, "(n=n.levels)")
    key.cols = eval(parse(text = color_brewer_Sel))
  }

  if (color_brewer == "grDevices") {
    if (palette_name == "rainbow") {
      color_brewer_Sel <-
        "grDevices::rainbow(n=n.levels, start = 0, end = 0.7)"
      key.cols <- rev(eval(parse(text = color_brewer_Sel)))
    }
    else if (palette_name == "heat.colors" |
             palette_name == "terrain.colors" |
             palette_name == "topo.colors" |
             palette_name == "cm.colors") {
      color_brewer_Sel <-
        paste("grDevices::",
              palette_name,
              "(n=n.levels, start = 0, end = 1)")
      key.cols <- rev(eval(parse(text = color_brewer_Sel)))
    }
    else{
      key.cols <-
        hcl.colors(
          n = n.levels,
          palette = palette_name,
          alpha = NULL,
          rev = FALSE,
          fixup = TRUE
        )
    }
  }

  if (color_brewer == "viridis") {
    color_brewer_Sel <-
      paste("viridis::", palette_name, "(n=n.levels,direction = -1)")
    key.cols = rev(eval(parse(text = color_brewer_Sel)))
  }


  maximum.level = max(wavelet$Power)
  power_max_mat.levels = quantile(wavelet$Power, probs = seq(from = 0,
                                                             to = 1, length.out = n.levels + 1))

  if (plot_horizontal == FALSE) {
    image(
      y = wavelet$x,
      x = wavelet$axis.2,
      z = (wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xlab = periodlab,
      ylab = depth_time_lab,
      #axes = FALSE,
      xaxt = "n" ,
      yaxt = yaxt ,
      main = main,
      ylim = xlim_vals,
      xlim = log2(c(lowerPeriod, upperPeriod))
    )

    if (plot.COI == TRUE) {
      polygon(
        y = wavelet$coi.1 ,
        x = wavelet$coi.2,
        border = NA,
        col = rgb(1, 1, 1, 0.5),
        ylim =  xlim_vals,
        xlim = log2(c(lowerPeriod, upperPeriod))
      )
    }

    periodtck = 0.02
    periodtcl = 0.5
    main = NULL
    lwd = 2
    lwd.axis = 1
    box(lwd = lwd.axis)

    period.tick = unique(trunc(wavelet$axis.2))
    period.tick <- period.tick[period.tick >= log2(ylim_vals[1])]
    period.tick <- period.tick[period.tick <= log2(ylim_vals[2])]
    nrs <-
      seq(period.tick[1], period.tick[length(period.tick)], by = period_ticks)
    nrs

    period.tick <- nrs
    period.tick[period.tick < log2(wavelet$Period[1])] = NA
    period.tick = na.omit(period.tick)
    period.tick.label = 2 ^ (period.tick)


    axis(
      1,
      lwd = lwd.axis,
      at = period.tick,
      labels = NA,
      tck = periodtck,
      tcl = periodtcl
    )

    mtext(
      period.tick.label,
      side = 1,
      at = period.tick,
      las = 2,
      line = par()$mgp[2] - 0.5,
      font = par()$font.axis,
      cex = par()$cex.axis
    )

    if (is.null(add_lines) != TRUE) {
      for (i  in 2:ncol(add_lines))
        lines(y = add_lines[, 1], x = log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        lines(y = add_points[, 1], x = log2(add_points[, i]))
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }

  }



  if (plot_horizontal == TRUE) {
    image(
      x = wavelet$x,
      y = wavelet$axis.2,
      z = t(wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xlab = depth_time_lab,
      ylab = periodlab,
      xaxt = "s",
      yaxt = "n" ,
      main = main,
      xlim = xlim_vals,
      ylim = log2(c(lowerPeriod, upperPeriod))
    )


    if (plot.COI == TRUE) {
      polygon(
        x = wavelet$coi.1 ,
        y = wavelet$coi.2,
        border = NA,
        col = rgb(1, 1, 1, 0.5),
        xlim =  xlim_vals,
        ylim = log2(c(lowerPeriod, upperPeriod))
      )
    }

    periodtck = 0.02
    periodtcl = 0.5
    main = NULL
    lwd = 2
    lwd.axis = 1
    box(lwd = lwd.axis)

    log2(ylim_vals[1])
    log2(ylim_vals[2])

    period.tick = unique(trunc(wavelet$axis.2))
    period.tick <- period.tick[period.tick >= log2(ylim_vals[1])]
    period.tick <- period.tick[period.tick <= log2(ylim_vals[2])]
    nrs <-
      seq(period.tick[1], period.tick[length(period.tick)], by = period_ticks)
    nrs

    period.tick <- nrs
    period.tick[period.tick < log2(wavelet$Period[1])] = NA
    period.tick = na.omit(period.tick)
    period.tick.label = 2 ^ (period.tick)


    axis(
      2,
      lwd = lwd.axis,
      at = period.tick,
      labels = NA,
      tck = periodtck,
      tcl = periodtcl
    )

    mtext(
      period.tick.label,
      side = 2,
      at = period.tick,
      las = 2,
      line = par()$mgp[2] - 0.5,
      font = par()$font.axis,
      cex = par()$cex.axis
    )

    if (is.null(add_lines) != TRUE) {
      for (i  in 2:ncol(add_lines))
        lines(x = add_lines[, 1], y = log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        lines(x = add_points[, 1], y = log2(add_points[, i]))
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = (add_abline_v))
    }
  }
}
