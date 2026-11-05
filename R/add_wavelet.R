#' @title Add a wavelet plot
#'
#' @description Generates a plot of a wavelet scalogram which can be integrated into
#' a larger composite plot
#'
#'@param wavelet Wavelet object (list) created with \code{\link{analyse_wavelet}} (or another function returning
#'the same structure, with elements \code{x}, \code{Period}, \code{Power}, \code{axis.2}, \code{coi.1} and
#'\code{coi.2}). Required (\code{Default=NULL} gives an error).
#'@param lowerPeriod Shortest period shown on the period axis, in the same unit as the depth/time axis of the
#'wavelet (numeric, > 0); \code{NULL} uses the shortest period of \code{wavelet}. \code{Default=NULL}.
#'@param upperPeriod Longest period shown on the period axis, in the same unit as the depth/time axis of the
#'wavelet (numeric, > \code{lowerPeriod}); \code{NULL} uses the longest period of \code{wavelet}.
#'\code{Default=NULL}.
#'@param lower_depth_time Lowest depth/time value shown on the depth/time axis (numeric, same unit as the
#'record); \code{NULL} uses the minimum depth/time of \code{wavelet}. \code{Default=NULL}.
#'@param upper_depth_time Highest depth/time value shown on the depth/time axis (numeric, same unit as the
#'record); \code{NULL} uses the maximum depth/time of \code{wavelet}. \code{Default=NULL}.
#'@param plot.COI Draw the cone of influence as a semi-transparent white overlay (logical, TRUE or FALSE).
#'\code{Default=TRUE}.
#'@param n.levels Number of colour levels of the power scale (integer, > 1); the colour breaks are quantiles of
#'the wavelet power, so each colour covers the same fraction of the scalogram. \code{Default=100}.
#' @param palette_name Name of the colour palette of the power (colour) scale, a character string. Which names are
#' valid depends on \code{color_brewer}:
#' \itemize{
#' \item \code{color_brewer = "grDevices"}: \code{"rainbow"}, \code{"heat.colors"}, \code{"terrain.colors"}, \code{"topo.colors"}, \code{"cm.colors"},
#' or any of the 115 \code{\link[grDevices]{hcl.colors}} palettes (list them with \code{grDevices::hcl.pals()}), for
#' example \code{"viridis"}, \code{"Plasma"}, \code{"Inferno"}, \code{"YlOrRd"}, \code{"Spectral"}, \code{"Blue-Red"}, \code{"Zissou 1"}.
#' \item \code{color_brewer = "viridis"}: \code{"viridis"}, \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"cividis"}, \code{"rocket"}, \code{"mako"}, \code{"turbo"}.
#' \item \code{color_brewer = "RColorBrewer"}: \code{"BrBG"}, \code{"PiYG"}, \code{"PRGn"}, \code{"PuOr"}, \code{"RdBu"}, \code{"RdGy"}, \code{"RdYlBu"}, \code{"RdYlGn"}, \code{"Spectral"}, \code{"Accent"}, \code{"Dark2"}, \code{"Paired"}, \code{"Pastel1"}, \code{"Pastel2"}, \code{"Set1"}, \code{"Set2"}, \code{"Set3"}, \code{"Blues"}, \code{"BuGn"}, \code{"BuPu"}, \code{"GnBu"}, \code{"Greens"}, \code{"Greys"}, \code{"Oranges"}, \code{"OrRd"}, \code{"PuBu"}, \code{"PuBuGn"}, \code{"PuRd"}, \code{"Purples"}, \code{"RdPu"}, \code{"Reds"}, \code{"YlGn"}, \code{"YlGnBu"}, \code{"YlOrBr"}, \code{"YlOrRd"} (also listed by \code{wr_brewer_palettes()}).
#' \item \code{color_brewer = "colorRamps"}: \code{"blue2green"}, \code{"blue2green2red"}, \code{"blue2red"}, \code{"blue2yellow"}, \code{"cyan2yellow"}, \code{"green2red"}, \code{"magenta2green"}, \code{"matlab.like"}, \code{"matlab.like2"}, \code{"ygobb"}.
#' \item \code{color_brewer = "scico"} (all except \code{"batlowK"} need the \pkg{scico} package): \code{"acton"}, \code{"bam"}, \code{"bamako"}, \code{"bamO"}, \code{"batlow"}, \code{"batlowK"}, \code{"batlowW"}, \code{"berlin"}, \code{"bilbao"}, \code{"broc"}, \code{"brocO"}, \code{"buda"}, \code{"bukavu"}, \code{"cork"}, \code{"corkO"}, \code{"davos"}, \code{"devon"}, \code{"fes"}, \code{"glasgow"}, \code{"grayC"}, \code{"hawaii"}, \code{"imola"}, \code{"lajolla"}, \code{"lapaz"}, \code{"lipari"}, \code{"lisbon"}, \code{"managua"}, \code{"navia"}, \code{"nuuk"}, \code{"oleron"}, \code{"oslo"}, \code{"roma"}, \code{"romaO"}, \code{"tofino"}, \code{"tokyo"}, \code{"turku"}, \code{"vanimo"}, \code{"vik"}, \code{"vikO"}.
#' }
#' Low power is drawn with the first colour of the palette and high power with the last. \code{Default="viridis"}.
#' @param color_brewer Palette family from which \code{palette_name} is taken, a character string: \code{"grDevices"},
#' \code{"viridis"}, \code{"RColorBrewer"}, \code{"colorRamps"} or \code{"scico"} (the last one, except for \code{"batlowK"}, needs the
#' \pkg{scico} package). See \code{palette_name} for the palette names of each family. \code{Default="viridis"}.
#'@param plot_dir Direction of the depth/time axis (logical, TRUE or FALSE). \code{TRUE} draws depth/time
#'increasing away from the origin (e.g. a borehole where age increases with depth); \code{FALSE} reverses the
#'depth/time axis (e.g. an outcrop log where 0 m is the base of the section). \code{Default=FALSE}.
#'@param add_lines Curves to draw on top of the scalogram (e.g. a tracked period): a matrix or data frame with
#'depth/time in column 1 and one or more columns of period values (same unit as the period axis) after it;
#'every period column is drawn as a line. \code{NULL} draws none. \code{Default=NULL}.
#'@param add_points Same input as \code{add_lines} (matrix or data frame, depth/time in column 1, period values
#'in the next columns); the period columns are currently also drawn as lines. \code{NULL} draws none.
#'\code{Default=NULL}.
#'@param add_abline_h Positions of horizontal reference lines (numeric vector, e.g. \code{c(2, 3, 5, 6)}). With
#'\code{plot_horizontal=TRUE} these are periods (log2-transformed internally to match the period axis); with
#'\code{plot_horizontal=FALSE} they are depth/time values. \code{NULL} draws none. \code{Default=NULL}.
#'@param add_abline_v Positions of vertical reference lines (numeric vector, e.g. \code{c(2, 3, 5, 6)}). With
#'\code{plot_horizontal=TRUE} these are depth/time values; with \code{plot_horizontal=FALSE} they are periods
#'(log2-transformed internally). \code{NULL} draws none. \code{Default=NULL}.
#'@param plot_horizontal Orientation of the plot (logical, TRUE or FALSE): \code{TRUE} puts depth/time on the
#'x axis and period on the y axis; \code{FALSE} puts period on the x axis and depth/time on the y axis (for
#'vertical log-style composite plots). \code{Default=TRUE}.
#'@param period_ticks Spacing of the period tick marks in octaves (numeric, > 0): a tick is placed every
#'\code{period_ticks} powers of 2, so 1 labels every power of 2, 2 every second one, and values below 1 add
#'intermediate ticks. \code{Default=1}.
#'@param periodlab Label of the period axis (character string). \code{Default="period (m)"}.
#'@param main Main title of the plot (character string); \code{NULL} gives no title. \code{Default=NULL}.
#'@param yaxt Draw the depth/time axis when \code{plot_horizontal=FALSE} (character): \code{"s"} draws it,
#'\code{"n"} suppresses it (see \code{\link[graphics]{par}}). Ignored when \code{plot_horizontal=TRUE}, where the
#'period axis is always drawn by the function itself. \code{Default="s"}.
#'@param xaxt Accepted for symmetry with \code{yaxt} but currently not used: the x axis is always drawn
#'(\code{plot_horizontal=TRUE}) or replaced by the period axis (\code{plot_horizontal=FALSE}) (character,
#'\code{"s"} or \code{"n"}). \code{Default="s"}.
#'@param depth_time_lab Label of the depth/time axis (character string). \code{Default="depth (m)"}.
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
#'oldpar <- par(no.readonly = TRUE)
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
#'  analyse_wavelet(
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
#'  plot_horizontal = FALSE,
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
#'# amplitude of the ~110-kyr cycle from the wavelet: the tracked 405-kyr curve scaled to 110 kyr
#'track_110 <- cbind(mag_track_solution[, 1], mag_track_solution[, 2] * 110 / 405)
#'mag_110_amp <- extract_amplitude(mag_110, cycle = track_110, band = c(0.75, 1.33))
#'mag_110_mean <- mean(mag_110[, 2])
#'
#'plot(mag_110[,2],mag_110[,1],type="l",
#'     yaxt="n", yaxs = "i",
#'     xlab="110-kyr ecc")
#'
#'lines(mag_110_mean + mag_110_amp[, 2], mag_110_amp[, 1], col = "red")
#'lines(mag_110_mean - mag_110_amp[, 2], mag_110_amp[, 1], col = "red")
#'par(oldpar)
#'}
#'@return
#'returns a plot of a wavelet scalogram
#' @param colour_brewer British spelling of \code{color_brewer}, with the same options; if given (not \code{NULL}) it overrides
#' \code{color_brewer}. \code{Default=NULL}.
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
#' @importFrom graphics abline
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors
#' @importFrom grDevices colorRampPalette
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
    color_brewer = "viridis",
    palette_name = "viridis",
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
    depth_time_lab = "depth (m)",
    colour_brewer = NULL){
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym
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


  key.cols <- .wr_key_cols(palette_name, color_brewer, n.levels)


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
      for (i  in 2:ncol(add_lines)) {
        lines(y = add_lines[, 1], x = log2(add_lines[, i]))
      }
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points)) {
        lines(y = add_points[, 1], x = log2(add_points[, i]))
      }
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
      for (i  in 2:ncol(add_lines)) {
        lines(x = add_lines[, 1], y = log2(add_lines[, i]))
      }
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points)) {
        lines(x = add_points[, 1], y = log2(add_points[, i]))
      }
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = (add_abline_v))
    }
  }
}
