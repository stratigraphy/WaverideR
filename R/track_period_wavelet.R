#' @title Track the period of a cycle in a wavelet spectra
#'
#' @description Interactively select points in a wavelet spectra to trace a period in a wavelet spectra.
#'The \code{\link{track_period_wavelet}} function plots a wavelet spectra in which spectral peaks can selected
#'allowing one to track a ridge hence one can track the a cycle with a changing period.
#'Tracking points can be selected in the Interactive interface and will be shown as white dots
#'when one wants to deselect a point the white dots can be re-clicked/re-selected and will turn red which
#'indicates that the previously selected point is deselected. Deselecting points can be quite tricky
#'due to the close spacing of  points and such the \code{\link{delpts_tracked_period_wt}} can be used to
#'delete points were previously selected using the \code{\link{track_period_wavelet}} function.
#' @param wavelet Wavelet object (class \code{analyze.wavelet}) created with \code{\link{analyse_wavelet}} in
#' which the cycle is tracked; the elements \code{x}, \code{y}, \code{Power}, \code{Period}, \code{axis.2},
#' \code{coi.1}, \code{coi.2} and \code{Power.avg} are used. Must be given (the \code{NULL} default fails).
#' \code{Default=NULL}.
#' @param astro_cycle Duration of the astronomical cycle that is tracked (numeric, > 0, in kyr, e.g. 405 for the long
#' eccentricity cycle). It is used to convert the tracked period into the sedimentation rate of the output
#' (period / \code{astro_cycle} * 100, i.e. cm/kyr for a record in metres). \code{Default=405}.
#' @param n.levels Number of colour levels of the power scale (integer, > 1); the breaks are quantiles of the power.
#' \code{Default=100}.
#' @param track_peaks Restrict the selection to spectral peaks (logical, TRUE or FALSE). With \code{TRUE} the local
#' power maxima along the period axis are drawn as black circles and a click selects the nearest one; with
#' \code{FALSE} any point of the wavelet spectrum can be clicked (the click is snapped to the nearest depth of the wavelet
#' and its period is used as clicked). In both modes clicking a selected point again removes it; finish with
#' right-click/Esc (or "Finish" in RStudio). \code{Default=TRUE}.
#' @param periodlab Label of the period axis (character). \code{Default="Period (metres)"}.
#' @param x_lab Label of the depth/time axis (character). \code{Default="depth (metres)"}.
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
#' @param plot_horizontal Layout of the plot (logical, TRUE or FALSE): \code{TRUE} depth/time on the x-axis and period
#' on the y-axis; \code{FALSE} period on the x-axis and depth/time on the y-axis. The output is the same in both layouts.
#' \code{Default=TRUE}.
#' @param lowerPeriod Shortest period shown on the period axis (numeric, in the depth/time unit of the record);
#' \code{NULL} uses the shortest period of \code{wavelet}. \code{Default=NULL}.
#' @param upperPeriod Longest period shown on the period axis (numeric, in the depth/time unit of the record);
#' \code{NULL} uses the longest period of \code{wavelet}. \code{Default=NULL}.
#' @param plot_dir Direction of the depth/time axis (logical, TRUE or FALSE): \code{TRUE} plots it from low to high
#' values (e.g. borehole data that get older with increasing depth); \code{FALSE} reverses the axis (e.g. outcrop logs
#' in which 0 m is the base of the section). \code{Default=TRUE}.
#' @param add_lines Optional lines to draw on the wavelet plot (matrix or data frame: depth/time in column 1, one or more
#' columns of periods after it; one line per period column). \code{Default=NULL}.
#' @param add_points Optional points to draw on the wavelet plot (matrix or data frame: depth/time in column 1, one or
#' more columns of periods after it). \code{Default=NULL}.
#' @param add_abline_h Optional horizontal reference lines (numeric vector, e.g. \code{c(2, 3, 5)}), drawn at these
#' depths/times; only used in the vertical layout (\code{plot_horizontal = FALSE}), the horizontal layout ignores
#' it. \code{Default=NULL}.
#' @param add_abline_v Optional vertical reference lines (numeric vector), drawn at these periods; only used in the
#' vertical layout (\code{plot_horizontal = FALSE}), the horizontal layout ignores it. \code{Default=NULL}.
#'@return Results of the tracking of a cycle in the wavelet spectra is a matrix with 3 columns.
#'The first column is depth/time
#'The second column is the period of the tracked cycle
#'The third column is the sedimentation rate based on the duration (in time) of the tracked cycle
#'
#' @author
#' The function is based/inspired on the \link[astrochron]{traceFreq}
#'function of the 'astrochron' R package
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#'@examples
#'\donttest{
#'if (interactive()) { # needs mouse input (identify())
#'#Track the 405kyr eccentricity cycle in the magnetic susceptibility record
#'# of the Sullivan core of Pas et al., (2018)
#'
#'mag_wt <- analyse_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#' mag_track <- track_period_wavelet(wavelet = mag_wt,
#'astro_cycle = 405,
#'n.levels = 100,
#'track_peaks = TRUE,
#'periodlab =  "Period (metres)",
#'x_lab = "depth (metres)",
#'palette_name = "rainbow",
#'color_brewer = "grDevices",
#'plot_horizontal = TRUE,
#'plot_dir = TRUE,
#'lowerPeriod = NULL,
#'upperPeriod = NULL,
#'add_lines = NULL,
#'add_points = NULL,
#'add_abline_h = NULL,
#'add_abline_v = NULL)
#'}}
#'
#' @param colour_brewer British spelling of \code{color_brewer}, with the same options; if given (not \code{NULL}) it overrides
#' \code{color_brewer}. \code{Default=NULL}.
#' @export
#' @importFrom stats quantile
#' @importFrom graphics par
#' @importFrom grDevices dev.new
#' @importFrom graphics image
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics text
#' @importFrom graphics box
#' @importFrom graphics polygon
#' @importFrom grDevices rgb
#' @importFrom graphics points
#' @importFrom stats aggregate
#' @importFrom stats na.omit
#' @importFrom astrochron traceFreq
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors

track_period_wavelet <- function(wavelet = NULL,
    astro_cycle = 405,
    n.levels = 100,
    track_peaks = TRUE,
    periodlab =  "Period (metres)",
    x_lab = "depth (metres)",
    palette_name = "viridis",
    color_brewer = "viridis",
    plot_horizontal = TRUE,
    plot_dir = TRUE,
    lowerPeriod = NULL,
    upperPeriod = NULL,
    add_lines = NULL,
    add_points = NULL,
    add_abline_h = NULL,
    add_abline_v = NULL,
    colour_brewer = NULL)
{
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))



  maximum.level = max(wavelet$Power)
  power_max_mat.levels = quantile(wavelet$Power, probs = seq(
      from = 0,
      to = 1,
      length.out = n.levels + 1
  ))


  key.cols <- .wr_key_cols(palette_name, color_brewer, n.levels)


  periodtck = 0.02
  periodtcl = 0.5
  main = NULL
  lwd = 2
  lwd.axis = 1

  legend.params = list(
      width = 1.2,
      shrink = 0.9,
      mar = 5.1,
      n.ticks = 6,
      label.digits = 3,
      label.format = "f",
      lab = NULL,
      lab.line = 2.5
  )

  key.marks = round(seq(
      from = 0,
      to = 1,
      length.out = legend.params$n.ticks
  ) *
      n.levels)

  key.labels = formatC(
      as.numeric(power_max_mat.levels),
      digits = legend.params$label.digits,
      format = legend.params$label.format
  )[key.marks + 1]








  y_axis <- as.numeric(unlist(wavelet$Period))
  pmax_avg_sel <- t(wavelet$Power)

  depth <-  wavelet$x
  y_axis <- wavelet$Period
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)




  if (plot_dir != TRUE) {
    xlim_vals = rev(c(min(wavelet$x), max(wavelet$x)))
  } else{
    xlim_vals = c(min(wavelet$x), max(wavelet$x))
  }


  if (is.null(lowerPeriod) == TRUE) {
    lowerPeriod  <- min(wavelet$Period)
  }
  if (is.null(upperPeriod) == TRUE) {
    upperPeriod  <- max(wavelet$Period)
  }

  ylim_vals = c(lowerPeriod, upperPeriod)

  if (plot_horizontal == TRUE) {
    wr_open_device(width = 15, height = 7, interactive_input = TRUE)


    layout.matrix <- matrix(c(1, 2, 4, 3),
        nrow = 2,
        ncol = 2 ,
        byrow = TRUE)
    graphics::layout(mat = layout.matrix,
        heights = c(0.25, 1),
        # Heights of the two rows
        widths = c(8, 2))

    power_max_mat.levels = quantile(pmax_avg_sel,
        probs = seq(
            from = 0,
            to = 1,
            length.out = n.levels + 1
        ))

    par(mar = c(0, 4, 2, 0))

    plot(
        y = wavelet$y,
        x = wavelet$x,
        type = "l",
        xaxs = "i",
        xlab = "",
        ylab = "proxy value",
        xaxt = "n",
        xlim = xlim_vals

    )

    add_abline_v = FALSE

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }


    par(new = FALSE,
        mar = c(3, 2, 2, 2),
        mgp = c(2, 1, 0))

    image(
        x = seq(from = 0, to = n.levels),
        y = 1,
        z = t(matrix(power_max_mat.levels,
            nrow = 1)),
        col = key.cols,
        breaks = power_max_mat.levels,
        useRaster = TRUE,
        yaxt = "n",
        xaxt = "n",
        xlab = "Power",
        ylab = ""
    )

    axis(
        1,
        lwd = lwd.axis,
        at = key.marks,
        labels = NA,
        tck = 0.02,
        tcl =  1.2
    )

    mtext(
        key.labels,
        side = 1,
        at = key.marks,
        line = 0.1,
        las = 2,
        cex = 0.75
    )
    box(lwd = lwd.axis)

    par(new = FALSE, mar = c(4, 0, 0, 0.5))

    plot(
        x = wavelet$Power.avg,
        y = wavelet$Period,
        log = "y",
        type = "l",
        yaxs = "i",
        yaxt = "n",
        xlab = "Wt. power",
        xaxs = "i",
        ylim = ylim_vals
    )


    add_abline_h <- FALSE

    if (is.null(add_abline_h) != TRUE) {
      abline(h = add_abline_h)
    }


    par(new = FALSE,
        mar = c(4, 4, 0, 0),
        xpd = FALSE)

    image(
        x = wavelet$x,
        y = wavelet$axis.2,
        z = t(wavelet$Power),
        col = key.cols,
        breaks =  power_max_mat.levels,
        useRaster = TRUE,
        ylab = periodlab,
        xlab = x_lab,
        axes = TRUE,
        yaxt = "n" ,
        main = main,
        xlim = xlim_vals,
        ylim = log2(ylim_vals)
    )


    polygon(
        x = wavelet$coi.1 ,
        y = wavelet$coi.2,
        border = NA,
        col = rgb(1, 1, 1, 0.5),
        ylim = xlim_vals,
        xlim = log2(ylim_vals)
    )

    box(lwd = lwd.axis)
    period.tick = unique(trunc(wavelet$axis.2))
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
    axis(
        4,
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
        lines(add_lines[, 1], log2(add_lines[, i]))
      }
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points)) {
        points(add_points[, 1], log2(add_points[, i]))
      }
    }


    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }



    if (track_peaks == TRUE) {
      Pwert <- wavelet$Power

      maxdetect <-
          matrix(nrow = (nrow(Pwert)), ncol = ncol(Pwert), 0)

      for (j in 1:ncol(Pwert)) {
        for (i in 2:(nrow(maxdetect) - 1)) {
          if ((Pwert[i, j] - Pwert[(i + 1), j] > 0) &
              (Pwert[i, j] - Pwert[(i - 1), j]  > 0))
          {
            maxdetect[i, j] <- 1
          }
        }
      }

      maxdetect2 <- wr_melt(maxdetect)



      depth <- rep(wavelet$x, each = length(wavelet$axis.2))
      period <- rep(wavelet$axis.2, times = length(wavelet$x))

      maxdetect2 <- as.data.frame(maxdetect2)
      maxdetect2[, 2] <- period
      maxdetect2[, 1] <- depth
      maxdetect2 <- maxdetect2[maxdetect2$value > 0, ]

      colnames(maxdetect2) <- c("y_val", "x_val", "ridge")

      points(
          y = maxdetect2$x_val,
          x = maxdetect2$y_val,
          type = "p",
          pch = 1,
          col = "black",
          lwd = "0.5"
      )

      n <- nrow(maxdetect2)
      y = maxdetect2$x_val
      x = maxdetect2$y_val
    }

    out <- .tpw_pick(wavelet, astro_cycle, horizontal = TRUE,
        candidates = if (isTRUE(track_peaks)) { maxdetect2 } else { NULL })
  }


  if (plot_horizontal == FALSE) {
    wr_open_device(width = 7, height = 10, interactive_input = TRUE)




    layout.matrix <- matrix(c(1, 2, 3, 4),
        nrow = 2,
        ncol = 2 ,
        byrow = TRUE)
    graphics::layout(mat = layout.matrix,
        heights = c(1, 4),
        # Heights of the two rows
        widths = c(1, 4))

    power_max_mat.levels = quantile(pmax_avg_sel,
        probs = seq(
            from = 0,
            to = 1,
            length.out = n.levels + 1
        ))

    par(mar = c(2, 0.5, 2, 6), xpd = NA)


    image(
        y = seq(from = 0, to = n.levels),
        x = 1,
        z = (matrix(power_max_mat.levels,
            nrow = 1)),
        col = key.cols,
        breaks = power_max_mat.levels,
        useRaster = TRUE,
        yaxt = "n",
        xaxt = "n",
        xlab = "",
        ylab = "",
    )


    axis(
        2,
        lwd = lwd.axis,
        at = key.marks,
        labels = NA,
        tck = 0.02,
        tcl = 1.24
    )

    mtext(
        key.labels,
        side = 2,
        at = key.marks,
        line = 0.5,
        las = 2,
        font = par()$font.axis,
        cex = par()$cex.axis
    )
    mtext(
        c("Power"),
        side = 1,
        at = 1,
        line = 0.5,
        font = par()$font.axis,
        cex = par()$cex.axis,
        las = 1,
        xpd = NA
    )
    #

    box(lwd = lwd.axis)



    par(mar = c(0, 0, 2, 2))


    plot(
        y = wavelet$Power.avg,
        x = wavelet$Period,
        log = "x",
        type = "l",
        xaxs = "i",
        xaxt = "n",
        ylab = "Wt. power",
        xlab = "",
        xaxs = "i",
        xlim = ylim_vals
    )

    title(ylab = "Wt. power", xpd = NA)


    par(mar = c(4, 4, 0, 0), xpd = TRUE)


    plot(
        x = wavelet$y,
        y = wavelet$x,
        type = "l",
        yaxs = "i",
        xlab = "proxy value",
        ylab = x_lab,
        #yaxt = "n",
        ylim = xlim_vals

    )


    par(new = FALSE,
        mar = c(4, 0, 0, 2),
        xpd = FALSE)

    image(
        y = wavelet$x,
        x = wavelet$axis.2,
        z = (wavelet$Power),
        col = key.cols,
        breaks = power_max_mat.levels,
        useRaster = TRUE,
        xlab = periodlab,
        ylab = "",
        #axes = FALSE,
        #yaxt = "n" ,
        xaxt = "n" ,
        yaxt = "n" ,
        main = main,
        ylim = xlim_vals,
        xlim = log2(ylim_vals)
    )

    polygon(
        y = wavelet$coi.1 ,
        x = wavelet$coi.2,
        border = NA,
        col = rgb(1, 1, 1, 0.5),
        ylim = xlim_vals,
        xlim = log2(ylim_vals)
    )


    box(lwd = lwd.axis)

    period.tick = unique(trunc(wavelet$axis.2))
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
        points(y = add_points[, 1], x = log2(add_points[, i]))
      }
    }



    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }



    #track_peaks <- TRUE

    if (track_peaks == TRUE) {
      Pwert <- wavelet$Power

      maxdetect <-
          matrix(nrow = (nrow(Pwert)), ncol = ncol(Pwert), 0)

      for (j in 1:ncol(Pwert)) {
        for (i in 2:(nrow(maxdetect) - 1)) {
          if ((Pwert[i, j] - Pwert[(i + 1), j] > 0) &
              (Pwert[i, j] - Pwert[(i - 1), j]  > 0))
          {
            maxdetect[i, j] <- 1
          }
        }
      }


      maxdetect2 <- wr_melt(maxdetect)



      depth <- rep(wavelet$x, each = length(wavelet$axis.2))
      period <- rep(wavelet$axis.2, times = length(wavelet$x))

      maxdetect2 <- as.data.frame(maxdetect2)
      maxdetect2[, 2] <- period
      maxdetect2[, 1] <- depth
      maxdetect2 <- maxdetect2[maxdetect2$value > 0, ]

      colnames(maxdetect2) <- c("y_val", "x_val", "ridge")

      points(
          x = maxdetect2$x_val,
          y = maxdetect2$y_val,
          type = "p",
          pch = 1,
          col = "black",
          lwd = "0.5"
      )

      n <- nrow(maxdetect2)
      x = maxdetect2$x_val
    }

    out <- .tpw_pick(wavelet, astro_cycle, horizontal = FALSE,
        candidates = if (isTRUE(track_peaks)) { maxdetect2 } else { NULL })
  }





  return(out)
}

# Interactive picking used by track_period_wavelet()
## Mouse input goes through these two wrappers so that the tests can replace them by scripted clicks.
.tpw_identify <- function(x, y) graphics::identify(x, y, n = 1, plot = FALSE, tolerance = 0.1)
.tpw_locator  <- function() graphics::locator(n = 1)

#' Pick points on the scalogram drawn by track_period_wavelet()
#'
#' Plot coordinates: horizontal layout x = depth, y = log2(period); vertical layout x = log2(period), y = depth.
#' With \code{candidates} (the local power maxima, columns 1 = depth, 2 = log2 period) a click selects the nearest
#' maximum (identify()); without, any point of the scalogram can be clicked (locator()): the click is snapped to the
#' nearest depth of the wavelet and its period is taken as clicked. Clicking a selected point again removes it.
#' Finish with right-click / Esc (or "Finish" in RStudio).
#' @return data.frame(depth, period, sedrate), sorted by depth, one row per depth.
#' @noRd
.tpw_pick <- function(wavelet, astro_cycle, horizontal = TRUE, candidates = NULL) {
  to_plot <- function(d, lp) if (horizontal) { list(x = d, y = lp) } else { list(x = lp, y = d) }
  usr <- graphics::par("usr")
  tolx <- 0.01 * abs(diff(usr[1:2]))
  toly <- 0.01 * abs(diff(usr[3:4]))
  if (!is.null(candidates)) {
    cd <- as.numeric(candidates[, 1])
    clp <- as.numeric(candidates[, 2])
    pp <- to_plot(cd, clp)
    sel <- logical(length(cd))
    repeat {
      ans <- .tpw_identify(pp$x, pp$y)
      if (!length(ans)) {
        break
      }
      sel[ans] <- !sel[ans]
      graphics::points(pp$x[ans], pp$y[ans], pch = 19, col = if (sel[ans]) { "white" } else { "red" })
    }
    d <- cd[sel]
    lp <- clp[sel]
  } else {
    d <- numeric(0)
    lp <- numeric(0)
    repeat {
      cl <- .tpw_locator()
      if (is.null(cl) || !length(cl$x)) {
        break
      }
      cd <- if (horizontal) { cl$x } else { cl$y }
      clp <- if (horizontal) { cl$y } else { cl$x }
      cd <- wavelet$x[which.min(abs(wavelet$x - cd))]            # snap to the depth grid of the wavelet
      pp <- to_plot(d, lp)
      hit <- which(abs(pp$x - (if (horizontal) { cd } else { clp })) < tolx &
          abs(pp$y - (if (horizontal) { clp } else { cd })) < toly)
      if (length(hit)) {                                           # clicked an existing pick: remove it
        q <- to_plot(d[hit], lp[hit])
        graphics::points(q$x, q$y, pch = 19, col = "red")
        d <- d[-hit]
        lp <- lp[-hit]
      } else {
        d <- c(d, cd)
        lp <- c(lp, clp)
        q <- to_plot(cd, clp)
        graphics::points(q$x, q$y, pch = 19, col = "white")
      }
    }
  }
  .tpw_output(d, lp, astro_cycle)
}

#' depth + log2(period) of the picks -> data.frame(depth, period, sedrate) for loess_auto() / curve2time()
#' @noRd
.tpw_output <- function(depth, log2period, astro_cycle) {
  ok <- is.finite(depth) & is.finite(log2period)
  if (!any(ok)) {
    return(data.frame(depth = numeric(0), period = numeric(0), sedrate = numeric(0)))
  }
  a <- stats::aggregate(log2period[ok], by = list(depth = depth[ok]), FUN = mean)   # sorted by depth
  out <- data.frame(depth = a$depth, period = 2^a$x)
  out$sedrate <- out$period / astro_cycle * 100
  out
}
