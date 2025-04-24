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
#' @param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#' @param astro_cycle Duration (in kyr) of the cycle which traced.
#' @param n.levels Number of color levels \code{Default=100}.
#' @param track_peaks Setting which indicates whether tracking is restricted
#' to spectral peaks (track_peaks=TRUE) or whether any point within the wavelet
#' spectra can be selected (track_peaks=FALSE) \code{Default=TRUE}.
#' @param periodlab label for the y-axis \code{Default="Period (metres)"}.
#' @param x_lab label for the x-axis \code{Default="depth (metres)"}.
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
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis is depth or y axis power  \code{Default=TRUE}
#'@param lowerPeriod Lowest period value which will be plotted
#'@param upperPeriod Highest period value which will be plotted
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
#'#Track the 405kyr eccentricity cycle in the magnetic susceptibility record
#'# of the Sullivan core of Pas et al., (2018)
#'
#'mag_wt <- analyze_wavelet(data = mag,
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
#'}
#'
#' @export
#' @importFrom reshape2 melt
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

track_period_wavelet <- function(wavelet = NULL,
                                 astro_cycle = 405,
                                 n.levels = 100,
                                 track_peaks = TRUE,
                                 periodlab =  "Period (metres)",
                                 x_lab = "depth (metres)",
                                 palette_name = "rainbow",
                                 color_brewer = "grDevices",
                                 plot_horizontal = TRUE,
                                 plot_dir = TRUE,
                                 lowerPeriod = NULL,
                                 upperPeriod = NULL,
                                 add_lines = NULL,
                                 add_points = NULL,
                                 add_abline_h = NULL,
                                 add_abline_v = NULL)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))



  maximum.level = max(wavelet$Power)
  power_max_mat.levels = quantile(wavelet$Power, probs = seq(
    from = 0,
    to = 1,
    length.out = n.levels + 1
  ))


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


  plot_horizontal <- TRUE






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
    dev.new(width = 15,
            height = 7,
            noRStudioGD = TRUE)


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
      for (i  in 2:ncol(add_lines))
        lines(add_lines[, 1], log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        points(add_points[, 1], log2(add_points[, i]))
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

      maxdetect2 <- melt(maxdetect)



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
    } else {
      x  <- rep(wavelet$x, each = length(wavelet$axis.2))
      y  <- rep(wavelet$axis.2, times = length(wavelet$x))
      n <- length(wavelet$x)
    }


    defaultW <- getOption("warn")
    options(warn = -1)
    xy <- xy.coords(x, y)
    y <- xy$y
    x <- xy$x
    sel <- cbind(rep(FALSE, length(y)), rep(FALSE, length(y)))

    while (sum(sel) < n) {
      ans <- identify(x,
                      y,
                      n = 1,
                      plot = F,
                      tolerance = 0.1)

      if (!length(ans))
        break


      if (sel[ans, 1] == FALSE) {
        sel[ans, 1] <- TRUE
        sel[ans, 2] <- FALSE
      } else{
        sel[ans, 1] <- FALSE
        sel[ans, 2] <- TRUE
      }

      points(x[sel[, 1]], y[sel[, 1]], pch = 19, col = "white")
      points(x[sel[, 2]], y[sel[, 2]], pch = 19, col = "red")
    }

    pts <- sel[, 1]
    if (track_peaks == TRUE) {
      out <- data.frame(cbind(maxdetect2[pts, 1], maxdetect2[pts, 2]))
      out <- na.omit(out)
    }


    if (track_peaks == FALSE) {
      out <- data.frame(y[sel[, 1]], x[sel[, 1]])
      out <- na.omit(out)
    }

    if (nrow(out) != 0) {
      out <- na.omit(out)
      out <- out[order(out[, 1]),]
      out <- na.omit(out)
      out <- aggregate(out,
                       by = list(name = out[, 1]),
                       data = out,
                       FUN = mean)
      out <- out[, c(2, 3)]
      out[, 2] <- 2 ^ out[, 2]
      out$sedrate <- out[, 2] / astro_cycle * 100
      colnames(out) <- c("depth", "period", "sedrate")
    }
  }


  if (plot_horizontal == FALSE) {
    dev.new(width = 7,
            height = 10,
            noRStudioGD = TRUE)




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
      for (i  in 2:ncol(add_lines))
        lines(y = add_lines[, 1], x = log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        points(y = add_points[, 1], x = log2(add_points[, i]))
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


      maxdetect2 <- melt(maxdetect)



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
      y = maxdetect2$y_val
    } else {
      y  <- rep(wavelet$x, each = length(wavelet$axis.2))
      x  <- rep(wavelet$axis.2, times = length(wavelet$x))
      n <- length(wavelet$x)
    }


    defaultW <- getOption("warn")
    options(warn = -1)
    xy <- xy.coords(x, y)
    y <- xy$y
    x <- xy$x
    sel <- cbind(rep(FALSE, length(y)), rep(FALSE, length(y)))

    while (sum(sel) < n) {
      ans <- identify(x,
                      y,
                      n = 1,
                      plot = F,
                      tolerance = 0.1)

      if (!length(ans))
        break


      if (sel[ans, 1] == FALSE) {
        sel[ans, 1] <- TRUE
        sel[ans, 2] <- FALSE
      } else{
        sel[ans, 1] <- FALSE
        sel[ans, 2] <- TRUE
      }

      points(x[sel[, 1]], y[sel[, 1]], pch = 19, col = "white")
      points(x[sel[, 2]], y[sel[, 2]], pch = 19, col = "red")
    }

    pts <- sel[, 1]
    if (track_peaks == TRUE) {
      out <- data.frame(cbind(maxdetect2[pts, 1], maxdetect2[pts, 2]))
      out <- na.omit(out)
    }


    if (track_peaks == FALSE) {
      out <- data.frame(y[sel[, 1]], x[sel[, 1]])
      out <- na.omit(out)
    }

    if (nrow(out) != 0) {
      out <- na.omit(out)
      out <- out[order(out[, 1]),]
      out <- na.omit(out)
      out <- aggregate(out,
                       by = list(name = out[, 1]),
                       data = out,
                       FUN = mean)
      out <- out[, c(2, 3)]
      out[, 2] <- 2 ^ out[, 2]
      out$sedrate <- out[, 2] / astro_cycle * 100
      colnames(out) <- c("depth", "period", "sedrate")
    }
  }





  return(out)
}
