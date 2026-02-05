#' @title Track the period of a cycle in a wavelet or superlet scalogram
#'
#' @description
#' Interactively select points in a wavelet or superlet power spectrum
#' to trace the evolution of a cycle with changing period.
#'
#' The \code{track_period} function plots a time-frequency spectrum
#' in which spectral peaks can be selected to track a ridge through time
#' or depth. This allows the user to follow a cycle whose period varies
#' along the record.
#'
#' Tracking points are selected interactively and displayed as white dots.
#' Previously selected points can be deselected by clicking them again,
#' after which they are shown as red dots. Because points may be closely
#' spaced, de-selection can be difficult. In such cases,
#' \code{\link{delpts_tracked_period_wt}} can be used to remove points that
#' were previously selected.
#'
#' @param scalogram
#' A wavelet or superlet object created using
#' \code{\link{analyze_wavelet}} or \code{\link{analyze_superlet}}.
#'
#' @param astro_cycle
#' Duration (in kyr) of the astronomical cycle that is being tracked.
#'
#' @param n.levels
#' Number of colour levels used for plotting.
#' Default is \code{100}.
#'
#' @param track_peaks
#' Logical flag indicating whether tracking is restricted to spectral
#' peaks (\code{TRUE}) or whether any point within the spectrum can be
#' selected (\code{FALSE}). Default is \code{TRUE}.
#'
#' @param periodlab
#' Label for the period axis.
#' Default is \code{"Period (metres)"}.
#'
#' @param x_lab
#' Label for the x-axis.
#' Default is \code{"depth (metres)"}.
#'
#' @param palette_name
#' Name of the colour palette used for plotting.
#'
#' @param color_brewer
#' Name of the R package from which the colour palette is selected.
#' Supported packages are \code{RColorBrewer}, \code{grDevices},
#' \code{ColorRamps}, and \code{viridis}.
#' Default is \code{"grDevices"}.
#'
#' @param plot_horizontal
#' Logical flag indicating whether the spectrum is plotted horizontally
#' or vertically. Default is \code{TRUE}.
#'
#' @param lowerPeriod
#' Lowest period value to be displayed.
#'
#' @param upperPeriod
#' Highest period value to be displayed.
#'
#' @param plot_dir
#' Logical flag defining the direction of the record.
#' If time increases with increasing depth (e.g. borehole data),
#' set to \code{TRUE}. If time decreases with increasing depth,
#' set to \code{FALSE}. Default is \code{TRUE}.
#'
#' @param add_lines
#' Optional matrix of additional lines to overlay on the spectrum.
#' The first column must be depth or time, and subsequent columns
#' contain period values.
#'
#' @param add_points
#' Optional matrix of additional points to overlay on the spectrum.
#' The first column must be depth or time, and subsequent columns
#' contain period values.
#'
#' @param add_abline_h
#' Optional numeric vector specifying horizontal reference lines.
#'
#' @param add_abline_v
#' Optional numeric vector specifying vertical reference lines.
#'
#' @return
#' A data frame with three columns:
#'   first - depth - Depth or time of the tracked points
#'   second -  period - Tracked period of the cycle
#'   third  - sedrate - Estimated sedimentation rate based on the cycle duration
#'
#' @author
#' The function is based on and inspired by the
#'  \link[astrochron]{traceFreq} function from the astrochron package.
#'
#'@references
#' Routines for astrochronologic testing, astronomical time-scale construction,
#' and time series analysis.
#' \doi{10.1016/j.earscirev.2018.11.015}
#'
#' @examples
#' \donttest{
#' ## Track the 405 kyr eccentricity cycle in a magnetic susceptibility record
#' mag_wt <- analyze_wavelet(
#'   data = mag,
#'   dj = 1/100,
#'   lowerPeriod = 0.1,
#'   upperPeriod = 254,
#'   verbose = FALSE,
#'   omega_nr = 10
#' )
#'
#' mag_track <- track_period(
#'   scalogram = mag_wt,
#'   astro_cycle = 405,
#'   track_peaks = TRUE
#' )
#' }
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

track_period <- function (scalogram = NULL, astro_cycle = 405, n.levels = 100,
                          track_peaks = TRUE, periodlab = "Period (metres)", x_lab = "depth (metres)",
                          palette_name = "rainbow", color_brewer = "grDevices", plot_horizontal = TRUE,
                          plot_dir = TRUE, lowerPeriod = NULL, upperPeriod = NULL,
                          add_lines = NULL, add_points = NULL, add_abline_h = NULL,
                          add_abline_v = NULL)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  maximum.level = max(scalogram$Power)
  power_max_mat.levels = quantile(scalogram$Power, probs = seq(from = 0,
                                                               to = 1, length.out = n.levels + 1))
  if (color_brewer == "RColorBrewer") {
    key.cols <- rev(colorRampPalette(brewer.pal(brewer.pal.info[palette_name,
                                                                1], palette_name))(n.levels))
  }
  if (color_brewer == "colorRamps") {
    color_brewer_Sel <- paste("colorRamps::", palette_name,
                              "(n=n.levels)")
    key.cols = eval(parse(text = color_brewer_Sel))
  }
  if (color_brewer == "grDevices") {
    if (palette_name == "rainbow") {
      color_brewer_Sel <- "grDevices::rainbow(n=n.levels, start = 0, end = 0.7)"
      key.cols <- rev(eval(parse(text = color_brewer_Sel)))
    }
    else if (palette_name == "heat.colors" | palette_name ==
             "terrain.colors" | palette_name == "topo.colors" |
             palette_name == "cm.colors") {
      color_brewer_Sel <- paste("grDevices::", palette_name,
                                "(n=n.levels, start = 0, end = 1)")
      key.cols <- rev(eval(parse(text = color_brewer_Sel)))
    }
    else {
      key.cols <- hcl.colors(n = n.levels, palette = palette_name,
                             alpha = NULL, rev = FALSE, fixup = TRUE)
    }
  }
  if (color_brewer == "viridis") {
    color_brewer_Sel <- paste("viridis::", palette_name,
                              "(n=n.levels,direction = -1)")
    key.cols = rev(eval(parse(text = color_brewer_Sel)))
  }
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
  plot_horizontal <- TRUE
  y_axis <- as.numeric(unlist(scalogram$Period))
  pmax_avg_sel <- t(scalogram$Power)
  depth <- scalogram$x
  y_axis <- scalogram$Period
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)
  if (plot_dir != TRUE) {
    xlim_vals = rev(c(min(scalogram$x), max(scalogram$x)))
  }  else {
    xlim_vals = c(min(scalogram$x), max(scalogram$x))
  }
  if (is.null(lowerPeriod) == TRUE) {
    lowerPeriod <- min(scalogram$Period)
  }
  if (is.null(upperPeriod) == TRUE) {
    upperPeriod <- max(scalogram$Period)
  }
  ylim_vals = c(lowerPeriod, upperPeriod)
  if (plot_horizontal == TRUE) {
    dev.new(width = 15, height = 7, noRStudioGD = TRUE)
    layout.matrix <- matrix(c(1, 2, 4, 3), nrow = 2, ncol = 2,
                            byrow = TRUE)
    graphics::layout(mat = layout.matrix, heights = c(0.25,
                                                      1), widths = c(8, 2))
    power_max_mat.levels = quantile(pmax_avg_sel, probs = seq(from = 0,
                                                              to = 1, length.out = n.levels + 1))
    par(mar = c(0, 4, 2, 0))
    plot(y = scalogram$y, x = scalogram$x, type = "l", xaxs = "i",
         xlab = "", ylab = "proxy value", xaxt = "n", xlim = xlim_vals)
    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }
    par(new = FALSE, mar = c(3, 2, 2, 2), mgp = c(2, 1,
                                                  0))
    image(x = seq(from = 0, to = n.levels), y = 1, z = t(matrix(power_max_mat.levels,
                                                                nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
          useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "Power",
          ylab = "")
    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.2)
    mtext(key.labels, side = 1, at = key.marks, line = 0.1,
          las = 2, cex = 0.75)
    box(lwd = lwd.axis)
    par(new = FALSE, mar = c(4, 0, 0, 0.5))
    plot(x = scalogram$Power.avg, y = scalogram$Period, log = "y",
         type = "l", yaxs = "i", yaxt = "n", xlab = "Wt. power",
         xaxs = "i", ylim = ylim_vals)
    if (is.null(add_abline_h) != TRUE) {
      abline(h = add_abline_h)
    }
    par(new = FALSE, mar = c(4, 4, 0, 0), xpd = FALSE)
    if (inherits(scalogram,"analyze.wavelet") == TRUE ){
      image(x = scalogram$x, y = scalogram$axis.2, z = t(scalogram$Power),
            col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
            ylab = periodlab, xlab = x_lab, axes = TRUE, yaxt = "n",
            main = main, xlim = xlim_vals, ylim = log2(ylim_vals))
      polygon(x = scalogram$coi.1, y = scalogram$coi.2, border = NA,
              col = rgb(1, 1, 1, 0.5), ylim = xlim_vals, xlim = log2(ylim_vals))
      box(lwd = lwd.axis)
      period.tick = unique(trunc(scalogram$axis.2))
      period.tick[period.tick < log2(scalogram$Period[1])] = NA
      period.tick = na.omit(period.tick)
      period.tick.label = 2^(period.tick)
      axis(2, lwd = lwd.axis, at = period.tick, labels = NA,
           tck = periodtck, tcl = periodtcl)
      axis(4, lwd = lwd.axis, at = period.tick, labels = NA,
           tck = periodtck, tcl = periodtcl)
      mtext(period.tick.label, side = 2, at = period.tick,
            las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
            cex = par()$cex.axis)
      if (is.null(add_lines) != TRUE) {
        for (i in 2:ncol(add_lines)) lines(add_lines[, 1],
                                           log2(add_lines[, i]))
      }
      if (is.null(add_points) != TRUE) {
        for (i in 2:ncol(add_points)) points(add_points[,
                                                        1], log2(add_points[, i]))
      }
      if (is.null(add_abline_h) != TRUE) {
        abline(h = log2(add_abline_h))
      }
      if (is.null(add_abline_v) != TRUE) {
        abline(v = add_abline_v)
      }
      if (track_peaks == TRUE) {
        Pwert <- scalogram$Power
        maxdetect <- matrix(nrow = (nrow(Pwert)), ncol = ncol(Pwert),
                            0)
        for (j in 1:ncol(Pwert)) {
          for (i in 2:(nrow(maxdetect) - 1)) {
            if ((Pwert[i, j] - Pwert[(i + 1), j] > 0) &
                (Pwert[i, j] - Pwert[(i - 1), j] > 0)) {
              maxdetect[i, j] <- 1
            }
          }
        }
        maxdetect2 <- melt(maxdetect)
        depth <- rep(scalogram$x, each = length(scalogram$axis.2))
        period <- rep(scalogram$axis.2, times = length(scalogram$x))
        maxdetect2 <- as.data.frame(maxdetect2)
        maxdetect2[, 2] <- period
        maxdetect2[, 1] <- depth
        maxdetect2 <- maxdetect2[maxdetect2$value > 0, ]
        colnames(maxdetect2) <- c("y_val", "x_val", "ridge")
        points(y = maxdetect2$x_val, x = maxdetect2$y_val,
               type = "p", pch = 1, col = "black", lwd = "0.5")
        n <- nrow(maxdetect2)
        y = maxdetect2$x_val
        x = maxdetect2$y_val
      }
      else {
        x <- rep(scalogram$x, each = length(scalogram$axis.2))
        y <- rep(scalogram$axis.2, times = length(scalogram$x))
        n <- length(scalogram$x)
      }
      defaultW <- getOption("warn")
      options(warn = -1)
      xy <- xy.coords(x, y)
      y <- xy$y
      x <- xy$x
      sel <- cbind(rep(FALSE, length(y)), rep(FALSE, length(y)))
      while (sum(sel) < n) {
        ans <- identify(x, y, n = 1, plot = F, tolerance = 0.1)
        if (!length(ans))
          break
        if (sel[ans, 1] == FALSE) {
          sel[ans, 1] <- TRUE
          sel[ans, 2] <- FALSE
        }
        else {
          sel[ans, 1] <- FALSE
          sel[ans, 2] <- TRUE
        }
        points(x[sel[, 1]], y[sel[, 1]], pch = 19, col = "white")
        points(x[sel[, 2]], y[sel[, 2]], pch = 19, col = "red")
      }
      pts <- sel[, 1]
      if (track_peaks == TRUE) {
        out <- data.frame(cbind(maxdetect2[pts, 1], maxdetect2[pts,
                                                               2]))
        out <- na.omit(out)
      }
      if (track_peaks == FALSE) {
        out <- data.frame(y[sel[, 1]], x[sel[, 1]])
        out <- na.omit(out)
      }
      if (nrow(out) != 0) {
        out <- na.omit(out)
        out <- out[order(out[, 1]), ]
        out <- na.omit(out)
        out <- aggregate(out, by = list(name = out[, 1]),
                         data = out, FUN = mean)
        out <- out[, c(2, 3)]
        out[, 2] <- 2^out[, 2]
        out$sedrate <- out[, 2]/astro_cycle * 100
        colnames(out) <- c("depth", "period", "sedrate")
      }
    }
  }

  if(inherits(scalogram,"analyze.superlet") == TRUE){
    image(x = scalogram$x, y = scalogram$axis.2, z = (scalogram$Power),
          col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
          ylab = periodlab, xlab = x_lab, axes = TRUE, yaxt = "n",
          main = main, xlim = xlim_vals, ylim = log2(ylim_vals))
    box(lwd = lwd.axis)
    period.tick = unique(trunc(scalogram$axis.2))
    period.tick <- period.tick[period.tick >= min(log2(ylim_vals))]
    period.tick <- period.tick[period.tick <= max(log2(ylim_vals))]
    period.tick.label = 2^(period.tick)
    axis(2, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    axis(4, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    mtext(period.tick.label, side = 2, at = period.tick,
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)
    if (is.null(add_lines) != TRUE) {
      for (i in 2:ncol(add_lines)) lines(add_lines[, 1],
                                         log2(add_lines[, i]))
    }
    if (is.null(add_points) != TRUE) {
      for (i in 2:ncol(add_points)) points(add_points[,
                                                      1], log2(add_points[, i]))
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }
    if (track_peaks == TRUE) {
      Pwert <- t(scalogram$Power)
      maxdetect <- matrix(nrow = (nrow(Pwert)), ncol = ncol(Pwert),
                          0)
      for (j in 1:ncol(Pwert)) {
        for (i in 2:(nrow(maxdetect) - 1)) {
          if ((Pwert[i, j] - Pwert[(i + 1), j] > 0) &
              (Pwert[i, j] - Pwert[(i - 1), j] > 0)) {
            maxdetect[i, j] <- 1
          }
        }
      }
      maxdetect2 <- melt(maxdetect)
      depth <- rep(scalogram$x, each = length(scalogram$axis.2))
      period <- rep(scalogram$axis.2, times = length(scalogram$x))
      maxdetect2 <- as.data.frame(maxdetect2)
      maxdetect2[, 2] <- period
      maxdetect2[, 1] <- depth
      maxdetect2 <- maxdetect2[maxdetect2$value > 0, ]
      colnames(maxdetect2) <- c("y_val", "x_val", "ridge")
      points(y = maxdetect2$x_val, x = maxdetect2$y_val,
             type = "p", pch = 1, col = "black", lwd = "0.5")
      n <- nrow(maxdetect2)
      y = maxdetect2$x_val
      x = maxdetect2$y_val
    }else {
      x <- rep(scalogram$x, each = length(scalogram$axis.2))
      y <- rep(scalogram$axis.2, times = length(scalogram$x))
      n <- length(scalogram$x)
    }
    defaultW <- getOption("warn")
    options(warn = -1)
    xy <- xy.coords(x, y)
    y <- xy$y
    x <- xy$x
    sel <- cbind(rep(FALSE, length(y)), rep(FALSE, length(y)))
    while (sum(sel) < n) {
      ans <- identify(x, y, n = 1, plot = F, tolerance = 0.1)
      if (!length(ans))
        break
      if (sel[ans, 1] == FALSE) {
        sel[ans, 1] <- TRUE
        sel[ans, 2] <- FALSE
      }
      else {
        sel[ans, 1] <- FALSE
        sel[ans, 2] <- TRUE
      }
      points(x[sel[, 1]], y[sel[, 1]], pch = 19, col = "white")
      points(x[sel[, 2]], y[sel[, 2]], pch = 19, col = "red")
    }
    pts <- sel[, 1]
    if (track_peaks == TRUE) {
      out <- data.frame(cbind(maxdetect2[pts, 1], maxdetect2[pts,
                                                             2]))
      out <- na.omit(out)
    }
    if (track_peaks == FALSE) {
      out <- data.frame(y[sel[, 1]], x[sel[, 1]])
      out <- na.omit(out)
    }
    if (nrow(out) != 0) {
      out <- na.omit(out)
      out <- out[order(out[, 1]), ]
      out <- na.omit(out)
      out <- aggregate(out, by = list(name = out[, 1]),
                       data = out, FUN = mean)
      out <- out[, c(2, 3)]
      out[, 2] <- 2^out[, 2]
      out$sedrate <- out[, 2]/astro_cycle * 100
      colnames(out) <- c("depth", "period", "sedrate")
    }
  }


  if (plot_horizontal == FALSE) {
    dev.new(width = 7, height = 10, noRStudioGD = TRUE)
    layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2,
                            byrow = TRUE)
    graphics::layout(mat = layout.matrix, heights = c(1,
                                                      4), widths = c(1, 4))
    power_max_mat.levels = quantile(pmax_avg_sel, probs = seq(from = 0,
                                                              to = 1, length.out = n.levels + 1))
    par(mar = c(2, 0.5, 2, 6), xpd = NA)
    image(y = seq(from = 0, to = n.levels), x = 1, z = (matrix(power_max_mat.levels,
                                                               nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
          useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "",
          ylab = "", )
    axis(2, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)
    mtext(key.labels, side = 2, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 1, at = 1, line = 0.5, font = par()$font.axis,
          cex = par()$cex.axis, las = 1, xpd = NA)
    box(lwd = lwd.axis)
    par(mar = c(0, 0, 2, 2))
    plot(y = scalogram$Power.avg, x = scalogram$Period, log = "x",
         type = "l", xaxs = "i", xaxt = "n", ylab = "Wt. power",
         xlab = "", xaxs = "i", xlim = ylim_vals)
    title(ylab = "Wt. power", xpd = NA)
    if (is.null(add_abline_v) != TRUE) {
      abline(v = (add_abline_v),xpd=FALSE)
    }
    par(mar = c(4, 4, 0, 0), xpd = TRUE)
    plot(x = scalogram$y, y = scalogram$x, type = "l", yaxs = "i",
         xlab = "proxy value", ylab = x_lab, ylim = xlim_vals)
    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h),xpd=FALSE)
    }

    par(new = FALSE, mar = c(4, 0, 0, 2), xpd = FALSE)
    if(inherits(scalogram, "analyze.wavelet") == TRUE){

      image(y = scalogram$x, x = scalogram$axis.2, z = (scalogram$Power),
            col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
            xlab = periodlab, ylab = "", xaxt = "n", yaxt = "n",
            main = main, ylim = xlim_vals, xlim = log2(ylim_vals))
      polygon(y = scalogram$coi.1, x = scalogram$coi.2, border = NA,
              col = rgb(1, 1, 1, 0.5), ylim = xlim_vals, xlim = log2(ylim_vals))
      box(lwd = lwd.axis)
      period.tick = unique(trunc(scalogram$axis.2))
      period.tick[period.tick < log2(scalogram$Period[1])] = NA
      period.tick = na.omit(period.tick)
      period.tick.label = 2^(period.tick)
      axis(1, lwd = lwd.axis, at = period.tick, labels = NA,
           tck = periodtck, tcl = periodtcl)
      mtext(period.tick.label, side = 1, at = period.tick,
            las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
            cex = par()$cex.axis)




    }
    if(inherits(scalogram, "analyze.superlet") == TRUE){
      image(
        y = scalogram$x,
        x = scalogram$axis.2,
        z = t(scalogram$Power),
        col = key.cols,
        breaks = power_max_mat.levels,
        useRaster = TRUE,
        xlab = periodlab,
        ylab = x_lab,
        xaxt = "n",
        yaxt="n",
        main = main,
        ylim = xlim_vals,
        xlim = log2(ylim_vals)
      )
      box(lwd = lwd.axis)
      period.tick = unique(trunc(scalogram$axis.2))
      period.tick <- period.tick[period.tick >= min(log2(ylim_vals))]
      period.tick <- period.tick[period.tick <= max(log2(ylim_vals))]
      period.tick.label = 2^(period.tick)
      axis(1, lwd = lwd.axis, at = period.tick, labels = NA,
           tck = periodtck, tcl = periodtcl)
      mtext(period.tick.label, side = 1, at = period.tick,
            las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
            cex = par()$cex.axis)

    }


    if (is.null(add_lines) != TRUE) {
      for (i in 2:ncol(add_lines)) lines(y = add_lines[,
                                                       1], x = log2(add_lines[, i]))
    }
    if (is.null(add_points) != TRUE) {
      for (i in 2:ncol(add_points)) points(y = add_points[,
                                                          1], x = log2(add_points[, i]))
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }
    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }
    if (track_peaks == TRUE) {
      if(inherits(scalogram,"analyze.superlet") == TRUE ){
        Pwert <- t(scalogram$Power)}
      if(inherits(scalogram,"analyze.wavelet") == TRUE){
        Pwert <- scalogram$Power


      }

      maxdetect <- matrix(nrow = (nrow(Pwert)), ncol = ncol(Pwert),
                          0)
      for (j in 1:ncol(Pwert)) {
        for (i in 2:(nrow(maxdetect) - 1)) {
          if ((Pwert[i, j] - Pwert[(i + 1), j] > 0) &
              (Pwert[i, j] - Pwert[(i - 1), j] > 0)) {
            maxdetect[i, j] <- 1
          }
        }
      }
      maxdetect2 <- melt(maxdetect)
      depth <- rep(scalogram$x, each = length(scalogram$axis.2))
      period <- rep(scalogram$axis.2, times = length(scalogram$x))
      maxdetect2 <- as.data.frame(maxdetect2)
      maxdetect2[, 2] <- period
      maxdetect2[, 1] <- depth
      maxdetect2 <- maxdetect2[maxdetect2$value > 0, ]
      colnames(maxdetect2) <- c("y_val", "x_val", "ridge")
      points(x = maxdetect2$x_val, y = maxdetect2$y_val,
             type = "p", pch = 1, col = "black", lwd = "0.5")
      n <- nrow(maxdetect2)
      x = maxdetect2$x_val
      y = maxdetect2$y_val
    }
    else {
      y <- rep(scalogram$x, each = length(scalogram$axis.2))
      x <- rep(scalogram$axis.2, times = length(scalogram$x))
      n <- length(scalogram$x)
    }
    defaultW <- getOption("warn")
    options(warn = -1)
    xy <- xy.coords(x, y)
    y <- xy$y
    x <- xy$x
    sel <- cbind(rep(FALSE, length(y)), rep(FALSE, length(y)))
    while (sum(sel) < n) {
      ans <- identify(x, y, n = 1, plot = F, tolerance = 0.1)
      if (!length(ans))
        break
      if (sel[ans, 1] == FALSE) {
        sel[ans, 1] <- TRUE
        sel[ans, 2] <- FALSE
      }
      else {
        sel[ans, 1] <- FALSE
        sel[ans, 2] <- TRUE
      }
      points(x[sel[, 1]], y[sel[, 1]], pch = 19, col = "white")
      points(x[sel[, 2]], y[sel[, 2]], pch = 19, col = "red")
    }
    pts <- sel[, 1]
    if (track_peaks == TRUE) {
      out <- data.frame(cbind(maxdetect2[pts, 1], maxdetect2[pts,
                                                             2]))
      out <- na.omit(out)
    }
    if (track_peaks == FALSE) {
      out <- data.frame(y[sel[, 1]], x[sel[, 1]])
      out <- na.omit(out)
    }
    if (nrow(out) != 0) {
      out <- na.omit(out)
      out <- out[order(out[, 1]), ]
      out <- na.omit(out)
      out <- aggregate(out, by = list(name = out[, 1]),
                       data = out, FUN = mean)
      out <- out[, c(2, 3)]
      out[, 2] <- 2^out[, 2]
      out$sedrate <- out[, 2]/astro_cycle * 100
      colnames(out) <- c("depth", "period", "sedrate")
    }
  }
  return(out)
}
