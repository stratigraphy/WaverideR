#' @title Extract a signal in between tracked boundaries in a wavelet scalogram
#'
#' @description Interactively select points in a wavelet scalogram to trace the upper and
#' lower period of an cycle.  The \code{\link{dynamic_extraction}} function plots a wavelet scalogram in which points peaks can selected
#'allowing one to track the lower and upper period of a cycle. First track the upper or lower period of the to
#'be extracted cycle and then track the other boundary. Tracking points can be selected in the Interactive interface and will be shown as white dots
#'connected by a black line. When one wants to deselect a point the white dots can be re-clicked/re-selected and will turn red which
#'indicates that the previously selected point is deselected. Deselecting points can be quite tricky.
#'After tracking the lower and upper boundaries of the cycle the \code{\link{dynamic_extraction}} function
#'will extract the signal in between the boundaries. the output can then used as input for the
#'\code{\link{minimal_tuning}} function to create an age model.
#' @param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#' @param n.levels Number of color levels \code{Default=100}.
#' @param add_peaks Setting which indicates whether spectral peaks should be
#' added to the tracking plot  \code{Default=FALSE}.
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
#'@param smooth smooth the tracked period using the "loess_auto" function
#'@param add_mean add the mean to the extracted signal
#'
#'@return Results of the tracking of a cycle in the wavelet spectra is a matrix with 3 columns.
#'The first column is depth/time
#'The second column is the extracted tracked cycle
#'The third column is upper tracked period
#'The fourth column is lower tracked period
#' @author
#' The function is based/inspired on the \link[astrochron]{traceFreq}
#'function of the 'astrochron' R package
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#'@examples
#'\dontrun{
#'#Track the 405kyr upper and lower periods of the eccentricity cycle in the
#'#magnetic susceptibility record of the Sullivan core of Pas et al., (2018)
#'
#'mag_wt <- analyze_wavelet(
#'  data = mag,
#'  dj = 1 / 100,
#'  lowerPeriod = 0.1,
#'  upperPeriod = 254,
#'  verbose = FALSE,
#'  omega_nr = 10
#')
#'
#'mag_ext <- dynamic_extraction(
#'  wavelet = mag_wt,
#'  n.levels = 100,
#'  add_peaks = FALSE,
#'  periodlab = "Period (metres)",
#'  x_lab = "depth (metres)",
#'  palette_name = "rainbow",
#'  color_brewer = "grDevices",
#'  plot_horizontal = TRUE,
#'  smooth = TRUE,
#'  add_mean = TRUE
#')
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
#' @importFrom DescTools Closest
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats rowMins

dynamic_extraction <- function(wavelet = NULL,
                               n.levels = 100,
                               add_peaks = FALSE,
                               periodlab = "Period (metres)",
                               x_lab = "depth (metres)",
                               palette_name = "rainbow",
                               color_brewer = "grDevices",
                               plot_horizontal = TRUE,
                               smooth = FALSE,
                               add_mean = TRUE) {
  plot.COI = TRUE

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

  useRaster = TRUE
  plot.legend = TRUE
  exponent = 1
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
  axis.1 <- wavelet$axis.1
  axis.2 <- wavelet$axis.2
  Power = wavelet$Power ^ exponent
  wavelet.levels = quantile(Power, probs = seq(
    from = 0,
    to = 1,
    length.out = n.levels + 1
  ))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  image.plt = par()$plt
  legend.plt = NULL


  if (plot_horizontal == TRUE) {
    dev.new(width = 15,
            height = 7,
            noRStudioGD = TRUE)

    if (plot.legend == T) {
      legend.plt = par()$plt
      char.size = par()$cin[1] / par()$din[1]
      hoffset = char.size * par()$mar[4]
      legend.width = char.size * legend.params$width
      legend.mar = char.size * legend.params$mar
      legend.plt[2] = 1 - legend.mar
      legend.plt[1] = legend.plt[2] - legend.width
      vmar = (legend.plt[4] - legend.plt[3]) * ((1 - legend.params$shrink) /
                                                  2)
      legend.plt[4] = legend.plt[4] - vmar
      legend.plt[3] = legend.plt[3] + vmar
      image.plt[2] = min(image.plt[2], legend.plt[1] - hoffset)
      par(plt = legend.plt)
      key.marks = round(seq(
        from = 0,
        to = 1,
        length.out = legend.params$n.ticks
      ) *
        n.levels)
      key.labels = formatC(
        as.numeric(wavelet.levels),
        digits = legend.params$label.digits,
        format = legend.params$label.format
      )[key.marks +
          1]
      image(
        1,
        seq(from = 0, to = n.levels),
        matrix(wavelet.levels,
               nrow = 1),
        col = key.cols,
        breaks = wavelet.levels,
        useRaster = useRaster,
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
          legend.params$width - 0.04
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
      text(
        x = par()$usr[2] + (1.5 + legend.params$lab.line) *
          par()$cxy[1],
        y = n.levels / 2,
        labels = legend.params$lab,
        xpd = NA,
        srt = 270,
        font = par()$font.lab,
        cex = par()$cex.lab
      )
      box(lwd = lwd.axis)
      par(new = TRUE, plt = image.plt)
    }

    par(mar = c(4, 4, 3, 5))

    image(
      x = wavelet$x,
      y = axis.2,
      z = t(Power),
      col = key.cols,
      breaks = wavelet.levels,
      useRaster = TRUE,
      ylab = periodlab,
      xlab = x_lab,
      axes = TRUE,
      yaxt = "n" ,
      main = main
    )


    if (plot.COI == T) {
      polygon(wavelet$coi.1 ,
              wavelet$coi.2,
              border = NA,
              col = rgb(1, 1, 1, 0.5))
    }



    box(lwd = lwd.axis)
    period.tick = unique(trunc(axis.2))
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
      las = 1,
      line = par()$mgp[2] - 0.5,
      font = par()$font.axis,
      cex = par()$cex.axis
    )


    if (add_peaks == TRUE) {
      Pwert <- wavelet$Power

      for (j in 1:ncol(Pwert)) {
        data  <- cbind(log2(wavelet$Period), Pwert[, j])
        data_dif <- data[1:(nrow(data) - 1), 2] - data[2:(nrow(data)), 2]
        data_dif_v2 <-
          data_dif[1:(length(data_dif) - 1)] - data_dif[2:(length(data_dif))]
        Pwert[, j] <-
          c(data_dif_v2, data_dif_v2[length(data_dif_v2)], data_dif_v2[length(data_dif_v2)]) *
          -1
      }

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
      maxdetect2[, 1] <- period
      maxdetect2[, 2] <- depth
      maxdetect2 <- maxdetect2[maxdetect2$value > 0, ]

      colnames(maxdetect2) <- c("y_val", "x_val", "ridge")

    }

    if (add_peaks == TRUE) {
      points(
        x = maxdetect2$x_val,
        y = maxdetect2$y_val,
        type = "p",
        pch = 1,
        col = "black",
        lwd = "0.5"
      )
    }




    x  <- rep(wavelet$x, each = length(wavelet$axis.2))
    y  <- rep(wavelet$axis.2, times = length(wavelet$x))
    n <- length(wavelet$x)


    defaultW <- getOption("warn")
    options(warn = -1)
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    sel <- cbind(rep(FALSE, length(x)), rep(FALSE, length(x)))


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

      image(
        x = wavelet$x,
        y = axis.2,
        z = t(Power),
        col = key.cols,
        breaks = wavelet.levels,
        useRaster = TRUE,
        ylab = periodlab,
        xlab = x_lab,
        axes = TRUE,
        yaxt = "n" ,
        main = main
      )


      if (plot.COI == T) {
        polygon(
          wavelet$coi.1 ,
          wavelet$coi.2,
          border = NA,
          col = rgb(1, 1, 1, 0.5)
        )
      }



      box(lwd = lwd.axis)
      period.tick = unique(trunc(axis.2))
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
        las = 1,
        line = par()$mgp[2] - 0.5,
        font = par()$font.axis,
        cex = par()$cex.axis
      )

      if (add_peaks == TRUE) {
        points(
          x = maxdetect2$x_val,
          y = maxdetect2$y_val,
          type = "p",
          pch = 1,
          col = "black",
          lwd = "0.5"
        )
      }

      points(x[sel[, 1]], y[sel[, 1]], pch = 19, col = "white")
      loc_sort <- data.frame(x[sel[, 1]], y[sel[, 1]])
      lines(loc_sort[order(loc_sort[, 1]), ], col = "black")
      points(x[sel[, 2]], y[sel[, 2]], pch = 19, col = "red")

    }

    out <- data.frame(x[sel[, 1]], y[sel[, 1]])


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
      colnames(out) <- c("depth", "period")
    }


    out_1 <- out

    if (plot.legend == T) {
      legend.plt = par()$plt
      char.size = par()$cin[1] / par()$din[1]
      hoffset = char.size * par()$mar[4]
      legend.width = char.size * legend.params$width
      legend.mar = char.size * legend.params$mar
      legend.plt[2] = 1 - legend.mar
      legend.plt[1] = legend.plt[2] - legend.width
      vmar = (legend.plt[4] - legend.plt[3]) * ((1 - legend.params$shrink) /
                                                  2)
      legend.plt[4] = legend.plt[4] - vmar
      legend.plt[3] = legend.plt[3] + vmar
      image.plt[2] = min(image.plt[2], legend.plt[1] - hoffset)
      par(plt = legend.plt)
      key.marks = round(seq(
        from = 0,
        to = 1,
        length.out = legend.params$n.ticks
      ) *
        n.levels)
      key.labels = formatC(
        as.numeric(wavelet.levels),
        digits = legend.params$label.digits,
        format = legend.params$label.format
      )[key.marks +
          1]
      image(
        1,
        seq(from = 0, to = n.levels),
        matrix(wavelet.levels,
               nrow = 1),
        col = key.cols,
        breaks = wavelet.levels,
        useRaster = useRaster,
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
          legend.params$width - 0.04
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
      text(
        x = par()$usr[2] + (1.5 + legend.params$lab.line) *
          par()$cxy[1],
        y = n.levels / 2,
        labels = legend.params$lab,
        xpd = NA,
        srt = 270,
        font = par()$font.lab,
        cex = par()$cex.lab
      )
      box(lwd = lwd.axis)
      par(new = TRUE, plt = image.plt)
    }

    par(mar = c(4, 4, 3, 5))

    image(
      x = wavelet$x,
      y = axis.2,
      z = t(Power),
      col = key.cols,
      breaks = wavelet.levels,
      useRaster = TRUE,
      ylab = periodlab,
      xlab = x_lab,
      axes = TRUE,
      yaxt = "n" ,
      main = main
    )

    lines(out_1[, 1],
          log2(out_1[, 2]),
          col = "black",
          lwd = 2)


    if (plot.COI == T) {
      polygon(wavelet$coi.1 ,
              wavelet$coi.2,
              border = NA,
              col = rgb(1, 1, 1, 0.5))
    }



    box(lwd = lwd.axis)
    period.tick = unique(trunc(axis.2))
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
      las = 1,
      line = par()$mgp[2] - 0.5,
      font = par()$font.axis,
      cex = par()$cex.axis
    )


    if (add_peaks == TRUE) {
      points(
        x = maxdetect2$x_val,
        y = maxdetect2$y_val,
        type = "p",
        pch = 1,
        col = "black",
        lwd = "0.5"
      )
    }

    x  <- rep(wavelet$x, each = length(wavelet$axis.2))
    y  <- rep(wavelet$axis.2, times = length(wavelet$x))
    n <- length(wavelet$x)


    defaultW <- getOption("warn")
    options(warn = -1)
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    sel <- cbind(rep(FALSE, length(x)), rep(FALSE, length(x)))


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

      image(
        x = wavelet$x,
        y = axis.2,
        z = t(Power),
        col = key.cols,
        breaks = wavelet.levels,
        useRaster = TRUE,
        ylab = periodlab,
        xlab = x_lab,
        axes = TRUE,
        yaxt = "n" ,
        main = main
      )

      lines(out_1[, 1],
            log2(out_1[, 2]),
            col = "black",
            lwd = 2)

      if (plot.COI == T) {
        polygon(
          wavelet$coi.1 ,
          wavelet$coi.2,
          border = NA,
          col = rgb(1, 1, 1, 0.5)
        )
      }



      box(lwd = lwd.axis)
      period.tick = unique(trunc(axis.2))
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
        las = 1,
        line = par()$mgp[2] - 0.5,
        font = par()$font.axis,
        cex = par()$cex.axis
      )

      if (add_peaks == TRUE) {
        points(
          x = maxdetect2$x_val,
          y = maxdetect2$y_val,
          type = "p",
          pch = 1,
          col = "black",
          lwd = "0.5"
        )
      }

      points(x[sel[, 1]], y[sel[, 1]], pch = 19, col = "white")
      loc_sort <- data.frame(x[sel[, 1]], y[sel[, 1]])
      lines(loc_sort[order(loc_sort[, 1]), ], col = "black")
      points(x[sel[, 2]], y[sel[, 2]], pch = 19, col = "red")

    }

    out <- data.frame(x[sel[, 1]], y[sel[, 1]])


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
      colnames(out) <- c("depth", "period")
    }




    res_list <- list(out, out_1)
  }

  if (plot_horizontal == FALSE) {
    dev.new(width = 7,
            height = 10,
            noRStudioGD = TRUE)

    if (plot.legend == T) {
      legend.plt = par()$plt
      char.size = par()$cin[1] / par()$din[1]
      hoffset = char.size * par()$mar[4]
      legend.width = char.size * legend.params$width
      legend.mar = char.size * legend.params$mar
      legend.plt[2] = 1 - legend.mar
      legend.plt[1] = legend.plt[2] - legend.width
      vmar = (legend.plt[4] - legend.plt[3]) * ((1 - legend.params$shrink) /
                                                  2)
      legend.plt[4] = legend.plt[4] - vmar
      legend.plt[3] = legend.plt[3] + vmar
      image.plt[2] = min(image.plt[2], legend.plt[1] - hoffset)
      par(plt = legend.plt)
      key.marks = round(seq(
        from = 0,
        to = 1,
        length.out = legend.params$n.ticks
      ) *
        n.levels)
      key.labels = formatC(
        as.numeric(wavelet.levels),
        digits = legend.params$label.digits,
        format = legend.params$label.format
      )[key.marks +
          1]
      image(
        1,
        seq(from = 0, to = n.levels),
        matrix(wavelet.levels,
               nrow = 1),
        col = key.cols,
        breaks = wavelet.levels,
        useRaster = useRaster,
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
          legend.params$width - 0.04
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
      text(
        x = par()$usr[2] + (1.5 + legend.params$lab.line) *
          par()$cxy[1],
        y = n.levels / 2,
        labels = legend.params$lab,
        xpd = NA,
        srt = 270,
        font = par()$font.lab,
        cex = par()$cex.lab
      )
      box(lwd = lwd.axis)
      par(new = TRUE, plt = image.plt)
    }

    par(mar = c(4, 4, 3, 5))

    image(
      y = wavelet$x,
      x = axis.2,
      z = (Power),
      col = key.cols,
      breaks = wavelet.levels,
      useRaster = TRUE,
      xlab = periodlab,
      ylab = x_lab,
      axes = TRUE,
      xaxt = "n" ,
      main = main
    )


    if (plot.COI == T) {
      polygon(wavelet$coi.2 ,
              wavelet$coi.1,
              border = NA,
              col = rgb(1, 1, 1, 0.5))
    }



    box(lwd = lwd.axis)
    period.tick = unique(trunc(axis.2))
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
      las = 1,
      line = par()$mgp[2] - 0.5,
      font = par()$font.axis,
      cex = par()$cex.axis
    )



    if (add_peaks == TRUE) {
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

    }

    if (add_peaks == TRUE) {
      points(
        x = maxdetect2$x_val,
        y = maxdetect2$y_val,
        type = "p",
        pch = 1,
        col = "black",
        lwd = "0.5"
      )
    }




    y  <- rep(wavelet$x, each = length(wavelet$axis.2))
    x  <- rep(wavelet$axis.2, times = length(wavelet$x))
    n <- length(wavelet$x)


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


      image(
        y = wavelet$x,
        x = axis.2,
        z = (Power),
        col = key.cols,
        breaks = wavelet.levels,
        useRaster = TRUE,
        xlab = periodlab,
        ylab = x_lab,
        axes = TRUE,
        xaxt = "n" ,
        main = main
      )


      if (plot.COI == T) {
        polygon(
          wavelet$coi.2 ,
          wavelet$coi.1,
          border = NA,
          col = rgb(1, 1, 1, 0.5)
        )
      }



      box(lwd = lwd.axis)
      period.tick = unique(trunc(axis.2))
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
        las = 1,
        line = par()$mgp[2] - 0.5,
        font = par()$font.axis,
        cex = par()$cex.axis
      )


      if (add_peaks == TRUE) {
        points(
          x = maxdetect2$x_val,
          y = maxdetect2$y_val,
          type = "p",
          pch = 1,
          col = "black",
          lwd = "0.5"
        )
      }

      points(x[sel[, 1]], y[sel[, 1]], pch = 19, col = "white")
      loc_sort <- data.frame(x[sel[, 1]], y[sel[, 1]])
      lines(loc_sort[order(loc_sort[, 2]), ], col = "black")
      points(x[sel[, 2]], y[sel[, 2]], pch = 19, col = "red")
    }


    out <- data.frame(y[sel[, 1]], x[sel[, 1]])

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
      colnames(out) <- c("depth", "period")
    }


    out_1 <- out

    if (plot.legend == T) {
      legend.plt = par()$plt
      char.size = par()$cin[1] / par()$din[1]
      hoffset = char.size * par()$mar[4]
      legend.width = char.size * legend.params$width
      legend.mar = char.size * legend.params$mar
      legend.plt[2] = 1 - legend.mar
      legend.plt[1] = legend.plt[2] - legend.width
      vmar = (legend.plt[4] - legend.plt[3]) * ((1 - legend.params$shrink) /
                                                  2)
      legend.plt[4] = legend.plt[4] - vmar
      legend.plt[3] = legend.plt[3] + vmar
      image.plt[2] = min(image.plt[2], legend.plt[1] - hoffset)
      par(plt = legend.plt)
      key.marks = round(seq(
        from = 0,
        to = 1,
        length.out = legend.params$n.ticks
      ) *
        n.levels)
      key.labels = formatC(
        as.numeric(wavelet.levels),
        digits = legend.params$label.digits,
        format = legend.params$label.format
      )[key.marks +
          1]
      image(
        1,
        seq(from = 0, to = n.levels),
        matrix(wavelet.levels,
               nrow = 1),
        col = key.cols,
        breaks = wavelet.levels,
        useRaster = useRaster,
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
          legend.params$width - 0.04
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
      text(
        x = par()$usr[2] + (1.5 + legend.params$lab.line) *
          par()$cxy[1],
        y = n.levels / 2,
        labels = legend.params$lab,
        xpd = NA,
        srt = 270,
        font = par()$font.lab,
        cex = par()$cex.lab
      )
      box(lwd = lwd.axis)
      par(new = TRUE, plt = image.plt)
    }

    par(mar = c(4, 4, 3, 5))

    image(
      y = wavelet$x,
      x = axis.2,
      z = (Power),
      col = key.cols,
      breaks = wavelet.levels,
      useRaster = TRUE,
      xlab = periodlab,
      ylab = x_lab,
      axes = TRUE,
      xaxt = "n" ,
      main = main
    )

    lines(
      y = out_1[, 1],
      x = log2(out_1[, 2]),
      col = "black",
      lwd = 2
    )

    if (plot.COI == T) {
      polygon(wavelet$coi.2 ,
              wavelet$coi.1,
              border = NA,
              col = rgb(1, 1, 1, 0.5))
    }



    box(lwd = lwd.axis)
    period.tick = unique(trunc(axis.2))
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
      las = 1,
      line = par()$mgp[2] - 0.5,
      font = par()$font.axis,
      cex = par()$cex.axis
    )



    if (add_peaks == TRUE) {
      points(
        x = maxdetect2$x_val,
        y = maxdetect2$y_val,
        type = "p",
        pch = 1,
        col = "black",
        lwd = "0.5"
      )
    }




    y  <- rep(wavelet$x, each = length(wavelet$axis.2))
    x  <- rep(wavelet$axis.2, times = length(wavelet$x))
    n <- length(wavelet$x)


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


      image(
        y = wavelet$x,
        x = axis.2,
        z = (Power),
        col = key.cols,
        breaks = wavelet.levels,
        useRaster = TRUE,
        xlab = periodlab,
        ylab = x_lab,
        axes = TRUE,
        xaxt = "n" ,
        main = main
      )

      lines(
        y = out_1[, 1],
        x = log2(out_1[, 2]),
        col = "black",
        lwd = 2
      )

      if (plot.COI == T) {
        polygon(
          wavelet$coi.2 ,
          wavelet$coi.1,
          border = NA,
          col = rgb(1, 1, 1, 0.5)
        )
      }



      box(lwd = lwd.axis)
      period.tick = unique(trunc(axis.2))
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
        las = 1,
        line = par()$mgp[2] - 0.5,
        font = par()$font.axis,
        cex = par()$cex.axis
      )


      if (add_peaks == TRUE) {
        points(
          x = maxdetect2$x_val,
          y = maxdetect2$y_val,
          type = "p",
          pch = 1,
          col = "black",
          lwd = "0.5"
        )
      }

      points(x[sel[, 1]], y[sel[, 1]], pch = 19, col = "white")
      loc_sort <- data.frame(x[sel[, 1]], y[sel[, 1]])
      lines(loc_sort[order(loc_sort[, 2]), ], col = "black")
      points(x[sel[, 2]], y[sel[, 2]], pch = 19, col = "red")
    }


    out <- data.frame(y[sel[, 1]], x[sel[, 1]])

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
      colnames(out) <- c("depth", "period")
    }



    res_list <- list(out, out_1)
  }

  res_1 <- res_list[[1]]
  x_axis <- wavelet$x

  yleft_comp <- res_1[1, 2]
  yright_com <- res_1[nrow(res_1), 2]

  seq <-
    seq(
      from = min(x_axis),
      to = max(x_axis),
      by = abs(x_axis[2] - x_axis[1])
    )
  app <-
    approx(
      x = res_1[, 1],
      y = res_1[, 2],
      xout = seq,
      method = "linear",
      yleft = yleft_comp,
      yright = yright_com
    )
  completed_series_1 <- as.data.frame(cbind(app$x, app$y))

  if (smooth == TRUE) {
    completed_series_1 <- loess_auto(
      time_series = completed_series_1,
      genplot = FALSE,
      print_span = FALSE,
      keep_editable = FALSE
    )
  }



  res_2 <- res_list[[2]]
  yleft_comp <- res_2[1, 2]
  yright_com <- res_2[nrow(res_2), 2]

  app <-
    approx(
      x = res_2[, 1],
      y = res_2[, 2],
      xout = seq,
      method = "linear",
      yleft = yleft_comp,
      yright = yright_com
    )
  completed_series_2 <- as.data.frame(cbind(app$x, app$y))

  if (smooth == TRUE) {
    completed_series_2 <- loess_auto(
      time_series = completed_series_2,
      genplot = FALSE,
      print_span = FALSE,
      keep_editable = FALSE
    )
  }



  completed_series <-
    cbind(completed_series_2[, c(1, 2)], completed_series_1[, 2])


  my.w <- wavelet
  my.data <- cbind(wavelet$x, wavelet$y)
  filtered_cycle <- my.data[, 1]
  filtered_cycle <- as.data.frame(filtered_cycle)
  filtered_cycle$value <- NA


  Wave = my.w$Wave
  Power = my.w$Power

  nc = my.w$nc
  nr = my.w$nr
  dt = my.w$dt
  dj = my.w$dj

  Scale = my.w$Scale
  Period = my.w$Period
  loess.span = my.w$loess.span
  rec.waves = matrix(0, nrow = nr, ncol = nc)


  for (s.ind in seq_len(nr)) {
    rec.waves[s.ind,] = (Re(Wave[s.ind,]) / sqrt(Scale[s.ind])) *
      dj * sqrt(dt) / (pi ^ (-1 / 4))
  }


  interpolated <- as.data.frame(completed_series)
  interpolated[, 2] <- rowMaxs(as.matrix(completed_series[, c(2, 3)]))
  interpolated[, 3] <- rowMins(as.matrix(completed_series[, c(2, 3)]))

  for (i in 1:nrow(filtered_cycle)) {
    row_nr_high <- Closest(Period[], interpolated[i, 2], which = TRUE)
    row_nr_low <- Closest(Period[], interpolated[i, 3], which = TRUE)

    row_nr_high <- row_nr_high[1]
    row_nr_low <- row_nr_low[1]

    value <- rec.waves[c(row_nr_low:row_nr_high), i]
    value  <- sum(value, na.rm = T)
    value <- as.numeric(value)
    filtered_cycle[i, 2] <- value
  }

  rec_value  <- colSums(rec.waves, na.rm = T)

  filtered_cycle[, 2] <-
    (filtered_cycle[, 2]) * sd(my.data[, 2]) / sd(rec_value)

  if (add_mean == TRUE) {
    filtered_cycle[, 2] <- filtered_cycle[, 2] + mean(my.data[, 2])
  }


  return(cbind(filtered_cycle, interpolated[, c(2, 3)]))
}





