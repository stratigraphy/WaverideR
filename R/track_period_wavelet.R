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
#'
#' @param astro_cycle Duration (in kyr) of the cycle which traced.
#' @param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
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
#'
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
#' mag_track <- track_period_wavelet(astro_cycle = 405,
#'                                   wavelet=mag_wt,
#'                                   n.levels = 100,
#'                                   track_peaks=TRUE,
#'                                   periodlab = "Period (metres)",
#'                                   x_lab = "depth (metres)",
#'                                  palette_name = "rainbow",
#'                                  color_brewer="grDevices",
#'                                  plot_horizontal=TRUE)
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


track_period_wavelet <- function(astro_cycle = 405,
                                 wavelet = NULL,
                                 n.levels = 100,
                                 track_peaks=TRUE,
                                 periodlab = "Period (metres)",
                                 x_lab = "depth (metres)",
                                 palette_name = "rainbow",
                                 color_brewer="grDevices",
                                 plot_horizontal=TRUE) {

  plot.COI = TRUE

  if (color_brewer== "RColorBrewer"){
    key.cols <-   rev(colorRampPalette(brewer.pal(brewer.pal.info[palette_name,1],palette_name))(n.levels))

  }


  if (color_brewer== "colorRamps"){
    color_brewer_Sel <- paste("colorRamps::",palette_name,"(n=n.levels)")
    key.cols = eval(parse(text = color_brewer_Sel))
  }


  if (color_brewer == "grDevices"){
    if (palette_name == "rainbow"){
      color_brewer_Sel <- "grDevices::rainbow(n=n.levels, start = 0, end = 0.7)"
      key.cols <- rev(eval(parse(text = color_brewer_Sel)))
    }
    else if (palette_name == "heat.colors"|
             palette_name == "terrain.colors"|
             palette_name == "topo.colors"|
             palette_name == "cm.colors"){
      color_brewer_Sel <- paste("grDevices::",palette_name,"(n=n.levels, start = 0, end = 1)")
      key.cols <- rev(eval(parse(text = color_brewer_Sel)))
    }
    else{
      key.cols <-  hcl.colors(n=n.levels, palette = palette_name, alpha = NULL, rev = FALSE, fixup = TRUE)}}



  if (color_brewer== "viridis"){
    color_brewer_Sel <- paste("viridis::",palette_name,"(n=n.levels,direction = -1)")
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


  if (plot_horizontal == TRUE){
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


    if (track_peaks==TRUE){
      Pwert <- wavelet$Power

      maxdetect <- matrix(nrow = (nrow(Pwert)), ncol = ncol(Pwert), 0)

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
      maxdetect2 <- maxdetect2[maxdetect2$value > 0,]

      colnames(maxdetect2) <- c("y_val", "x_val", "ridge")

      points(
        x = maxdetect2$x_val,
        y = maxdetect2$y_val,
        type = "p",
        pch = 1,
        col = "black",
        lwd = "0.5"
      )


      x = maxdetect2$x_val
      y = maxdetect2$y_val
    }else {
      x  <- rep(wavelet$x, each = length(wavelet$axis.2))
      y  <- rep(wavelet$axis.2, times = length(wavelet$x))}




    defaultW <- getOption("warn")
    options(warn = -1)
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    sel <- cbind(rep(FALSE, length(x)), rep(FALSE, length(x)))
    n <- nrow(maxdetect2)

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

    out <- data.frame(cbind(maxdetect2[pts, 2], maxdetect2[pts, 1]))
    out <- na.omit(out)

    if (nrow(out) != 0) {
      out <- na.omit(out)
      out <- out[order(out[, 1]), ]
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

    return(out)
  }




  if (plot_horizontal == FALSE){
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



    if (track_peaks==TRUE){
      Pwert <- wavelet$Power

      maxdetect <- matrix(nrow = (nrow(Pwert)), ncol = ncol(Pwert), 0)

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
      maxdetect2 <- maxdetect2[maxdetect2$value > 0,]

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
    }else {
      y  <- rep(wavelet$x, each = length(wavelet$axis.2))
      x  <- rep(wavelet$axis.2, times = length(wavelet$x))
      n <- length(wavelet$x)}


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

    out <- data.frame(cbind(maxdetect2[pts, 1], maxdetect2[pts, 2]))
    out <- na.omit(out)

    if (nrow(out) != 0) {
      out <- na.omit(out)
      out <- out[order(out[, 1]), ]
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

    return(out)
  }


}

