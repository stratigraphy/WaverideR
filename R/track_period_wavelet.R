#' @title Track the period of a cycle in a wavelet spectra
#'
#' @description Interactively select points in a wavelet spectra to trace a period in a wavelet spectra.
#'The \code{\link{track_period_wavelet}} function plots a wavelet spectra in which spectral peaks can selected
#'allowing one to track a ridge hence one can track the a cycle with a changing period.
#'
#' @param astro_cycle Duration (in kyr) of the cycle which traced.
#' @param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#' @param n.levels Number of color levels \code{Default=100}.
#' @param periodlab label for the y-axis \code{Default="Period (metres)"}.
#' @param x_lab label for the x-axis \code{Default="depth (metres)"}.
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
#'time series analysis \doi{<doi:10.1016/j.earscirev.2018.11.015>}
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
#'                                   periodlab = "Period (metres)",
#'                                   x_lab = "depth (metres)")
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

track_period_wavelet <- function(astro_cycle = 405,
                                 wavelet = NULL,
                                 n.levels = 100,
                                 periodlab = "Period (metres)",
                                 x_lab = "depth (metres)") {
  plot.COI = TRUE
  color.palette = "rainbow(n.levels, start = 0, end = .7)"
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
  key.cols = rev(eval(parse(text = color.palette)))


  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  image.plt = par()$plt
  legend.plt = NULL
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


  pts <- ID_points(x = maxdetect2$x_val,
                   y = maxdetect2$y_val)



  out <- data.frame(cbind(maxdetect2[pts, 2], maxdetect2[pts, 1]))
  out <- na.omit(out)

  if(nrow(out)!= 0){
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
    colnames(out) <- c("depth", "period", "sedrate")}

  return(out)

}
