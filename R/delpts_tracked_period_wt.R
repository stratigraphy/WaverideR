#' @title Remove tracking points which were tracked in a wavelet spectra
#'
#' @description Interactively select points for deletion
#'With the  \code{\link{track_period_wavelet}} function it is possible to track points in a wavelet spectra,
#'however errors can be made and as such it is possible to delete these points with the \code{\link{delpts_tracked_period_wt}} function.
#'This function allows one to select points for deletion.
#'#'
#' @param tracking_pts Points tracked using the \code{\link{track_period_wavelet}} function.
#' @param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#' @param n.levels Number of color levels \code{Default=100}.
#' @param periodlab label for the y-axis \code{Default="Period (metres)"}.
#' @param x_lab label for the x-axis \code{Default="depth (metres)"}.
#'
#'@return The results of the deletion of the tracking points is a matrix with 3 columns.
#'The first column is depth/time
#'The second column is the period of the tracked cycle
#'The third column is the sedimentation rate based on the duration (in time) of the tracked cycle
#'
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
#'#mag_track <- track_period_wavelet(astro_cycle = 405,
#'#                                   wavelet=mag_wt,
#'#                                   n.levels = 100,
#'#                                   periodlab = "Period (metres)",
#'#                                   x_lab = "depth (metres)")
#'
#'#load the mag_track_solution data set to get an example data set from which
#'#data points can be deleted
#'
#'
#' mag_track_corr <- delpts_tracked_period_wt(tracking_pts = mag_track_solution,
#'                                     wavelet = mag_wt,
#'                                     n.levels = 100,
#'                                     periodlab = "Period (metres)",
#'                                     x_lab = "depth (metres)")
#'}
#'
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
#' @importFrom stats na.omit


delpts_tracked_period_wt <- function(tracking_pts = NULL,
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


  points(
    x = tracking_pts[, 1],
    y = log2(tracking_pts[, 2]),
    type = "p",
    pch = 19,
    col = "black",
    lwd = "0.5"
  )

  x <-  tracking_pts[, 1]
  y <-  log2(tracking_pts[, 2])

  defaultW <- getOption("warn")
  options(warn = -1)
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y

  sel <- rep(FALSE, length(x))
  res <- integer(0)
  n <- length(x)

  while (sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n = 1, plot = F,tolerance = 0.15)
    if (!length(ans))
      break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = 19, col = "white")
    sel[ans] <- TRUE
  }

  out <- tracking_pts[!sel,]


  return(out)

}
