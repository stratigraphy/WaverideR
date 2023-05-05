#' @title Plots a wavelet power spectra
#'
#' @description Plot wavelet spectra using the outcome of the \code{\link{analyze_wavelet}} function.
#'
#' @param wavelet wavelet object created using the \code{\link{analyze_wavelet}} function.
#' @param plot.COI Option to plot the cone of influence \code{Default=TRUE}.
#' @param n.levels  Number of color levels \code{Default=100}.
#' @param color.palette Definition of color palette \code{Default="rainbow(n.levels, start = 0, end = 0.7)"}.
#' @param useRaster plot as a raster or vector image \code{Default=TRUE}.
#' WARNING plotting as a vector image is computationally intensive.
#' @param periodlab label for the y-axis \code{Default="Period (metres)"}.
#' @param x_lab label for the x-axis \code{Default="depth (metres)"}.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#' @return
#' The output is a plot of a wavelet spectra.
#'
#' @author
#' Code based on  ased on the \link[WaveletComp]{analyze.wavelet} function of the WaveletComp R package
#' and \link[biwavelet]{wt} function of the biwavelet R package which are based on the
#' wavelet MATLAB code written by Christopher Torrence and Gibert P. Compo.
#'
#' @references
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
#' \doi{<doi:10.1190/1.1441328>}
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#'  Geophysics 1982 47 (2): 222–236. \doi{<doi:10.1190/1.1441329>}
#'
#' @examples
#' \donttest{
#'#Example 1. A plot of a wavelet spectra using the Total Solar Irradiance
#'# data set of Steinhilver et al., (2012)
#'TSI_wt <-
#'  analyze_wavelet(
#'    data = TSI,
#'    dj = 1/200,
#'    lowerPeriod = 16,
#'    upperPeriod = 8192,
#'    verbose = TRUE,
#'    omega_nr = 6
#'  )
#'
#'plot_wavelet(
#'  wavelet = TSI_wt,
#' plot.COI = TRUE,
#'  n.levels = 100,
#'  color.palette = "rainbow(n.levels, start = 0, end = 0.7)",
#'  useRaster = TRUE,
#'  periodlab = "Period (years)",
#' x_lab = "years (before present)",
#'  keep_editable=FALSE
#')
#'
#'#Example 2. A plot of a wavelet spectra using the magnetic susceptibility
#'#data set of Pas et al., (2018)
#'mag_wt <-
#'analyze_wavelet(
#'data = mag,
#'dj = 1/100,
#'lowerPeriod = 0.1,
#'upperPeriod = 254,
#'verbose = TRUE,
#'omega_nr = 10
#')
#'plot_wavelet(
#'  wavelet = mag_wt,
#'  plot.COI = TRUE,
#'  n.levels = 100,
#'  color.palette = "rainbow(n.levels, start = 0, end = 0.7)",
#'  useRaster = TRUE,
#'  periodlab = "Period (metres)",
#'  x_lab = "depth (metres)",
#'   keep_editable=FALSE
#')
#'
#'#Example 3. A plot of a wavelet spectra using the greyscale
#'# data set of Zeeden et al., (2013)
#'grey_wt <-
#'  analyze_wavelet(
#'    data = grey,
#'    dj = 1/200,
#'    lowerPeriod = 0.02,
#'    upperPeriod = 256,
#'    verbose = TRUE,
#'    omega_nr = 8
#'  )
#'plot_wavelet(
#'  wavelet = grey_wt,
#'  plot.COI = TRUE,
#'  n.levels = 100,
#'  color.palette = "rainbow(n.levels, start = 0, end = 0.7)",
#'  useRaster = TRUE,
#'  periodlab = "Period (metres)",
#'  x_lab = "depth (metres)",
#'   keep_editable=FALSE
#')
#'
#'
#'
#'}
#' @export
#' @importFrom stats quantile
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics text
#' @importFrom graphics box
#' @importFrom graphics polygon
#' @importFrom grDevices rgb
#' @importFrom WaveletComp analyze.wavelet
#' @importFrom biwavelet wt


plot_wavelet <- function(wavelet = NULL,
                         plot.COI = TRUE,
                         n.levels = 100,
                         color.palette = "rainbow(n.levels, start = 0, end = 0.7)",
                         useRaster = TRUE,
                         periodlab = "Period (metres)",
                         x_lab = "depth (metres)",
                         keep_editable = FALSE) {
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

  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }
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

}
