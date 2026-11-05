#' @title Remove tracking points which were tracked in a wavelet spectra
#'
#' @description Interactively select points for deletion
#'With the  \code{\link{track_period_wavelet}} function it is possible to track points in a wavelet spectra,
#'however errors can be made and as such it is possible to delete these points with the \code{\link{delpts_tracked_period_wt}} function.
#'This function allows one to select points for deletion.
#'#'
#' @param tracking_pts Tracked points to clean (matrix or data frame, at least two columns): depth in column 1 and
#' the tracked period in column 2 (same units as the wavelet), e.g. the output of \code{\link{track_period_wavelet}}.
#' The points are drawn on the scalogram and the ones clicked are removed. Required (\code{Default=NULL} gives an
#' error).
#' @param wavelet Scalogram on which the points are shown: the output of \code{\link{analyse_wavelet}} (class
#' \code{analyze.wavelet}). Required (\code{Default=NULL} gives an error).
#' @param n.levels Number of colour levels of the power scale (integer, > 1). \code{Default=100}.
#' @param periodlab Label of the y-axis (character). \code{Default="Period (metres)"}.
#' @param x_lab Label of the x-axis (character). \code{Default="depth (metres)"}.
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
#'@return The results of the deletion of the tracking points is a matrix with 3 columns.
#'The first column is depth/time
#'The second column is the period of the tracked cycle
#'The third column is the sedimentation rate based on the duration (in time) of the tracked cycle
#'
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
#'#mag_track <- track_period_wavelet(astro_cycle = 405,
#'#                                   wavelet=mag_wt,
#'#                                   n.levels = 100,
#'#                                   periodlab = "Period (metres)",
#'#                                   x_lab = "depth (metres)"
#'#                                    palette_name ="rainbow",
#'#                                    color_brewer ="grDevices)
#'
#'#load the mag_track_solution data set to get an example data set from which
#'#data points can be deleted
#'
#'
#' mag_track_corr <- delpts_tracked_period_wt(tracking_pts = mag_track_solution,
#'                                     wavelet = mag_wt,
#'                                     n.levels = 100,
#'                                     periodlab = "Period (metres)",
#'                                     x_lab = "depth (metres)",
#'                                    palette_name ="rainbow",
#'                                    color_brewer ="grDevices")
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
#' @importFrom stats na.omit


delpts_tracked_period_wt <- function(tracking_pts = NULL,
    wavelet = NULL,
    n.levels = 100,
    periodlab = "Period (metres)",
    x_lab = "depth (metres)",
    palette_name = "viridis",
    color_brewer = "viridis",
    colour_brewer = NULL){
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym
  plot.COI = TRUE

  key.cols <- .wr_key_cols(palette_name, color_brewer, n.levels)


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
  on.exit(par(oldpar), add = TRUE)

  image.plt = par()$plt
  legend.plt = NULL
  wr_open_device(width = 15, height = 7, interactive_input = TRUE)

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

  oldop <- options(warn = -1)
  on.exit(options(oldop), add = TRUE)
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y

  sel <- rep(FALSE, length(x))
  res <- integer(0)
  n <- length(x)

  while (sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n = 1, plot = F,tolerance = 0.15)
    if (!length(ans)) {
      break
    }
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = 19, col = "white")
    sel[ans] <- TRUE
  }

  out <- tracking_pts[!sel,]


  return(out)

}
