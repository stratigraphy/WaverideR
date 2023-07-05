#' @title Plots a wavelet power spectra
#'
#' @description Plot wavelet spectra using the outcome of the \code{\link{analyze_wavelet}} function.
#'
#'@param wavelet wavelet object created using the \code{\link{analyze_wavelet}} function.
#'@param lowerPeriod Lowest period value which will be plotted
#'@param upperPeriod Highest period value which will be plotted
#'@param plot.COI Option to plot the cone of influence \code{Default=TRUE}.
#'@param n.levels  Number of color levels \code{Default=100}.
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
#'@param useRaster Plot as a raster or vector image \code{Default=TRUE}.
#'WARNING plotting as a vector image is computationally intensive.
#'@param periodlab Label for the y-axis \code{Default="Period (metres)"}.
#'@param x_lab Label for the x-axis \code{Default="depth (metres)"}.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'@param dev_new Opens a new plotting window to plot the plot, this guarantees a "nice" looking plot however when plotting in an R markdown
#'document the plot might not plot  \code{Default=TRUE}
#'@param time_dir The direction of the proxy record which is assumed for tuning if time increases with increasing depth/time values
#'(e.g. bore hole data which gets older with increasing depth ) then time_dir should be set to TRUE
#'if time decreases with depth/time values (eg stratospheric logs where 0m is the bottom of the section)
#'then time_dir should be set to FALSE \code{time_dir=TRUE}
#'@param add_lines Add  lines to the wavelet plot input should be matrix with first axis being depth/time the columns after that
#'should be period values  \code{Default=NULL}
#'@param add_points Add points to the wavelet plot input should be matrix with first axis being depth/time and columns after that
#'should be period values \code{Default=NULL}
#'@param add_abline_h Add horizontal lines to the plot. Specify the lines as a vector e.g. c(2,3,5,6)  \code{Default=NULL}
#'@param add_abline_v Add vertical lines to the plot. Specify the lines as a vector e.g. c(2,3,5,6)  \code{Default=NULL}
#'@param add_MTM_peaks Add the MTM peak periods as horizontal lines \code{Default=FALSE}
#'@param add_data Plot the data on top of the wavelet \code{Default=TRUE}
#'@param add_avg Plot the average wavelet spectral power to the side of the wavelet \code{Default=FALSE}
#'@param add_MTM Add the MTM  plot next to the wavelet plot \code{Default=FALSE}
#'@param siglvl Set the significance level for the MTM analysis (0-1) \code{Default=0.95}
#'@param demean_mtm Remove mean from data before conducting the MTM analysis \code{Default=TRUE}
#'@param detrend_mtm Remove mean from data before conducting the MTM analysis \code{Default=TRUE}
#'@param padfac_mtm Pad factor for the MTM analysis \code{Default=5}
#'@param tbw_mtm time bandwidth product of the MTM analysis  \code{Default=3}
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis is depth or y axis power \code{Default=TRUE}
#'
#'
#' @return
#' The output is a plot of a wavelet spectra.
#' if add_MTM_peaks = TRUE then the output of the MTM analysis will given as matrix
#'
#' @author
#' Code based on the \link[WaveletComp]{analyze.wavelet} and \link[WaveletComp]{wt.image} functions of the 'WaveletComp' R package
#' and \link[biwavelet]{wt} function of the 'biwavelet' R package which are based on the
#' wavelet MATLAB code written by Christopher Torrence and Gibert P. Compo (1998).
#' The MTM analysis is from the astrochron R package of Meyers et al., (2012)
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
#'  Geophysics 1982 47 (2): 222–236. <\doi{doi:10.1190/1.1441329}>
#'
#'S.R. Meyers, 2012, Seeing Red in Cyclic Stratigraphy: Spectral Noise Estimation for
#'Astrochronology: Paleoceanography, 27, PA3228, <\doi{doi:10.1029/2012PA002307}>
#'
#' @examples
#' \donttest{
#'#Example 1. A plot of a wavelet spectra using the Total Solar Irradiance
#'# data set of Steinhilber et al., (2012)
#'
#'TSI_wt <-
#'  analyze_wavelet(
#'    data = TSI,
#'    dj = 1/200,
#'    lowerPeriod = 16,
#'    upperPeriod = 8192,
#'    verbose = FALSE,
#'    omega_nr = 6
#'  )
#'
#'plot_wavelet(
#'  wavelet = TSI_wt,
#'  lowerPeriod = NULL,
#'  upperPeriod = NULL,
#'  plot.COI = TRUE,
#'  n.levels = 100,
#'  palette_name = "rainbow",
#' color_brewer= "grDevices",
#'  useRaster = TRUE,
#'  periodlab = "Period (metres)",
#'  x_lab = "depth (metres)",
#'  keep_editable = FALSE,
#'  dev_new=TRUE,
#'  time_dir = TRUE,
#'  add_lines = NULL,
#'  add_points= NULL,
#'  add_abline_h = NULL,
#'  add_abline_v = NULL,
#'  add_MTM_peaks = FALSE,
#'  add_data = TRUE,
#'  add_avg = TRUE,
#'  add_MTM = FALSE,
#'  siglvl = 0.95,
#'  demean_mtm = TRUE,
#'  detrend_mtm = TRUE,
#'  padfac_mtm = 5,
#'  tbw_mtm = 3,
#'  plot_horizontal=TRUE)
#'
#'#Example 2. A plot of a wavelet spectra using the magnetic susceptibility
#'#data set of Pas et al., (2018)
#'mag_wt <-
#'analyze_wavelet(
#'data = mag,
#'dj = 1/100,
#'lowerPeriod = 0.1,
#'upperPeriod = 254,
#'verbose = FALSE,
#'omega_nr = 10
#')
#'
#'plot_wavelet(
#'wavelet = mag_wt,
#'lowerPeriod = NULL,
#'upperPeriod = NULL,
#'plot.COI = TRUE,
#'n.levels = 100,
#' palette_name = "rainbow",
#' color_brewer= "grDevices",
#'useRaster = TRUE,
#'periodlab = "Period (metres)",
#'x_lab = "depth (metres)",
#'keep_editable = FALSE,
#'dev_new=TRUE,
#'time_dir = TRUE,
#'add_lines= NULL,
#'add_points= NULL,
#'add_abline_h = NULL,
#'add_abline_v = NULL,
#'add_MTM_peaks = FALSE,
#'add_data = TRUE,
#'add_avg = TRUE,
#'add_MTM = FALSE,
#'siglvl = 0.95,
#'demean_mtm = TRUE,
#'detrend_mtm = TRUE,
#'padfac_mtm = 5,
#'tbw_mtm = 3,
#'plot_horizontal=TRUE)
#'
#'
#'#Example 3. A plot of a wavelet spectra using the greyscale
#'# data set of Zeeden et al., (2013)
#'grey_wt <-
#'  analyze_wavelet(
#'    data = grey,
#'    dj = 1/200,
#'    lowerPeriod = 0.02,
#'    upperPeriod = 256,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'plot_wavelet(
#'wavelet = grey_wt,
#'lowerPeriod = NULL,
#'upperPeriod = NULL,
#'plot.COI = TRUE,
#'n.levels = 100,
#' palette_name = "rainbow",
#' color_brewer= "grDevices",
#'useRaster = TRUE,
#'periodlab = "Period (metres)",
#'x_lab = "depth (metres)",
#'keep_editable = FALSE,
#'dev_new=TRUE,
#'time_dir = TRUE,
#'add_lines = NULL,
#'add_points= NULL,
#'add_abline_h = NULL,
#'add_abline_v = NULL,
#'add_MTM_peaks = FALSE,
#'add_data = TRUE,
#'add_avg = TRUE,
#'add_MTM = FALSE,
#'siglvl = 0.95,
#'demean_mtm = TRUE,
#'detrend_mtm = TRUE,
#'padfac_mtm = 5,
#'tbw_mtm = 3,
#'plot_horizontal=TRUE)
#'
#'}
#'
#' @export
#' @importFrom stats quantile
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics text
#' @importFrom graphics box
#' @importFrom graphics polygon
#' @importFrom graphics title
#' @importFrom grDevices rgb
#' @importFrom WaveletComp analyze.wavelet
#' @importFrom WaveletComp wt.image
#' @importFrom biwavelet wt
#' @importFrom astrochron mtm
#' @importFrom DescTools Closest
#' @importFrom graphics abline
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

plot_wavelet <- function(wavelet = NULL,
                         lowerPeriod = NULL,
                         upperPeriod = NULL,
                         plot.COI = TRUE,
                         n.levels = 100,
                         palette_name = "rainbow",
                         color_brewer= "grDevices",
                         useRaster = TRUE,
                         periodlab = "Period (metres)",
                         x_lab = "depth (metres)",
                         keep_editable = FALSE,
                         dev_new=TRUE,
                         time_dir = TRUE,
                         add_lines = NULL,
                         add_points= NULL,
                         add_abline_h = NULL,
                         add_abline_v = NULL,
                         add_MTM_peaks = FALSE,
                         add_data = TRUE,
                         add_avg = FALSE,
                         add_MTM = FALSE,
                         siglvl = 0.95,
                         demean_mtm = TRUE,
                         detrend_mtm = TRUE,
                         padfac_mtm = 5,
                         tbw_mtm = 3,
                         plot_horizontal = TRUE) {

    if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }


  maximum.level = max(wavelet$Power)
  power_max_mat.levels = quantile(wavelet$Power, probs = seq(from = 0,
                                                             to = 1, length.out = n.levels + 1))


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

  key.marks = round(seq(from = 0, to = 1, length.out = legend.params$n.ticks) *
                      n.levels)

  key.labels = formatC(as.numeric(power_max_mat.levels), digits = legend.params$label.digits,
                       format = legend.params$label.format)[key.marks + 1]

  if(dev_new==TRUE & plot_horizontal==TRUE){
    dev.new(width = 15,
            height = 7,
            noRStudioGD = TRUE)}

  if(dev_new==TRUE & plot_horizontal==FALSE){
    dev.new(width = 7,
            height = 10,
            noRStudioGD = TRUE)}



  y_axis <- as.numeric(unlist(wavelet$Period))
  pmax_avg_sel <- t(wavelet$Power)

  depth <-  wavelet$x
  y_axis <- wavelet$Period
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)

  #library(astrochron)
  #library(DescTools)

  if (add_MTM_peaks == TRUE) {
    MTM_res_1 <- mtm(
      cbind(wavelet$x, wavelet$y),
      tbw = tbw_mtm,
      padfac = padfac_mtm,
      output = 1,
      siglevel = siglvl,
      genplot = FALSE,
      verbose = FALSE,
      demean = demean_mtm,
      detrend = detrend_mtm
    )
    MTM_res_2 <-
      mtm(
        cbind(wavelet$x, wavelet$y),
        tbw = tbw_mtm,
        padfac = padfac_mtm,
        output = 3,
        siglevel = siglvl,
        genplot = FALSE,
        verbose = FALSE,
        demean = demean_mtm,
        detrend = detrend_mtm
      )
    MTM_res_1$period <- 1 / MTM_res_1[, 1]
    MTM_res_1$Peak <- 0

    for (i  in 1:nrow(MTM_res_2)) {
      row_nr <- Closest(MTM_res_1[, 1], MTM_res_2[i, 1], which = TRUE)
      row_nr <- row_nr[1]
      MTM_res_1[row_nr, 10] <- 1
    }
    mtm_res <- MTM_res_1

  }

  if (time_dir != TRUE) {
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





  if (add_data == TRUE & add_avg == FALSE & add_MTM == FALSE & plot_horizontal==FALSE) {

    layout.matrix <- matrix(c(1,  3, 0,2),
                            nrow = 2,
                            ncol = 2 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(1, 0.25),
           # Heights of the two rows
           widths = c(1, 4))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(4, 4, 2, 0))


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

    if (is.null(add_abline_h) != TRUE) {
      abline(h = add_abline_h)
    }


    par(new = FALSE, mar = c(3, 0, 2, 2))

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
      xlab = "",
      ylab = "",
    )


    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 1, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 2, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)


    par(new = FALSE, mar = c(4, 0, 2, 2))

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
      yaxt = "n" ,
      xaxt = "n" ,
      main = main,
      ylim = xlim_vals,
      xlim = log2(ylim_vals)
    )


    polygon(
      y=wavelet$coi.1 ,
      x=wavelet$coi.2,
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


    # axis(
    #   4,
    #   lwd = lwd.axis,
    #   at = period.tick,
    #   labels = NA,
    #   tck = periodtck,
    #   tcl = periodtcl
    # )

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
        lines(y=add_lines[, 1], x=log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        lines(y=add_points[, 1], x=log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(v = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }


  }

  if (add_data == FALSE & add_avg == FALSE & add_MTM == FALSE & plot_horizontal==FALSE) {

    layout.matrix <- matrix(c(1,  2),
                            nrow = 2,
                            ncol = 1 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(1, 4),
           # Heights of the two rows
           widths = c(1))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(2, 4, 2, 3))


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
      xlab = "",
      ylab = "",
    )


    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 1, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 2, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)


    par(new = FALSE, mar = c(4, 4, 2, 3))

    image(
      y = wavelet$x,
      x = wavelet$axis.2,
      z = (wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xlab = periodlab,
      ylab = x_lab,
      #axes = FALSE,
      #yaxt = "n" ,
      xaxt = "n" ,
      main = main,
      ylim = xlim_vals,
      xlim = log2(ylim_vals)
    )


    polygon(
      y=wavelet$coi.1 ,
      x=wavelet$coi.2,
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


    # axis(
    #   4,
    #   lwd = lwd.axis,
    #   at = period.tick,
    #   labels = NA,
    #   tck = periodtck,
    #   tcl = periodtcl
    # )

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
        lines(y=add_lines[, 1], x=log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        lines(y=add_points[, 1], x=log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(v = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }


  }

  if (add_data == FALSE & add_avg == FALSE & add_MTM == TRUE & plot_horizontal ==FALSE) {


    layout.matrix <- matrix(c(2,1,3,0,4,0,5,0),
                            nrow = 4,
                            ncol = 2 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(1, 1,1,4),
           # Heights of the two rows
           widths = c(4, 1))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(2, 2, 2, 3))


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
      xlab = "",
      ylab = "",
    )


    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 1, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 2, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)



    par(mar = c(0, 4, 2, 2))

    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 2],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "MTM power",
      xlab = "",
      log = "x",
      xlim = ylim_vals
    )

    for (i in 5:8) {
      lines(
        x = 1 / MTM_res_1[, 1],
        y = MTM_res_1[, i],
        xlim = ylim_vals,
        lty = 3,
        col = "grey",
        lwd = 2
      )
    }

    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    par(mar = c(0, 4, 0, 2))

    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 4],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "Ar. conf lvl",
      xlab = "",
      log = "x",
      ylim = c(80, 101),
      xaxs = "i",
      xlim = ylim_vals
    )

    abline(
      h = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )

    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 3],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "Har. conf lvl",
      xlab = "",
      log = "xy",
      ylim = c(80, 101),
      yaxs = "i",
      xlim = ylim_vals
    )

    abline(
      h = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )
    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }


    par(new = FALSE, mar = c(4, 4, 0, 2))

    image(
      y = wavelet$x,
      x = wavelet$axis.2,
      z = (wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xlab = periodlab,
      ylab = x_lab,
      #axes = FALSE,
      #yaxt = "n" ,
      xaxt = "n" ,
      main = main,
      ylim = xlim_vals,
      xlim = log2(ylim_vals)
    )


    polygon(
      y=wavelet$coi.1 ,
      x=wavelet$coi.2,
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


    # axis(
    #   4,
    #   lwd = lwd.axis,
    #   at = period.tick,
    #   labels = NA,
    #   tck = periodtck,
    #   tcl = periodtcl
    # )

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
        lines(y=add_lines[, 1], x=log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        lines(y=add_points[, 1], x=log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(v = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }


  }

  if (add_data == FALSE & add_avg == TRUE & add_MTM == TRUE & plot_horizontal ==FALSE) {

    layout.matrix <- matrix(c(2,1,3,0,4,0,5,0,6,0),
                            nrow = 5,
                            ncol = 2 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(1, 1,1,1,4),
           # Heights of the two rows
           widths = c(4, 1))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(2, 2, 2, 3))


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
      xlab = "",
      ylab = "",
    )


    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 1, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 2, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)



    par(mar = c(0, 4, 2, 2))

    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 2],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "MTM power",
      xlab = "",
      log = "x",
      xlim = ylim_vals
    )

    for (i in 5:8) {
      lines(
        x = 1 / MTM_res_1[, 1],
        y = MTM_res_1[, i],
        xlim = ylim_vals,
        lty = 3,
        col = "grey",
        lwd = 2
      )
    }

    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    par(mar = c(0, 4, 0, 2))

    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 4],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "Ar. conf lvl",
      xlab = "",
      log = "x",
      ylim = c(80, 101),
      xaxs = "i",
      xlim = ylim_vals
    )

    abline(
      h = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )

    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 3],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "Har. conf lvl",
      xlab = "",
      log = "xy",
      ylim = c(80, 101),
      yaxs = "i",
      xlim = ylim_vals
    )

    abline(
      h = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )
    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      y = wavelet$Power.avg,
      x = wavelet$Period,
      log = "x",
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "Wt. power",
      xaxs = "i",
      xlim = ylim_vals
    )





    par(new = FALSE, mar = c(4, 4, 0, 2))

    image(
      y = wavelet$x,
      x = wavelet$axis.2,
      z = (wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xlab = periodlab,
      ylab = x_lab,
      #axes = FALSE,
      #yaxt = "n" ,
      xaxt = "n" ,
      main = main,
      ylim = xlim_vals,
      xlim = log2(ylim_vals)
    )


    polygon(
      y=wavelet$coi.1 ,
      x=wavelet$coi.2,
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


    # axis(
    #   4,
    #   lwd = lwd.axis,
    #   at = period.tick,
    #   labels = NA,
    #   tck = periodtck,
    #   tcl = periodtcl
    # )

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
        lines(y=add_lines[, 1], x=log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        lines(y=add_points[, 1], x=log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(v = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }


  }

  if (add_data == TRUE & add_avg == TRUE & add_MTM == TRUE & plot_horizontal ==FALSE) {


    layout.matrix <- matrix(c(1,2,0,3,0,4,0,5,6,7),
                            nrow = 5,
                            ncol = 2 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(1, 1,1,1,4),
           # Heights of the two rows
           widths = c(1, 4))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(0, 1, 2, 6),xpd=NA)


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
      xlab = "",
      ylab = "",
    )


    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 1, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 1, at = mean(key.marks), line = 4,
          font = par()$font.axis, cex = par()$cex.axis,las=1)
    box(lwd = lwd.axis)



    par(mar = c(0,0, 2, 2),xpd=FALSE)

    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 2],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "MTM power",
      xlab = "",
      log = "x",
      xlim = ylim_vals
    )


    for (i in 5:8) {
      lines(
        x = 1 / MTM_res_1[, 1],
        y = MTM_res_1[, i],
        xlim = ylim_vals,
        lty = 3,
        col = "grey",
        lwd = 2
      )
    }

    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    par(mar = c(0, 0, 0, 2),xpd=FALSE)

    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 4],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "",
      xlab = "",
      log = "x",
      ylim = c(80, 101),
      xaxs = "i",
      xlim = ylim_vals
    )

    abline(
      h = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )

    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    title(ylab="Ar. conf lvl",xpd=NA)


    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 3],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "",
      xlab = "",
      log = "xy",
      ylim = c(80, 101),
      yaxs = "i",
      xlim = ylim_vals
    )
    title(ylab="Har. conf lvl",xpd=NA)

    abline(
      h = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )
    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      y = wavelet$Power.avg,
      x = wavelet$Period,
      log = "x",
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "Wt. power",
      xaxs = "i",
      xlim = ylim_vals
    )

    title(ylab="Wt. power",xpd=NA)


    par( mar = c(4, 4, 0, 0),xpd=FALSE)


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


    par(new = FALSE, mar = c(4, 0, 0, 2))

    image(
      y = wavelet$x,
      x = wavelet$axis.2,
      z = (wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xlab = periodlab,
      #ylab = x_lab,
      #axes = FALSE,
      #yaxt = "n" ,
      xaxt = "n" ,
      yaxt = "n" ,
      main = main,
      ylim = xlim_vals,
      xlim = log2(ylim_vals)
    )


    polygon(
      y=wavelet$coi.1 ,
      x=wavelet$coi.2,
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


    # axis(
    #   4,
    #   lwd = lwd.axis,
    #   at = period.tick,
    #   labels = NA,
    #   tck = periodtck,
    #   tcl = periodtcl
    # )

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
        lines(y=add_lines[, 1], x=log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        lines(y=add_points[, 1], x=log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(v = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }


  }

  if (add_data == TRUE & add_avg == TRUE & add_MTM == FALSE & plot_horizontal ==FALSE) {


    layout.matrix <- matrix(c(1,2,3,4),
                            nrow = 2,
                            ncol = 2 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(1, 4),
           # Heights of the two rows
           widths = c(1, 4))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(2, 0.5, 2, 6),xpd=NA)


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


    axis(2, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 2, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 1, at = 1, line = 0.5,
          font = par()$font.axis, cex = par()$cex.axis,las=1,xpd=NA)
    #

    box(lwd = lwd.axis)



    par(mar = c(0,0, 2, 2))


    plot(
      y = wavelet$Power.avg,
      x = wavelet$Period,
      log = "x",
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "Wt. power",
      xlab="",
      xaxs = "i",
      xlim = ylim_vals
    )

    title(ylab="Wt. power",xpd=NA)


    par( mar = c(4, 4, 0, 0))


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


    par(new = FALSE, mar = c(4, 0, 0, 2))

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
      y=wavelet$coi.1 ,
      x=wavelet$coi.2,
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


    # axis(
    #   4,
    #   lwd = lwd.axis,
    #   at = period.tick,
    #   labels = NA,
    #   tck = periodtck,
    #   tcl = periodtcl
    # )

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
        lines(y=add_lines[, 1], x=log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        lines(y=add_points[, 1], x=log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(v = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }


  }

  if (add_data == TRUE & add_avg == FALSE & add_MTM == TRUE & plot_horizontal ==FALSE) {


    layout.matrix <- matrix(c(1,2,0,3,0,4,5,6),
                            nrow = 4,
                            ncol = 2 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(1, 1,1,4),
           # Heights of the two rows
           widths = c(1, 4))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(0, 1, 2, 6),xpd=NA)


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
      xlab = "",
      ylab = "",
    )


    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 1, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 1, at = mean(key.marks), line = 4,
          font = par()$font.axis, cex = par()$cex.axis,las=1)
    box(lwd = lwd.axis)



    par(mar = c(0,0, 2, 2))

    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 2],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "MTM power",
      xlab = "",
      log = "x",
      xlim = ylim_vals
    )


    for (i in 5:8) {
      lines(
        x = 1 / MTM_res_1[, 1],
        y = MTM_res_1[, i],
        xlim = ylim_vals,
        lty = 3,
        col = "grey",
        lwd = 2
      )
    }

    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    par(mar = c(0, 0, 0, 2),xpd=FALSE)

    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 4],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "",
      xlab = "",
      log = "x",
      ylim = c(80, 101),
      xaxs = "i",
      xlim = ylim_vals
    )

    abline(
      h = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )

    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    title(ylab="Ar. conf lvl",xpd=NA)


    plot(
      x = 1 / MTM_res_1[, 1],
      y = MTM_res_1[, 3],
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "",
      xlab = "",
      log = "xy",
      ylim = c(80, 101),
      yaxs = "i",
      xlim = ylim_vals
    )
    title(ylab="Har. conf lvl",xpd=NA)

    abline(
      h = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )
    if (add_MTM_peaks == TRUE) {
      abline(v = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }



    par( mar = c(4, 4, 0, 0))


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


    par(new = FALSE, mar = c(4, 0, 0, 2))

    image(
      y = wavelet$x,
      x = wavelet$axis.2,
      z = (wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xlab = periodlab,
      #ylab = x_lab,
      #axes = FALSE,
      #yaxt = "n" ,
      xaxt = "n" ,
      yaxt = "n" ,
      main = main,
      ylim = xlim_vals,
      xlim = log2(ylim_vals)
    )


    polygon(
      y=wavelet$coi.1 ,
      x=wavelet$coi.2,
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


    # axis(
    #   4,
    #   lwd = lwd.axis,
    #   at = period.tick,
    #   labels = NA,
    #   tck = periodtck,
    #   tcl = periodtcl
    # )

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
        lines(y=add_lines[, 1], x=log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        lines(y=add_points[, 1], x=log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(v = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }


  }

  if (add_data == FALSE & add_avg == TRUE & add_MTM == FALSE & plot_horizontal ==FALSE) {

    layout.matrix <- matrix(c(1,2,0,3),
                            nrow = 2,
                            ncol = 2 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(1,4),
           # Heights of the two rows
           widths = c(1, 4))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(2, 2, 2, 3))


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
      xlab = "",
      ylab = "",
    )


    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 1, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 1, at = median(key.marks), line = 3,
          las = 2, font = par()$font.axis, cex = par()$cex.axis,xpd=NA)
    box(lwd = lwd.axis)


    par(mar = c(0, 4, 2, 2))

    plot(
      y = wavelet$Power.avg,
      x = wavelet$Period,
      log = "x",
      type = "l",
      xaxs = "i",
      xaxt = "n",
      ylab = "Wt. power",
      xaxs = "i",
      xlim = ylim_vals
    )


    if (is.null(add_abline_v) != TRUE) {
      abline(v = (add_abline_v))
    }


    par(new = FALSE, mar = c(4, 4, 0, 2))

    image(
      y = wavelet$x,
      x = wavelet$axis.2,
      z = (wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xlab = periodlab,
      ylab = x_lab,
      #axes = FALSE,
      #yaxt = "n" ,
      xaxt = "n" ,
      main = main,
      ylim = xlim_vals,
      xlim = log2(ylim_vals)
    )


    polygon(
      y=wavelet$coi.1 ,
      x=wavelet$coi.2,
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


    # axis(
    #   4,
    #   lwd = lwd.axis,
    #   at = period.tick,
    #   labels = NA,
    #   tck = periodtck,
    #   tcl = periodtcl
    # )

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
        lines(y=add_lines[, 1], x=log2(add_lines[, i]))
    }

    if (is.null(add_points) != TRUE) {
      for (i  in 2:ncol(add_points))
        lines(y=add_points[, 1], x=log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(v = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }


  }


  if (add_data == TRUE & add_avg == FALSE & add_MTM == TRUE & plot_horizontal==TRUE) {
    layout.matrix <- matrix(
      c(1, 1, 1, 1, 1, 1, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 4, 5, 6),
      nrow = 2,
      ncol = 9 ,
      byrow = TRUE
    )
    layout(mat = layout.matrix,
           heights = c(0.25, 1),
           # Heights of the two rows
           widths = c(1))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0 ,
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


    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }

    par(new = FALSE,mar = c(3, 2, 2, 2),
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

    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl =  1.2)

    mtext(key.labels, side = 1, at = key.marks, line = 0.1,
          las = 2, cex=0.75)
    box(lwd = lwd.axis)


    par(new = FALSE,
        mar = c(4, 4, 0, 0),
        mgp = c(2, 1, 0))


    image(
      x = wavelet$x,
      y = wavelet$axis.2,
      z = t(wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
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
      wavelet$coi.1 ,
      wavelet$coi.2,
      border = NA,
      col = rgb(1, 1, 1, 0.5),
      xlim = xlim_vals,
      ylim = log2(ylim_vals)
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
        lines(add_points[, 1], log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(h = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }


    par(new = FALSE, mar = c(4, 0, 0, 0.5))

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 2],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "MTM power",
      ylab = "",
      log = "y",
      ylim = ylim_vals
    )

    for (i in 5:8) {
      lines(
        y = 1 / MTM_res_1[, 1],
        x = MTM_res_1[, i],
        ylim = ylim_vals,
        lty = 3,
        col = "grey",
        lwd = 2
      )
    }
    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 4],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "Ar. conf lvl",
      ylab = "",
      log = "y",
      xlim = c(80, 101),
      xaxs = "i",
      ylim = ylim_vals
    )

    abline(
      v = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )

    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 3],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "Har. conf lvl",
      ylab = "",
      log = "xy",
      xlim = c(80, 101),
      xaxs = "i",
      ylim = ylim_vals
    )

    abline(
      v = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )
    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }
  }

  if (add_data == TRUE & add_avg == TRUE & add_MTM == FALSE & plot_horizontal==TRUE) {
    layout.matrix <- matrix(c(1, 2, 4, 3),
                            nrow = 2,
                            ncol = 2 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
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

    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl =  1.2 )

    mtext(key.labels, side = 1, at = key.marks, line = 0.1,
          las = 2, cex=0.75)
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

    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = add_abline_h)
    }




    par(new = FALSE, mar = c(4, 4, 0, 0))

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
      wavelet$coi.1 ,
      wavelet$coi.2,
      border = NA,
      col = rgb(1, 1, 1, 0.5),
      xlim = xlim_vals,
      ylim = log2(ylim_vals)
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
        lines(add_points[, 1], log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(h = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }

  }



  if (add_data == TRUE & add_avg == FALSE & add_MTM == FALSE & plot_horizontal ==TRUE ) {
    layout.matrix <- matrix(c(1, 0, 3, 2),
                            nrow = 2,
                            ncol = 2 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(0.25, 1),
           # Heights of the two rows
           widths = c(8, 2))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(0, 4, 2, 2))


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

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }


    par(new = FALSE, mar = c(4, 0, 2, 5))

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

    axis(4, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 4, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 4, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)


    par(new = FALSE, mar = c(4, 4, 0, 2))

    image(
      x = wavelet$x,
      y = wavelet$axis.2,
      z = t(wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
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
      wavelet$coi.1 ,
      wavelet$coi.2,
      border = NA,
      col = rgb(1, 1, 1, 0.5),
      xlim = xlim_vals,
      ylim = log2(ylim_vals)
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
        lines(add_points[, 1], log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(h = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }


  }

  if (add_data == TRUE & add_avg == TRUE & add_MTM == TRUE & plot_horizontal ==TRUE) {
    layout.matrix <- matrix(
      c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 4, 5, 6, 7),
      nrow = 2,
      ncol = 10 ,
      byrow = TRUE
    )
    layout(mat = layout.matrix,
           heights = c(0.25, 1),
           # Heights of the two rows
           widths = c(1))

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

    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl =  1.2 )

    mtext(key.labels, side = 1, at = key.marks, line = 0.1,
          las = 2, cex=0.75)
    box(lwd = lwd.axis)


    par(new = FALSE,
        mar = c(4, 4, 0, 0),
        mgp = c(2, 1, 0))


    image(
      x = wavelet$x,
      y = wavelet$axis.2,
      z = t(wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
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
      wavelet$coi.1 ,
      wavelet$coi.2,
      border = NA,
      col = rgb(1, 1, 1, 0.5),
      xlim = xlim_vals,
      ylim = log2(ylim_vals)
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
        lines(add_points[, 1], log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(h = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }



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

    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "black",
             lty = 3)
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = add_abline_h)
    }

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 2],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "MTM power",
      ylab = "",
      log = "y",
      ylim = ylim_vals
    )

    for (i in 5:8) {
      lines(
        y = 1 / MTM_res_1[, 1],
        x = MTM_res_1[, i],
        ylim = ylim_vals,
        lty = 3,
        col = "grey",
        lwd = 2
      )
    }
    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 4],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "Ar. conf lvl",
      ylab = "",
      log = "y",
      xlim = c(80, 101),
      xaxs = "i",
      ylim = ylim_vals
    )

    abline(
      v = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )

    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 3],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "Har. conf lvl",
      ylab = "",
      log = "xy",
      xlim = c(80, 101),
      xaxs = "i",
      ylim = ylim_vals
    )

    abline(
      v = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )
    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }
  }

  if (add_data == FALSE & add_avg == TRUE & add_MTM == FALSE & plot_horizontal ==TRUE) {
    layout.matrix <- matrix(c(3, 2, 1),
                            nrow = 1,
                            ncol = 3 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(1),
           # Heights of the two rows
           widths = c(6, 1, 1))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))
    par(mar = c(4, 0, 2, 5))


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

    axis(4, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 4, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 4, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)


    par(new = FALSE, mar = c(4, 0, 2, 0.5))

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

    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = add_abline_h)
    }



    par(new = FALSE,
        mar = c(4, 4, 2, 0),
        mgp = c(2, 1, 0))

    image(
      x = wavelet$x,
      y = wavelet$axis.2,
      z = t(wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
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
      wavelet$coi.1 ,
      wavelet$coi.2,
      border = NA,
      col = rgb(1, 1, 1, 0.5),
      xlim = xlim_vals,
      ylim = log2(ylim_vals)
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
        lines(add_points[, 1], log2(add_points[, i]))
    }

    if (add_MTM_peaks == TRUE) {
      abline(h = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }







  }

  if (add_data == FALSE & add_avg == TRUE & add_MTM == TRUE & plot_horizontal ==TRUE) {
    layout.matrix <-
      matrix(c(6, 5, 4, 3, 2, 1),
             nrow = 1,
             ncol = 6 ,
             byrow = TRUE)
    layout(
      mat = layout.matrix,
      heights = c(1),
      # Heights of the two rows
      widths = c(6, 1, 1, 1, 1, 1)
    )

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(4, 0, 2, 5))

    image(
      y = seq(from = 0, to = n.levels),
      x = 1,
      z = (matrix(power_max_mat.levels,
                  nrow = 1)),
      col =key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      yaxt = "n",
      xaxt = "n",
      xlab = "",
      ylab = "",
    )

    axis(4, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 4, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 4, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)

    par(new = FALSE, mar = c(4, 0, 2, 0.5))

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 3],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "Har. conf lvl",
      ylab = "",
      log = "xy",
      xlim = c(80, 101),
      xaxs = "i",
      ylim = ylim_vals
    )

    abline(
      v = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )
    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 4],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "Ar. conf lvl",
      ylab = "",
      log = "y",
      xlim = c(80, 101),
      xaxs = "i",
      ylim = ylim_vals
    )

    abline(
      v = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )

    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 2],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "MTM power",
      ylab = "",
      log = "y",
      ylim = ylim_vals
    )

    for (i in 5:8) {
      lines(
        y = 1 / MTM_res_1[, 1],
        x = MTM_res_1[, i],
        ylim = c(min(y_axis), max(y_axis)),
        lty = 3,
        col = "grey",
        lwd = 2
      )
    }
    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

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

    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = add_abline_h)
    }


    par(new = FALSE,
        mar = c(4, 4, 2, 0),
        mgp = c(2, 1, 0))

    image(
      x = wavelet$x,
      y = wavelet$axis.2,
      z = t(wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
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
      wavelet$coi.1 ,
      wavelet$coi.2,
      border = NA,
      col = rgb(1, 1, 1, 0.5),
      xlim = xlim_vals,
      ylim = log2(ylim_vals)
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
        lines(add_points[, 1], log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(h = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }


  }

  if (add_data == FALSE & add_avg == FALSE & add_MTM == TRUE & plot_horizontal==TRUE) {
    layout.matrix <-
      matrix(c(5, 4, 3, 2, 1),
             nrow = 1,
             ncol = 5 ,
             byrow = TRUE)
    layout(
      mat = layout.matrix,
      heights = c(1),
      # Heights of the two rows
      widths = c(7, 1, 1, 1, 1)
    )

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))
    par(mar = c(4, 0, 2, 5))

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

    axis(4, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 4, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 4, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)


    par(new = FALSE, mar = c(4, 0, 2, 0.5))

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 3],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "Har. conf lvl",
      ylab = "",
      log = "xy",
      xlim = c(80, 101),
      xaxs = "i",
      ylim = ylim_vals
    )

    abline(
      v = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )
    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 4],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "Ar. conf lvl",
      ylab = "",
      log = "y",
      xlim = c(80, 101),
      xaxs = "i",
      ylim = ylim_vals
    )

    abline(
      v = c(90, 95, 99),
      lty = 3,
      col = "grey",
      lwd = 2
    )

    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    plot(
      y = 1 / MTM_res_1[, 1],
      x = MTM_res_1[, 2],
      type = "l",
      yaxs = "i",
      yaxt = "n",
      xlab = "MTM power",
      ylab = "",
      log = "y",
      ylim = ylim_vals
    )

    for (i in 5:8) {
      lines(
        y = 1 / MTM_res_1[, 1],
        x = MTM_res_1[, i],
        ylim = ylim_vals,
        lty = 3,
        col = "grey",
        lwd = 2
      )
    }
    if (add_MTM_peaks == TRUE) {
      abline(h = 1 / MTM_res_2[, 1],
             col = "red",
             lty = 3)
    }

    par(new = FALSE,
        mar = c(4, 4, 2, 0),
        mgp = c(2, 1, 0))

    image(
      x = wavelet$x,
      y = wavelet$axis.2,
      z = t(wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
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
      wavelet$coi.1 ,
      wavelet$coi.2,
      border = NA,
      col = rgb(1, 1, 1, 0.5),
      xlim = xlim_vals,
      ylim = log2(ylim_vals)
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
        lines(add_points[, 1], log2(add_points[, i]))
    }


    if (add_MTM_peaks == TRUE) {
      abline(h = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }

  }

  if (add_data == FALSE & add_avg == FALSE & add_MTM == FALSE & plot_horizontal==TRUE) {
    layout.matrix <- matrix(c(2, 1),
                            nrow = 1,
                            ncol = 2 ,
                            byrow = TRUE)
    layout(mat = layout.matrix,
           heights = c(1),
           # Heights of the two rows
           widths = c(10, 2.25))

    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = 0,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))

    par(mar = c(4, 0, 2, 5))

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

    axis(4, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)

    mtext(key.labels, side = 4, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 4, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)



    par(new = FALSE, mar = c(4, 4, 2, 0.5))

    image(
      x = wavelet$x,
      y = wavelet$axis.2,
      z = t(wavelet$Power),
      col = key.cols,
      breaks = power_max_mat.levels,
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
      wavelet$coi.1 ,
      wavelet$coi.2,
      border = NA,
      col = rgb(1, 1, 1, 0.5),
      xlim = xlim_vals,
      ylim = log2(ylim_vals)
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
        lines(add_points[, 1], log2(add_points[, i]))
    }

    if (add_MTM_peaks == TRUE) {
      abline(h = log2(1 / MTM_res_2[, 1]),
             col = "black",
             lty = 3)
    }

    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }




  }

  if (add_MTM_peaks == TRUE) {
    return(invisible(mtm_res))
  }
}
