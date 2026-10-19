#' @title Plots a superlet scalogram
#'
#' @description Plot superlet scalogram using the outcome of the \code{\link{analyse_superlet}} function.
#'
#'@param superlet superlet object created using the \code{\link{analyse_superlet}} function.
#'@param lowerPeriod Lowest period value which will be plotted
#'@param upperPeriod Highest period value which will be plotted
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
#'the wr_brewer_palettes() function
#'The R package 'colorRamps' has the color palette options:"blue2green",
#'"blue2green2red", "blue2red",	"blue2yellow", "colorRamps",	"cyan2yellow",
#'"green2red", "magenta2green", "matlab.like", "matlab.like2" and	"ygobb"
#'The R package 'grDevices' has the built in  palette options:"rainbow",
#'"heat.colors", "terrain.colors","topo.colors" and "cm.colors"
#'To see even more color palette options of the The R pacakge 'grDevices' run
#'the grDevices::hcl.pals() function
#'The R package 'scico' has the color palette options: “acton”, “bam”,“bamako”,
#' “bamO”, “batlow”, “batlowK”,“batlowW”,“berlin”,“bilbao”,”broc”,”brocO”,
#' ”buda”,”bukavu”,”cork”,”CorkO”,”davos”,”devon”,”fes”,”Glasgow”,”grayC”,
#' “hawaii”,”imola”,”lajolla”,”lapaz”,”lipari”,”lisbon”,”manague”,”navia”,
#' ”nuuk”,”oleron”,”oslo”,”roma”,”romaO”,”Tofino”,”Tokyo”,”turku”,”Vanimo”,
#' ”vik”,”vikO”
#'The R package 'Viridis' has the color palette options: “magma”, “plasma”,
#'“inferno”, “viridis”, “mako”, and “rocket”  and “turbo”
#'@param color_brewer Name of the R package from which the color palette is chosen from.
#'The included R packages from which palettes can be chosen
#'are; the RColorBrewer, grDevices, ColorRamps,scico and Viridis R packages.
#'There are many options to choose from so please
#'read the documentation of these packages. "\code{Default=grDevices}
#'@param useRaster Plot as a raster or vector image \code{Default=TRUE}.
#'WARNING plotting as a vector image is computationally intensive.
#'@param periodlab Label for the y-axis \code{Default="Period (metres)"}.
#'@param x_lab Label for the x-axis \code{Default="depth (metres)"}.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'@param dev_new Open a new plotting window for the plot. The wavelet plots contain a raster image that needs more room
#'than the small RStudio plot pane, which gives "figure margins too large" errors or clipped text, so a separate window is
#'opened by default in the RStudio GUI. In every other front end (plain R, ESS, Positron, knitr, Rscript) the plot is drawn
#'on the current device and no window is opened. Set TRUE or FALSE to force either behaviour, or set
#'\code{options(waverider.dev_new = TRUE/FALSE)}. Default: TRUE in the RStudio GUI, FALSE elsewhere.
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
#'@param add_data Plot the data on top of the wavelet \code{Default=TRUE}
#'@param add_avg Plot the average wavelet spectral power to the side of the wavelet \code{Default=FALSE}
#' @param add_signif Optional significance object created with \code{\link{wavelet_significance}}. The significant
#' cells are outlined on the scalogram and the critical level of the global spectrum is drawn as a red dashed line on
#' the average power panel \code{Default=NULL}.
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis is depth or y axis power \code{Default=TRUE}
#'
#'
#' @return
#' The output is a plot of a superlet scalogram.
#'
#' @author
#' Code based on the "wt.image" functions of the 'WaveletComp' R package
#' Whereas the "analyse_superlet" that generates the input for the plotting function
#' is based on the matlab code in Moca et al. (2021)
#'
#' @references
#' Roesch, A., & Schmidbauer, H. (2018). WaveletComp: Computational
#' Wavelet Analysis. R package version 1.1.
#' \url{https://CRAN.R-project.org/package=WaveletComp}
#'
#' Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., & Mureșan, R. C. (2021).
#' Time-frequency super-resolution with superlets.
#' Nature Communications, 12(1), 337.
#' \doi{10.1038/s41467-020-20539-9}
#'
#'@examples
#' \donttest{
#'#Example 1. A plot of a wavelet spectra using the Total Solar Irradiance
#'# data set of Steinhilber et al., (2012)
#'
#'TSI_super <-
#'  analyse_superlet(
#'    data = TSI,
#'    Nf = 128,
#'    lowerPeriod = 16,
#'    upperPeriod = 8192,
#'    verbose = FALSE,
#'    c1=1,
#'    o = c(1,5),
#'    mult = TRUE)
#'
#'plot_superlet(
#'  superlet = TSI_super,
#'  lowerPeriod = 16,
#'  upperPeriod = 8192,
#'  n.levels = 100,
#'  palette_name = "rainbow",
#' color_brewer= "grDevices",
#'  useRaster = TRUE,
#'  periodlab = "Period (metres)",
#'  x_lab = "depth (metres)",
#'  keep_editable = FALSE,
#'  dev_new=TRUE,
#'  plot_dir = TRUE,
#'  add_lines = NULL,
#'  add_points= NULL,
#'  add_abline_h = NULL,
#'  add_abline_v = NULL,
#'  add_data = TRUE,
#'  add_avg = TRUE,
#'  plot_horizontal = TRUE)
#'
#'#Example 2. A plot of a wavelet spectra using the magnetic susceptibility
#'#data set of Pas et al., (2018)
#'mag_super <-
#'analyse_superlet(
#'    data = mag,
#'    Nf = 128,
#'    lowerPeriod = 0.1,
#'    upperPeriod = 254,
#'    verbose = FALSE,
#'    c1=1,
#'    o = c(1,5),
#'    mult = TRUE)
#'
#'
#'plot_superlet(
#'  superlet = mag_super,
#'  lowerPeriod = 0.1,
#'  upperPeriod = 254,
#'  n.levels = 100,
#'  palette_name = "rainbow",
#' color_brewer= "grDevices",
#'  useRaster = TRUE,
#'  periodlab = "Period (metres)",
#'  x_lab = "depth (metres)",
#'  keep_editable = FALSE,
#'  dev_new=TRUE,
#'  plot_dir = TRUE,
#'  add_lines = NULL,
#'  add_points= NULL,
#'  add_abline_h = NULL,
#'  add_abline_v = NULL,
#'  add_data = TRUE,
#'  add_avg = TRUE,
#'  plot_horizontal = TRUE)
#'
#'
#'
#'
#'#Example 3. A plot of a wavelet spectra using the greyscale
#'# data set of Zeeden et al., (2013)
#'grey_super <-
#'analyse_superlet(
#'    data = grey,
#'    Nf = 128,
#'    lowerPeriod = 0.02,
#'    upperPeriod = 256,
#'    verbose = FALSE,
#'    c1=1,
#'    o = c(1,5),
#'    mult = TRUE)
#'
#'plot_superlet(
#'  superlet = grey_super,
#'  lowerPeriod = 0.02,
#'  upperPeriod = 256,
#'  n.levels = 100,
#'  palette_name = "rainbow",
#' color_brewer= "grDevices",
#'  useRaster = TRUE,
#'  periodlab = "Period (metres)",
#'  x_lab = "depth (metres)",
#'  keep_editable = FALSE,
#'  dev_new=TRUE,
#'  plot_dir = TRUE,
#'  add_lines = NULL,
#'  add_points= NULL,
#'  add_abline_h = NULL,
#'  add_abline_v = NULL,
#'  add_data = TRUE,
#'  add_avg = TRUE,
#'  plot_horizontal = TRUE)
#'}
#'
#' @param scale Colour scale as in \code{\link{plot_wavelet}} \code{Default="auto"}.
#' @param colour_brewer British spelling of \code{color_brewer}; if given it overrides \code{color_brewer}.
#' @export
#' @importFrom stats quantile
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics text
#' @importFrom graphics box
#' @importFrom graphics polygon
#' @importFrom graphics layout
#' @importFrom graphics title
#' @importFrom grDevices rgb
#' @importFrom graphics abline
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors

plot_superlet <- function (superlet = NULL,
                           lowerPeriod = NULL,
                           upperPeriod = NULL,
                           n.levels = 100,
                           palette_name = "rainbow",
                           color_brewer = "grDevices",
                           useRaster = TRUE,
                           periodlab = "Period (metres)",
                           x_lab = "depth (metres)",
                           keep_editable = FALSE,
                           dev_new = wr_dev_new_default(),
                           plot_dir = TRUE,
                           add_lines = NULL,
                           add_points = NULL,
                           add_abline_h = NULL,
                           add_abline_v = NULL,
                           add_data = TRUE,
                           add_avg = FALSE,
                           plot_horizontal = TRUE,
                           add_signif = NULL,
  colour_brewer = NULL, scale = "auto")
{
  if (!is.null(colour_brewer)) color_brewer <- colour_brewer  # British spelling synonym
  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }
  maximum.level = max(superlet$Power,na.rm=TRUE)
  power_max_mat.levels = quantile(superlet$Power, probs = seq(
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
  if (dev_new == TRUE & plot_horizontal == TRUE) {
    wr_open_device(width = 15, height = 7)
  }
  if (dev_new == TRUE & plot_horizontal == FALSE) {
    wr_open_device(width = 7, height = 10)
  }
  depth <- superlet$x
  y_axis <- superlet$Period
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)

  if (plot_dir != TRUE) {
    xlim_vals = rev(c(min(superlet$x), max(superlet$x)))
  }else {
    xlim_vals = c(min(superlet$x), max(superlet$x))
  }
  if (is.null(lowerPeriod) == TRUE) {
    lowerPeriod <- min(superlet$Period)
  }
  if (is.null(upperPeriod) == TRUE) {
    upperPeriod <- max(superlet$Period)
  }
  ylim_vals = (c(lowerPeriod, upperPeriod))

  p <- list(wavelet = superlet, data_x = superlet$x, data_y = superlet$y, depth = superlet$x, ax2 = superlet$axis.2,
            z = t(superlet$Power), avg = superlet$Power.avg, per = superlet$Period,
            tick_min = min(log2(ylim_vals)), tick_max = max(log2(ylim_vals)), overlay = TRUE,
            gws_abline_all = TRUE, levels = power_max_mat.levels,
            n.levels = n.levels, key.cols = key.cols, key.marks = key.marks,
            key.labels = key.labels, lwd.axis = lwd.axis, pmax_avg_sel = superlet$Power, xlim_vals = xlim_vals,
            ylim_vals = sort(ylim_vals), periodlab = periodlab, x_lab = x_lab, main = main, periodtck = periodtck,
            periodtcl = periodtcl, add_lines = add_lines, add_points = add_points, add_abline_h = add_abline_h,
            add_abline_v = add_abline_v, add_MTM_peaks = FALSE, add_pval = FALSE,
            pval_abline = NULL, pval_cutoff = NULL, add_signif = add_signif,
            MTM_res_1 = NULL, MTM_res_2 = matrix(data = NA, ncol = 2, nrow = 1),
            plot_horizontal = plot_horizontal, add_data = add_data, add_avg = add_avg, add_MTM = FALSE)
  p$scale <- scale
  .pw_draw(p)

}
