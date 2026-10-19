#' @title Plots a cross superlet scalogram
#'
#' @description Plot cross superlet scalogram using the outcome of the \code{\link{analyse_Xsuperlet}} function.
#'
#'@param Xsuperlet Xsuperlet object created using the \code{\link{analyse_Xsuperlet}} function.
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
#'@param plot_arrows add the phase arrows of the cross superlet transform \code{Default=TRUE}
#'@param pwr_quant Power quantile above which the arrows are plotted  \code{Default= 0.85}
#'@param n_arrows number arrows width wise plotted on the time axis \code{Default= 50}
#'@param T_unit length of the arrows in depth/time units \code{Default=NULL}: automatic, a fixed visual length of 0.3 inch whatever the record length
#'@param P_unit length of the arrows in log2 period units \code{Default=NULL}: automatic (see \code{T_unit})
#'@param arrow_length length of the slanted part of the arrow  \code{Default= 0.1}
#'@param arrow_lwd thickness of the arrow \code{Default= 2}
#'@param arrow_angle angle of the slanted part of the arrow \code{Default= 25}
#'@param add_lines Add  lines to the wavelet plot input should be matrix with first axis being depth/time the columns after that
#'should be period values  \code{Default=NULL}
#'@param add_points Add points to the wavelet plot input should be matrix with first axis being depth/time and columns after that
#'should be period values \code{Default=NULL}
#'@param add_abline_h Add horizontal lines to the plot. Specify the lines as a vector e.g. c(2,3,5,6)  \code{Default=NULL}
#'@param add_abline_v Add vertical lines to the plot. Specify the lines as a vector e.g. c(2,3,5,6)  \code{Default=NULL}
#'@param add_data Plot the data on top of the wavelet \code{Default=TRUE}
#'@param add_avg Plot the average wavelet spectral power to the side of the wavelet \code{Default=FALSE}
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis is depth or y axis power \code{Default=TRUE}
#'
#'
#' @return
#' The output is a plot of a cross superlet scalogram.
#'
#' @author
#' plotting code based on the "wt.image" function of the 'WaveletComp' R package
#' Whereas the "analyse_Xsuperlet" that generates the input for the plotting function
#' is based on the matlab code in Moca et al. (2021) and the the "analyze.coherency" function of the  'WaveletComp' R package
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
#'#Example 1. A cross superlet plot of two etp solutions with noise overprint
#'
#'
#'etp_1 <- astrochron::etp(
#'  tmin = 0,
#'  tmax = 500,
#'  dt = 1,
#'  eWt = 1.5,
#'  oWt = 0.75,
#'  pWt = 1,
#'  esinw = TRUE,
#'  standardize = TRUE,
#'  genplot = FALSE,
#'  verbose = FALSE
#')
#'
#'etp_2 <- astrochron::etp(
#'  tmin = 0,
#'  tmax = 500,
#'  dt = 1,
#'  eWt = 1,
#'  oWt = 0.5,
#'  pWt = 1.5,
#'  esinw = TRUE,
#'  standardize = TRUE,
#'  genplot = FALSE,
#' verbose = FALSE
#')
#'
#'etp_1[, 2] <- etp_1[, 2] + wr_coloured_noise(
#'  nrow(etp_1),
#'  sd = sd(etp_1[, 2]) / 1.5,
#'  mean = mean(etp_1[, 2]),
#'  phi = 0.9
#')
#'etp_2[, 2] <- etp_2[, 2] + wr_coloured_noise(
#'  nrow(etp_2),
#'  sd = sd(etp_2[, 2]) / 1.5,
#'  mean = mean(etp_2[, 2]),
#'  phi = 0.9
#')
#'
#'X_super_etp <- analyse_Xsuperlet(
#'  data_1 = etp_1,
#'  data_2  = etp_2,
#'  upperPeriod = 1024,
#'  lowerPeriod = 2,
#'  Nf = 32,
#'  c1 = 3,
#'  o = c(1, 10),
#'  mult = TRUE,
#'  verbose = FALSE
#')
#'plot_Xsuperlet(Xsuperlet = X_super_etp, lowerPeriod = 2, upperPeriod = 1024,
#'               n.levels = 100, palette_name = "rainbow", color_brewer = "grDevices",
#'               useRaster = TRUE, periodlab = "Period (metres)", x_lab = "depth (metres)",
#'               keep_editable = FALSE, dev_new = TRUE, plot_dir = TRUE,
#'               plot_arrows = TRUE,
#'               pwr_quant = 0.85,
#'               n_arrows = 50,
#'               T_unit = NULL,
#'               P_unit  = NULL,
#'               arrow_length = 0.1,
#'               arrow_lwd = 2,
#'               arrow_angle = 25,add_lines = NULL, add_points = NULL, add_abline_h = NULL,
#'               add_abline_v = NULL, add_data = TRUE, add_avg = FALSE, plot_horizontal = TRUE)
#'}
#'
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
#' @importFrom graphics arrows


plot_Xsuperlet <- function (Xsuperlet = NULL, lowerPeriod = NULL, upperPeriod = NULL,
                            n.levels = 100, palette_name = "rainbow", color_brewer = "grDevices",
                            useRaster = TRUE, periodlab = "Period (metres)", x_lab = "depth (metres)",
                            keep_editable = FALSE, dev_new = wr_dev_new_default(), plot_dir = TRUE,
                            plot_arrows = TRUE,
                            pwr_quant = 0.85,
                            n_arrows = 75,
                            T_unit = NULL,
                            P_unit  = NULL,
                            arrow_length = 0.075,
                            arrow_lwd = 1.9,
                            arrow_angle = 25,add_lines = NULL, add_points = NULL, add_abline_h = NULL,
                            add_abline_v = NULL, add_data = TRUE, add_avg = FALSE, plot_horizontal = TRUE,
  colour_brewer = NULL)
{
  if (!is.null(colour_brewer)) color_brewer <- colour_brewer  # British spelling synonym

  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  Power <- t(Xsuperlet$Power)
  Phase <- t(Xsuperlet$Phase)
  maximum.level = max(Xsuperlet$Power, na.rm = TRUE)
  power_max_mat.levels = quantile(Xsuperlet$Power, probs = seq(from = 0,
                                                               to = 1, length.out = n.levels + 1))
  key.cols <- .wr_key_cols(palette_name, color_brewer, n.levels)
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
  if (dev_new == TRUE & plot_horizontal == TRUE) {
    wr_open_device(width = 15, height = 7)
  }
  if (dev_new == TRUE & plot_horizontal == FALSE) {
    wr_open_device(width = 7, height = 10)
  }
  depth <- Xsuperlet$x1
  y_axis <- Xsuperlet$Period
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)
  if (plot_dir != TRUE) {
    xlim_vals = rev(c(min(Xsuperlet$x1), max(Xsuperlet$x1)))
  }
  else {
    xlim_vals = c(min(Xsuperlet$x1), max(Xsuperlet$x1))
  }
  if (is.null(lowerPeriod) == TRUE) {
    lowerPeriod <- min(Xsuperlet$Period)
  }
  if (is.null(upperPeriod) == TRUE) {
    upperPeriod <- max(Xsuperlet$Period)
  }
  ylim_vals = (c(lowerPeriod, upperPeriod))
  p <- list(wavelet = Xsuperlet, data_x = Xsuperlet$x1, data_y = Xsuperlet$y1, data2 = list(x = Xsuperlet$x2, y = Xsuperlet$y2),
            depth = Xsuperlet$x1, ax2 = Xsuperlet$axis.2, z = Power, avg = Xsuperlet$Power.avg, per = Xsuperlet$Period,
            tick_min = min(log2(ylim_vals)), tick_max = max(log2(ylim_vals)), overlay = TRUE,
            gws_abline_all = TRUE, levels = power_max_mat.levels,
            n.levels = n.levels, key.cols = key.cols, key.marks = key.marks,
            key.labels = key.labels, lwd.axis = lwd.axis, xlim_vals = xlim_vals,
            ylim_vals = sort(ylim_vals), periodlab = periodlab, x_lab = x_lab, main = main, periodtck = periodtck,
            periodtcl = periodtcl, add_lines = add_lines, add_points = add_points, add_abline_h = add_abline_h,
            add_abline_v = add_abline_v, add_MTM_peaks = FALSE, add_pval = FALSE,
            pval_abline = NULL, pval_cutoff = NULL, add_signif = NULL,
            MTM_res_1 = NULL, MTM_res_2 = matrix(data = NA, ncol = 2, nrow = 1),
            plot_horizontal = plot_horizontal, add_data = add_data, add_avg = add_avg, add_MTM = FALSE,
            arrows = if (plot_arrows == TRUE) list(Power = Power, Phase = Phase, nr = Xsuperlet$nr, nc = Xsuperlet$nc,
                                                   p_thresh = quantile(Xsuperlet$Power, pwr_quant), n_arrows = n_arrows,
                                                   depth = Xsuperlet$x1, ax2 = Xsuperlet$axis.2, T_unit = T_unit, P_unit = P_unit,
                                                   length = arrow_length, lwd = arrow_lwd, angle = arrow_angle))
  .pw_draw(p)
}
