#' @title Plots an log2 spaced EHA spectra
#'
#' @description Plot wavelet scalogram using the outcome of the \code{\link{eha_log2}} function.
#'
#'@param eha_log2 eha_log2 object created using the \code{\link{eha_log2}} function.
#'@param plot_opt plot options are  "Power", "Amplitude", "Probability" or "F_test"
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
#'@param add_MTM_peaks Add the MTM peak periods as horizontal lines \code{Default=FALSE}
#'@param add_data Plot the data on top of the wavelet \code{Default=TRUE}
#'@param add_avg Plot the average wavelet spectral power to the side of the wavelet \code{Default=FALSE}
#'@param add_MTM Add the MTM  plot next to the wavelet plot \code{Default=FALSE}
#'@param mtm_siglvl select the significance level  (0-1) for the MTM spectrum \code{Default=0.95}
#'@param demean_mtm Remove mean from data before conducting the MTM analysis \code{Default=TRUE}
#'@param detrend_mtm Remove mean from data before conducting the MTM analysis \code{Default=TRUE}
#'@param padfac_mtm Pad factor for the MTM analysis \code{Default=5}
#'@param tbw_mtm time bandwidth product of the MTM analysis  \code{Default=3}
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis is depth or y axis power \code{Default=TRUE}
#'@param pval_abline p value of the abline
#'@param pval_cutoff p value cutoff below which no data is displayed
#'
#' @return
#' The output is a plot of a EHA spectra.
#' if add_MTM_peaks = TRUE then the output of the MTM analysis will given as matrix
#'
#' @author
#' Code based on the wt.image" functions of the 'WaveletComp' R package
#' The EHA and MTM analysis parts are from the astrochron R package of Meyers et al., (2012)
#'
#' @references
#'Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
#'Wavelet Analysis. R package version 1.1.
#'\url{https://CRAN.R-project.org/package=WaveletComp}
#'
#'S.R. Meyers, 2012, Seeing Red in Cyclic Stratigraphy: Spectral Noise Estimation for
#'Astrochronology: Paleoceanography, 27, PA3228, <doi:10.1029/2012PA002307>
#'
#'S.R. Meyers, 2019 Cyclostratigraphy and the problem of astrochronologic testing, Earth-Science Reviews,
#'Volume 190, <doi.org/10.1016/j.earscirev.2018.11.015.>
#'
#'@examples
#' \donttest{
#'#Example 1. A plot of a wavelet spectra using the Total Solar Irradiance
#'# data set of Steinhilber et al., (2012)
#'
#'TSI_eha_log2 <- eha_log2(data = TSI,
#'win =  8192,
#'tbw = 4,
#'demean = TRUE,
#'detrend = TRUE,
#'upperPeriod = 8192,
#'lowerPeriod = 50,
#'pad = NULL,
#'padding = "noise")
#'
#'plot_eha_log2(
#'eha_log2 = TSI_eha_log2,
#'plot_opt = "Amplitude",
#'lowerPeriod = 50,
#'upperPeriod = 8192,
#'n.levels = 100,
#'palette_name = "rainbow",
#'color_brewer = "grDevices",
#'useRaster = TRUE,
#'periodlab = "Period (metres)",
#'x_lab = "depth (metres)",
#'keep_editable = FALSE,
#'dev_new = TRUE,
#'plot_dir = TRUE,
#'add_lines = NULL,
#'add_points = NULL,
#'add_abline_h = NULL,
#'add_abline_v = NULL,
#'add_MTM_peaks = FALSE,
#'add_data = TRUE,
#'add_avg = FALSE,
#'add_MTM = FALSE,
#'mtm_siglvl = 0.95,
#'demean_mtm = TRUE,
#'detrend_mtm = TRUE,
#'padfac_mtm = 5,
#'tbw_mtm = 3,
#'plot_horizontal = TRUE)
#'
#'#Example 2. A plot of a wavelet spectra using the magnetic susceptibility
#'#data set of Pas et al., (2018)
#'mag_eha_log2 <- eha_log2(data = mag,
#'win =  50,
#'tbw = 4,
#'demean = TRUE,
#'detrend = TRUE,
#'upperPeriod = 50,
#'lowerPeriod = 1,
#'pad = NULL,
#'padding = "noise")
#'
#'plot_eha_log2(
#'eha_log2 = mag_eha_log2,
#'plot_opt = "Amplitude",
#'lowerPeriod = 1,
#'upperPeriod = 50,
#'n.levels = 100,
#'palette_name = "rainbow",
#'color_brewer = "grDevices",
#'useRaster = TRUE,
#'periodlab = "Period (metres)",
#'x_lab = "depth (metres)",
#'keep_editable = FALSE,
#'dev_new = TRUE,
#'plot_dir = TRUE,
#'add_lines = NULL,
#'add_points = NULL,
#'add_abline_h = NULL,
#'add_abline_v = NULL,
#'add_MTM_peaks = FALSE,
#'add_data = TRUE,
#'add_avg = FALSE,
#'add_MTM = FALSE,
#'mtm_siglvl = 0.95,
#'demean_mtm = TRUE,
#'detrend_mtm = TRUE,
#'padfac_mtm = 5,
#'tbw_mtm = 3,
#'plot_horizontal = TRUE)
#'
#'
#'#Example 3. A plot of a wavelet spectra using the greyscale
#'# data set of Zeeden et al., (2013)
#'grey_eha_log2 <- eha_log2(data = grey,
#'win =  20,
#'tbw = 4,
#'demean = TRUE,
#'detrend = TRUE,
#'upperPeriod = 20,
#'lowerPeriod = 1,
#'pad = NULL,
#'padding = "noise")
#'
#'
#'plot_eha_log2(
#'eha_log2 = grey_eha_log2,
#'plot_opt = "Amplitude",
#'lowerPeriod = 1,
#'upperPeriod = 20,
#'n.levels = 100,
#'palette_name = "rainbow",
#'color_brewer = "grDevices",
#'useRaster = TRUE,
#'periodlab = "Period (metres)",
#'x_lab = "depth (metres)",
#'keep_editable = FALSE,
#'dev_new = TRUE,
#'plot_dir = TRUE,
#'add_lines = NULL,
#'add_points = NULL,
#'add_abline_h = NULL,
#'add_abline_v = NULL,
#'add_MTM_peaks = FALSE,
#'add_data = TRUE,
#'add_avg = FALSE,
#'add_MTM = FALSE,
#'mtm_siglvl = 0.95,
#'demean_mtm = TRUE,
#'detrend_mtm = TRUE,
#'padfac_mtm = 5,
#'tbw_mtm = 3,
#'plot_horizontal = TRUE)
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
#' @importFrom astrochron mtm
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

plot_eha_log2 <- function (eha_log2 = NULL, plot_opt = "Amplitude", lowerPeriod = NULL,
                           upperPeriod = NULL, n.levels = 100, palette_name = "rainbow",
                           color_brewer = "grDevices", useRaster = TRUE, periodlab = "Period (metres)",
                           x_lab = "depth (metres)", keep_editable = FALSE, dev_new = wr_dev_new_default(),
                           plot_dir = TRUE, add_lines = NULL, add_points = NULL, add_abline_h = NULL,
                           add_abline_v = NULL, add_MTM_peaks = FALSE, add_data = TRUE,
                           add_avg = FALSE, pval_abline = c(0.1, 0.05), pval_cutoff = c(0.1),
                           add_MTM = FALSE, mtm_siglvl = 0.95, demean_mtm = TRUE, detrend_mtm = TRUE,
                           padfac_mtm = 5, tbw_mtm = 3, plot_horizontal = TRUE,
  colour_brewer = NULL, scale = "auto")
{
  if (!is.null(colour_brewer)) color_brewer <- colour_brewer  # British spelling synonym

  if (plot_opt == "Amplitude") {
    mat <- eha_log2$amp
    average <- eha_log2$amp.avg
  }
  else if (plot_opt == "Power") {
    mat <- eha_log2$pwr
    average <- eha_log2$Power.avg
  }
  else if (plot_opt == "Probability") {
    mat <- eha_log2$prob
    average <- eha_log2$prob.avg
  }
  else if (plot_opt == "F_test") {
    mat <- eha_log2$f_test
    average <- eha_log2$f_test.avg
  }
  else {
    mat <- eha_log2$amp
    average <- eha_log2$amp.avg
  }
  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }
  maximum.level = max(mat, na.rm = TRUE)
  power_max_mat.levels = quantile(mat, na.rm = TRUE, probs = seq(from = 0,
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
  depth <- eha_log2$depth
  y_axis <- eha_log2$log2_period
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)
  MTM_res_2 <- matrix(data = NA, ncol = 2, nrow = 1)
  if (add_MTM == TRUE | add_MTM_peaks == TRUE) {
    MTM_res_1 <- mtm(cbind(eha_log2$x, eha_log2$y), tbw = tbw_mtm,
                     padfac = padfac_mtm, output = 1, siglevel = mtm_siglvl,
                     genplot = FALSE, verbose = FALSE, demean = demean_mtm,
                     detrend = detrend_mtm)
    MTM_res_2 <- mtm(cbind(eha_log2$x, eha_log2$y), tbw = tbw_mtm,
                     padfac = padfac_mtm, output = 3, siglevel = mtm_siglvl,
                     genplot = FALSE, verbose = FALSE, demean = demean_mtm,
                     detrend = detrend_mtm)
    MTM_res_1$period <- 1/MTM_res_1[, 1]
    MTM_res_1$Peak <- 0
    if (is.na(MTM_res_2[1, 2]) == FALSE) {
      for (i in 1:nrow(MTM_res_2)) {
        row_nr <- wr_closest(MTM_res_1[, 1],
                                     MTM_res_2[i, 1], which = TRUE)
        row_nr <- row_nr[1]
        MTM_res_1[row_nr, 10] <- 1
      }
    }
    mtm_res <- MTM_res_1
  }
  if (plot_dir != TRUE) {
    xlim_vals = rev(c(min(eha_log2$x), max(eha_log2$x)))
  }
  else {
    xlim_vals = c(min(eha_log2$x), max(eha_log2$x))
  }
  if (is.null(lowerPeriod) == TRUE) {
    lowerPeriod <- min(2^eha_log2$log2_period)
  }
  if (is.null(upperPeriod) == TRUE) {
    upperPeriod <- max(2^eha_log2$log2_period)
  }
  ylim_vals = c(lowerPeriod, upperPeriod)
  p <- list(wavelet = NULL, data_x = eha_log2$x, data_y = eha_log2$y, depth = eha_log2$depth, ax2 = eha_log2$log2_period,
            z = mat, avg = average, per = 2^eha_log2$log2_period, tick_min = -Inf,
            overlay = FALSE, n.levels = n.levels, key.cols = key.cols, key.marks = key.marks,
            key.labels = key.labels, lwd.axis = lwd.axis, pmax_avg_sel = mat, xlim_vals = xlim_vals,
            ylim_vals = ylim_vals, periodlab = periodlab, x_lab = x_lab, main = main, periodtck = periodtck,
            periodtcl = periodtcl, add_lines = add_lines, add_points = add_points, add_abline_h = add_abline_h,
            add_abline_v = add_abline_v, add_MTM_peaks = add_MTM_peaks, add_pval = FALSE,
            pval_abline = NULL, pval_cutoff = NULL, add_signif = NULL,
            MTM_res_1 = if (add_MTM == TRUE | add_MTM_peaks == TRUE) MTM_res_1 else NULL,
            MTM_res_2 = MTM_res_2, plot_horizontal = plot_horizontal, add_data = add_data, add_avg = add_avg,
            add_MTM = add_MTM)
  p$scale <- scale
  .pw_draw(p)
  if (add_MTM_peaks == TRUE & is.na(MTM_res_2[1, 2]) == FALSE) {
    return(invisible(mtm_res))
  }
}

