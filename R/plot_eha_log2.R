#' @title Plots an log2 spaced EHA spectra
#'
#' @description Plot wavelet scalogram using the outcome of the \code{\link{eha_log2}} function.
#'
#' @param eha_log2 Evolutive harmonic analysis result to plot: the list returned by \code{\link{eha_log2}}. Required
#' (\code{Default=NULL} gives an error).
#' @param plot_opt Quantity shown in the spectrogram and in the average panel (character), one of
#' \code{"Amplitude"}, \code{"Power"}, \code{"Probability"} (harmonic F-test probability) or \code{"F_test"}
#' (F-test value); any other value plots the amplitude. \code{Default="Amplitude"}.
#' @param lowerPeriod Shortest period shown on the period axis, in the depth/time unit of the input record (numeric,
#' > 0). \code{NULL} uses the shortest period of \code{eha_log2}. \code{Default=NULL}.
#' @param upperPeriod Longest period shown on the period axis, in the depth/time unit of the input record (numeric,
#' > \code{lowerPeriod}). \code{NULL} uses the longest period of \code{eha_log2}. \code{Default=NULL}.
#' @param n.levels Number of colour levels of the colour scale (integer, >= 2). \code{Default=100}.
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
#' @param useRaster Draw the spectrogram as a raster image (\code{TRUE}) or as vector polygons (\code{FALSE}) (logical).
#' WARNING: plotting as a vector image is computationally intensive and gives large files. \code{Default=TRUE}.
#' @param periodlab Label of the period axis (character). \code{Default="Period (metres)"}.
#' @param x_lab Label of the depth/time axis (character). \code{Default="depth (metres)"}.
#' @param keep_editable Keep the plot layout and graphical parameters (\code{par}) after plotting, so that extra
#' elements can be added to the spectrogram panel with \code{lines}, \code{points} etc. (logical). With \code{FALSE} the
#' previous \code{par} settings are restored on exit. \code{Default=FALSE}.
#' @param dev_new Open a new plotting window for the plot (logical, TRUE or FALSE). The wavelet plots contain a raster
#' image that needs more room than the small RStudio plot pane, which gives "figure margins too large" errors or clipped
#' text, so a separate window is opened by default in the RStudio GUI. In every other front end (plain R, ESS,
#' Positron, knitr, Rscript) the plot is drawn on the current device and no window is opened. Set TRUE or FALSE to
#' force either behaviour, or set \code{options(waverider.dev_new = TRUE/FALSE)}.
#' \code{Default=wr_dev_new_default()}: TRUE in the RStudio GUI, FALSE elsewhere.
#' @param plot_dir Direction of the depth/time axis (logical). \code{TRUE} plots increasing depth/time values from left
#' to right (or bottom to top when \code{plot_horizontal=FALSE}); \code{FALSE} reverses the axis. \code{Default=TRUE}.
#' @param add_lines Lines to draw on the spectrogram: a matrix or data frame with depth/time in column 1 and, in each
#' further column, period values (same unit as the period axis) of one line; \code{NULL} draws none.
#' \code{Default=NULL}.
#' @param add_points Points to draw on the spectrogram: a matrix or data frame with depth/time in column 1 and, in each
#' further column, period values of one series of points; \code{NULL} draws none. \code{Default=NULL}.
#' @param add_abline_h Numeric vector of positions of horizontal lines, e.g. \code{c(2, 3, 5, 6)}: periods when
#' \code{plot_horizontal=TRUE}, depth/time values when \code{plot_horizontal=FALSE}; \code{NULL} draws none.
#' \code{Default=NULL}.
#' @param add_abline_v Numeric vector of positions of vertical lines, e.g. \code{c(2, 3, 5, 6)}: depth/time values
#' when \code{plot_horizontal=TRUE} (also drawn on the data panel), periods when \code{plot_horizontal=FALSE};
#' \code{NULL} draws none. \code{Default=NULL}.
#' @param add_MTM_peaks Run a multitaper (MTM) analysis of the whole record (\code{\link[astrochron]{mtm}}) and draw
#' the periods of its significant peaks (at \code{mtm_siglvl}) as dotted lines on the spectrogram (and in red on the
#' MTM panels); the MTM result with a \code{Peak} column (1 = significant peak) is then returned invisibly
#' (logical). \code{Default=FALSE}.
#' @param add_data Plot the input record in a panel beside the spectrogram (logical). \code{Default=TRUE}.
#' @param add_avg Plot the average of the plotted quantity (\code{plot_opt}) over the record in a panel beside the
#' spectrogram (logical). \code{Default=FALSE}.
#' @param pval_abline Not used; kept for compatibility with older versions (numeric vector of p-values).
#' \code{Default=c(0.1, 0.05)}.
#' @param pval_cutoff Not used; kept for compatibility with older versions (numeric p-value).
#' \code{Default=c(0.1)}.
#' @param add_MTM Add three panels beside the spectrogram with the MTM power spectrum of the whole record (with its
#' AR(1) confidence levels), the AR(1) confidence level and the harmonic F-test confidence level (in \%), computed with
#' \code{\link[astrochron]{mtm}} (logical). \code{Default=FALSE}.
#' @param mtm_siglvl Significance level used to select the MTM peaks (numeric, between 0 and 1; \code{siglevel} of
#' \code{\link[astrochron]{mtm}}). \code{Default=0.95}.
#' @param demean_mtm Remove the mean from the record before the MTM analysis (logical; \code{demean} of
#' \code{\link[astrochron]{mtm}}). \code{Default=TRUE}.
#' @param detrend_mtm Remove a linear trend from the record before the MTM analysis (logical; \code{detrend} of
#' \code{\link[astrochron]{mtm}}). \code{Default=TRUE}.
#' @param padfac_mtm Zero-padding factor of the MTM analysis (numeric, >= 1; \code{padfac} of
#' \code{\link[astrochron]{mtm}}). \code{Default=5}.
#' @param tbw_mtm Time-bandwidth product of the MTM analysis (numeric, e.g. 2, 3 or 4; \code{tbw} of
#' \code{\link[astrochron]{mtm}}); higher values give smoother spectra with lower frequency resolution.
#' \code{Default=3}.
#' @param plot_horizontal Orientation of the plot (logical): \code{TRUE} puts depth/time on the x axis and period on
#' the y axis; \code{FALSE} puts depth/time on the y axis and period on the x axis. \code{Default=TRUE}.
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
#' @param scale Spacing of the colour breaks of the colour scale, as in all image plots of the package (character):
#' \code{"quantile"} (histogram-equalised), \code{"log"}, \code{"sqrt"}, \code{"linear"} or \code{"auto"}
#' (logarithmic for mostly empty maps or values spanning more than three orders of magnitude, quantiles otherwise);
#' see \code{\link{plot_wavelet}}. \code{Default="auto"}.
#' @param colour_brewer British spelling of \code{color_brewer}, with the same options; if given (not \code{NULL}) it overrides
#' \code{color_brewer}. \code{Default=NULL}.
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
    upperPeriod = NULL, n.levels = 100, palette_name = "viridis",
    color_brewer = "viridis", useRaster = TRUE, periodlab = "Period (metres)",
    x_lab = "depth (metres)", keep_editable = FALSE, dev_new = wr_dev_new_default(),
    plot_dir = TRUE, add_lines = NULL, add_points = NULL, add_abline_h = NULL,
    add_abline_v = NULL, add_MTM_peaks = FALSE, add_data = TRUE,
    add_avg = FALSE, pval_abline = c(0.1, 0.05), pval_cutoff = c(0.1),
    add_MTM = FALSE, mtm_siglvl = 0.95, demean_mtm = TRUE, detrend_mtm = TRUE,
    padfac_mtm = 5, tbw_mtm = 3, plot_horizontal = TRUE,
    colour_brewer = NULL, scale = "auto")
{
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym

  if (plot_opt == "Amplitude") {
    mat <- eha_log2$amp
    average <- eha_log2$amp.avg
  } else if (plot_opt == "Power") {
    mat <- eha_log2$pwr
    average <- eha_log2$Power.avg
  } else if (plot_opt == "Probability") {
    mat <- eha_log2$prob
    average <- eha_log2$prob.avg
  } else if (plot_opt == "F_test") {
    mat <- eha_log2$f_test
    average <- eha_log2$f_test.avg
  } else {
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
  } else {
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
      MTM_res_1 = if (add_MTM == TRUE | add_MTM_peaks == TRUE) { MTM_res_1 } else { NULL },
      MTM_res_2 = MTM_res_2, plot_horizontal = plot_horizontal, add_data = add_data, add_avg = add_avg,
      add_MTM = add_MTM)
  p$scale <- scale
  .pw_draw(p)
  if (add_MTM_peaks == TRUE & is.na(MTM_res_2[1, 2]) == FALSE) {
    return(invisible(mtm_res))
  }
}
