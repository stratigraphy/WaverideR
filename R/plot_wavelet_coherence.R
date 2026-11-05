#' @title Plots a cross wavelet scalogram
#'
#' @description Plot cross wavelet scalogram using the outcome of the \code{\link{analyse_wavelet_coherence}} function.
#'
#'@param wavelet_coh Wavelet coherence object to plot (required): the list returned by
#'\code{\link{analyse_wavelet_coherence}}. Its \code{Coherence}, \code{Phase}, \code{Period}, \code{x1}/\code{y1} and
#'\code{x2}/\code{y2} elements are used. \code{Default=NULL}.
#'@param lowerPeriod Shortest period shown on the period axis, in the depth/time unit of the records (numeric, > 0).
#'\code{NULL} uses the shortest period of \code{wavelet_coh}. \code{Default=NULL}.
#'@param upperPeriod Longest period shown on the period axis, in the depth/time unit of the records (numeric,
#'> \code{lowerPeriod}). \code{NULL} uses the longest period of \code{wavelet_coh}. \code{Default=NULL}.
#'@param n.levels Number of colour levels of the coherence scale (integer, > 1); the breaks are quantiles of the
#'coherence. \code{Default=100}.
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
#'@param useRaster Draw the scalogram as a raster image (\code{TRUE}) or as a vector image of individual cells
#'(\code{FALSE}) (logical). Vector images are slow to draw and give large files. \code{Default=TRUE}.
#'@param periodlab Label of the period axis (character). \code{Default="Period (metres)"}.
#'@param x_lab Label of the depth/time axis (character). \code{Default="depth (metres)"}.
#'@param keep_editable Keep the graphical parameters (\code{par}, layout) of the plot after the function returns, so that
#'extra elements can be added to the plot afterwards (logical). With \code{FALSE} the previous graphical parameters
#'are restored on exit. \code{Default=FALSE}.
#'@param dev_new Open a new plotting window for the plot. The wavelet plots contain a raster image that needs more room
#'than the small RStudio plot pane, which gives "figure margins too large" errors or clipped text, so a separate window is
#'opened by default in the RStudio GUI. In every other front end (plain R, ESS, Positron, knitr, Rscript) the plot is drawn
#'on the current device and no window is opened. Set TRUE or FALSE to force either behaviour, or set
#'\code{options(waverider.dev_new = TRUE/FALSE)}. Logical.
#'\code{Default=wr_dev_new_default()}: TRUE in the RStudio GUI, FALSE elsewhere.
#'@param plot_dir Direction of the depth/time axis (logical): \code{TRUE} draws it with values increasing from left to
#'right (or top to bottom in a vertical plot), \code{FALSE} reverses the axis. \code{Default=TRUE}.
#'@param plot_arrows Draw phase arrows on the scalogram (logical). Arrows are drawn at the local coherence maxima along
#'the period axis whose coherence is above the \code{pwr_quant} quantile; their direction gives the phase difference
#'between the two records (pointing right: in phase, left: anti-phase). \code{Default=TRUE}.
#'@param pwr_quant Quantile of the coherence (numeric, between 0 and 1) above which phase arrows are drawn; higher values
#'give fewer arrows. \code{Default=0.85}.
#'@param n_arrows Approximate number of arrow positions along the depth/time axis (integer, > 0); the columns of the
#'scalogram are thinned to about this number. \code{Default=75}.
#'@param T_unit Horizontal size of the arrows in depth/time units (numeric). \code{NULL} (or \code{P_unit = NULL}) uses an
#'automatic, fixed visual length of 0.3 inch whatever the record length. \code{Default=NULL}.
#'@param P_unit Vertical size of the arrows in log2 period units (numeric). \code{NULL} (or \code{T_unit = NULL}) uses the
#'automatic size (see \code{T_unit}). \code{Default=NULL}.
#'@param arrow_length Length of the edges of the arrow head, in inches (numeric; \code{length} of
#'\code{\link[graphics]{arrows}}). \code{Default=0.075}.
#'@param arrow_lwd Line width of the arrows (numeric). \code{Default=1.9}.
#'@param arrow_angle Angle between the shaft and the edges of the arrow head, in degrees (numeric; \code{angle} of
#'\code{\link[graphics]{arrows}}). \code{Default=25}.
#'@param add_lines Lines to draw on the scalogram, e.g. tracked periods (matrix or data frame, or \code{NULL} for none):
#'depth/time in column 1 and one or more columns of period values after it; each period column is drawn as a separate
#'line. \code{Default=NULL}.
#'@param add_points Points to draw on the scalogram (matrix or data frame, or \code{NULL} for none): depth/time in
#'column 1 and one or more columns of period values after it. \code{Default=NULL}.
#'@param add_abline_h Positions of extra horizontal lines (numeric vector, e.g. \code{c(2, 3, 5)}, or \code{NULL} for none).
#'In a horizontal plot these are periods; in a vertical plot (\code{plot_horizontal = FALSE}) they are depth/time
#'values. \code{Default=NULL}.
#'@param add_abline_v Positions of extra vertical lines (numeric vector or \code{NULL} for none). In a horizontal plot
#'these are depth/time values; in a vertical plot they are periods. \code{Default=NULL}.
#'@param add_data Draw the two data series in a panel next to the scalogram, the second one in red (logical).
#'\code{Default=TRUE}.
#'@param add_avg Draw the average coherence per period (\code{Coh.avg}) in a panel to the side of the scalogram
#'(logical). \code{Default=FALSE}.
#'@param plot_horizontal Orientation of the plot (logical): \code{TRUE} puts depth/time on the x-axis and period on the
#'y-axis; \code{FALSE} draws a vertical plot with depth/time on the y-axis. \code{Default=TRUE}.
#'
#'
#' @return
#' The output is a plot of a cross wavelet scalogram.
#'
#' @author
#' plotting code based on the "wt.image" and "analyze.coherency" functions
#'  of the 'WaveletComp' R package
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
#'  tmax = 1500,
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
#'  tmax = 1500,
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
#'coh_etp <- analyse_wavelet_coherence(
#'  data_1 = etp_1,
#'  data_2  = etp_2,
#'  upperPeriod = 1024,
#'  lowerPeriod = 2
#')
#'plot_wavelet_coherence(wavelet_coh = coh_etp, lowerPeriod = 2, upperPeriod = 1024,
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

plot_wavelet_coherence<- function (wavelet_coh = NULL, lowerPeriod = NULL, upperPeriod = NULL,
    n.levels = 100, palette_name = "viridis", color_brewer = "viridis",
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
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym

  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  Coherence <- t(wavelet_coh$Coherence)
  Power <-  Coherence
  Phase <- t(wavelet_coh$Phase)
  maximum.level = max(wavelet_coh$Coherence, na.rm = TRUE)
  power_max_mat.levels = quantile(wavelet_coh$Coherence, probs = seq(from = 0,
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
  depth <- wavelet_coh$x1
  y_axis <- wavelet_coh$Period
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)
  if (plot_dir != TRUE) {
    xlim_vals = rev(c(min(wavelet_coh$x1), max(wavelet_coh$x1)))
  } else {
    xlim_vals = c(min(wavelet_coh$x1), max(wavelet_coh$x1))
  }
  if (is.null(lowerPeriod) == TRUE) {
    lowerPeriod <- min(wavelet_coh$Period)
  }
  if (is.null(upperPeriod) == TRUE) {
    upperPeriod <- max(wavelet_coh$Period)
  }
  ylim_vals = (c(lowerPeriod, upperPeriod))
  p <- list(wavelet = wavelet_coh, data_x = wavelet_coh$x1, data_y = wavelet_coh$y1, data2 = list(x = wavelet_coh$x2, y = wavelet_coh$y2),
      depth = wavelet_coh$x1, ax2 = wavelet_coh$axis.2, z = Power, avg = wavelet_coh$Coh.avg, per = wavelet_coh$Period,
      tick_min = min(log2(ylim_vals)), tick_max = max(log2(ylim_vals)), overlay = FALSE,
      gws_abline_all = TRUE, levels = power_max_mat.levels,
      n.levels = n.levels, key.cols = key.cols, key.marks = key.marks,
      key.labels = key.labels, lwd.axis = lwd.axis, xlim_vals = xlim_vals,
      ylim_vals = sort(ylim_vals), periodlab = periodlab, x_lab = x_lab, main = main, periodtck = periodtck,
      periodtcl = periodtcl, add_lines = add_lines, add_points = add_points, add_abline_h = add_abline_h,
      add_abline_v = add_abline_v, add_MTM_peaks = FALSE, add_pval = FALSE,
      pval_abline = NULL, pval_cutoff = NULL, add_signif = NULL,
      MTM_res_1 = NULL, MTM_res_2 = matrix(data = NA, ncol = 2, nrow = 1),
      plot_horizontal = plot_horizontal, add_data = add_data, add_avg = add_avg, add_MTM = FALSE,
      arrows = if (plot_arrows == TRUE) {
        list(Power = Power, Phase = Phase, nr = wavelet_coh$nr, nc = wavelet_coh$nc,
            p_thresh = quantile(wavelet_coh$Coherence, pwr_quant), n_arrows = n_arrows,
            depth = wavelet_coh$x1, ax2 = wavelet_coh$axis.2, T_unit = T_unit, P_unit = P_unit,
            length = arrow_length, lwd = arrow_lwd, angle = arrow_angle)
      })
  .pw_draw(p)
}
