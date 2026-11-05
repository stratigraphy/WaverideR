#' @title Plots a cross superlet scalogram
#'
#' @description Plot cross superlet scalogram using the outcome of the \code{\link{analyse_Xsuperlet}} function.
#'
#' @param Xsuperlet Cross superlet result to plot: the list returned by \code{\link{analyse_Xsuperlet}}
#' (computed from the two records \code{data_1} and \code{data_2}). Required (\code{Default=NULL} gives an error).
#' @param lowerPeriod Shortest period shown on the period axis, in the depth/time unit of the input records (numeric,
#' > 0). \code{NULL} uses the shortest period of \code{Xsuperlet}. \code{Default=NULL}.
#' @param upperPeriod Longest period shown on the period axis, in the depth/time unit of the input records (numeric,
#' > \code{lowerPeriod}). \code{NULL} uses the longest period of \code{Xsuperlet}. \code{Default=NULL}.
#' @param n.levels Number of colour levels of the power scale (integer, >= 2); the level breaks are the quantiles of
#' the cross power, so each colour covers the same share of the scalogram. \code{Default=100}.
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
#' @param useRaster Draw the scalogram as a raster image (\code{TRUE}) or as vector polygons (\code{FALSE}) (logical).
#' WARNING: plotting as a vector image is computationally intensive and gives large files. \code{Default=TRUE}.
#' @param periodlab Label of the period axis (character). \code{Default="Period (metres)"}.
#' @param x_lab Label of the depth/time axis (character). \code{Default="depth (metres)"}.
#' @param keep_editable Keep the plot layout and graphical parameters (\code{par}) after plotting, so that extra
#' elements can be added to the scalogram panel with \code{lines}, \code{points} etc. (logical). With \code{FALSE} the
#' previous \code{par} settings are restored on exit. \code{Default=FALSE}.
#' @param dev_new Open a new plotting window for the plot (logical, TRUE or FALSE). The wavelet plots contain a raster
#' image that needs more room than the small RStudio plot pane, which gives "figure margins too large" errors or clipped
#' text, so a separate window is opened by default in the RStudio GUI. In every other front end (plain R, ESS,
#' Positron, knitr, Rscript) the plot is drawn on the current device and no window is opened. Set TRUE or FALSE to
#' force either behaviour, or set \code{options(waverider.dev_new = TRUE/FALSE)}.
#' \code{Default=wr_dev_new_default()}: TRUE in the RStudio GUI, FALSE elsewhere.
#' @param plot_dir Direction of the depth/time axis (logical). \code{TRUE} plots increasing depth/time values from left
#' to right (or bottom to top when \code{plot_horizontal=FALSE}); \code{FALSE} reverses the axis. \code{Default=TRUE}.
#' @param plot_arrows Draw the phase arrows of the cross superlet transform on the scalogram (logical). An arrow is
#' drawn at each local maximum of cross power along the period axis whose power exceeds the \code{pwr_quant} quantile;
#' its direction gives the phase difference between the two records (radians, measured anticlockwise from the
#' depth/time axis). \code{Default=TRUE}.
#' @param pwr_quant Quantile of the cross power (numeric, between 0 and 1) above which phase arrows are drawn; higher
#' values keep only the arrows on the strongest power ridges. \code{Default=0.85}.
#' @param n_arrows Approximate number of arrow positions along the depth/time axis (integer, >= 1); every
#' \code{round(n/n_arrows)}-th depth/time step of the \code{n} steps gets arrows. \code{Default=75}.
#' @param T_unit Length of the arrows along the depth/time axis, in depth/time units (numeric, > 0). \code{NULL} (or
#' \code{P_unit=NULL}) gives arrows of a fixed visual length of 0.3 inch whatever the record length.
#' \code{Default=NULL}.
#' @param P_unit Length of the arrows along the period axis, in log2(period) units (numeric, > 0). \code{NULL} (or
#' \code{T_unit=NULL}) uses the automatic arrow length (see \code{T_unit}). \code{Default=NULL}.
#' @param arrow_length Length of the arrow head edges in inches (numeric, > 0; \code{length} argument of
#' \code{\link[graphics]{arrows}}). \code{Default=0.075}.
#' @param arrow_lwd Line width of the arrows (numeric, > 0). \code{Default=1.9}.
#' @param arrow_angle Angle between the arrow shaft and the edges of the arrow head, in degrees (numeric, 0-90).
#' \code{Default=25}.
#' @param add_lines Lines to draw on the scalogram: a matrix or data frame with depth/time in column 1 and, in each
#' further column, period values (same unit as the period axis) of one line; \code{NULL} draws none.
#' \code{Default=NULL}.
#' @param add_points Points to draw on the scalogram: a matrix or data frame with depth/time in column 1 and, in each
#' further column, period values of one series of points; \code{NULL} draws none. \code{Default=NULL}.
#' @param add_abline_h Numeric vector of positions of horizontal lines, e.g. \code{c(2, 3, 5, 6)}: periods when
#' \code{plot_horizontal=TRUE}, depth/time values when \code{plot_horizontal=FALSE}; \code{NULL} draws none.
#' \code{Default=NULL}.
#' @param add_abline_v Numeric vector of positions of vertical lines, e.g. \code{c(2, 3, 5, 6)}: depth/time values
#' when \code{plot_horizontal=TRUE} (also drawn on the data panel), periods when \code{plot_horizontal=FALSE};
#' \code{NULL} draws none. \code{Default=NULL}.
#' @param add_data Plot the two input records (record 1 in black, record 2 in red) in a panel beside the scalogram
#' (logical). \code{Default=TRUE}.
#' @param add_avg Plot the time-averaged cross power spectrum in a panel beside the scalogram (logical).
#' \code{Default=FALSE}.
#' @param plot_horizontal Orientation of the plot (logical): \code{TRUE} puts depth/time on the x axis and period on
#' the y axis; \code{FALSE} puts depth/time on the y axis and period on the x axis. \code{Default=TRUE}.
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


plot_Xsuperlet <- function (Xsuperlet = NULL, lowerPeriod = NULL, upperPeriod = NULL,
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
  } else {
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
      arrows = if (plot_arrows == TRUE) {
        list(Power = Power, Phase = Phase, nr = Xsuperlet$nr, nc = Xsuperlet$nc,
            p_thresh = quantile(Xsuperlet$Power, pwr_quant), n_arrows = n_arrows,
            depth = Xsuperlet$x1, ax2 = Xsuperlet$axis.2, T_unit = T_unit, P_unit = P_unit,
            length = arrow_length, lwd = arrow_lwd, angle = arrow_angle)
      })
  .pw_draw(p)
}
