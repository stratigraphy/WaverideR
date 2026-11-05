#' @title Plots a adaptive wavelet scalogram
#'
#' @description Plot adaptive wavelet scalogram using the outcome of the \code{\link{analyse_Awavelet}} function.
#'
#' @param Awavelet Adaptive wavelet object, the output of \code{\link{analyse_Awavelet}}. Required
#' (\code{Default=NULL}).
#' @param lowerPeriod Shortest period shown on the period axis, in the depth/time unit of the record (numeric, > 0).
#' \code{NULL} uses the shortest period of \code{Awavelet}. \code{Default=NULL}.
#' @param upperPeriod Longest period shown on the period axis, in the depth/time unit of the record (numeric,
#' > \code{lowerPeriod}). \code{NULL} uses the longest period of \code{Awavelet}. \code{Default=NULL}.
#' @param n.levels Number of colour levels of the power scale (integer, > 1); the level boundaries are quantiles of
#' the power. \code{Default=100}.
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
#' @param useRaster Draw the scalogram as a raster image (logical, TRUE or FALSE). Kept for compatibility: the
#' current plotting code always draws a raster image, whatever the value. \code{Default=TRUE}.
#' @param periodlab Label of the period axis (character). \code{Default="Period (metres)"}.
#' @param x_lab Label of the depth/time axis (character). \code{Default="depth (metres)"}.
#' @param keep_editable Keep the graphical parameters (\code{par}) of the plot after the function returns, so that
#' features can be added to it (logical, TRUE or FALSE); \code{FALSE} restores the previous settings.
#' \code{Default=FALSE}.
#' @param dev_new Open a new plotting window for the plot (logical, TRUE or FALSE). The wavelet plots contain a raster
#' image that needs more room than the small RStudio plot pane, which gives "figure margins too large" errors or
#' clipped text, so a separate window is opened by default in the RStudio GUI. In every other front end (plain R, ESS,
#' Positron, knitr, Rscript) the plot is drawn on the current device and no window is opened. Set TRUE or FALSE to
#' force either behaviour, or set \code{options(waverider.dev_new = TRUE/FALSE)}.
#' \code{Default=wr_dev_new_default()} (TRUE in the RStudio GUI, FALSE elsewhere).
#' @param plot_dir Direction of the depth/time axis (logical, TRUE or FALSE): \code{TRUE} plots increasing
#' depth/time values from left to right (e.g. bore hole data that get older with depth); \code{FALSE} reverses the
#' axis (e.g. outcrop logs where 0 m is the base of the section). \code{Default=TRUE}.
#' @param add_lines Lines drawn on the scalogram (numeric matrix or data frame, or \code{NULL}): depth/time in
#' column 1 and one or more columns of period values (in the unit of the period axis), one line per column.
#' \code{Default=NULL}.
#' @param add_points Points drawn on the scalogram (numeric matrix or data frame, or \code{NULL}): depth/time in
#' column 1 and one or more columns of period values, one set of points per column. \code{Default=NULL}.
#' @param add_abline_h Positions of horizontal lines drawn on the scalogram (numeric vector, or \code{NULL}), e.g.
#' \code{c(2, 3, 5, 6)}: periods when \code{plot_horizontal = TRUE}, depths/times when \code{plot_horizontal = FALSE}.
#' \code{Default=NULL}.
#' @param add_abline_v Positions of vertical lines drawn on the scalogram (numeric vector, or \code{NULL}), e.g.
#' \code{c(2, 3, 5, 6)}: depths/times when \code{plot_horizontal = TRUE}, periods when \code{plot_horizontal = FALSE}.
#' \code{Default=NULL}.
#' @param add_MTM_peaks Draw the periods of the significant MTM spectral peaks (see \code{mtm_siglvl}) as dotted
#' lines across the scalogram, and return these peaks invisibly (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param add_data Plot the proxy record in an extra panel along the scalogram (logical, TRUE or FALSE).
#' \code{Default=TRUE}.
#' @param add_avg Plot the average wavelet spectral power beside the scalogram (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param add_MTM Plot the MTM spectrum (\code{\link[astrochron]{mtm}}) of the record beside the scalogram (logical,
#' TRUE or FALSE). \code{Default=FALSE}.
#' @param mtm_siglvl Significance level of the MTM spectrum (numeric, between 0 and 1) above which peaks are
#' significant, passed to \code{astrochron::mtm} as \code{siglevel}. \code{Default=0.95}.
#' @param demean_mtm Remove the mean of the record before the MTM analysis (logical, TRUE or FALSE).
#' \code{Default=TRUE}.
#' @param detrend_mtm Remove a linear trend from the record before the MTM analysis (logical, TRUE or FALSE).
#' \code{Default=TRUE}.
#' @param padfac_mtm Zero-padding factor of the MTM analysis, passed to \code{astrochron::mtm} as \code{padfac}
#' (numeric, >= 1). \code{Default=5}.
#' @param tbw_mtm Time-bandwidth product of the MTM analysis, passed to \code{astrochron::mtm} as \code{tbw}
#' (numeric, e.g. 2 or 3). \code{Default=3}.
#' @param plot_horizontal Orientation of the plot (logical, TRUE or FALSE): \code{TRUE} puts depth/time on the
#' horizontal axis and period on the vertical axis; \code{FALSE} puts depth/time on the vertical axis.
#' \code{Default=TRUE}.
#'
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
#'#Example 1. A plot of a wavelet spectra using the Total Solar Irradiance
#'# data set of Steinhilber et al., (2012)
#'
#'TSI_wt <-
#'  analyse_Awavelet(
#'    data = TSI,
#'    dj = 1/200,
#'    lowerPeriod = 16,
#'    upperPeriod = 8192,
#'    verbose = FALSE,
#' omega_min = 6,
#' omega_max = 12,
#' scaling = "log2",
#' alpha = 1
#'  )
#'
#'plot_Awavelet(
#'  Awavelet = TSI_wt,
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
#'  add_MTM_peaks = FALSE,
#'  add_data = TRUE,
#'  add_avg = TRUE,
#'  add_MTM = FALSE,
#'  mtm_siglvl = 0.95,
#'  demean_mtm = TRUE,
#'  detrend_mtm = TRUE,
#'  padfac_mtm = 5,
#'  tbw_mtm = 3,
#'  plot_horizontal=TRUE)
#'
#'#Example 2. A plot of a wavelet spectra using the magnetic susceptibility
#'#data set of Pas et al., (2018)
#'mag_wt <-
#'analyse_Awavelet(
#'data = mag,
#'dj = 1/100,
#'lowerPeriod = 1,
#'upperPeriod = 254,
#'verbose = FALSE,
#' omega_min = 6,
#' omega_max = 12,
#' scaling = "log2",
#' alpha = 1
#')
#'
#'plot_Awavelet(
#'Awavelet = mag_wt,
#'lowerPeriod = 1,
#'upperPeriod = 254,
#'n.levels = 100,
#' palette_name = "rainbow",
#' color_brewer= "grDevices",
#'useRaster = TRUE,
#'periodlab = "Period (metres)",
#'x_lab = "depth (metres)",
#'keep_editable = FALSE,
#'dev_new=TRUE,
#'plot_dir = TRUE,
#'add_lines= NULL,
#'add_points= NULL,
#'add_abline_h = NULL,
#'add_abline_v = NULL,
#'add_MTM_peaks = FALSE,
#'add_data = TRUE,
#'add_avg = TRUE,
#'add_MTM = FALSE,
#'mtm_siglvl = 0.95,
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
#'  analyse_Awavelet(
#'    data = grey,
#'    dj = 1/200,
#'    lowerPeriod = 2,
#'    upperPeriod = 25,
#'    verbose = FALSE,
#' omega_min = 6,
#' omega_max = 12,
#' scaling = "log2",
#' alpha = 1)
#'
#'plot_Awavelet(
#'Awavelet = grey_wt,
#'lowerPeriod = 2,
#'upperPeriod = 25,
#'n.levels = 100,
#' palette_name = "rainbow",
#' color_brewer= "grDevices",
#'useRaster = TRUE,
#'periodlab = "Period (metres)",
#'x_lab = "depth (metres)",
#'keep_editable = FALSE,
#'dev_new=TRUE,
#'plot_dir = TRUE,
#'add_lines = NULL,
#'add_points= NULL,
#'add_abline_h = NULL,
#'add_abline_v = NULL,
#'add_MTM_peaks = FALSE,
#'add_data = TRUE,
#'add_avg = TRUE,
#'add_MTM = FALSE,
#'mtm_siglvl = 0.95,
#'demean_mtm = TRUE,
#'detrend_mtm = TRUE,
#'padfac_mtm = 5,
#'tbw_mtm = 3,
#'plot_horizontal=TRUE)
#'
#'}
#'
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

plot_Awavelet <- function(Awavelet = NULL,
    lowerPeriod = NULL,
    upperPeriod = NULL,
    n.levels = 100,
    palette_name = "viridis",
    color_brewer = "viridis",
    useRaster = TRUE,
    periodlab = "Period (metres)",
    x_lab = "depth (metres)",
    keep_editable = FALSE,
    dev_new = wr_dev_new_default(),
    plot_dir = TRUE,
    add_lines = NULL,
    add_points= NULL,
    add_abline_h = NULL,
    add_abline_v = NULL,
    add_MTM_peaks = FALSE,
    add_data = TRUE,
    add_avg = FALSE,
    add_MTM = FALSE,
    mtm_siglvl = 0.95,
    demean_mtm = TRUE,
    detrend_mtm = TRUE,
    padfac_mtm = 5,
    tbw_mtm = 3,
    plot_horizontal = TRUE,
    colour_brewer = NULL) {
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym


  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }


  maximum.level = max(Awavelet$Power)
  power_max_mat.levels = quantile(Awavelet$Power, probs = seq(from = 0,
      to = 1, length.out = n.levels + 1))


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

  key.marks = round(seq(from = 0, to = 1, length.out = legend.params$n.ticks) *
      n.levels)

  key.labels = formatC(as.numeric(power_max_mat.levels), digits = legend.params$label.digits,
      format = legend.params$label.format)[key.marks + 1]

  if(dev_new==TRUE & plot_horizontal==TRUE){
    wr_open_device(width = 15, height = 7)
  }

  if(dev_new==TRUE & plot_horizontal==FALSE){
    wr_open_device(width = 7, height = 10)
  }



  # y_axis <- as.numeric(unlist(Awavelet$Period))
  pmax_avg_sel <- t(Awavelet$Power)

  depth <-  Awavelet$x
  y_axis <- Awavelet$Period
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)


  MTM_res_2 <- matrix(data=NA,ncol=2,nrow=1)

  if (add_MTM == TRUE |   add_MTM_peaks == TRUE) {
    MTM_res_1 <- mtm(
        cbind(Awavelet$x, Awavelet$y),
        tbw = tbw_mtm,
        padfac = padfac_mtm,
        output = 1,
        siglevel = mtm_siglvl,
        genplot = FALSE,
        verbose = FALSE,
        demean = demean_mtm,
        detrend = detrend_mtm
    )
    MTM_res_2 <-
        mtm(
            cbind(Awavelet$x, Awavelet$y),
            tbw = tbw_mtm,
            padfac = padfac_mtm,
            output = 3,
            siglevel = mtm_siglvl,
            genplot = FALSE,
            verbose = FALSE,
            demean = demean_mtm,
            detrend = detrend_mtm
        )
    MTM_res_1$period <- 1 / MTM_res_1[, 1]
    MTM_res_1$Peak <- 0

    if (is.na(MTM_res_2[1,2]) == FALSE){
      for (i  in 1:nrow(MTM_res_2)) {
        row_nr <- wr_closest(MTM_res_1[, 1], MTM_res_2[i, 1], which = TRUE)
        row_nr <- row_nr[1]
        MTM_res_1[row_nr, 10] <- 1
      }
    }


    mtm_res <- MTM_res_1

  }


  if (plot_dir != TRUE) {
    xlim_vals = rev(c(min(Awavelet$x), max(Awavelet$x)))
  } else{
    xlim_vals = c(min(Awavelet$x), max(Awavelet$x))
  }

  if (is.null(lowerPeriod) == TRUE) {
    lowerPeriod  <- min(Awavelet$Period)
  }
  if (is.null(upperPeriod) == TRUE) {
    upperPeriod  <- max(Awavelet$Period)
  }

  ylim_vals = c(lowerPeriod, upperPeriod)





  p <- list(wavelet = NULL, data_x = Awavelet$x, data_y = Awavelet$y, depth = Awavelet$x, ax2 = Awavelet$axis.2,
      z = Awavelet$Power, avg = Awavelet$Power.avg, per = Awavelet$Period, tick_min = log2(Awavelet$Period[1]),
      overlay = FALSE, n.levels = n.levels, key.cols = key.cols, key.marks = key.marks,
      key.labels = key.labels, lwd.axis = lwd.axis, pmax_avg_sel = pmax_avg_sel, xlim_vals = xlim_vals,
      ylim_vals = ylim_vals, periodlab = periodlab, x_lab = x_lab, main = main, periodtck = periodtck,
      periodtcl = periodtcl, add_lines = add_lines, add_points = add_points, add_abline_h = add_abline_h,
      add_abline_v = add_abline_v, add_MTM_peaks = add_MTM_peaks, add_pval = FALSE,
      pval_abline = NULL, pval_cutoff = NULL, add_signif = NULL,
      MTM_res_1 = if (add_MTM == TRUE | add_MTM_peaks == TRUE) { MTM_res_1 } else { NULL },
      MTM_res_2 = MTM_res_2, plot_horizontal = plot_horizontal, add_data = add_data, add_avg = add_avg,
      add_MTM = add_MTM)
  .pw_draw(p)

  if (add_MTM_peaks == TRUE & is.na(MTM_res_2[1,2]) == FALSE) {
    return(invisible(mtm_res))
  }
}

