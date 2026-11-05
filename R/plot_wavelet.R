#' @title Plots a wavelet scalogram
#'
#' @description Plot wavelet scalogram using the outcome of the \code{\link{analyse_wavelet}} function.
#'
#'@param wavelet Wavelet object to plot (required): the list returned by \code{\link{analyse_wavelet}} (or by one of
#'the other transforms of the package with the same structure, e.g. a synchrosqueezed transform). Its
#'\code{Power}, \code{Period}, \code{x} and \code{y} elements are used. \code{Default=NULL}.
#'@param lowerPeriod Shortest period shown on the period axis, in the depth/time unit of the analysed record
#'(numeric, > 0). \code{NULL} uses the shortest period of \code{wavelet}. \code{Default=NULL}.
#'@param upperPeriod Longest period shown on the period axis, in the depth/time unit of the analysed record
#'(numeric, > \code{lowerPeriod}). \code{NULL} uses the longest period of \code{wavelet}. \code{Default=NULL}.
#'@param n.levels Number of colour levels of the power scale (integer, > 1); more levels give a smoother colour
#'gradient. \code{Default=100}.
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
#'right (or top to bottom in a vertical plot), \code{FALSE} reverses the axis (e.g. for a section where 0 m is the base
#'and the record gets younger upwards). \code{Default=TRUE}.
#'@param add_lines Lines to draw on the scalogram, e.g. tracked periods (matrix or data frame, or \code{NULL} for none):
#'depth/time in column 1 and one or more columns of period values (same units as the wavelet) after it; each period
#'column is drawn as a separate line. \code{Default=NULL}.
#'@param add_points Points to draw on the scalogram (matrix or data frame, or \code{NULL} for none): depth/time in
#'column 1 and one or more columns of period values after it; each period column is drawn as a separate set of
#'points. \code{Default=NULL}.
#'@param add_abline_h Positions of extra horizontal lines (numeric vector, e.g. \code{c(2, 3, 5)}, or \code{NULL} for none).
#'In a horizontal plot these are periods; in a vertical plot (\code{plot_horizontal = FALSE}) they are depth/time
#'values. \code{Default=NULL}.
#'@param add_abline_v Positions of extra vertical lines (numeric vector or \code{NULL} for none). In a horizontal plot
#'these are depth/time values; in a vertical plot they are periods. \code{Default=NULL}.
#'@param add_MTM_peaks Run an MTM analysis (\code{astrochron::mtm}) of the data and draw the periods of its significant
#'peaks as dotted lines on the scalogram and the side panels (logical); the MTM result with a \code{Peak} column is
#'then returned invisibly. \code{Default=FALSE}.
#'@param add_data Draw the data series (\code{wavelet$x}, \code{wavelet$y}) in a panel next to the scalogram (logical).
#'\code{Default=TRUE}.
#'@param add_avg Draw the average (global) wavelet power spectrum in a panel to the side of the scalogram (logical).
#'\code{Default=FALSE}.
#'@param add_pval Use the Monte Carlo p-values stored in \code{wavelet} (\code{Power.pval}, \code{Power.avg.pval}; only
#'present when \code{\link{analyse_wavelet}} was run with \code{pval = TRUE}) (logical). The cells of the scalogram
#'with a p-value \eqn{\ge} \code{pval_cutoff} are covered with a semi-transparent white veil, and the p-value curve of the
#'average power is drawn in red on the average power panel (needs \code{add_avg = TRUE}). Without p-values in
#'\code{wavelet} the cone of influence is drawn instead. The p-values come from surrogate series that combine red noise
#'(with a windowed lag-1 autocorrelation, window half the record length) and a shuffled copy of the data; see
#'\code{\link{analyse_wavelet}}. \code{Default=FALSE}.
#'@param pval_abline P-values at which grey dashed reference lines are drawn on the p-value axis of the average power
#'panel (numeric vector, values between 0 and 1); used only with \code{add_pval = TRUE}. \code{Default=c(0.1,0.05)}.
#'@param pval_cutoff P-value below which a cell of the scalogram counts as significant and is left unveiled (numeric,
#'between 0 and 1); used only with \code{add_pval = TRUE}. \code{Default=c(0.1)}.
#'@param add_MTM Draw a multitaper (MTM) power spectrum of the data, computed with \code{astrochron::mtm}, in a panel next
#'to the scalogram (logical). \code{Default=FALSE}.
#'@param mtm_siglvl Significance level (confidence, numeric between 0 and 1) used by \code{astrochron::mtm}
#'(\code{siglevel}) to identify the significant MTM peaks. \code{Default=0.95}.
#'@param demean_mtm Remove the mean of the data before the MTM analysis (logical; \code{demean} of \code{astrochron::mtm}).
#'\code{Default=TRUE}.
#'@param detrend_mtm Remove a linear trend from the data before the MTM analysis (logical; \code{detrend} of
#'\code{astrochron::mtm}). \code{Default=TRUE}.
#'@param padfac_mtm Zero-padding factor of the MTM analysis (numeric, >= 1; \code{padfac} of \code{astrochron::mtm});
#'larger values give a finer frequency grid. \code{Default=5}.
#'@param tbw_mtm Time-bandwidth product of the MTM analysis (numeric, typically 2 to 4; \code{tbw} of
#'\code{astrochron::mtm}); larger values give smoother spectra with lower frequency resolution. \code{Default=3}.
#' @param add_signif Significance object created with \code{\link{wavelet_significance}} (list, or \code{NULL} for none).
#' The significant cells are outlined on the scalogram and the critical level of the global spectrum is drawn as a red
#' dashed line on the average power panel. \code{Default=NULL}.
#'@param plot_horizontal Orientation of the plot (logical): \code{TRUE} puts depth/time on the x-axis and period on the
#'y-axis; \code{FALSE} draws a vertical plot with depth/time on the y-axis and period on the x-axis. \code{Default=TRUE}.
#'
#'
#' @return
#' The output is a plot of a wavelet spectra.
#' if add_MTM_peaks = TRUE then the output of the MTM analysis will given as matrix
#'
#' @author
#' Code based on the "analyze.wavelet" and "wt.image" functions of the 'WaveletComp' R package
#' and "wt" function of the 'biwavelet' R package which are based on the
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
#'Wave propagation and sampling theory—Part I: Complex signal and scattering in multilayered media.
#'Geophysics 47, no. 2 (1982): 203-221.
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#'  Geophysics 1982 47 (2): 222–236.
#'
#'S.R. Meyers, 2012, Seeing Red in Cyclic Stratigraphy: Spectral Noise Estimation for
#'Astrochronology: Paleoceanography, 27, PA3228, <doi:10.1029/2012PA002307>
#'
#'@examples
#' \donttest{
#'#Example 1. A plot of a wavelet spectra using the Total Solar Irradiance
#'# data set of Steinhilber et al., (2012)
#'
#'TSI_wt <-
#'  analyse_wavelet(
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
#'  add_pval = FALSE,
#'  pval_abline = c(0.1,0.05),
#'  pval_cutoff = c(0.1),
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
#'analyse_wavelet(
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
#'add_pval = FALSE,
#'pval_abline = c(0.1,0.05),
#'pval_cutoff = c(0.1),
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
#'  analyse_wavelet(
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
#'add_pval = FALSE,
#'pval_abline = c(0.1,0.05),
#'pval_cutoff = c(0.1),
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
#' @param colour_brewer British spelling of \code{color_brewer}, with the same options; if given (not \code{NULL}) it overrides
#' \code{color_brewer}. \code{Default=NULL}.
#' @param scale Spacing of the colour breaks of the power scale (character or \code{NULL}): \code{"quantile"}
#' (histogram-equalised quantiles), \code{"log"} (logarithmic; falls back to quantiles for signed values),
#' \code{"sqrt"}, \code{"linear"} or \code{"auto"} (logarithmic for mostly empty maps or power spanning more than three
#' orders of magnitude, quantiles otherwise). \code{NULL} uses \code{"quantile"} for a standard CWT and \code{"auto"} for a
#' synchrosqueezed transform. \code{Default=NULL}.
#' @param whiten Divide the power of every period by its mean along the record before plotting (logical).
#' Removes the red-noise background, so that weak ridges at short periods (precession in a noisy record) are shown with
#' the same contrast as the strong long periods; on CIP Case 3 (very red noise) the oscillating precession and obliquity
#' ridges only become visible this way. The average panel is not informative with \code{whiten = TRUE}.
#' \code{Default=FALSE}.
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

plot_wavelet <- function(wavelet = NULL,
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
    add_pval = FALSE,
    pval_abline = c(0.1,0.05),
    pval_cutoff = c(0.1),
    add_MTM = FALSE,
    mtm_siglvl = 0.95,
    demean_mtm = TRUE,
    detrend_mtm = TRUE,
    padfac_mtm = 5,
    tbw_mtm = 3,
    plot_horizontal = TRUE,
    add_signif = NULL,
    colour_brewer = NULL,
    whiten = FALSE, scale = NULL) {
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym
  if (isTRUE(whiten)) {
    # divide every period row by its mean power: removes the red-noise background so that ridges at short and long
    # periods are shown with comparable contrast (the average panel is then flat and is not informative)
    rm <- rowMeans(wavelet$Power, na.rm = TRUE)
    rm[!is.finite(rm) | rm <= 0] <- 1
    wavelet$Power <- wavelet$Power / rm
    wavelet$Power.avg <- rowMeans(wavelet$Power, na.rm = TRUE)
    if (!is.null(wavelet$Ampl)) {
      wavelet$Ampl <- wavelet$Ampl / sqrt(rm)
    }
  }


  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }


  maximum.level = max(wavelet$Power)
  power_max_mat.levels = quantile(wavelet$Power, probs = seq(from = 0,
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



  # y_axis <- as.numeric(unlist(wavelet$Period))
  pmax_avg_sel <- t(wavelet$Power)

  depth <-  wavelet$x
  y_axis <- wavelet$Period
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)


  MTM_res_2 <- matrix(data=NA,ncol=2,nrow=1)

  if (add_MTM == TRUE |   add_MTM_peaks == TRUE) {
    MTM_res_1 <- mtm(
        cbind(wavelet$x, wavelet$y),
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
            cbind(wavelet$x, wavelet$y),
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





  # the columns of the map are at axis.1; for windowed maps (chirp EHA) that differs from the data positions x
  depth_ax <- if (!is.null(wavelet$axis.1) && length(wavelet$axis.1) == ncol(wavelet$Power)) { wavelet$axis.1 } else { wavelet$x }
  p <- list(wavelet = wavelet, data_x = wavelet$x, data_y = wavelet$y, depth = depth_ax, ax2 = wavelet$axis.2,
      z = wavelet$Power, avg = wavelet$Power.avg, per = wavelet$Period, tick_min = log2(wavelet$Period[1]),
      overlay = TRUE, n.levels = n.levels, key.cols = key.cols, key.marks = key.marks,
      key.labels = key.labels, lwd.axis = lwd.axis, pmax_avg_sel = pmax_avg_sel, xlim_vals = xlim_vals,
      ylim_vals = ylim_vals, periodlab = periodlab, x_lab = x_lab, main = main, periodtck = periodtck,
      periodtcl = periodtcl, add_lines = add_lines, add_points = add_points, add_abline_h = add_abline_h,
      add_abline_v = add_abline_v, add_MTM_peaks = add_MTM_peaks, add_pval = add_pval,
      pval_abline = pval_abline, pval_cutoff = pval_cutoff, add_signif = add_signif, MTM_res_1 = if (add_MTM == TRUE | add_MTM_peaks == TRUE) { MTM_res_1 } else { NULL },
      MTM_res_2 = MTM_res_2, plot_horizontal = plot_horizontal, add_data = add_data, add_avg = add_avg,
      add_MTM = add_MTM, scale = scale)
  .pw_draw(p)

  if (add_MTM_peaks == TRUE & is.na(MTM_res_2[1,2]) == FALSE) {
    return(invisible(mtm_res))
  }
}

