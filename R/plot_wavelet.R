#' @title Plots a wavelet scalogram
#'
#' @description Plot wavelet scalogram using the outcome of the \code{\link{analyse_wavelet}} function.
#'
#'@param wavelet wavelet object created using the \code{\link{analyse_wavelet}} function.
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
#'@param add_pval add an transparent overlay on the wavelet scalogram based on the p-value and add the p-value curve to the
#'average spectral power curve. The p-value is based on a Monte Carlo simulation of the \code{\link{analyse_wavelet}} function.
#'The p-value is based on Monte Carlo modelling runs on surrogate data generated based on autocorrelated noise (red noise)
#'the calculated using a windowed (the window is half the size of the data set) temporal autocorrelation
#' and on shuffling the data set resulting in a random data sets which has similar spectral characteristics
#' to the original data set.The shuffling of the data set creates white noise which ensures that high amplitude high frequency/short
#' period cycles do not result in statistical significant peaks. The part of the data generated using the  autocorrelated noise (red noise)
#' based on the windowed  (the window is half the size of the data set) temporal autocorrelation represent a spectral signature similar to
#' to that of the original data. The original data might include spectral peaks which are the result of astronomical
#' forcing. The result is that the spectral power profile is biased towards rejecting the 0-hypothesis (e.g. no astronomical forcing).
#' By combining the shuffling of the data set with autocorrelated noise a surrogate data set is created which rejects
#' high amplitude high frequency/short period cycles and a reduced biased towards towards rejecting the 0-hypothesis if the data was
#' solely the result of autocorrelated noise. \code{Default=FALSE}
#'@param pval_abline Add straight lines to the average spectral power plot which indicate certain p-values \code{Default=c(0.1,0.5)}
#'@param pval_cutoff Cutoff p-value to be used in the transparent overlay of the wavelet scalogram \code{Default=c(0.1)}
#'@param add_MTM Add the MTM  plot next to the wavelet plot \code{Default=FALSE}
#'@param mtm_siglvl select the significance level  (0-1) for the MTM spectrum \code{Default=0.95}
#'@param demean_mtm Remove mean from data before conducting the MTM analysis \code{Default=TRUE}
#'@param detrend_mtm Remove mean from data before conducting the MTM analysis \code{Default=TRUE}
#'@param padfac_mtm Pad factor for the MTM analysis \code{Default=5}
#'@param tbw_mtm time bandwidth product of the MTM analysis  \code{Default=3}
#' @param add_signif Optional significance object created with \code{\link{wavelet_significance}}. The significant
#' cells are outlined on the scalogram and the critical level of the global spectrum is drawn as a red dashed line on
#' the average power panel \code{Default=NULL}.
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis is depth or y axis power \code{Default=TRUE}
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
#' @param colour_brewer British spelling of \code{color_brewer}; if given it overrides \code{color_brewer}.
#' @param scale Colour scale, the same for every image plot of the package: \code{"auto"} (default: logarithmic for
#' synchrosqueezed, superlet-SST, chirplet-SST and other mostly empty maps, histogram-equalised quantiles otherwise),
#' \code{"log"}, \code{"quantile"}, \code{"sqrt"} or \code{"linear"}.
#' @param whiten Divide the power of every period by its mean along the record before plotting \code{Default=FALSE}.
#' Removes the red-noise background, so that weak ridges at short periods (precession in a noisy record) are shown with
#' the same contrast as the strong long periods; on CIP Case 3 (very red noise) the oscillating precession and obliquity
#' ridges only become visible this way. The average panel is not informative with \code{whiten = TRUE}.
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
                         palette_name = "rainbow",
                         color_brewer= "grDevices",
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
  whiten = FALSE, scale = "auto") {
  if (!is.null(colour_brewer)) color_brewer <- colour_brewer  # British spelling synonym
  if (isTRUE(whiten)) {
    # divide every period row by its mean power: removes the red-noise background so that ridges at short and long
    # periods are shown with comparable contrast (the average panel is then flat and is not informative)
    rm <- rowMeans(wavelet$Power, na.rm = TRUE); rm[!is.finite(rm) | rm <= 0] <- 1
    wavelet$Power <- wavelet$Power / rm
    wavelet$Power.avg <- rowMeans(wavelet$Power, na.rm = TRUE)
    if (!is.null(wavelet$Ampl)) wavelet$Ampl <- wavelet$Ampl / sqrt(rm)
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
    wr_open_device(width = 15, height = 7)}

  if(dev_new==TRUE & plot_horizontal==FALSE){
     wr_open_device(width = 7, height = 10)}



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
  depth_ax <- if (!is.null(wavelet$axis.1) && length(wavelet$axis.1) == ncol(wavelet$Power)) wavelet$axis.1 else wavelet$x
  p <- list(wavelet = wavelet, data_x = wavelet$x, data_y = wavelet$y, depth = depth_ax, ax2 = wavelet$axis.2,
            z = wavelet$Power, avg = wavelet$Power.avg, per = wavelet$Period, tick_min = log2(wavelet$Period[1]),
            overlay = TRUE, n.levels = n.levels, key.cols = key.cols, key.marks = key.marks,
            key.labels = key.labels, lwd.axis = lwd.axis, pmax_avg_sel = pmax_avg_sel, xlim_vals = xlim_vals,
            ylim_vals = ylim_vals, periodlab = periodlab, x_lab = x_lab, main = main, periodtck = periodtck,
            periodtcl = periodtcl, add_lines = add_lines, add_points = add_points, add_abline_h = add_abline_h,
            add_abline_v = add_abline_v, add_MTM_peaks = add_MTM_peaks, add_pval = add_pval,
            pval_abline = pval_abline, pval_cutoff = pval_cutoff, add_signif = add_signif, MTM_res_1 = if (add_MTM == TRUE | add_MTM_peaks == TRUE) MTM_res_1 else NULL,
            MTM_res_2 = MTM_res_2, plot_horizontal = plot_horizontal, add_data = add_data, add_avg = add_avg,
            add_MTM = add_MTM, scale = scale)
  .pw_draw(p)

  if (add_MTM_peaks == TRUE & is.na(MTM_res_2[1,2]) == FALSE) {
    return(invisible(mtm_res))
  }
}



