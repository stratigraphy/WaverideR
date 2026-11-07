#' @title Calculate the uncertainty associated with the wavelet analysis based on the Gabor uncertainty principle
#'
#' @description The \code{\link{wavelet_uncertainty}} function is used to calculate uncertainties associated
#' with the wavelet analysis based on the Gabor uncertainty principle applied to the
#' continuous wavelet transform using a Morlet wavelet. The calculated uncertainty is the underlying
#' analytical uncertainty which is the result of applying the Gabor uncertainty principle to the
#' continuous wavelet transform using a Morlet wavelet.
#'
#'
#' @param tracked_cycle Curve of the tracked cycle (data frame with depth or time in column 1 and the period of the
#' tracked cycle, in the same unit, in column 2; further columns are ignored), e.g. from
#' \code{\link{track_period_wavelet}} completed with \code{\link{completed_series}} and smoothed with
#' \code{\link{loess_auto}}. Must be given (the \code{NULL} default fails). \code{Default=NULL}.
#' @param period_of_tracked_cycle Duration of the tracked cycle in time (numeric, > 0, in kyr, e.g. 405 for the long
#' eccentricity cycle), passed to \code{\link{curve2time}} to convert the period curves into time. Must be given.
#' \code{Default=NULL}.
#' @param wavelet Wavelet object created with \code{\link{analyse_wavelet}} (class \code{analyze.wavelet}); its
#' \code{omega_nr} sets the width of the Gabor uncertainty, and it is plotted with \code{genplot_uncertainty_wt}. Must
#' be given. \code{Default=NULL}.
#' @param multi Multiple of the standard deviation of the Gabor uncertainty (in frequency) used for the upper and lower
#' curves (numeric, > 0; e.g. 2 for about 95 percent). \code{Default=1}.
#' @param verbose Print the total duration, the uncertainty of the duration and the relative uncertainty per cycle
#' (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param genplot_time Plot the depth-time curves (black: tracked period, red and blue: period plus and minus
#' \code{multi} standard deviations of the Gabor uncertainty of the Morlet wavelet) (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param genplot_uncertainty Plot the tracked period against depth with the upper (red) and lower (blue) period
#' curves of \code{multi} standard deviations (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param genplot_uncertainty_wt Plot the wavelet spectrum (\code{\link{plot_wavelet}}) with the tracked period
#' (black), the period plus (red) and minus (blue) the uncertainty of \code{multi} standard deviations, and the area in
#' between shaded in grey (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param keep_editable Keep the graphical parameters after plotting, so that extra features can be added to the last
#' plot (logical, TRUE or FALSE); \code{FALSE} restores the previous settings. \code{Default=FALSE}.

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
#'@return Results pertaining to the uncertainty calculated based on the Gabor uncertainty principle.\cr
#'If the genplot_time is TRUE then a depth time plot will be plotted with 3 lines, the mean age,age plus
#' x times the standard deviation and age minus x times the standard deviation . \cr
#'If the genplot_uncertainty is TRUE then a curve will be plotted with the mean period, the tracked period plus
#'x times the standard deviation and the tracked period minus x times the standard deviation. \cr
#'If the genplot_uncertainty_wt is TRUE a wavelet spectra will be plotted with the tracked period, the tracked period plus
#'x times the standard deviation,the tracked period minus x times the standard deviation and the area in between will be shaded in grey.\cr
#'
#'Returns a matrix with 8 columns.\cr
#'The first column is called "depth" eg. depth \cr
#'The second column is "period" of the originally tracked period. \cr
#'The third column is "frequency" of the originally tracked period. \cr
#'The fourth column "uncertainty in frequency FWHM" is the uncertainty in frequency based on the Gabor uncertainty principle defined as
#' (FWHM) full width at half maximum. \cr
#'The fifth column "uncertainty in frequency x_times SD" is the uncertainty in frequency based on the Gabor uncertainty principle defined as
#' times x standard deviations. \cr
#'The sixth column "time mean" is the mean time based on the tracked period. \cr
#'The seventh column "time plus x_times sd" is the time based on the tracked period plus x times the standard deviation. \cr
#'The eight column "time min x_times sd" is the time based on the tracked period min x times the standard deviation. \cr
#'
#' @author
#' Code based on the "analyze.wavelet" function of the 'WaveletComp' R package
#' and "wt" function of the 'biwavelet' R package which are based on the
#' wavelet 'MATLAB' code written by Christopher Torrence and Gibert P. Compo (1998).
#' The assignment of the standard deviation of the uncertainty of the wavelet
#' is based on the work of Gabor (1946) and Russell et al., (2016)
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
#'Gabor, Dennis. Theory of communication. Part 1: The analysis of information.
#' Journal of the Institution of Electrical Engineers-part III: radio and
#' communication engineering 93, no. 26 (1946): 429-441. \doi{10.1049/ji-3-2.1946.0074}
#'
#'Russell, Brian, and Jiajun Han.Jean Morlet and the continuous wavelet transform.
#'CREWES Res. Rep 28 (2016): 115. \url{https://www.crewes.org/Documents/ResearchReports/2016/CRR201668.pdf}
#'
#'
#'Morlet, Jean, Georges Arens, Eliane Fourgeau, and Dominique Glard.
#'Wave propagation and sampling theory—Part I: Complex signal and scattering in multilayered media.
#'Geophysics 47, no. 2 (1982): 203-221.
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#'  Geophysics 1982 47 (2): 222–236.
#'
#'
#'@examples
#'\donttest{
#'#calculate the Gabor uncertainty derived mathematical uncertainty of the
#'#magnetic susceptibility record of the Sullivan core of Pas et al., (2018)
#'
#'mag_wt <- analyse_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'#Track the 405 kyr eccentricity cycle in a wavelet spectra
#'
#'#mag_track <- track_period_wavelet(astro_cycle = 405,
#'#                                   wavelet=mag_wt,
#'#                                   n.levels = 100,
#'#                                   periodlab = "Period (metres)",
#'#                                   x_lab = "depth (metres)",
#'#                                palette_name="rainbow",
#'#                                color_brewer= "grDevices")
#'
#'#Instead of tracking, the tracked solution data set mag_track_solution is used
#'mag_track <- mag_track_solution
#'
#'mag_track_complete <- completed_series(
#'   wavelet = mag_wt,
#'   tracked_curve = mag_track,
#'   period_up = 1.2,
#'   period_down = 0.8,
#'   extrapolate = FALSE,
#'   genplot = FALSE,
#'   keep_editable=FALSE
#' )
#'
#' mag_track_complete <- loess_auto(time_series = mag_track_complete,
#'  genplot = FALSE, print_span = FALSE,keep_editable=FALSE)
#'
#'uncertainty <- wavelet_uncertainty(
#'   tracked_cycle = mag_track_complete,
#'   period_of_tracked_cycle = 405,
#'   wavelet = mag_wt,
#'   multi=1,
#'   verbose = FALSE,
#'   genplot_time = FALSE,
#'   genplot_uncertainty = FALSE,
#'   genplot_uncertainty_wt = FALSE,
#'   keep_editable=FALSE,
#'   palette_name="rainbow",
#'   color_brewer= "grDevices"
#' )
#'}
#' @param colour_brewer British spelling of \code{color_brewer}, with the same options; if given (not \code{NULL}) it overrides
#' \code{color_brewer}. \code{Default=NULL}.
#' @export
#' @importFrom grDevices dev.new
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors

wavelet_uncertainty <- function(tracked_cycle = NULL,
    period_of_tracked_cycle = NULL,
    wavelet = NULL,
    multi = 1,
    verbose = FALSE,
    genplot_time = FALSE,
    genplot_uncertainty = FALSE,
    genplot_uncertainty_wt = FALSE,
    keep_editable = FALSE,
    palette_name = "viridis",
    color_brewer = "viridis",
    colour_brewer = NULL){
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym
  # a matrix or a data frame: the columns are added with $ below
  data <- as.data.frame(tracked_cycle)[, c(1, 2)]
  ncycles <- wavelet$omega_nr
  b <- (2 * sqrt(2 * log(2)))
  a <- ((8 * log(2) / (2 * pi)))
  k <- (ncycles / (8 * log(2))) * 2

  data$f0 <- (1 / data[, 2])
  data$df <- (a * data$f0) / k
  data$sd_morlet <- data$df / b

  sedrate_min <- 1 / (data$f0 - data$sd_morlet * multi)
  sedrate_plus <-  1 / (data$f0 + data$sd_morlet * multi)




  sedrate_min <- cbind(data[, 1], sedrate_min)
  sedrate_plus <- cbind(data[, 1], sedrate_plus)
  sedrate <-  cbind(data[, 1], 1 / (data$f0))

  time_min <- curve2time(
      tracked_cycle_curve = sedrate_min,
      tracked_cycle_period = period_of_tracked_cycle,
      genplot = FALSE
  )

  time_plus <- curve2time(
      tracked_cycle_curve = sedrate_plus,
      tracked_cycle_period = period_of_tracked_cycle,
      genplot = FALSE
  )

  time <- curve2time(
      tracked_cycle_curve = sedrate,
      tracked_cycle_period = period_of_tracked_cycle,
      genplot = FALSE
  )


  if (verbose == TRUE) {
    cat("total duration is: ", time[nrow(time), 2])
    cat("uncertainty is: ", time_plus[nrow(time_plus), 2] - time[nrow(time), 2])
    cat("uncertainty is defined as ", multi, " standard deviations ")
    cat("uncertainty defined as standard deviation is: ",
        ((a / k) / b),
        " cycle per cycle")

  }


  if (genplot_time == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    wr_open_device(width = 7, height = 7)
    plot(
        time_min,
        type = "l",
        col = "red",
        ylim = c(0, time_plus[nrow(time_plus), 2])
    )
    lines(time_plus, col = "blue")
    lines(time, col = "black", lwd = 3)
  }



  if (genplot_uncertainty == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    wr_open_device(width = 7, height = 7)
    max_y <- max(1 / (data$f0 - data$sd_morlet * multi))
    min_y <- min(1 / (data$f0 + data$sd_morlet * multi))

    plot(
        data[, 1],
        1 / data$f0,
        type = "l",
        ylim = c(min_y, max_y),
        lwd = 2,
        xlab = "depth",
        ylab = "period (m)"
    )
    lines(data[, 1],
        1 / (data$f0 - data$sd_morlet * multi),
        col = "red",
        lwd = 2)
    lines(data[, 1],
        1 / (data$f0 + data$sd_morlet * multi),
        col = "blue",
        lwd = 2)
  }


  if (genplot_uncertainty_wt == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }

    res <- plot_wavelet(
        wavelet = wavelet,
        n.levels = 100,
        palette_name = palette_name,
        color_brewer= color_brewer,
        useRaster = TRUE,
        periodlab = "Period (metres)",
        x_lab = "depth (metres)",
        keep_editable = TRUE,
        add_lines = NULL,
        add_points= NULL,
        add_abline_h = NULL,
        add_abline_v = NULL,
        add_MTM_peaks = FALSE,
        add_data = TRUE,
        add_avg = FALSE,
        add_MTM = FALSE,
        demean_mtm = TRUE,
        detrend_mtm = TRUE,
        padfac_mtm = 5,
        tbw_mtm = 3)



    lines(data[, 1], log2(1 / data$f0), lwd = 2)
    lines(data[, 1], log2(1 / (data$f0 - data$sd_morlet * multi)), col =
        "red", lwd = 2)
    lines(data[, 1], log2(1 / (data$f0 + data$sd_morlet * multi)), col =
        "blue", lwd = 2)

    combined_sedrate <-
        cbind(data[, 1], data$f0, data$sd_morlet * multi)

    xcords <-
        c(combined_sedrate[, 1],
            sort(combined_sedrate[, 1], decreasing = TRUE))
    xcords
    data_sort1 <-
        combined_sedrate[order(combined_sedrate[, 1], decreasing = TRUE), ]
    ycords <-
        c(1 / (combined_sedrate[, 2] + combined_sedrate[, 3]),
            1 / (data_sort1[, 2] - data_sort1[, 3]))

    polygon(
        x = xcords,
        y = log2(ycords),
        col = rgb(0.5, 0.5, 0.5, 0.5),
        border = "black"
    )

  }

  data$time_mean <- time[, 2]
  data$time_min_sd <- time_min[, 2]
  data$time_plus_sd <- time_plus[, 2]
  colnames(data) <- c(
      "depth",
      "period",
      "frequency",
      "uncertainty in frequency FWHM",
      "uncertainty in frequency x_times SD",
      "time mean",
      "time plus x_times sd",
      "time min x_times sd"
  )
  return(data)
}
