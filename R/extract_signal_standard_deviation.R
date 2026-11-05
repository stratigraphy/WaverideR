#' @title Extract a signal using standard deviation
#'
#' @description Extract signal from a
#' wavelet spectra in the depth domain using a the standard deviation of the omega (number of cycles)
#' as boundaries. The uncertainty is based on the Gabor uncertainty principle applied to the
#' continuous wavelet transform using a Morlet wavelet. The calculated uncertainty is the underlying
#' analytical uncertainty which is the result of applying the Gabor uncertainty principle to the
#' continuous wavelet transform using a Morlet wavelet.
#'
#'
#' @param wavelet Wavelet object (list) created with \code{\link{analyse_wavelet}}; the signal is reconstructed from
#' its coefficients and its number of cycles \code{omega_nr} sets the Gabor (Morlet) uncertainty. Required
#' (\code{Default=NULL}).
#' @param tracked_cycle_curve Tracked period of the cycle along the record (matrix or data frame: depth/time in
#' column 1, period in column 2, in the depth/time unit of the wavelet), e.g. the output of
#' \code{\link{track_period_wavelet}} completed with \code{\link{completed_series}}. Rows with \code{NA} are dropped
#' and the curve is linearly interpolated (and extrapolated) to the depth/time axis of the wavelet. Required
#' (\code{Default=NULL}).
#' @param multi Multiple of the standard deviation of the analytical frequency uncertainty of the Morlet wavelet
#' that sets the band boundaries, \code{1 / (f0 -/+ multi * sd)} around the frequency \code{f0} of the (rescaled)
#' tracked period (numeric, > 0); larger values give a wider band. \code{Default=1}.
#' @param extract_cycle Duration of the cycle to extract (numeric, > 0), in the same unit as
#' \code{tracked_cycle_period}; use the same value to extract the tracked cycle itself. The tracked curve is
#' multiplied by \code{extract_cycle / tracked_cycle_period}. Required (\code{Default=NULL}).
#' @param tracked_cycle_period Duration in time of the tracked cycle, e.g. 405 (kyr) for the long eccentricity cycle
#' (numeric, > 0), in the same unit as \code{extract_cycle}; also used for the tuning when \code{tune=TRUE}.
#' Required (\code{Default=NULL}).
#' @param add_mean Add the mean of the original data series to the extracted cycle (logical, TRUE or FALSE);
#' \code{FALSE} returns it centred on zero. \code{Default=TRUE}.
#' @param tune Convert the extracted cycle (and the plots) from the depth to the time domain with
#' \code{\link{curve2tune}}, using \code{tracked_cycle_curve} and \code{tracked_cycle_period} (logical, TRUE or
#' FALSE). \code{Default=FALSE}.
#' @param genplot_uncertainty_wt Plot the wavelet scalogram (with \code{\link{plot_wavelet}}) with the tracked curve
#' (black) and the upper (red) and lower (blue) boundaries set by \code{multi} times the standard deviation of the
#' uncertainty superimposed (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param genplot_extracted Plot the data with the extracted cycle (blue) on top of it (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param keep_editable Keep the plot settings active after the call so that extra features can be added to the
#' plots (logical, TRUE or FALSE); \code{FALSE} restores the previous graphical parameters. \code{Default=FALSE}.
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
#' @author
#' Code based on the reconstruct function of the 'WaveletComp' R package
#' which is based on the wavelet 'MATLAB' code written by Christopher Torrence and Gibert P. Compo (1998).
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
#'Gabor, Dennis. "Theory of communication. Part 1: The analysis of information."
#' Journal of the Institution of Electrical Engineers-part III: radio and
#' communication engineering 93, no. 26 (1946): 429-441. \doi{10.1049/ji-3-2.1946.0074}
#'
#'Russell, Brian, and Jiajun Han. "Jean Morlet and the continuous wavelet transform.
#'" CREWES Res. Rep 28 (2016): 115. \url{https://www.crewes.org/Documents/ResearchReports/2016/CRR201668.pdf}
#'
#'
#'Morlet, Jean, Georges Arens, Eliane Fourgeau, and Dominique Glard.
#'"Wave propagation and sampling theory—Part I: Complex signal and scattering in multilayered media.
#'" Geophysics 47, no. 2 (1982): 203-221.
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#'Geophysics 1982 47 (2): 222–236.
#'
#'
#'@examples
#'\donttest{
#'#Extract the 405 kyr eccentricity cycle from the magnetic susceptibility
#'#record of the Sullivan core of Pas et al., (2018) and use the Gabor
#'# uncertainty principle to define the mathematical uncertainty of the
#'# analysis and use a factor of that standard deviation to define
#'# boundaries
#'
#'# perform the CWT
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
#'#                                   palette_name="rainbow",
#'#                                   color_brewer="grDevices")
#'
#'#Instead of tracking, the tracked solution data set mag_track_solution is used
#'mag_track <- mag_track_solution
#'
#' mag_track_complete <- completed_series(
#'   wavelet = mag_wt,
#'   tracked_curve = mag_track,
#'   period_up = 1.2,
#'   period_down = 0.8,
#'   extrapolate = TRUE,
#'   genplot = FALSE
#' )
#'
#'# smooth the tracking of the 405 kyr eccentricity cycle
#' mag_track_complete <- loess_auto(time_series = mag_track_complete,
#' genplot = FALSE, print_span = FALSE)
#'
#'# extract the 405 kyr eccentricity cycle from the wavelet spectrum and use
#'# the Gabor uncertainty principle to define the mathematical uncertainty of
#'# the analysis and use a multiple of the derived standard deviation to define boundaries
#'
#'mag_405_ecc <- extract_signal_standard_deviation(
#'wavelet = mag_wt,
#'tracked_cycle_curve = mag_track_complete,
#'multi = 1,
#'extract_cycle = 405,
#'tracked_cycle_period = 405,
#'add_mean = TRUE,
#'tune = FALSE,
#'genplot_uncertainty_wt = FALSE,
#'genplot_extracted = FALSE,
#'keep_editable=FALSE,
#'palette_name="rainbow",
#'color_brewer="grDevices"
#')
#'}
#' @return Signal extracted from the wavelet spectra.
#' Output is a matrix with the first column being depth/time
#'and the second column is the astronomical cycle extracted from the proxy record
#'
#'If \code{genplot_uncertainty_wt=TRUE} then a wavelet spectra will be plotted
#'with the uncertainty superimposed on top of it. In the plot the red curve and
#' blue curves are the upper and lower bounds
#'based on the \code{multi} parameter.The black curve is the \code{Default=tracked_cycle_curve} curve.
#'If \code{genplot_extracted=TRUE} plot with the data set and
#'the extracted cycle on top of it will be plotted.
#'
#' @param colour_brewer British spelling of \code{color_brewer}, with the same options; if given (not \code{NULL}) it overrides
#' \code{color_brewer}. \code{Default=NULL}.
#' @export
#' @importFrom stats na.omit
#' @importFrom graphics par
#' @importFrom graphics hist
#' @importFrom graphics lines
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors



extract_signal_standard_deviation <- function(wavelet = NULL,
    tracked_cycle_curve = NULL,
    multi = 1,
    extract_cycle = NULL,
    tracked_cycle_period = NULL,
    add_mean = TRUE,
    tune = FALSE,
    genplot_uncertainty_wt = FALSE,
    genplot_extracted = FALSE,
    keep_editable = FALSE,
    palette_name = "viridis",
    color_brewer = "viridis",
    colour_brewer = NULL) {
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym
  my.w <- wavelet
  my.data <- cbind(wavelet$x, wavelet$y)
  filtered_cycle <- my.data[, 1]
  filtered_cycle <- as.data.frame(filtered_cycle)
  filtered_cycle$value <- NA

  completed_series <- na.omit(tracked_cycle_curve)
  completed_series[, 2] <-
      completed_series[, 2] * (extract_cycle / tracked_cycle_period)
  app <-
      wr_approx_extrap(
          x = completed_series[, 1],
          y = completed_series[, 2],
          xout = my.data[, 1],
          method = "linear"
      )
  interpolated <- cbind(app$x, app$y)

  Wave = my.w$Wave
  Power = my.w$Power

  nc = my.w$nc
  nr = my.w$nr
  dt = my.w$dt
  dj = my.w$dj

  Scale = my.w$Scale
  Period = my.w$Period
  loess.span = my.w$loess.span
  rec.waves = matrix(0, nrow = nr, ncol = nc)


  for (s.ind in seq_len(nr)) {
    rec.waves[s.ind, ] = (Re(Wave[s.ind, ]) / sqrt(Scale[s.ind])) *
        dj * sqrt(dt) / .wr_psi0(my.w)
  }


  interpolated <- as.data.frame(interpolated)

  ncycles <- my.w$omega_nr
  b <- (2 * sqrt(2 * log(2)))
  a <- ((8 * log(2) / (2 * pi)))
  k <- (ncycles / (8 * log(2))) * 2



  interpolated$f0 <- (1 / (interpolated[, 2]))
  interpolated$df <- (a * interpolated$f0) / k
  interpolated$sd_morlet <- interpolated$df / b


  fact_high <-
      1 / (interpolated$f0 - (interpolated$sd_morlet * multi))
  fact_low <-
      1 / (interpolated$f0 + (interpolated$sd_morlet * multi))

  fact_high[fact_high > max(my.w$Period)] <- max(my.w$Period)
  fact_low[fact_low < min(my.w$Period)] <- min(my.w$Period)


  interpolated[, 3] <- fact_high
  interpolated[, 4] <- fact_low

  for (i in 1:nrow(filtered_cycle)) {
    row_nr_high <- wr_closest(Period[], interpolated[i, 3], which = TRUE)
    row_nr_low <-
        wr_closest(Period[], interpolated[i, 4], which = TRUE)

    row_nr_high <- row_nr_high[1]
    row_nr_low <- row_nr_low[1]

    value <- rec.waves[c(row_nr_low:row_nr_high), i]
    value  <- sum(value, na.rm = T)
    value <- as.numeric(value)
    filtered_cycle[i, 2] <- value
  }

  rec_value  <- colSums(rec.waves, na.rm = T)

  filtered_cycle[, 2] <-
      (filtered_cycle[, 2]) * sd(my.data[, 2]) / sd(rec_value)

  if (add_mean == TRUE) {
    filtered_cycle[, 2] <- filtered_cycle[, 2] + mean(my.data[, 2])
  }


  if (genplot_uncertainty_wt == TRUE & tune == FALSE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    plot_wavelet(
        wavelet = wavelet,
        n.levels = 100,
        palette_name = palette_name,
        color_brewer= color_brewer,
        useRaster = TRUE,
        periodlab = "Period (metres)",
        x_lab = "depth (metres)",
        keep_editable = TRUE
    )



    combined_sedrate <-
        cbind(interpolated[, 1],
            1 / interpolated[, 2],
            interpolated$sd_morlet * multi)

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

    lines(interpolated[, 1], log2(interpolated[, 2]), lwd = 2)
    lines(interpolated[, 1],
        log2(interpolated[, 3]),
        col = "red",
        lwd = 2)
    lines(interpolated[, 1],
        log2(interpolated[, 4]),
        col = "blue",
        lwd = 2)

  }


  if (genplot_uncertainty_wt == TRUE & tune == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    data_set <- cbind(wavelet$x, wavelet$y)
    data_set_time <- curve2tune(
        data = data_set,
        tracked_cycle_curve = tracked_cycle_curve,
        tracked_cycle_period = tracked_cycle_period,
        genplot = FALSE
    )

    dat <- as.matrix(data_set_time)
    dat <- na.omit(dat)
    dat <- dat[order(dat[, 1], na.last = NA, decreasing = F),]
    npts <- length(dat[, 1])
    start <- dat[1, 1]
    end <- dat[length(dat[, 1]), 1]
    x1 <- dat[1:(npts - 1), 1]
    x2 <- dat[2:(npts), 1]
    dx = x2 - x1
    dt = median(dx)
    steps_size <- dt
    xout <- seq(start, end, by = dt)
    npts <- length(xout)
    interp <- approx(dat[, 1], dat[, 2], xout, method = "linear",
        n = npts)
    data_set_time <- as.data.frame(interp)

    data_set_time_wt <-
        analyse_wavelet(
            data = data_set_time,
            dj = 1 / 200,
            lowerPeriod = steps_size,
            upperPeriod =  data_set_time[nrow(data_set_time), 1] - data_set_time[1, 1],
            verbose = FALSE,
            omega_nr = wavelet$omega_nr
        )

    plot_wavelet(
        wavelet = data_set_time_wt,
        n.levels = 100,
        palette_name = "viridis",
        color_brewer = "viridis",
        useRaster = TRUE,
        periodlab = "Period (kyr)",
        x_lab = "depth (metres)",
        keep_editable = TRUE
    )


    interpolated_time <- curve2tune(
        data = interpolated[, c(1, 2)],
        tracked_cycle_curve = tracked_cycle_curve,
        tracked_cycle_period = tracked_cycle_period,
        genplot = FALSE
    )


    combined_sedrate <-
        cbind(
            interpolated_time[, 1],
            1 / (
                interpolated[, 2] / (tracked_cycle_curve[, 2] * tracked_cycle_period)
            )
            ,
            interpolated$sd_morlet * multi * (tracked_cycle_curve[, 2] * tracked_cycle_period)
        )
    xcords <-
        c(combined_sedrate[, 1],
            sort(combined_sedrate[, 1], decreasing = TRUE))
    data_sort1 <-
        combined_sedrate[order(combined_sedrate[, 1], decreasing = TRUE), ]
    ycords <-
        c(1 / (combined_sedrate[, 2] + combined_sedrate[, 3]),
            1 / (data_sort1[, 2] - data_sort1[, 3]))
    polygon(
        x = xcords,
        y = log2(1 / ycords),
        col = rgb(0.5, 0.5, 0.5, 0.5),
        border = "black"
    )

    lines(interpolated_time[, 1], log2(interpolated[, 2]), lwd = 2)
    lines(interpolated_time[, 1],
        log2((
            combined_sedrate[, 2] + combined_sedrate[, 3]
        )),
        col = "red",
        lwd = 2)
    lines(interpolated_time[, 1],
        log2((
            combined_sedrate[, 2] - combined_sedrate[, 3]
        )),
        col = "blue",
        lwd = 2)

  }



  if (tune == TRUE) {
    filtered_cycle <- curve2tune(
        data = filtered_cycle,
        tracked_cycle_curve = tracked_cycle_curve,
        tracked_cycle_period = tracked_cycle_period,
        genplot = FALSE
    )
  }


  if (genplot_extracted == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    wr_open_device(width = 7, height = 4)
    if (tune == TRUE) {
      my.data_time <- curve2tune(
          data = my.data,
          tracked_cycle_curve = tracked_cycle_curve,
          tracked_cycle_period = tracked_cycle_period,
          genplot = FALSE
      )
      plot(my.data_time, type = "l")
    } else{
      plot(my.data[, 1], my.data[, 2], type = "l")
    }

    if (add_mean == TRUE) {
      lines(filtered_cycle[, 1], filtered_cycle[, 2], col = "blue")
    } else{
      lines(filtered_cycle[, 1],
          filtered_cycle[, 2] + mean(my.data[, 2]),
          col = "blue")
    }


  }

  return(filtered_cycle)
}
