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
#' @param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#' @param tracked_cycle_curve Curve of the cycle tracked using the
#' \code{\link{track_period_wavelet}} function. Any input (matrix or data frame)
#'  in which the first column is depth or time and the second column is period should work.
#' @param multi multiple of the standard deviation to be used as boundaries for the cycle extraction
#'  \code{Default=1}.
#' @param extract_cycle Period of the cycle to be extracted.
#' @param tracked_cycle_period Period of the tracked cycle.
#' @param add_mean Add mean to the extracted cycle \code{Default=TRUE}.
#' @param tune Tune data set using the \code{Default=tracked_cycle_curve} curve \code{Default=FALSE}.
#' @param genplot_uncertainty_wt Generate a wavelet spectra plot with the tracked curve and its
#' analytical uncertainty based the Gabor uncertainty principle applied
#' continuous wavelet transform using a Morlet wavelet on superimposed on top of it.
#' In the plot the red curve and blue curves are the upper and lower bounds
#'based on the \code{multi} parameter which x-times the standard deviation of uncertainty.
#'The black curve is the \code{Default=FALSE} curve.
#' @param genplot_extracted Generates a plot with the data set and
#' the extracted cycle on top \code{Default=FALSE} of it.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'@param palette_name Name of the color palette which is used for plotting.
#'The color palettes than can be chosen depends on which the R package is specified in
#'the color_brewer parameter. The included R packages from which palettes can be chosen
#'from are; the 'RColorBrewer', 'grDevices', 'ColorRamps' and 'Viridis' R packages.
#'There are many options to choose from so please
#'read the documentation of these packages \code{Default=rainbow}.
#'The R package 'viridis' has the color palette options: “magma”, “plasma”,
#'“inferno”, “viridis”, “mako”, and “rocket”  and “turbo”
#'To see the color palette options of the The R pacakge 'RColorBrewer' run
#'the RColorBrewer::brewer.pal.info() function
#'The R package 'colorRamps' has the color palette options:"blue2green",
#'"blue2green2red", "blue2red",	"blue2yellow", "colorRamps",	"cyan2yellow",
#'"green2red", "magenta2green", "matlab.like", "matlab.like2" and	"ygobb"
#'The R package 'grDevices' has the built in  palette options:"rainbow",
#'"heat.colors", "terrain.colors","topo.colors" and "cm.colors"
#'To see even more color palette options of the The R pacakge 'grDevices' run
#'the grDevices::hcl.pals() function
#'@param color_brewer Name of the R package from which the color palette is chosen from.
#'The included R packages from which palettes can be chosen
#'are; the RColorBrewer, grDevices, ColorRamps and Viridis R packages.
#'There are many options to choose from so please
#'read the documentation of these packages. "\code{Default=grDevices}
#' @author
#' Code based on the \link[WaveletComp]{reconstruct} function of the 'WaveletComp' R package
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
#' communication engineering 93, no. 26 (1946): 429-441.\url{http://genesis.eecg.toronto.edu/gabor1946.pdf}
#'
#'Russell, Brian, and Jiajun Han. "Jean Morlet and the continuous wavelet transform.
#'" CREWES Res. Rep 28 (2016): 115. \url{https://www.crewes.org/Documents/ResearchReports/2016/CRR201668.pdf}
#'
#'
#'Morlet, Jean, Georges Arens, Eliane Fourgeau, and Dominique Glard.
#'"Wave propagation and sampling theory—Part I: Complex signal and scattering in multilayered media.
#'" Geophysics 47, no. 2 (1982): 203-221.
#' \url{https://pubs.geoscienceworld.org/geophysics/article/47/2/203/68601/Wave-propagation-and-sampling-theory-Part-I}
#'
#'J. Morlet, G. Arens, E. Fourgeau, D. Giard;
#' Wave propagation and sampling theory; Part II, Sampling theory and complex waves.
#'Geophysics 1982 47 (2): 222–236. \url{https://pubs.geoscienceworld.org/geophysics/article/47/2/222/68604/Wave-propagation-and-sampling-theory-Part-II}
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
#'mag_wt <- analyze_wavelet(data = mag,
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
#' @export
#' @importFrom Hmisc approxExtrap
#' @importFrom stats na.omit
#' @importFrom graphics par
#' @importFrom graphics hist
#' @importFrom graphics lines
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom DescTools Closest
#' @importFrom WaveletComp reconstruct
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom colorRamps blue2green
#' @importFrom colorRamps blue2green2red
#' @importFrom colorRamps blue2red
#' @importFrom colorRamps blue2yellow
#' @importFrom colorRamps cyan2yellow
#' @importFrom colorRamps green2red
#' @importFrom colorRamps magenta2green
#' @importFrom colorRamps matlab.like
#' @importFrom colorRamps matlab.like2
#' @importFrom colorRamps ygobb
#' @importFrom viridis viridis
#' @importFrom viridis magma
#' @importFrom viridis plasma
#' @importFrom viridis inferno
#' @importFrom viridis cividis
#' @importFrom viridis mako
#' @importFrom viridis rocket
#' @importFrom viridis turbo
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
                                              palette_name="rainbow",
                                              color_brewer="grDevices") {
  my.w <- wavelet
  my.data <- cbind(wavelet$x, wavelet$y)
  filtered_cycle <- my.data[, 1]
  filtered_cycle <- as.data.frame(filtered_cycle)
  filtered_cycle$value <- NA

  completed_series <- na.omit(tracked_cycle_curve)
  completed_series[, 2] <-
    completed_series[, 2] * (extract_cycle / tracked_cycle_period)
  app <-
    approxExtrap(
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
      dj * sqrt(dt) / (pi ^ (-1 / 4))
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
    row_nr_high <- Closest(Period[], interpolated[i, 3], which = TRUE)
    row_nr_low <-
      Closest(Period[], interpolated[i, 4], which = TRUE)

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
      plot.COI = TRUE,
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
      analyze_wavelet(
        data = data_set_time,
        dj = 1 / 200,
        lowerPeriod = steps_size,
        upperPeriod =  data_set_time[nrow(data_set_time), 1] - data_set_time[1, 1],
        verbose = FALSE,
        omega_nr = wavelet$omega_nr
      )

    plot_wavelet(
      wavelet = data_set_time_wt,
      plot.COI = TRUE,
      n.levels = 100,
      palette_name = "rainbow",
      color_brewer= "grDevices",
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
    dev.new(width = 7,
            height = 4,
            noRStudioGD = TRUE)
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
