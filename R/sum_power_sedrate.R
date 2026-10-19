#' @title Calculate sum of maximum spectral power for sedimentation rates
#' for a wavelet spectra
#'
#' @description The \code{\link{sum_power_sedrate}} function is used calculate the sum of
#' maximum spectral power for a list of astronomical cycles from a wavelet spectra.
#' The data is first normalised using the average spectral power curves
#'  for a given percentile based on results of the \code{\link{model_red_noise_wt}} function
#'
#'@param red_noise Red noise curves generated using the \code{\link{model_red_noise_wt}} function
#'@param wavelet Wavelet object created using the \code{\link{analyse_wavelet}} function
#'@param percentile Percentile value (0-1) of the rednoise runs which is used to normalize the data for.
#'To account for the distribution/distortion of the spectral power distribution based on the analytical technique and
#'random red-noise the data is normalised against a percentile based red-noise curve which is the results of the
#''\code{\link{model_red_noise_wt}} modelling runs.
#'@param sedrate_low Minimum sedimentation rate (cm/kyr)for which the sum of maximum spectral power is calculated for.
#'@param sedrate_high Maximum sedimentation rate (cm/kyr) for which the sum of maximum spectral power is calculated  for.
#'@param spacing Spacing (cm/kyr) between sedimentation rates
#'@param cycles Astronomical cycles (in kyr) for which the combined sum of maximum spectral power is calculated for
#'@param x_lab label for the y-axis \code{Default="depth"}
#'@param y_lab label for the y-axis \code{Default="sedrate"}
#'@param run_multicore run simulation using multiple cores \code{Default=FALSE}
#'the simulation is run at x-2 cores to allow the 2 remaining processes to run background processes
#' @param n_cores Number of cores used when `run_multicore = TRUE` (`Default = NULL`: an
#'   automatic, safe number: physical cores minus one, at most 4, and limited by the free memory).
#'   Set `options(waverider.n_cores = ...)` to change the default. A value > 1 also switches on
#'   parallel computing. Results do not depend on the number of cores.
#'@param genplot Generate plot \code{Default="FALSE"}
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
#'@param color_brewer Name of the R package from which the color palette is chosen from.
#'The included R packages from which palettes can be chosen
#'are; the RColorBrewer, grDevices, ColorRamps and Viridis R packages.
#'There are many options to choose from so please
#'read the documentation of these packages. "\code{Default=grDevices}
#'@param plot_res plot options are 1: sum max power or 2: nr of components \code{Default=2}
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#' @param verbose Print text \code{Default=FALSE}.
#'
#' @author
#'Based on the \link[astrochron]{asm} and \link[astrochron]{eAsm}
#'functions of the 'astrochron' R package
#'and the 'eCOCO' and 'COCO' functions of the 'Acycle' software
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#'Acycle: Time-series analysis software for paleoclimate research and education,
#'Mingsong Li, Linda Hinnov, Lee Kump,
#'Computers & Geosciences,Volume 127,2019,Pages 12-22,ISSN 0098-3004,
#'<doi:10.1016/j.cageo.2019.02.011>
#'
#'Tracking variable sedimentation rates and astronomical forcing in Phanerozoic paleoclimate proxy series with evolutionary correlation coefficients and hypothesis testing,
#'Mingsong Li, Lee R. Kump, Linda A. Hinnov, Michael E. Mann,
#'Earth and Planetary Science Letters,Volume 501,
#'T2018,Pages 165-179,ISSN 0012-821X,<doi:10.1016/j.epsl.2018.08.041>
#'
#'@examples
#'\donttest{
#'#estimate sedimentation rate for the the magnetic susceptibility record
#'# of the Sullivan core of Pas et al., (2018).
#'
#'mag_wt <- analyse_wavelet(data = mag[1:200, ],
#' dj = 1/10,
#' lowerPeriod = 0.1,
#' upperPeriod = 32,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'#increase n_simulations to better define the red noise spectral power curve
#'mag_wt_red_noise <- model_red_noise_wt(wavelet=mag_wt,
#'n_simulations = 3,
#'run_multicore=FALSE,
#'verbose=FALSE)
#'
#'sedrates <- sum_power_sedrate(red_noise=mag_wt_red_noise,
#'wavelet=mag_wt,
#'percentile=0.75,
#'sedrate_low = 0.5,
#'sedrate_high = 4,
#'spacing = 1,
#'cycles = c(2376,1600,1180,696,406,110),
#'x_lab="depth",
#'y_lab="sedrate",
#'run_multicore=FALSE,
#'genplot = FALSE,
#'plot_res=1,
#'keep_editable=FALSE,
#'palette_name = "rainbow",
#'color_brewer="grDevices",
#'verbose=FALSE)
#'}
#'
#'
#' @return
#'Returns a list which contains 4 elements
#'element 1: sum of maximum spectral power
#'element 2: number of cycles used in the sum of maximum spectral power
#'element 3: y-axis values of the matrices which is sedimentation rate
#'element 4: x-axis values of the matrices which is depth
#'
#'If \code{Default="TRUE"} a plot is created with 3 subplots.
#'Subplot 1 is plot in which the the sum of maximum spectral power for a
#'given sedimentation rate or nr of cycles is plotted for each depth given depth.
#'Subplot 2 is a plot in which the average sum of maximum spectral power is plotted fro each sedimentation
#'Subplot 3 is a color scale for subplot 1.
#'
#'
#' @param colour_brewer British spelling of \code{color_brewer}; if given it overrides \code{color_brewer}.
#' @export
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics text
#' @importFrom graphics box
#' @importFrom graphics polygon
#' @importFrom grDevices rgb
#' @importFrom graphics layout
#' @importFrom astrochron asm
#' @importFrom astrochron eAsm
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors




sum_power_sedrate  <- function(red_noise = NULL,
                               wavelet = NULL,
                               percentile = NULL,
                               sedrate_low = NULL,
                               sedrate_high = NULL,
                               spacing = NULL,
                               cycles = c(NULL),
                               x_lab = "depth",
                               y_lab = "sedrate",
                               run_multicore = FALSE,
                               n_cores = NULL,
                               genplot = FALSE,
                               plot_res = 1,
                               keep_editable = FALSE,
                               palette_name = "rainbow",
                               color_brewer="grDevices",
                               verbose=FALSE,
  colour_brewer = NULL) {
  if (!is.null(colour_brewer)) color_brewer <- colour_brewer  # British spelling synonym
  my.data <- cbind(wavelet$x, wavelet$y)
  my.w <- wavelet

  # Wavelet power: rows = periods, columns = positions along the record.
  Powert <- matrix(as.numeric(my.w$Power), nrow = my.w$nr, ncol = my.w$nc)
  periods <- as.numeric(my.w$Period)
  nper <- length(periods)
  red_noise <- as.matrix(red_noise)
  if (is.unsorted(periods)) {
    ord <- order(periods)
    periods <- periods[ord]
    Powert <- Powert[ord, , drop = FALSE]
    red_noise <- red_noise[ord, , drop = FALSE]
  }

  # Local maxima of the power along the period axis (same rule for every position).
  maxdetect <- matrix(FALSE, nrow = nper, ncol = ncol(Powert))
  if (nper > 6) {
    i <- 3:(nper - 3)
    maxdetect[i, ] <- (Powert[i, , drop = FALSE] - Powert[i + 3, , drop = FALSE] > 0) &
      (Powert[i, , drop = FALSE] - Powert[i - 2, , drop = FALSE] > 0)
  }

  # Red-noise threshold curve (percentile of the red-noise power at every period).
  meanz <- rowMeans(red_noise)
  sdz <- wr_row_sds(red_noise)
  probabillity <- qnorm(p = percentile, mean = meanz, sd = sdz)
  probabillity[probabillity < 0] <- 0

  testsedrates <- seq(from = sedrate_low, to = sedrate_high, by = spacing)
  n_sed <- length(testsedrates)
  n_cyc <- length(cycles)
  simulations <- ncol(Powert)

  # Expected cycle periods (in depth) and the Morlet band around them for every sedimentation rate
  # and astronomical cycle. These do not depend on the position, so they are computed once.
  multi <- 2 * sqrt(2 * log(2))
  b <- 2 * sqrt(2 * log(2))
  a <- 8 * log(2) / (2 * pi)
  k <- (my.w$omega_nr / (8 * log(2))) * 2
  centre <- outer(testsedrates, cycles) / 100            # n_sed x n_cyc
  f0 <- 1 / centre
  sd_morlet <- ((a * f0) / k) / b
  long_p <- 1 / (f0 - sd_morlet * multi)                 # upper period limit of the band
  short_p <- 1 / (f0 + sd_morlet * multi)                # lower period limit of the band
  # Index range of the periods inside every band (periods are sorted).
  lo_idx <- matrix(findInterval(short_p, periods, left.open = TRUE) + 1L, nrow = n_sed)
  hi_idx <- matrix(findInterval(long_p, periods), nrow = n_sed)
  hi_idx[is.na(long_p)] <- 0L
  lo_idx[is.na(short_p)] <- 1L
  win_len <- pmax(hi_idx - lo_idx + 1L, 1L)
  win_lev <- floor(log2(win_len))
  tol <- sqrt(.Machine$double.eps)

  fit <-
    wr_lapply(seq_len(simulations), function(ijk) {
      pk <- periods[maxdetect[, ijk]]
      fits <- matrix(0, nrow = n_sed, ncol = 2)
      if (length(pk) == 0L || n_cyc == 0L) {
        return(fits)
      }
      Ps <- pmax(Powert[, ijk], 0) * probabillity
      rng <- range(Ps)

      # Sparse table for window maxima of the thresholded power.
      nlev <- floor(log2(nper)) + 1L
      st <- vector("list", nlev)
      st[[1]] <- Ps
      if (nlev > 1L) {
        for (lv in 2:nlev) {
          prev <- st[[lv - 1L]]
          len <- nper - 2^(lv - 1L) + 1L
          st[[lv]] <- pmax(prev[1:len], prev[(1:len) + 2^(lv - 2L)])
        }
      }

      # Nearest local maximum for every expected cycle (the first one in case of ties).
      pos <- findInterval(centre, pk)                    # last pk <= centre
      lower <- pmax(pos, 1L)
      upper <- pmin(pos + 1L, length(pk))
      dl <- abs(pk[lower] - centre)
      du <- abs(pk[upper] - centre)
      near <- matrix(ifelse(dl <= du + tol, lower, upper), nrow = n_sed)
      spectra <- matrix(pk[near], nrow = n_sed)
      diffs <- abs(spectra - centre)

      for (ij in seq_len(n_sed)) {
        nr_i <- near[ij, ]
        d_i <- diffs[ij, ]
        # Several cycles can share the same peak: keep the closest (the last one in case of ties).
        keep <- rep(FALSE, n_cyc)
        for (u in unique(nr_i)) {
          idx <- which(nr_i == u)
          dm <- d_i[idx]
          keep[idx[max(which(dm == min(dm)))]] <- TRUE
        }
        sp_i <- spectra[ij, ]
        keep <- keep & sp_i > short_p[ij, ] & sp_i < long_p[ij, ]
        keep[is.na(keep)] <- FALSE
        if (!any(keep)) {
          next
        }
        sel <- which(keep)
        sel <- sel[order(-sp_i[sel])]
        mx <- numeric(length(sel))
        for (s in seq_along(sel)) {
          cc <- sel[s]
          lv <- win_lev[ij, cc]
          t1 <- st[[lv + 1L]]
          mx[s] <- max(t1[lo_idx[ij, cc]], t1[hi_idx[ij, cc] - 2^lv + 1L])
        }
        n_comp <- sum(mx > 0)
        calc <- ((sum(mx) - rng[1]) / (rng[2] - rng[1])) * (n_comp / n_cyc)
        if (!is.finite(calc)) {
          calc <- 0
        }
        fits[ij, 1] <- calc
        fits[ij, 2] <- n_comp
      }
      fits
    }, n_cores = n_cores, run_multicore = run_multicore, combine = "list", progress = verbose)

  power_max_mat <- matrix(NA_real_, ncol = length(my.w$x), nrow = n_sed)
  nr_compenents_mat <- matrix(NA_real_, ncol = length(my.w$x), nrow = n_sed)
  for (kk in seq_along(fit)) {
    power_max_mat[, kk] <- fit[[kk]][, 1]
    nr_compenents_mat[, kk] <- fit[[kk]][, 2]
  }
  testsedrates <- as.data.frame(testsedrates)

  depth <- (my.data[, 1])
  y_axis <- unlist(testsedrates)

  results <- list(
    power_max_mat = power_max_mat,
    nr_compenents_ma = nr_compenents_mat,
    depth = my.data[, 1],
    y_axis = as.numeric(unlist(testsedrates))
  )

  if (genplot == TRUE) {
    wr_open_device(width = 14, height = 7)
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }

    layout.matrix <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
    layout(
      mat = layout.matrix,
      heights = c(1, 1, 1),
      # Heights of the two rows
      widths = c(8, 2, 2)
    )

    par(mar = c(4, 4, 2, 0))

    pmax_avg <- t(results[[plot_res]])


    n.levels = 100
    power_max_mat.levels = quantile(pmax_avg, probs = seq(
      from = 0,
      to = 1,
      length.out = n.levels + 1
    ))


    image.plt = par()$plt
    key.cols <- .wr_key_cols(palette_name, color_brewer, n.levels)

    depth <- (my.data[, 1])
    y_axis <- unlist(testsedrates)
    depth <- as.numeric(depth)
    y_axis <- as.numeric(y_axis)

    image(
      x = depth,
      y = y_axis,
      z = (pmax_avg),
      col = key.cols,
      breaks = power_max_mat.levels,
      xlab = x_lab,
      ylab = y_lab,
      useRaster = TRUE
    )

    r_sum <- colMeans(pmax_avg)

    par(new=FALSE, mar = c(4, 0, 2, 0.5))

    plot(
      y = y_axis,
      x = r_sum,
      type = "l",
      ylim = c(min(y_axis), max(y_axis)),
      yaxs = "i",
      yaxt = "n",
      xlab = "normalized power",
      ylab = y_lab
    )

    par(new=FALSE, mar = c(4, 0.5, 2, 5))

    lwd.axis = 1
    n.ticks = 6
    label.digits = 3
    label.format = "f"
    width = 1.2
    lab.line = 2.5
    lab = NULL

    key.marks = round(seq(
      from = 0,
      to = 1,
      length.out = n.ticks
    ) *
      n.levels)
    key.labels = formatC(as.numeric(power_max_mat.levels),
                         digits = label.digits,
                         format = label.format)[key.marks +
                                                  1]


    image(
      1,
      seq(from = 0, to = n.levels),
      matrix(power_max_mat.levels,
             nrow = 1),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = ""
    )


    axis(
      4,
      lwd = lwd.axis,
      at = key.marks,
      labels = NA,
      tck = 0.02,
      tcl = (par()$usr[2] - par()$usr[1]) *
        width - 0.04
    )
    mtext(
      key.labels,
      side = 4,
      at = key.marks,
      line = 0.5,
      las = 2,
      font = par()$font.axis,
      cex = par()$cex.axis
    )

    box(lwd = lwd.axis)
  }
  return(results)
}

