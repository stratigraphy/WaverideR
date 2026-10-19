#' @title Fit linear models to spectral peaks extracted from the wavelet spectra to astronomical cycles multiplied by sedimentation rate x
#'
#' @description The \code{\link{flmw}} function is used calculate the linear correlation
#' for a list of astronomical cycles transformed using a range of sedimentation rates and then compared
#' to spectral peaks of a wavelet spectra
#'
#'@param wavelet Wavelet object created using the \code{\link{analyse_wavelet}} function
#'@param sedrate_low Minimum sedimentation rate (cm/kyr)for which the sum of maximum spectral power is calculated for.
#'@param sedrate_high Maximum sedimentation rate (cm/kyr) for which the sum of maximum spectral power is calculated  for.
#'@param spacing Spacing (cm/kyr) between sedimentation rates
#'@param cycles Astronomical cycles (in kyr) for which the combined sum of maximum spectral power is calculated for
#'@param x_lab label for the y-axis \code{Default="depth"}
#'@param y_lab label for the y-axis \code{Default="sedrate"}
#'@param run_random run multiple simulation to calculate percentile against the 0 hypothesis
#'@param rand_simulations nr of simulations to calculate percentile against the 0 hypothesis
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
#'@param plot_res options 1-8 (only used with \code{method = "classic"}; the other methods always plot the score)
#'option 1: slope coefficient, option 2: r squared,
#'option 3: nr of components, option 4: difference to the  origin , option 5: slope coefficient percentile
#'option 6: r squared percentile, option 7: nr of components percentile,
#'option 8: difference to the origin percentile \code{Default=2}
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#' @param verbose Print text \code{Default=FALSE}.
#'@param method Method used to compare the cycles with the spectrum. \code{"classic"} (default) is the original
#'peak based method (output identical to earlier versions); \code{"weighted"} matches predicted cycles one-to-one to
#'spectral peaks, weighting by peak power and using a relative (log period) misfit; \code{"template"} correlates the
#'spectrum with Gaussian peaks at the predicted periods (no peak picking, fastest and usually most robust).
#'@param tolerance Relative period tolerance (standard deviation in log period, 0.1 is about 10 percent) used by the
#'\code{"weighted"} and \code{"template"} methods \code{Default=0.1}.
#'@param whiten Divide the wavelet power of every period by its average over depth before comparing (removes the red
#'background) for the \code{"weighted"} and \code{"template"} methods \code{Default=TRUE}.
#'
#' @author
#'Based on the \link[astrochron]{eAsm} function of the 'astrochron' R package and the 'eCOCO' and 'COCO' function of the 'Acycle' software
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
#'#estimate sedimentation rate for the magnetic susceptibility record
#'# of the Sullivan core of Pas et al., (2018).
#'
#'mag_wt <- analyse_wavelet(data = mag[1:200, ],
#' dj = 1/10,
#' lowerPeriod = 0.1,
#' upperPeriod = 32,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'sedrates <- flmw(wavelet = mag_wt,
#'     sedrate_low = 0.5,
#'     sedrate_high = 4,
#'     spacing = 1,
#'     cycles = c(2376,1600,1180,696,406,110),
#'     x_lab = "depth",
#'     y_lab = "sedrate",
#'     run_random = FALSE,
#'     rand_simulations = 50, # increase to get better constrainted resutls
#'     run_multicore = FALSE,
#'     genplot = FALSE,
#'     palette_name = "rainbow",
#'     color_brewer = "grDevices",
#'     plot_res = 2,
#'     keep_editable=FALSE,
#'     verbose=FALSE)
#'
#'# template method (no peak picking) followed by a smooth sedimentation rate curve
#'sed_t <- flmw(wavelet = mag_wt, sedrate_low = 0.5, sedrate_high = 4, spacing = 0.5,
#'     cycles = c(406, 110), method = "template")
#'sed_path <- flmw_path(sed_t)
#'}
#'
#' @section Choosing a method:
#' \code{"classic"} keeps the earlier results reproducible. In a synthetic test with strong noise and a sedimentation rate
#' that changes with depth, \code{"weighted"} and \code{"template"} followed by \code{\link{flmw_path}} recovered the
#' true curve far better than the classic score (see the NEWS file). \code{"template"} needs no peak picking and is the
#' fastest.
#'
#' @return
#'Returns a list. For \code{method = "classic"} it has the 10 elements below plus
#'\code{score_mat} (a heuristic combined score between 0 and 1), \code{best_sedrate} (sedimentation rate with the highest
#'score per depth) and \code{method}. For the other methods the list has \code{score_mat}, \code{best_sedrate}, \code{depth},
#'\code{y_axis} and \code{method}. Use \code{\link{flmw_path}} to turn \code{score_mat} into a smooth sedimentation
#'rate curve with uncertainty.
#'
#'Elements of the classic result:
#'element 1: slope coefficient
#'element 2: r squared
#'element 3: nr of components
#'element 4: difference to the origin
#'element 5: slope coefficient percentile
#'element 6: r squared percentile
#'element 7: nr of components percentile,
#'element 8: difference to the origin percentile
#'element 9: y-axis values of the matrices which is sedimentation rate
#'element 10: x-axis values of the matrices which is depth
#'
#' @param colour_brewer British spelling of \code{color_brewer}; if given it overrides \code{color_brewer}.
#' @export
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats ecdf
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

flmw  <- function(wavelet = NULL,
                  sedrate_low = NULL,
                  sedrate_high = NULL,
                  spacing = NULL,
                  cycles = c(NULL),
                  x_lab = "depth",
                  y_lab = "sedrate",
                  run_random = FALSE,
                  rand_simulations = 1000,
                  run_multicore = FALSE,
                  n_cores = NULL,
                  genplot = FALSE,
                  palette_name = "rainbow",
                  color_brewer = "grDevices",
                  plot_res = 2,
                  keep_editable = FALSE,
                  verbose = FALSE,
                  method = c("classic", "weighted", "template"),
                  tolerance = 0.1,
                  whiten = TRUE,
  colour_brewer = NULL) {
  if (!is.null(colour_brewer)) color_brewer <- colour_brewer  # British spelling synonym
  method <- match.arg(method)
  my.w <- wavelet
  if (inherits(my.w, "analyze.superlet")) {                    # superlet power is stored as depth x period
    Powert <- t(my.w$Power)
  } else {
    Powert <- matrix(as.numeric(unlist(as.data.frame(my.w$Power))), nrow = my.w$nr, ncol = my.w$nc)
  }
  period <- my.w$Period
  depth <- as.numeric(my.w$x)
  testsedrates <- seq(from = sedrate_low, to = sedrate_high, by = spacing)
  nsr <- length(testsedrates)
  ndepth <- ncol(Powert)

  # power relative to the average spectrum (removes the red background)
  if (method != "classic" && whiten) {
    bg <- rowMeans(Powert, na.rm = TRUE)
    bg[!is.finite(bg) | bg <= 0] <- 1
    Powert <- Powert / bg
  }

  if (method == "classic") {
    if (run_random == TRUE) {
      randomize_mat <-
        matrix(
          data = wr_rtruncnorm(
            n = (nrow(Powert) / 2) * rand_simulations,
            a = 1 / 10,
            b = 10,
            mean = 1,
            sd = 1
          ),
          ncol = rand_simulations
        )
    } else {
      randomize_mat <- NULL
    }
    fit <- wr_lapply(seq_len(ndepth), function(ijk) {
      .flmw_depth_fit(period, Powert[, ijk], cycles, testsedrates, run_random, rand_simulations, randomize_mat)
    }, n_cores = n_cores, run_multicore = run_multicore, combine = "list", progress = verbose)
    col_mat <- function(k) matrix(vapply(fit, function(m) m[, k], numeric(nsr)), nrow = nsr, ncol = ndepth)
    results <- list(
      slope_coeff_mat = col_mat(1),
      corr_coeff_mat = col_mat(2),
      nr_compenents_mat = col_mat(3),
      orig_mat = col_mat(4),
      slope_coeff_mat_per = col_mat(5),
      corr_coeff_mat_per = col_mat(6),
      nr_compenents_mat_per = col_mat(7),
      orig_mat_per = col_mat(8),
      depth = depth,
      y_axis = as.numeric(testsedrates)
    )
    # heuristic combined score: fit quality x slope agreement x share of matched cycles
    sc <- results$corr_coeff_mat * (1 - pmin(results$slope_coeff_mat, 1)) * (results$nr_compenents_mat / length(cycles))
    sc[!is.finite(sc)] <- 0
    sc[results$nr_compenents_mat < 2] <- 0
    results$score_mat <- sc
    plot_mat <- results[[plot_res]]
  } else {
    if (method == "weighted") {
      sc <- wr_lapply(seq_len(ndepth), function(ijk) {
        .flmw_weighted_depth(period, Powert[, ijk], cycles, testsedrates, tolerance)
      }, n_cores = n_cores, run_multicore = run_multicore, combine = "list", progress = verbose)
      sc <- matrix(unlist(sc), nrow = nsr, ncol = ndepth)
    } else {
      sc <- .flmw_template_scores(period, Powert, cycles, testsedrates, tolerance)
    }
    results <- list(score_mat = sc, depth = depth, y_axis = as.numeric(testsedrates))
    plot_mat <- sc
  }
  best <- max.col(t(results$score_mat), ties.method = "first")
  results$best_sedrate <- as.numeric(testsedrates)[best]
  results$method <- method
  if (genplot == TRUE) {
    .flmw_plot(plot_mat, depth, as.numeric(testsedrates), x_lab, y_lab, palette_name, color_brewer, keep_editable)
  }
  return(results)
}

# Plot of a (sedimentation rate x depth) matrix with the sum of the matrix next to it and a colour key
.flmw_plot <- function(mat, depth, y_axis, x_lab, y_lab, palette_name, color_brewer, keep_editable) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }


    wr_open_device(width = 14, height = 7)
    layout.matrix <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
    layout(
      mat = layout.matrix,
      heights = c(1, 1, 1),
      # Heights of the two rows
      widths = c(8, 2, 2)
    )

    par(mar = c(4, 4, 2, 4))
    pmax_avg <- t(as.matrix(mat))

    n.levels = 100


    # common colour rule of the package: quantiles, or log for a mostly empty map (e.g. method = "template")
    power_max_mat.levels = as.numeric(.wr_colour_breaks(pmax_avg, "auto", n.levels))


    image.plt = par()$plt
    key.cols <- .wr_key_cols(palette_name, color_brewer, n.levels)


    depth <- (depth)
    y_axis <- y_axis
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
    plot(
      y = y_axis,
      x = r_sum,
      type = "l",
      ylim = c(min(y_axis), max(y_axis)),
      yaxs = "i",
      xlab = "normalized power",
      ylab = y_lab
    )

    lwd.axis = 1
    n.ticks = 6
    label.digits = 6
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
