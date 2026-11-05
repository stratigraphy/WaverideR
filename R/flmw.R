#' @title Fit linear models to spectral peaks extracted from the wavelet spectra to astronomical cycles multiplied by sedimentation rate x
#'
#' @description The \code{\link{flmw}} function is used calculate the linear correlation
#' for a list of astronomical cycles transformed using a range of sedimentation rates and then compared
#' to spectral peaks of a wavelet spectra
#'
#' @param wavelet Wavelet object of a depth record (depth in metres), the output of \code{\link{analyse_wavelet}}
#' (superlet output of class \code{"analyze.superlet"} is also accepted). Required (\code{Default=NULL}).
#' @param sedrate_low Lowest sedimentation rate to test, in cm/kyr (numeric, > 0). Required (\code{Default=NULL}).
#' @param sedrate_high Highest sedimentation rate to test, in cm/kyr (numeric, > \code{sedrate_low}). Required
#' (\code{Default=NULL}).
#' @param spacing Step between the tested sedimentation rates, in cm/kyr (numeric, > 0); the trial rates are
#' \code{seq(sedrate_low, sedrate_high, by = spacing)}. Required (\code{Default=NULL}).
#' @param cycles Periods of the astronomical cycles in kyr (numeric vector, e.g. \code{c(405, 125, 95, 37, 20)}). For
#' each trial sedimentation rate they are converted to periods in metres (\code{cycles * sedrate / 100}) and compared
#' with the wavelet spectrum at every depth. Required (\code{Default=c(NULL)}).
#' @param x_lab Label of the x-axis (depth) of the plot (character). \code{Default="depth"}.
#' @param y_lab Label of the y-axis (sedimentation rate) of the plot (character). \code{Default="sedrate"}.
#' @param run_random Only for \code{method = "classic"}: also run Monte Carlo simulations in which the spectral peak
#' periods are randomly perturbed, to give the percentile of each fit statistic against this null hypothesis (logical,
#' TRUE or FALSE); without it the percentile elements (5-8) of the output are \code{NA}. \code{Default=FALSE}.
#' @param rand_simulations Number of random simulations used when \code{run_random = TRUE} (integer, > 0; more
#' simulations give more stable percentiles but take longer). \code{Default=1000}.
#' @param run_multicore Fit the depth levels in parallel on several cores (logical, TRUE or FALSE); not used by
#' \code{method = "template"}. Results do not depend on it. \code{Default=FALSE}.
#' @param n_cores Number of cores used when \code{run_multicore = TRUE} (integer or \code{NULL}). \code{NULL} gives an
#' automatic, safe number (physical cores minus one, at most 4, limited by the free memory); change this default with
#' \code{options(waverider.n_cores = ...)}. A value > 1 also switches on parallel computing. Results do not depend
#' on the number of cores. \code{Default=NULL}.
#' @param genplot Plot the result matrix (sedimentation rate against depth) with its depth-averaged profile and a
#' colour key (logical, TRUE or FALSE). \code{Default=FALSE}.
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
#' @param plot_res Which result matrix is plotted when \code{genplot = TRUE} (integer 1-8; only used with
#' \code{method = "classic"}, the other methods always plot the score): 1 slope coefficient, 2 r squared,
#' 3 number of matched components, 4 difference to the origin, 5 slope coefficient percentile, 6 r squared
#' percentile, 7 number of components percentile, 8 difference to the origin percentile (5-8 need
#' \code{run_random = TRUE}). \code{Default=2}.
#' @param keep_editable Keep the graphical parameters (\code{par}) of the plot after the function returns, so that
#' features can be added to it (logical, TRUE or FALSE); \code{FALSE} restores the previous settings.
#' \code{Default=FALSE}.
#' @param verbose Show a progress bar over the depth levels (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param method Method used to compare the cycles with the spectrum (character), one of:
#' \itemize{
#' \item \code{"classic"}: the original peak-based linear-fit method (output identical to earlier versions);
#' \item \code{"weighted"}: matches the predicted cycles one-to-one to spectral peaks, weighting by peak power and
#' using a relative (log period) misfit;
#' \item \code{"template"}: correlates the spectrum with Gaussian peaks at the predicted periods (no peak picking,
#' fastest and usually most robust).
#' }
#' \code{Default=c("classic", "weighted", "template")}, i.e. \code{"classic"}.
#' @param tolerance Relative period tolerance used by the \code{"weighted"} and \code{"template"} methods: the standard
#' deviation of the Gaussian misfit in natural-log period (numeric, > 0; 0.1 is about 10 percent). \code{Default=0.1}.
#' @param whiten For the \code{"weighted"} and \code{"template"} methods only: divide the wavelet power of every period
#' by its average over depth before the comparison, which removes the red-noise background (logical, TRUE or FALSE).
#' \code{Default=TRUE}.
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
#' @param colour_brewer British spelling of \code{color_brewer}, with the same options; if given (not \code{NULL}) it overrides
#' \code{color_brewer}. \code{Default=NULL}.
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
    palette_name = "viridis",
    color_brewer = "viridis",
    plot_res = 2,
    keep_editable = FALSE,
    verbose = FALSE,
    method = c("classic", "weighted", "template"),
    tolerance = 0.1,
    whiten = TRUE,
    colour_brewer = NULL) {
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym
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
    col_mat <- function(k) {
      out <- matrix(NA_real_, nrow = nsr, ncol = ndepth)
      for (jf in seq_along(fit)) {
        out[, jf] <- fit[[jf]][, k]
      }
      out
    }
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
