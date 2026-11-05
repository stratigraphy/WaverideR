#' @title Plot sedimentation modelling results
#'
#' @description The \code{\link{plot_sed_model}} function is used plot/re-plot the results from the
#' \code{\link{flmw}} and  \code{\link{sum_power_sedrate}} functions
#'
#' @param model_results Sedimentation-rate modelling result to plot: the list returned by \code{\link{flmw}} or by
#' \code{\link{sum_power_sedrate}} (a list of depth x sedimentation-rate matrices plus \code{depth} and
#' \code{y_axis}). Required (\code{Default=NULL} gives an error).
#' @param plot_res Which matrix of \code{model_results} is plotted (integer index of the list element, or its name).
#' For \code{\link{flmw}} with \code{method = "classic"}:
#' \itemize{
#' \item 1: slope coefficient, 2: r squared, 3: number of components, 4: difference to origin;
#' \item 5: slope coefficient percentile, 6: r squared percentile, 7: number of components percentile,
#' 8: difference to origin percentile;
#' \item 11 (or \code{"score_mat"}): the combined score.
#' }
#' For \code{\link{flmw}} with \code{method = "weighted"} or \code{"template"}: 1 (the score matrix). For
#' \code{\link{sum_power_sedrate}}: 1, sum of the maximum spectral power; 2, number of components (cycles).
#' \code{Default=1}.
#' @param x_lab Label of the depth axis (character). \code{Default="depth (m)"}.
#' @param y_lab Label of the sedimentation-rate axis (character). \code{Default="sed rate cm/kyr"}.
#' @param keep_editable Keep the plot layout and graphical parameters (\code{par}) after plotting, so that extra
#' elements can be added to the image (logical). With \code{FALSE} the previous \code{par} settings are restored on
#' exit. \code{Default=FALSE}.
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
#'n_simulations = 3, # increase for a better constrained result
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
#'palette_name = "rainbow",
#'color_brewer= "grDevices",
#'verbose=FALSE)
#'
#'plot_sed_model(model_results=sedrates,
#'plot_res=1,
#'x_lab = "depth (m)",
#'y_lab = "sed rate cm/kyr",
#'keep_editable=FALSE,
#'palette_name = "rainbow",
#'color_brewer= "grDevices")
#'}
#'
#'
#' @return
#'Returns a plot of sedimentation rates vs depth and a value which was generated using
#'the \code{\link{flmw}} or \code{\link{sum_power_sedrate}} functions
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
#' @importFrom grDevices rgb
#' @importFrom graphics layout
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors

plot_sed_model  <- function(model_results = NULL,
    plot_res = 1,
    x_lab = "depth (m)",
    y_lab = "sed rate cm/kyr",
    keep_editable = FALSE,
    palette_name = "viridis",
    color_brewer = "viridis",
    colour_brewer = NULL) {
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym
  depth <-   as.numeric(model_results$depth)
  y_axis <- as.numeric(model_results$y_axis)

  wr_open_device(width = 14, height = 7)
  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  layout.matrix <- matrix(c(3, 1, 2), nrow = 1, ncol = 3)
  layout(mat = layout.matrix,
      heights = c(1, 1, 1),
      # Heights of the two rows
      widths = c(8, 2, 2))

  par(mar = c(4, 4, 2, 4))
  pmax_avg <- t(as.matrix(model_results[[plot_res]]))

  n.levels = 100

  power_max_mat.levels = quantile(pmax_avg, probs = seq(
      from = 0,
      to = 1,
      length.out = n.levels + 1
  ))

  image.plt = par()$plt
  key.cols <- .wr_key_cols(palette_name, color_brewer, n.levels)




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
}
