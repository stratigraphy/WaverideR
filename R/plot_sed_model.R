#' @title Plot sedimentation modelling results
#'
#' @description The \code{\link{plot_sed_model}} function is used plot/re-plot the results from the
#' \code{\link{flmw}} and  \code{\link{sum_power_sedrate}} functions
#'
#'@param model_results Wavelet object created using the \code{\link{analyze_wavelet}} function
#'@param plot_res Numbers to be used as input form the \code{\link{flmw}}output
#' options 1-8 option 1: slope coefficient, option 2: r squared,
#'option 3: nr of components, option 4: difference to origin, option 5: slope coefficient percentile
#'option 6: r squared percentile, option 7: nr of components percentile,
#'option 8: difference to origin percentile. If the output of the  \code{\link{sum_power_sedrate}} function is used
#'then input should be option 1: sum max power option 2: nr of components
#'@param x_lab Label for x-axis \code{Default="depth (m)"}
#'@param y_lab Label for y-axis \code{Default=""sed rate cm/kyr""}
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
#'
#'@examples
#'\donttest{
#'#estimate sedimentation rate for the the magnetic susceptibility record
#'# of the Sullivan core of Pas et al., (2018).
#'
#'mag_wt <- analyze_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'#increase n_simulations to better define the red noise spectral power curve
#'mag_wt_red_noise <- model_red_noise_wt(wavelet=mag_wt,
#'n_simulations=100,
#'run_multicore=FALSE,
#'verbose=FALSE)
#'
#'sedrates <- sum_power_sedrate(red_noise=mag_wt_red_noise,
#'wavelet=mag_wt,
#'percentile=0.75,
#'sedrate_low = 0.5,
#'sedrate_high = 4,
#'spacing = 0.05,
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

plot_sed_model  <- function(model_results = NULL,
                            plot_res = 1,
                            x_lab = "depth (m)",
                            y_lab = "sed rate cm/kyr",
                            keep_editable = FALSE,
                            palette_name = "rainbow",
                            color_brewer= "grDevices" ) {
  depth <-   as.numeric(model_results$depth)
  y_axis <- as.numeric(model_results$y_axis)

  dev.new(width = 14, height = 7)
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
  if (color_brewer== "RColorBrewer"){
    key.cols <-   rev(colorRampPalette(brewer.pal(brewer.pal.info[palette_name,1],palette_name))(n.levels))

  }


  if (color_brewer== "colorRamps"){
    color_brewer_Sel <- paste("colorRamps::",palette_name,"(n=n.levels)")
    key.cols = eval(parse(text = color_brewer_Sel))
  }


  if (color_brewer == "grDevices"){
    if (palette_name == "rainbow"){
      color_brewer_Sel <- "grDevices::rainbow(n=n.levels, start = 0, end = 0.7)"
      key.cols <- rev(eval(parse(text = color_brewer_Sel)))
    }
    else if (palette_name == "heat.colors"|
             palette_name == "terrain.colors"|
             palette_name == "topo.colors"|
             palette_name == "cm.colors"){
      color_brewer_Sel <- paste("grDevices::",palette_name,"(n=n.levels, start = 0, end = 1)")
      key.cols <- rev(eval(parse(text = color_brewer_Sel)))
    }
    else{
      key.cols <-  hcl.colors(n=n.levels, palette = palette_name, alpha = NULL, rev = FALSE, fixup = TRUE)}}



  if (color_brewer== "viridis"){
    color_brewer_Sel <- paste("viridis::",palette_name,"(n=n.levels,direction = -1)")
    key.cols = rev(eval(parse(text = color_brewer_Sel)))
  }




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
