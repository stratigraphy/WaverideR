#' @title Plot windowed fft based spectral analysis results
#'
#' @description The \code{\link{plot_win_fft}} function allows for the (re)plotting of the results of the \code{\link{win_fft}}
#'
#'@param win_fft list which is the results of the \code{\link{win_fft}}
#'@param x_lab label for the y-axis \code{Default="depth"}
#'@param y_lab label for the y-axis \code{Default="sedrate"}
#'@param plot_res plot 1 of 8 options option 1: Amplitude matrix,
#'option 2: Power matrix,
#'option 3: Phase matrix,
#'option 4: AR1_CL matrix,
#'option 5: AR1_Fit matrix ,
#'option 6: AR1_90_power matrix,
#'option 7: AR1_95_power matrix,
#'option 8: AR1_99_power matrix, \code{Default=1}
#'@param perc_vis Cutoff percentile when plotting \code{Default=0}
#'@param freq_max Maximum frequency to plot
#'@param freq_min Minimum frequency to plot
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
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis is depth or y axis power  \code{Default=TRUE}
#'@param dev_new Open a new plotting window for the plot. The wavelet plots contain a raster image that needs more room
#'than the small RStudio plot pane, which gives "figure margins too large" errors or clipped text, so a separate window is
#'opened by default in the RStudio GUI. In every other front end (plain R, ESS, Positron, knitr, Rscript) the plot is drawn
#'on the current device and no window is opened. Set TRUE or FALSE to force either behaviour, or set
#'\code{options(waverider.dev_new = TRUE/FALSE)}. Default: TRUE in the RStudio GUI, FALSE elsewhere.
#'
#'@examples
#'\donttest{
#'#Conduct a windowed fft on the magnetic susceptibility record \cr
#'# of the Sullivan core of Pas et al., (2018).
#'
#'mag_win_fft <- win_fft(data= mag,
#'                    padfac = 5,
#'                    window_size = 12.5,
#'                    run_multicore = FALSE,
#'                    genplot = FALSE,
#'                    palette_name = "rainbow",
#'                    color_brewer="grDevices",
#'                    x_lab = c("depth (m)"),
#'                    y_lab = c("frequency cycle/meter"),
#'                    plot_res = 1,
#'                    perc_vis = 0.5,
#'                    freq_max = 5,
#'                    freq_min = 0.001,
#'                    keep_editable=FALSE,
#'                    verbose=FALSE)
#'
#'# Plot the amplitude spectra
#'plot_win_fft(win_fft= mag_win_fft,
#'x_lab = c("depth (m)"),
#'y_lab = c("frequency cycle/meter"),
#'plot_res = 1,
#'perc_vis = 0.5,
#'freq_max = 5,
#'freq_min = 0.001,
#'keep_editable=FALSE,
#'palette_name = "rainbow",
#'color_brewer="grDevices",
#'plot_horizontal=TRUE,
#'dev_new=TRUE)
#'
#'}
#'
#' @return
#'Returns a plot of, which plot 1 of 8 options,
#'option 1: Amplitude matrix
#'option 2: Power matrix
#'option 3: Phase matrix
#'option 4: AR1_CL matrix
#'option 5: AR1_Fit matrix
#'option 6: AR1_90_power matrix
#'option 7: AR1_95_power matrix
#'option 8: AR1_99_power matrix
#'
#' @param colour_brewer British spelling of \code{color_brewer}; if given it overrides \code{color_brewer}.
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






plot_win_fft <- function(win_fft = NULL,
                         x_lab = c("depth (m)"),
                         y_lab = c("frequency cycle/metre"),
                         plot_res = 1,
                         perc_vis = 0,
                         freq_max = NULL,
                         freq_min = NULL,
                         keep_editable = FALSE,
                         palette_name = "rainbow",
                         color_brewer="grDevices",
                         plot_horizontal=TRUE,
                         dev_new = wr_dev_new_default(),
  colour_brewer = NULL) {
  if (!is.null(colour_brewer)) color_brewer <- colour_brewer  # British spelling synonym

  results <- win_fft
  n.levels = 100
  y_axis <- as.numeric(unlist(results$y_axis))
  if (is.null(freq_max)) freq_max <- max(y_axis) + 1e-9
  if (is.null(freq_min)) freq_min <- min(y_axis) - 1e-9
  sel_cols_up <- max(which(y_axis < freq_max))
  sel_cols_down <- min(which(y_axis > freq_min))

  bottom_perc <- perc_vis

  pmax_avg_sel <- t(results[[plot_res]])
  pmax_avg_sel <- pmax_avg_sel[, sel_cols_down:sel_cols_up]


  if (dev_new==TRUE & plot_horizontal==TRUE){
    wr_open_device(width = 14, height = 7)

  }

  if (dev_new==TRUE & plot_horizontal==FALSE){
    wr_open_device(width = 7, height = 10)

  }


  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  key.cols <- .wr_key_cols(palette_name, color_brewer, n.levels)




  if (plot_horizontal==TRUE){



  layout.matrix <- matrix(c(3, 1, 2), nrow = 1, ncol = 3)
  layout(mat = layout.matrix,
         heights = c(1, 1, 1),
         # Heights of the two rows
         widths = c(8, 2, 2))

  par(mar = c(4, 4, 2, 3))
  power_max_mat.levels = quantile(pmax_avg_sel,
                                  probs = seq(
                                    from = bottom_perc,
                                    to = 1,
                                    length.out = n.levels + 1
                                  ))


  image.plt = par()$plt
  depth <-  results$depth
  y_axis <- results$y_axis
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)
  y_axis <- y_axis[sel_cols_down:sel_cols_up]


  r_sum <- colMeans(pmax_avg_sel)
  plot(
    y = y_axis,
    x = r_sum,
    type = "l",
    ylim = c(min(y_axis), max(y_axis)),
    yaxs = "i",
    xlab = "mean valure",
    ylab = y_lab
  )

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

  image(
    x = depth,
    y = y_axis,
    z = (pmax_avg_sel),
    col = key.cols,
    breaks = power_max_mat.levels,
    xlab = x_lab,
    ylab = y_lab,
    useRaster = TRUE
  )
}

  if (plot_horizontal==FALSE){

    layout.matrix <- matrix(c(1, 2,0,3), nrow = 2, ncol = 2,byrow=TRUE)
    layout(mat = layout.matrix,
           heights = c(1, 4),
           # Heights of the two rows
           widths = c(1, 3))

    par(mar = c(0, 1, 2, 6))
    power_max_mat.levels = quantile(pmax_avg_sel,
                                    probs = seq(
                                      from = bottom_perc,
                                      to = 1,
                                      length.out = n.levels + 1
                                    ))


    image.plt = par()$plt

    depth <-  results$depth
    y_axis <- results$y_axis
    depth <- as.numeric(depth)
    y_axis <- as.numeric(y_axis)
    y_axis <- y_axis[sel_cols_down:sel_cols_up]


    r_sum <- colMeans(pmax_avg_sel)

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



    image(y=1,
      x=seq(from = 0, to = n.levels),
      t(matrix(power_max_mat.levels,
             nrow = 1)),
      col = key.cols,
      breaks = power_max_mat.levels,
      useRaster = TRUE,
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = ""
    )

    axis(
      1,
      lwd = lwd.axis,
      at = key.marks,
      labels = NA,
      tck = 0.02,
      tcl = 1
    )


    mtext(
      key.labels,
      side = 1,
      at = key.marks,
      line = 0.5,
      las = 2,
      font = par()$font.axis,
      cex = 0.75
    )

    title(xlab="Power",xpd=NA)

    box(lwd = lwd.axis)


    par(mar = c(0, 0, 4, 2))

    plot(
      x = y_axis,
      y = r_sum,
      type = "l",
      xlim = c(min(y_axis), max(y_axis)),
      xaxs = "i",
      ylab = "",
      xlab = y_lab
    )

    title(ylab="mean power",xpd=NA)


    par(mar = c(4, 0, 0, 2))



    image(
      y = depth,
      x = y_axis,
      z = t(pmax_avg_sel),
      col = key.cols,
      breaks = power_max_mat.levels,
      ylab = "",
      xlab = y_lab,
      useRaster = TRUE
    )

    title(ylab=x_lab,xpd=NA)


  }

  }
