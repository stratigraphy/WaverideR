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
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis is depth or y axis power  \code{Default=TRUE}
#'@param dev_new Opens a new plotting window to plot the plot, this guarantees  a "nice" looking plot however when plotting in an R markdown
#'document the plot might not plot  \code{Default=TRUE}
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
#' @export
#' @importFrom Matrix rowMeans
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
                         dev_new=TRUE) {

  results <- win_fft
  n.levels = 100
  y_axis <- as.numeric(unlist(results$y_axis))
  sel_cols_up <- max(which(y_axis < freq_max))
  sel_cols_down <- min(which(y_axis > freq_min))

  bottom_perc <- perc_vis

  pmax_avg_sel <- t(results[[plot_res]])
  pmax_avg_sel <- pmax_avg_sel[, sel_cols_down:sel_cols_up]


  if (dev_new==TRUE & plot_horizontal==TRUE){
    dev.new(width = 14, height = 7)

  }

  if (dev_new==TRUE & plot_horizontal==FALSE){
    dev.new(width = 7, height = 10)

  }


  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

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
