#' @title Plot windowed fft based spectral analysis results
#'
#' @description The \code{\link{plot_win_fft}} function allows for the (re)plotting of the results of the \code{\link{win_fft}}
#'
#'@param win_fft Result of the windowed FFT (required): the list returned by \code{\link{win_fft}}. \code{Default=NULL}.
#'@param x_lab Label of the depth/time axis (character). \code{Default="depth (m)"}.
#'@param y_lab Label of the frequency axis (character). \code{Default="frequency cycle/metre"}.
#'@param plot_res Which matrix of \code{win_fft} to plot, given as its position in the list (integer, 1 to 9):
#'\itemize{
#'\item 1: amplitude (\code{Amplitude_mat}),
#'\item 2: power (\code{Power_mat}),
#'\item 3: phase (\code{Phase_mat}),
#'\item 4: frequency uncertainty (\code{uncer_mat}),
#'\item 5: AR1 confidence level in \% (\code{AR1_CL_mat}),
#'\item 6: AR1 red-noise fit (\code{AR1_Fit_mat}),
#'\item 7: AR1 90\% confidence power (\code{AR1_90_power_mat}),
#'\item 8: AR1 95\% confidence power (\code{AR1_95_power_mat}),
#'\item 9: AR1 99\% confidence power (\code{AR1_99_power_mat}).
#'}
#'\code{Default=1}.
#'@param perc_vis Lower quantile of the colour scale (numeric, between 0 and 1): the colour breaks run from this
#'quantile of the plotted values to the maximum, and values below it are left blank, which hides the background.
#'\code{0} colours all values. \code{Default=0}.
#'@param freq_max Highest frequency plotted, in cycles per depth/time unit (numeric, or \code{NULL} for the highest
#'frequency of \code{win_fft}). \code{Default=NULL}.
#'@param freq_min Lowest frequency plotted, in cycles per depth/time unit (numeric, or \code{NULL} for the lowest
#'frequency of \code{win_fft}). \code{Default=NULL}.
#'@param keep_editable Keep the graphical parameters (\code{par}, layout) of the plot after the function returns, so that
#'extra elements can be added to the plot afterwards (logical). With \code{FALSE} the previous graphical parameters
#'are restored on exit. \code{Default=FALSE}.
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
#'@param plot_horizontal Orientation of the plot (logical): \code{TRUE} puts depth/time on the x-axis and frequency on
#'the y-axis, with the mean spectrum to the left; \code{FALSE} draws depth/time on the y-axis and frequency on the
#'x-axis, with the mean spectrum on top. \code{Default=TRUE}.
#'@param dev_new Open a new plotting window for the plot. The wavelet plots contain a raster image that needs more room
#'than the small RStudio plot pane, which gives "figure margins too large" errors or clipped text, so a separate window is
#'opened by default in the RStudio GUI. In every other front end (plain R, ESS, Positron, knitr, Rscript) the plot is drawn
#'on the current device and no window is opened. Set TRUE or FALSE to force either behaviour, or set
#'\code{options(waverider.dev_new = TRUE/FALSE)}. Logical.
#'\code{Default=wr_dev_new_default()}: TRUE in the RStudio GUI, FALSE elsewhere.
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






plot_win_fft <- function(win_fft = NULL,
    x_lab = c("depth (m)"),
    y_lab = c("frequency cycle/metre"),
    plot_res = 1,
    perc_vis = 0,
    freq_max = NULL,
    freq_min = NULL,
    keep_editable = FALSE,
    palette_name = "viridis",
    color_brewer = "viridis",
    plot_horizontal=TRUE,
    dev_new = wr_dev_new_default(),
    colour_brewer = NULL) {
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym

  results <- win_fft
  n.levels = 100
  y_axis <- as.numeric(unlist(results$y_axis))
  if (is.null(freq_max)) {
    freq_max <- max(y_axis) + 1e-9
  }
  if (is.null(freq_min)) {
    freq_min <- min(y_axis) - 1e-9
  }
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
