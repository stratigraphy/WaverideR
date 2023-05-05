#' @title plot windowed fft based spectral analysis results
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
#'
#'@examples
#'\donttest{
#'#Conduct a windowed ftt on the magnetic susceptibility record \cr
#'# of the Sullivan core of Pas et al., (2018).
#'
#'mag_win_fft <- win_fft(data= mag,
#'                    padfac = 5,
#'                    window_size = 12.5,
#'                    run_multicore = FALSE,
#'                    genplot = FALSE,
#'                    x_lab = c("depth (m)"),
#'                    y_lab = c("frequency cycle/metre"),
#'                    plot_res = 1,
#'                    perc_vis = 0.5,
#'                    freq_max = 5,
#'                    freq_min = 0.001,
#'                    keep_editable=FALSE)
#'
#'# Plot the mplitude spectra
#'plot_win_fft(win_fft= mag_win_fft,
#'x_lab = c("depth (m)"),
#'y_lab = c("frequency cycle/metre"),
#'plot_res = 1,
#'perc_vis = 0.5,
#'freq_max = 5,
#'freq_min = 0.001,
#'keep_editable=FALSE)
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


plot_win_fft <- function(win_fft = NULL,
                         x_lab = c("depth (m)"),
                         y_lab = c("frequency cycle/metre"),
                         plot_res = 1,
                         perc_vis = 0,
                         freq_max = NULL,
                         freq_min = NULL,
                         keep_editable = FALSE) {
  results <- win_fft

  y_axis <- as.numeric(unlist(results$y_axis))
  sel_cols_up <- max(which(y_axis < freq_max))
  sel_cols_down <- min(which(y_axis > freq_min))

  bottom_perc <- perc_vis

  pmax_avg_sel <- t(results[[plot_res]])
  pmax_avg_sel <- pmax_avg_sel[, sel_cols_down:sel_cols_up]

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

  par(mar = c(4, 4, 2, 3))
  n.levels = 100
  power_max_mat.levels = quantile(pmax_avg_sel,
                                  probs = seq(
                                    from = bottom_perc,
                                    to = 1,
                                    length.out = n.levels + 1
                                  ))


  image.plt = par()$plt
  color.palette = "rainbow(n.levels, start = 0, end = 0.7)"
  key.cols = rev(eval(parse(text = color.palette)))


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
