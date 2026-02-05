#' @title Plots a cross superlet scalogram
#'
#' @description Plot cross superlet scalogram using the outcome of the \code{\link{analyze_Xsuperlet}} function.
#'
#'@param Xsuperlet Xsuperlet object created using the \code{\link{analyze_Xsuperlet}} function.
#'@param lowerPeriod Lowest period value which will be plotted
#'@param upperPeriod Highest period value which will be plotted
#'@param n.levels  Number of color levels \code{Default=100}.
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
#'@param useRaster Plot as a raster or vector image \code{Default=TRUE}.
#'WARNING plotting as a vector image is computationally intensive.
#'@param periodlab Label for the y-axis \code{Default="Period (metres)"}.
#'@param x_lab Label for the x-axis \code{Default="depth (metres)"}.
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'@param dev_new Opens a new plotting window to plot the plot, this guarantees a "nice" looking plot however when plotting in an R markdown
#'document the plot might not plot  \code{Default=TRUE}
#'@param plot_dir The direction of the proxy record which is assumed for tuning if time increases with increasing depth/time values
#'(e.g. bore hole data which gets older with increasing depth ) then plot_dir should be set to TRUE
#'if time decreases with depth/time values (eg stratospheric logs where 0m is the bottom of the section)
#'then plot_dir should be set to FALSE \code{plot_dir=TRUE}
#'@param plot_arrows add the phase arrows of the cross superlet transform \code{Default=TRUE}
#'@param pwr_quant Power quantile above which the arrows are plotted  \code{Default= 0.85}
#'@param n_arrows number arrows width wise plotted on the time axis \code{Default= 50}
#'@param T_unit length of the arrows in time units \code{Default= 15}
#'@param P_unit length of the arrows in period units \code{Default = 0.5}
#'@param arrow_length length of the slanted part of the arrow  \code{Default= 0.1}
#'@param arrow_lwd thickness of the arrow \code{Default= 2}
#'@param arrow_angle angle of the slanted part of the arrow \code{Default= 25}
#'@param add_lines Add  lines to the wavelet plot input should be matrix with first axis being depth/time the columns after that
#'should be period values  \code{Default=NULL}
#'@param add_points Add points to the wavelet plot input should be matrix with first axis being depth/time and columns after that
#'should be period values \code{Default=NULL}
#'@param add_abline_h Add horizontal lines to the plot. Specify the lines as a vector e.g. c(2,3,5,6)  \code{Default=NULL}
#'@param add_abline_v Add vertical lines to the plot. Specify the lines as a vector e.g. c(2,3,5,6)  \code{Default=NULL}
#'@param add_data Plot the data on top of the wavelet \code{Default=TRUE}
#'@param add_avg Plot the average wavelet spectral power to the side of the wavelet \code{Default=FALSE}
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis is depth or y axis power \code{Default=TRUE}
#'
#'
#' @return
#' The output is a plot of a cross superlet scalogram.
#'
#' @author
#' plotting code based on the "wt.image" function of the 'WaveletComp' R package
#' Whereas the "analyze_Xsuperlet" that generates the input for the plotting function
#' is based on the matlab code in Moca et al. (2021) and the the "analyze.coherency" function of the  'WaveletComp' R package
#'
#' @references
#'Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
#'Wavelet Analysis. R package version 1.1.
#'\url{https://CRAN.R-project.org/package=WaveletComp}
#'
#'Moca, Vasile V., Harald Bârzan, Adriana Nagy-Dăbâcan, and Raul C. Mureșan.
#'Time-frequency super-resolution with superlets.
#'Nature communications 12, no. 1 (2021): 337.
#'\url{https://doi.org/10.1038/s41467-020-20539-9}
#'
#'@examples
#' \donttest{
#'#Example 1. A cross superlet plot of two etp solutions with noise overprint
#'
#'
#'etp_1 <- etp(
#'  tmin = 0,
#'  tmax = 1500,
#'  dt = 1,
#'  eWt = 1.5,
#'  oWt = 0.75,
#'  pWt = 1,
#'  esinw = TRUE,
#'  standardize = TRUE,
#'  genplot = FALSE,
#'  verbose = FALSE
#')
#'
#'etp_2 <- etp(
#'  tmin = 0,
#'  tmax = 1500,
#'  dt = 1,
#'  eWt = 1,
#'  oWt = 0.5,
#'  pWt = 1.5,
#'  esinw = TRUE,
#'  standardize = TRUE,
#'  genplot = FALSE,
#' verbose = FALSE
#')
#'
#'etp_1[, 2] <- etp_1[, 2] + colorednoise::colored_noise(
#'  nrow(etp_1),
#'  sd = sd(etp_1[, 2]) / 1.5,
#'  mean = mean(etp_1[, 2]),
#'  phi = 0.9
#')
#'etp_2[, 2] <- etp_2[, 2] + colorednoise::colored_noise(
#'  nrow(etp_2),
#'  sd = sd(etp_2[, 2]) / 1.5,
#'  mean = mean(etp_2[, 2]),
#'  phi = 0.9
#')
#'
#'X_super_etp <- analyze_Xsuperlet(
#'  data_1 = etp_1,
#'  data_2  = etp_2,
#'  upperPeriod = 1024,
#'  lowerPeriod = 2,
#'  Nf = 128,
#'  c1 = 3,
#'  o = c(1, 10),
#'  mult = TRUE,
#'  verbose = FALSE
#')
#'plot_Xsuperlet(Xsuperlet = X_super_etp, lowerPeriod = 2, upperPeriod = 1024,
#'               n.levels = 100, palette_name = "rainbow", color_brewer = "grDevices",
#'               useRaster = TRUE, periodlab = "Period (metres)", x_lab = "depth (metres)",
#'               keep_editable = FALSE, dev_new = TRUE, plot_dir = TRUE,
#'               plot_arrows = TRUE,
#'               pwr_quant = 0.85,
#'               n_arrows = 50,
#'               T_unit = 15,
#'               P_unit  = 0.5,
#'               arrow_length = 0.1,
#'               arrow_lwd = 2,
#'               arrow_angle = 25,add_lines = NULL, add_points = NULL, add_abline_h = NULL,
#'               add_abline_v = NULL, add_data = TRUE, add_avg = FALSE, plot_horizontal = TRUE)
#'}
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
#' @importFrom graphics layout
#' @importFrom graphics title
#' @importFrom grDevices rgb
#' @importFrom DescTools Closest
#' @importFrom graphics abline
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

plot_Xsuperlet <- function (Xsuperlet = NULL, lowerPeriod = NULL, upperPeriod = NULL,
                            n.levels = 100, palette_name = "rainbow", color_brewer = "grDevices",
                            useRaster = TRUE, periodlab = "Period (metres)", x_lab = "depth (metres)",
                            keep_editable = FALSE, dev_new = TRUE, plot_dir = TRUE,
                            plot_arrows = TRUE,
                            pwr_quant = 0.85,
                            n_arrows = 75,
                            T_unit = 15,
                            P_unit  = 0.5,
                            arrow_length = 0.075,
                            arrow_lwd = 1.9,
                            arrow_angle = 25,add_lines = NULL, add_points = NULL, add_abline_h = NULL,
                            add_abline_v = NULL, add_data = TRUE, add_avg = FALSE, plot_horizontal = TRUE)
{

  if (keep_editable == FALSE) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  Power <- t(Xsuperlet$Power)
  Phase <- t(Xsuperlet$Phase)
  maximum.level = max(Xsuperlet$Power, na.rm = TRUE)
  power_max_mat.levels = quantile(Xsuperlet$Power, probs = seq(from = 0,
                                                               to = 1, length.out = n.levels + 1))
  if (color_brewer == "RColorBrewer") {
    key.cols <- rev(colorRampPalette(brewer.pal(brewer.pal.info[palette_name,
                                                                1], palette_name))(n.levels))
  }
  if (color_brewer == "colorRamps") {
    color_brewer_Sel <- paste("colorRamps::", palette_name,
                              "(n=n.levels)")
    key.cols = eval(parse(text = color_brewer_Sel))
  }
  if (color_brewer == "grDevices") {
    if (palette_name == "rainbow") {
      color_brewer_Sel <- "grDevices::rainbow(n=n.levels, start = 0, end = 0.7)"
      key.cols <- rev(eval(parse(text = color_brewer_Sel)))
    }
    else if (palette_name == "heat.colors" | palette_name ==
             "terrain.colors" | palette_name == "topo.colors" |
             palette_name == "cm.colors") {
      color_brewer_Sel <- paste("grDevices::", palette_name,
                                "(n=n.levels, start = 0, end = 1)")
      key.cols <- rev(eval(parse(text = color_brewer_Sel)))
    }
    else {
      key.cols <- hcl.colors(n = n.levels, palette = palette_name,
                             alpha = NULL, rev = FALSE, fixup = TRUE)
    }
  }
  if (color_brewer == "viridis") {
    color_brewer_Sel <- paste("viridis::", palette_name,
                              "(n=n.levels,direction = -1)")
    key.cols = rev(eval(parse(text = color_brewer_Sel)))
  }
  periodtck = 0.02
  periodtcl = 0.5
  main = NULL
  lwd = 2
  lwd.axis = 1
  legend.params = list(width = 1.2, shrink = 0.9, mar = 5.1,
                       n.ticks = 6, label.digits = 3, label.format = "f", lab = NULL,
                       lab.line = 2.5)
  key.marks = round(seq(from = 0, to = 1, length.out = legend.params$n.ticks) *
                      n.levels)
  key.labels = formatC(as.numeric(power_max_mat.levels), digits = legend.params$label.digits,
                       format = legend.params$label.format)[key.marks + 1]
  if (dev_new == TRUE & plot_horizontal == TRUE) {
    dev.new(width = 15, height = 7, noRStudioGD = TRUE)
  }
  if (dev_new == TRUE & plot_horizontal == FALSE) {
    dev.new(width = 7, height = 10, noRStudioGD = TRUE)
  }
  depth <- Xsuperlet$x1
  y_axis <- Xsuperlet$Period
  depth <- as.numeric(depth)
  y_axis <- as.numeric(y_axis)
  if (plot_dir != TRUE) {
    xlim_vals = rev(c(min(Xsuperlet$x1), max(Xsuperlet$x1)))
  }
  else {
    xlim_vals = c(min(Xsuperlet$x1), max(Xsuperlet$x1))
  }
  if (is.null(lowerPeriod) == TRUE) {
    lowerPeriod <- min(Xsuperlet$Period)
  }
  if (is.null(upperPeriod) == TRUE) {
    upperPeriod <- max(Xsuperlet$Period)
  }
  ylim_vals = (c(lowerPeriod, upperPeriod))
  if (add_data == TRUE & add_avg == FALSE & plot_horizontal ==
      FALSE) {
    layout.matrix <- matrix(c(1, 3, 0, 2), nrow = 2, ncol = 2,
                            byrow = TRUE)
    graphics::layout(mat = layout.matrix, heights = c(1,
                                                      0.25), widths = c(1, 4))
    par(mar = c(4, 4, 2, 0))
    plot(
      x = Xsuperlet$y1,
      y = Xsuperlet$x1,
      type = "l",
      col = "black",
      yaxs = "i",
      ylim = xlim_vals,
      xlab = "proxy_1",
      ylab = x_lab,
      xaxt = "s",
      yaxt = "s"
    )

    par(new = TRUE)
    plot(
      x = Xsuperlet$y2,
      y = Xsuperlet$x2,
      type = "l",
      col = "red",
      yaxs = "i",
      ylim = xlim_vals,
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n"
    )

    axis(
      side = 3,
      col = "red",
      col.axis = "red"
    )

    mtext(
      "proxy_2",
      side = 3,
      line = 2,
      col = "red"
    )

    if (is.null(add_abline_h) != TRUE) {
      abline(h = add_abline_h)
    }
    par(new = FALSE, mar = c(3, 0, 2, 2))
    image(x = seq(from = 0, to = n.levels), y = 1, z = t(matrix(power_max_mat.levels,
                                                                nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
          useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "",
          ylab = "")
    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)
    mtext(key.labels, side = 1, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 2, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)
    par(new = FALSE, mar = c(4, 0, 2, 2))
    image(y = Xsuperlet$x1, x = Xsuperlet$axis.2, z = t(Xsuperlet$Power),
          col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
          xlab = periodlab, ylab = "", yaxt = "n", xaxt = "n",
          main = main, ylim = xlim_vals, xlim = log2(sort(ylim_vals)))
    box(lwd = lwd.axis)
    period.tick = unique(trunc(Xsuperlet$axis.2))
    period.tick <- period.tick[period.tick >= min(log2(ylim_vals))]
    period.tick <- period.tick[period.tick <= max(log2(ylim_vals))]
    period.tick.label = 2^(period.tick)
    axis(1, lwd = lwd.axis, at = as.vector(period.tick),
         labels = NA, tck = periodtck, tcl = periodtcl)
    mtext(period.tick.label, side = 1, at = period.tick,
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)
    if (is.null(add_lines) != TRUE) {
      for (i in 2:ncol(add_lines)) lines(y = add_lines[,
                                                       1], x = log2(add_lines[, i]))
    }
    if (is.null(add_points) != TRUE) {
      for (i in 2:ncol(add_points)) points(y = add_points[,
                                                          1], x = log2(add_points[, i]))
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }
    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }


    if (plot_arrows == TRUE) {
      maxdetect <- matrix(0, nrow = Xsuperlet$nr, ncol = Xsuperlet$nc)
      for (j in 1:Xsuperlet$nc) {
        for (i in 2:(Xsuperlet$nr - 1)) {
          # Vertical check along the Period axis
          if ((Power[i, j] > Power[i + 1, j]) &
              (Power[i, j] > Power[i - 1, j])) {
            maxdetect[i, j] <- 1
          }
        }
      }


      # --- 4. PHASE ARROWS ---
      # Filter by top 20% power and use the ridge detector
      p_thresh = quantile(Xsuperlet$Power, pwr_quant)
      # Subsample time (t_idx) to avoid overlapping arrows in high-res Superlets
      t_idx = seq(1, Xsuperlet$nc, by = round(Xsuperlet$nc / n_arrows))
      s_idx = 1:Xsuperlet$nr

      for (i in t_idx) {
        for (j in s_idx) {
          if (Power[j, i] >= p_thresh && maxdetect[j, i] == 1) {
            x0 <- Xsuperlet$x1[i]
            y0 <- Xsuperlet$axis.2[j]
            angle <- Phase[j, i]

            # Draw the phase arrow
            # Angle 0: In-phase (Right), Pi: Anti-phase (Left)
            x1 <- x0 + cos(angle) * T_unit  # Time-unit length
            y1 <- y0 + sin(angle) * P_unit # Log2-unit length
            if (plot_horizontal == TRUE) {
              arrows(
                x0,
                y0,
                x1,
                y1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )
            } else{
              arrows(
                y0,
                x0,
                y1,
                x1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )


            }
          }
        }
      }
    }
  }
  if (add_data == FALSE & add_avg == FALSE & plot_horizontal ==
      FALSE) {
    layout.matrix <- matrix(c(1, 2), nrow = 2, ncol = 1,
                            byrow = TRUE)
    graphics::layout(mat = layout.matrix, heights = c(1,
                                                      4), widths = c(1))
    par(mar = c(2, 4, 2, 3))
    image(x = seq(from = 0, to = n.levels), y = 1, z = t(matrix(power_max_mat.levels,
                                                                nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
          useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "",
          ylab = "", )
    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)
    mtext(key.labels, side = 1, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 2, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)
    par(new = FALSE, mar = c(4, 4, 2, 3))
    image(y = Xsuperlet$x1, x = Xsuperlet$axis.2, z = t(Xsuperlet$Power),
          col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
          xlab = periodlab, ylab = x_lab, xaxt = "n", main = main,
          ylim = xlim_vals, xlim = log2(sort(ylim_vals)))
    box(lwd = lwd.axis)
    period.tick = unique(trunc(Xsuperlet$axis.2))
    period.tick <- period.tick[period.tick >= min(log2(ylim_vals))]
    period.tick <- period.tick[period.tick <= max(log2(ylim_vals))]
    period.tick.label = 2^(period.tick)
    axis(1, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    mtext(period.tick.label, side = 1, at = period.tick,
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)
    if (is.null(add_lines) != TRUE) {
      for (i in 2:ncol(add_lines)) lines(y = add_lines[,
                                                       1], x = log2(add_lines[, i]))
    }
    if (is.null(add_points) != TRUE) {
      for (i in 2:ncol(add_points)) points(y = add_points[,
                                                          1], x = log2(add_points[, i]))
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }
    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }
    if (plot_arrows == TRUE) {
      maxdetect <- matrix(0, nrow = Xsuperlet$nr, ncol = Xsuperlet$nc)
      for (j in 1:Xsuperlet$nc) {
        for (i in 2:(Xsuperlet$nr - 1)) {
          # Vertical check along the Period axis
          if ((Power[i, j] > Power[i + 1, j]) &
              (Power[i, j] > Power[i - 1, j])) {
            maxdetect[i, j] <- 1
          }
        }
      }


      # --- 4. PHASE ARROWS ---
      # Filter by top 20% power and use the ridge detector
      p_thresh = quantile(Xsuperlet$Power, pwr_quant)
      # Subsample time (t_idx) to avoid overlapping arrows in high-res Superlets
      t_idx = seq(1, Xsuperlet$nc, by = round(Xsuperlet$nc / n_arrows))
      s_idx = 1:Xsuperlet$nr

      for (i in t_idx) {
        for (j in s_idx) {
          if (Power[j, i] >= p_thresh && maxdetect[j, i] == 1) {
            x0 <- Xsuperlet$x1[i]
            y0 <- Xsuperlet$axis.2[j]
            angle <- Phase[j, i]

            # Draw the phase arrow
            # Angle 0: In-phase (Right), Pi: Anti-phase (Left)
            x1 <- x0 + cos(angle) * T_unit  # Time-unit length
            y1 <- y0 + sin(angle) * P_unit # Log2-unit length
            if (plot_horizontal == TRUE) {
              arrows(
                x0,
                y0,
                x1,
                y1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )
            } else{
              arrows(
                y0,
                x0,
                y1,
                x1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )


            }
          }
        }
      }
    }
  }
  if (add_data == FALSE & add_avg == TRUE & plot_horizontal ==
      FALSE) {
    layout.matrix <- matrix(c(1, 2, 0, 3), nrow = 2, ncol = 2,
                            byrow = TRUE)
    graphics::layout(mat = layout.matrix, heights = c(1,
                                                      4), widths = c(1, 4))
    par(mar = c(2, 2, 2, 3))
    image(x = seq(from = 0, to = n.levels), y = 1, z = t(matrix(power_max_mat.levels,
                                                                nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
          useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "",
          ylab = "", )
    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)
    mtext(key.labels, side = 1, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 1, at = median(key.marks),
          line = 3, las = 2, font = par()$font.axis, cex = par()$cex.axis,
          xpd = NA)
    box(lwd = lwd.axis)
    par(mar = c(0, 4, 2, 2))
    plot(y = Xsuperlet$Power.avg, x = Xsuperlet$Period, log = "x",
         type = "l", xaxs = "i", xaxt = "n", ylab = "Wt. power",
         xaxs = "i", xlim = sort(ylim_vals))
    if (is.null(add_abline_v) != TRUE) {
      abline(v = (add_abline_v))
    }
    par(new = FALSE, mar = c(4, 4, 0, 2))
    image(y = Xsuperlet$x1, x = Xsuperlet$axis.2, z = t(Xsuperlet$Power),
          col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
          xlab = periodlab, ylab = x_lab, xaxt = "n", main = main,
          ylim = xlim_vals, xlim = log2(sort(ylim_vals)))
    box(lwd = lwd.axis)
    period.tick = unique(trunc(Xsuperlet$axis.2))
    period.tick <- period.tick[period.tick >= min(log2(ylim_vals))]
    period.tick <- period.tick[period.tick <= max(log2(ylim_vals))]
    period.tick.label = 2^(period.tick)
    axis(1, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    mtext(period.tick.label, side = 1, at = period.tick,
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)
    if (is.null(add_lines) != TRUE) {
      for (i in 2:ncol(add_lines)) lines(y = add_lines[,
                                                       1], x = log2(add_lines[, i]))
    }
    if (is.null(add_points) != TRUE) {
      for (i in 2:ncol(add_points)) points(y = add_points[,
                                                          1], x = log2(add_points[, i]))
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }
    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }
    if (plot_arrows == TRUE) {
      maxdetect <- matrix(0, nrow = Xsuperlet$nr, ncol = Xsuperlet$nc)
      for (j in 1:Xsuperlet$nc) {
        for (i in 2:(Xsuperlet$nr - 1)) {
          # Vertical check along the Period axis
          if ((Power[i, j] > Power[i + 1, j]) &
              (Power[i, j] > Power[i - 1, j])) {
            maxdetect[i, j] <- 1
          }
        }
      }


      # --- 4. PHASE ARROWS ---
      # Filter by top 20% power and use the ridge detector
      p_thresh = quantile(Xsuperlet$Power, pwr_quant)
      # Subsample time (t_idx) to avoid overlapping arrows in high-res Superlets
      t_idx = seq(1, Xsuperlet$nc, by = round(Xsuperlet$nc / n_arrows))
      s_idx = 1:Xsuperlet$nr

      for (i in t_idx) {
        for (j in s_idx) {
          if (Power[j, i] >= p_thresh && maxdetect[j, i] == 1) {
            x0 <- Xsuperlet$x1[i]
            y0 <- Xsuperlet$axis.2[j]
            angle <- Phase[j, i]

            # Draw the phase arrow
            # Angle 0: In-phase (Right), Pi: Anti-phase (Left)
            x1 <- x0 + cos(angle) * T_unit  # Time-unit length
            y1 <- y0 + sin(angle) * P_unit # Log2-unit length
            if (plot_horizontal == TRUE) {
              arrows(
                x0,
                y0,
                x1,
                y1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )
            } else{
              arrows(
                y0,
                x0,
                y1,
                x1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )


            }
          }
        }
      }
    }
  }
  if (add_data == TRUE & add_avg == TRUE & plot_horizontal ==
      FALSE) {
    layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2,
                            byrow = TRUE)
    graphics::layout(mat = layout.matrix, heights = c(1,
                                                      4), widths = c(1, 4))
    par(mar = c(2, 0.5, 2, 6), xpd = NA)
    image(y = seq(from = 0, to = n.levels), x = 1, z = (matrix(power_max_mat.levels,
                                                               nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
          useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "",
          ylab = "", )
    axis(2, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)
    mtext(key.labels, side = 2, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 1, at = 1, line = 0.5, font = par()$font.axis,
          cex = par()$cex.axis, las = 1, xpd = NA)
    box(lwd = lwd.axis)
    par(mar = c(0, 0, 2, 2))
    plot(y = Xsuperlet$Power.avg, x = Xsuperlet$Period, log = "x",
         type = "l", xaxs = "i", xaxt = "n", ylab = "Wt. power",
         xlab = "", xaxs = "i", xlim = sort(ylim_vals))
    if (is.null(add_abline_v) != TRUE) {
      abline(v = (add_abline_v), xpd = FALSE)
    }
    title(ylab = "Wt. power", xpd = NA)
    par(mar = c(4, 4, 0, 0), xpd = TRUE)
    plot(
      x = Xsuperlet$y1,
      y = Xsuperlet$x1,
      type = "l",
      col = "black",
      yaxs = "i",
      ylim = xlim_vals,
      xlab = "proxy_1",
      ylab = x_lab,
      xaxt = "s",
      yaxt = "s"
    )

    par(new = TRUE)
    plot(
      x = Xsuperlet$y2,
      y = Xsuperlet$x2,
      type = "l",
      col = "red",
      yaxs = "i",
      ylim = xlim_vals,
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n"
    )

    axis(
      side = 3,
      col = "red",
      col.axis = "red"
    )

    mtext(
      "proxy_2",
      side = 3,
      line = 2,
      col = "red"
    )


    par(new = FALSE, mar = c(4, 0, 0, 2), xpd = FALSE)
    image(y = Xsuperlet$x1, x = Xsuperlet$axis.2, z = t(Xsuperlet$Power),
          col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
          xlab = periodlab, ylab = "", xaxt = "n", yaxt = "n",
          main = main, ylim = xlim_vals, xlim = log2(sort(ylim_vals)))
    box(lwd = lwd.axis)
    period.tick = unique(trunc(Xsuperlet$axis.2))
    period.tick <- period.tick[period.tick >= min(log2(ylim_vals))]
    period.tick <- period.tick[period.tick <= max(log2(ylim_vals))]
    period.tick.label = 2^(period.tick)
    axis(1, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    mtext(period.tick.label, side = 1, at = period.tick,
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)
    if (is.null(add_lines) != TRUE) {
      for (i in 2:ncol(add_lines)) lines(y = add_lines[,
                                                       1], x = log2(add_lines[, i]))
    }
    if (is.null(add_points) != TRUE) {
      for (i in 2:ncol(add_points)) points(y = add_points[,
                                                          1], x = log2(add_points[, i]))
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = (add_abline_h))
    }
    if (is.null(add_abline_v) != TRUE) {
      abline(v = log2(add_abline_v))
    }
    if (plot_arrows == TRUE) {
      maxdetect <- matrix(0, nrow = Xsuperlet$nr, ncol = Xsuperlet$nc)
      for (j in 1:Xsuperlet$nc) {
        for (i in 2:(Xsuperlet$nr - 1)) {
          # Vertical check along the Period axis
          if ((Power[i, j] > Power[i + 1, j]) &
              (Power[i, j] > Power[i - 1, j])) {
            maxdetect[i, j] <- 1
          }
        }
      }


      # --- 4. PHASE ARROWS ---
      # Filter by top 20% power and use the ridge detector
      p_thresh = quantile(Xsuperlet$Power, pwr_quant)
      # Subsample time (t_idx) to avoid overlapping arrows in high-res Superlets
      t_idx = seq(1, Xsuperlet$nc, by = round(Xsuperlet$nc / n_arrows))
      s_idx = 1:Xsuperlet$nr

      for (i in t_idx) {
        for (j in s_idx) {
          if (Power[j, i] >= p_thresh && maxdetect[j, i] == 1) {
            x0 <- Xsuperlet$x1[i]
            y0 <- Xsuperlet$axis.2[j]
            angle <- Phase[j, i]

            # Draw the phase arrow
            # Angle 0: In-phase (Right), Pi: Anti-phase (Left)
            x1 <- x0 + cos(angle) * T_unit  # Time-unit length
            y1 <- y0 + sin(angle) * P_unit # Log2-unit length
            if (plot_horizontal == TRUE) {
              arrows(
                x0,
                y0,
                x1,
                y1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )
            } else{
              arrows(
                y0,
                x0,
                y1,
                x1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )


            }
          }
        }
      }
    }
  }
  if (add_data == TRUE & add_avg == TRUE & plot_horizontal ==
      TRUE) {
    layout.matrix <- matrix(c(1, 2, 4, 3), nrow = 2, ncol = 2,
                            byrow = TRUE)
    graphics::layout(mat = layout.matrix, heights = c(0.25,
                                                      1), widths = c(8, 2))
    par(mar = c(0, 4, 2, 0))

    plot(
      y = Xsuperlet$y1,
      x = Xsuperlet$x1,
      type = "l",
      col = "black",
      xaxs = "i",
      xlim = xlim_vals,
      ylab = "proxy_1",
      xlab = x_lab,
      yaxt = "s",
      xaxt = "s"
    )

    par(new = TRUE)
    plot(
      y = Xsuperlet$y2,
      x = Xsuperlet$x2,
      type = "l",
      col = "red",
      xaxs = "i",
      xlim = xlim_vals,
      ylab = "",
      xlab = "",
      yaxt = "n",
      xaxt = "n"
    )

    axis(
      side = 4,
      col = "red",
      col.axis = "red"
    )

    mtext(
      "proxy_2",
      side = 4,
      line = 2,
      col = "red"
    )

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }
    par(new = FALSE, mar = c(3, 2, 2, 2), mgp = c(2, 1,
                                                  0))
    image(x = seq(from = 0, to = n.levels), y = 1, z = t(matrix(power_max_mat.levels,
                                                                nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
          useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "Power",
          ylab = "")
    axis(1, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.2)
    mtext(key.labels, side = 1, at = key.marks, line = 0.1,
          las = 2, cex = 0.75)
    box(lwd = lwd.axis)
    par(new = FALSE, mar = c(4, 0, 0, 0.5))
    plot(x = Xsuperlet$Power.avg, y = Xsuperlet$Period, log = "y",
         type = "l", yaxs = "i", xlab = "Wt. power", xaxs = "i",
         ylim = sort(ylim_vals))
    if (is.null(add_abline_h) != TRUE) {
      abline(h = add_abline_h)
    }
    par(new = FALSE, mar = c(4, 4, 0, 0), xpd = FALSE)
    image(x = Xsuperlet$x1, y = Xsuperlet$axis.2, z = (Xsuperlet$Power),
          col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
          ylab = periodlab, xlab = x_lab, axes = TRUE, yaxt = "n",
          main = main, xlim = xlim_vals, ylim = log2(sort(ylim_vals)))
    box(lwd = lwd.axis)
    period.tick = unique(trunc(Xsuperlet$axis.2))
    period.tick <- period.tick[period.tick >= min(log2(ylim_vals))]
    period.tick <- period.tick[period.tick <= max(log2(ylim_vals))]
    period.tick.label = 2^(period.tick)
    axis(2, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    axis(4, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    mtext(period.tick.label, side = 2, at = period.tick,
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)
    if (is.null(add_lines) != TRUE) {
      for (i in 2:ncol(add_lines)) lines(add_lines[, 1],
                                         log2(add_lines[, i]))
    }
    if (is.null(add_points) != TRUE) {
      for (i in 2:ncol(add_points)) points(add_points[,
                                                      1], log2(add_points[, i]))
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }
    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }
    if (plot_arrows == TRUE) {
      maxdetect <- matrix(0, nrow = Xsuperlet$nr, ncol = Xsuperlet$nc)
      for (j in 1:Xsuperlet$nc) {
        for (i in 2:(Xsuperlet$nr - 1)) {
          # Vertical check along the Period axis
          if ((Power[i, j] > Power[i + 1, j]) &
              (Power[i, j] > Power[i - 1, j])) {
            maxdetect[i, j] <- 1
          }
        }
      }


      # --- 4. PHASE ARROWS ---
      # Filter by top 20% power and use the ridge detector
      p_thresh = quantile(Xsuperlet$Power, pwr_quant)
      # Subsample time (t_idx) to avoid overlapping arrows in high-res Superlets
      t_idx = seq(1, Xsuperlet$nc, by = round(Xsuperlet$nc / n_arrows))
      s_idx = 1:Xsuperlet$nr

      for (i in t_idx) {
        for (j in s_idx) {
          if (Power[j, i] >= p_thresh && maxdetect[j, i] == 1) {
            x0 <- Xsuperlet$x1[i]
            y0 <- Xsuperlet$axis.2[j]
            angle <- Phase[j, i]

            # Draw the phase arrow
            # Angle 0: In-phase (Right), Pi: Anti-phase (Left)
            x1 <- x0 + cos(angle) * T_unit  # Time-unit length
            y1 <- y0 + sin(angle) * P_unit # Log2-unit length
            if (plot_horizontal == TRUE) {
              arrows(
                x0,
                y0,
                x1,
                y1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )
            } else{
              arrows(
                y0,
                x0,
                y1,
                x1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )


            }
          }
        }
      }
    }
  }
  if (add_data == TRUE & add_avg == FALSE & plot_horizontal ==
      TRUE) {
    layout.matrix <- matrix(c(1, 0, 3, 2), nrow = 2, ncol = 2,
                            byrow = TRUE)
    graphics::layout(mat = layout.matrix, heights = c(0.25,
                                                      1), widths = c(8, 2))
    par(mar = c(0, 4, 2, 2))

    plot(
      y = Xsuperlet$y1,
      x = Xsuperlet$x1,
      type = "l",
      col = "black",
      xaxs = "i",
      xlim = xlim_vals,
      ylab = "proxy_1",
      xlab = x_lab,
      yaxt = "s",
      xaxt = "s"
    )

    par(new = TRUE)
    plot(
      y = Xsuperlet$y2,
      x = Xsuperlet$x2,
      type = "l",
      col = "red",
      xaxs = "i",
      xlim = xlim_vals,
      ylab = "",
      xlab = "",
      yaxt = "n",
      xaxt = "n"
    )

    axis(
      side = 4,
      col = "red",
      col.axis = "red"
    )

    mtext(
      "proxy_2",
      side = 4,
      line = 2,
      col = "red"
    )

    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }
    par(new = FALSE, mar = c(4, 0, 2, 5))
    image(y = seq(from = 0, to = n.levels), x = 1, z = (matrix(power_max_mat.levels,
                                                               nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
          useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "",
          ylab = "", )
    axis(4, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)
    mtext(key.labels, side = 4, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 4, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)
    par(new = FALSE, mar = c(4, 4, 0, 2))
    image(x = Xsuperlet$x1, y = Xsuperlet$axis.2, z = (Xsuperlet$Power),
          col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
          ylab = periodlab, xlab = x_lab, axes = TRUE, yaxt = "n",
          main = main, xlim = xlim_vals, ylim = log2(sort(ylim_vals)))
    box(lwd = lwd.axis)
    period.tick = unique(trunc(Xsuperlet$axis.2))
    period.tick <- period.tick[period.tick >= min(log2(ylim_vals))]
    period.tick <- period.tick[period.tick <= max(log2(ylim_vals))]
    period.tick.label = 2^(period.tick)
    axis(2, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    axis(4, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    mtext(period.tick.label, side = 2, at = period.tick,
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)
    if (is.null(add_lines) != TRUE) {
      for (i in 2:ncol(add_lines)) lines(add_lines[, 1],
                                         log2(add_lines[, i]))
    }
    if (is.null(add_points) != TRUE) {
      for (i in 2:ncol(add_points)) points(add_points[,
                                                      1], log2(add_points[, i]))
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }
    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }
    if (plot_arrows == TRUE) {
      maxdetect <- matrix(0, nrow = Xsuperlet$nr, ncol = Xsuperlet$nc)
      for (j in 1:Xsuperlet$nc) {
        for (i in 2:(Xsuperlet$nr - 1)) {
          # Vertical check along the Period axis
          if ((Power[i, j] > Power[i + 1, j]) &
              (Power[i, j] > Power[i - 1, j])) {
            maxdetect[i, j] <- 1
          }
        }
      }


      # --- 4. PHASE ARROWS ---
      # Filter by top 20% power and use the ridge detector
      p_thresh = quantile(Xsuperlet$Power, pwr_quant)
      # Subsample time (t_idx) to avoid overlapping arrows in high-res Superlets
      t_idx = seq(1, Xsuperlet$nc, by = round(Xsuperlet$nc / n_arrows))
      s_idx = 1:Xsuperlet$nr

      for (i in t_idx) {
        for (j in s_idx) {
          if (Power[j, i] >= p_thresh && maxdetect[j, i] == 1) {
            x0 <- Xsuperlet$x1[i]
            y0 <- Xsuperlet$axis.2[j]
            angle <- Phase[j, i]

            # Draw the phase arrow
            # Angle 0: In-phase (Right), Pi: Anti-phase (Left)
            x1 <- x0 + cos(angle) * T_unit  # Time-unit length
            y1 <- y0 + sin(angle) * P_unit # Log2-unit length
            if (plot_horizontal == TRUE) {
              arrows(
                x0,
                y0,
                x1,
                y1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )
            } else{
              arrows(
                y0,
                x0,
                y1,
                x1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )


            }
          }
        }
      }
    }
  }
  if (add_data == FALSE & add_avg == TRUE & plot_horizontal ==
      TRUE) {
    layout.matrix <- matrix(c(3, 2, 1), nrow = 1, ncol = 3,
                            byrow = TRUE)
    graphics::layout(mat = layout.matrix, heights = c(1),
                     widths = c(6, 1, 1))
    par(mar = c(4, 0, 2, 5))
    image(y = seq(from = 0, to = n.levels), x = 1, z = (matrix(power_max_mat.levels,
                                                               nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
          useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "",
          ylab = "", )
    axis(4, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)
    mtext(key.labels, side = 4, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 4, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)
    par(new = FALSE, mar = c(4, 0, 2, 0.5))
    plot(x = Xsuperlet$Power.avg, y = Xsuperlet$Period, log = "y",
         type = "l", yaxs = "i", yaxt = "n", xlab = "Wt. power",
         xaxs = "i", ylim = sort(ylim_vals))
    if (is.null(add_abline_h) != TRUE) {
      abline(h = add_abline_h)
    }
    par(new = FALSE, mar = c(4, 4, 2, 0), mgp = c(2, 1,
                                                  0))
    image(x = Xsuperlet$x1, y = Xsuperlet$axis.2, z = (Xsuperlet$Power),
          col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
          ylab = periodlab, xlab = x_lab, axes = TRUE, yaxt = "n",
          main = main, xlim = xlim_vals, ylim = log2(sort(ylim_vals)))
    box(lwd = lwd.axis)
    period.tick = unique(trunc(Xsuperlet$axis.2))
    period.tick <- period.tick[period.tick >= min(log2(ylim_vals))]
    period.tick <- period.tick[period.tick <= max(log2(ylim_vals))]
    period.tick.label = 2^(period.tick)
    axis(2, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    axis(4, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    mtext(period.tick.label, side = 2, at = period.tick,
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)
    if (is.null(add_lines) != TRUE) {
      for (i in 2:ncol(add_lines)) lines(add_lines[, 1],
                                         log2(add_lines[, i]))
    }
    if (is.null(add_points) != TRUE) {
      for (i in 2:ncol(add_points)) points(add_points[,
                                                      1], log2(add_points[, i]))
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }
    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }
    if (plot_arrows == TRUE) {
      maxdetect <- matrix(0, nrow = Xsuperlet$nr, ncol = Xsuperlet$nc)
      for (j in 1:Xsuperlet$nc) {
        for (i in 2:(Xsuperlet$nr - 1)) {
          # Vertical check along the Period axis
          if ((Power[i, j] > Power[i + 1, j]) &
              (Power[i, j] > Power[i - 1, j])) {
            maxdetect[i, j] <- 1
          }
        }
      }


      # --- 4. PHASE ARROWS ---
      # Filter by top 20% power and use the ridge detector
      p_thresh = quantile(Xsuperlet$Power, pwr_quant)
      # Subsample time (t_idx) to avoid overlapping arrows in high-res Superlets
      t_idx = seq(1, Xsuperlet$nc, by = round(Xsuperlet$nc / n_arrows))
      s_idx = 1:Xsuperlet$nr

      for (i in t_idx) {
        for (j in s_idx) {
          if (Power[j, i] >= p_thresh && maxdetect[j, i] == 1) {
            x0 <- Xsuperlet$x1[i]
            y0 <- Xsuperlet$axis.2[j]
            angle <- Phase[j, i]

            # Draw the phase arrow
            # Angle 0: In-phase (Right), Pi: Anti-phase (Left)
            x1 <- x0 + cos(angle) * T_unit  # Time-unit length
            y1 <- y0 + sin(angle) * P_unit # Log2-unit length
            if (plot_horizontal == TRUE) {
              arrows(
                x0,
                y0,
                x1,
                y1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )
            } else{
              arrows(
                y0,
                x0,
                y1,
                x1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )


            }
          }
        }
      }
    }
  }
  if (add_data == FALSE & add_avg == FALSE & plot_horizontal ==
      TRUE) {
    layout.matrix <- matrix(c(2, 1), nrow = 1, ncol = 2,
                            byrow = TRUE)
    graphics::layout(mat = layout.matrix, heights = c(1),
                     widths = c(10, 2.25))
    par(mar = c(4, 0, 2, 5))
    image(y = seq(from = 0, to = n.levels), x = 1, z = (matrix(power_max_mat.levels,
                                                               nrow = 1)), col = key.cols, breaks = power_max_mat.levels,
          useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "",
          ylab = "", )
    axis(4, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = 1.24)
    mtext(key.labels, side = 4, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    mtext(c("Power"), side = 4, at = mean(key.marks), line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    box(lwd = lwd.axis)
    par(new = FALSE, mar = c(4, 4, 2, 0.5))
    image(x = Xsuperlet$x1, y = Xsuperlet$axis.2, z = (Xsuperlet$Power),
          col = key.cols, breaks = power_max_mat.levels, useRaster = TRUE,
          ylab = periodlab, xlab = x_lab, axes = TRUE, yaxt = "n",
          main = main, xlim = xlim_vals, ylim = log2(sort(ylim_vals)))
    box(lwd = lwd.axis)
    period.tick = unique(trunc(Xsuperlet$axis.2))
    period.tick <- period.tick[period.tick >= min(log2(ylim_vals))]
    period.tick <- period.tick[period.tick <= max(log2(ylim_vals))]
    period.tick.label = 2^(period.tick)
    axis(2, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    axis(4, lwd = lwd.axis, at = period.tick, labels = NA,
         tck = periodtck, tcl = periodtcl)
    mtext(period.tick.label, side = 2, at = period.tick,
          las = 2, line = par()$mgp[2] - 0.5, font = par()$font.axis,
          cex = par()$cex.axis)
    if (is.null(add_lines) != TRUE) {
      for (i in 2:ncol(add_lines)) lines(add_lines[, 1],
                                         log2(add_lines[, i]))
    }
    if (is.null(add_points) != TRUE) {
      for (i in 2:ncol(add_points)) points(add_points[,
                                                      1], log2(add_points[, i]))
    }
    if (is.null(add_abline_h) != TRUE) {
      abline(h = log2(add_abline_h))
    }
    if (is.null(add_abline_v) != TRUE) {
      abline(v = add_abline_v)
    }
    if (plot_arrows == TRUE) {
      maxdetect <- matrix(0, nrow = Xsuperlet$nr, ncol = Xsuperlet$nc)
      for (j in 1:Xsuperlet$nc) {
        for (i in 2:(Xsuperlet$nr - 1)) {
          # Vertical check along the Period axis
          if ((Power[i, j] > Power[i + 1, j]) &
              (Power[i, j] > Power[i - 1, j])) {
            maxdetect[i, j] <- 1
          }
        }
      }


      # --- 4. PHASE ARROWS ---
      # Filter by top 20% power and use the ridge detector
      p_thresh = quantile(Xsuperlet$Power, pwr_quant)
      # Subsample time (t_idx) to avoid overlapping arrows in high-res Superlets
      t_idx = seq(1, Xsuperlet$nc, by = round(Xsuperlet$nc / n_arrows))
      s_idx = 1:Xsuperlet$nr

      for (i in t_idx) {
        for (j in s_idx) {
          if (Power[j, i] >= p_thresh && maxdetect[j, i] == 1) {
            x0 <- Xsuperlet$x1[i]
            y0 <- Xsuperlet$axis.2[j]
            angle <- Phase[j, i]

            # Draw the phase arrow
            # Angle 0: In-phase (Right), Pi: Anti-phase (Left)
            x1 <- x0 + cos(angle) * T_unit  # Time-unit length
            y1 <- y0 + sin(angle) * P_unit # Log2-unit length
            if (plot_horizontal == TRUE) {
              arrows(
                x0,
                y0,
                x1,
                y1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )
            } else{
              arrows(
                y0,
                x0,
                y1,
                x1,
                length = arrow_length,
                lwd = arrow_lwd ,
                angle = arrow_angle,
                col = "black"
              )


            }
          }
        }
      }
    }
  }
}
