#' @title Windowed fft based spectral analysis
#'
#' @description The \code{\link{win_fft}} function conducts a windowed spectral analysis based using the fft
#'
#' @param data Record to analyse (matrix or data frame with two columns: depth/time in column 1, proxy value in
#' column 2), evenly sampled (the sampling interval is taken from the first two rows; rows with missing values are
#' removed). Must be given (the \code{NULL} default fails). \code{Default=NULL}.
#' @param padfac Zero-padding factor of each window (numeric, >= 1): every window is padded with zeros to
#' \code{padfac} times its length before the FFT, which gives a finer (interpolated) frequency grid and a smoother
#' spectrum; 1 or less means no padding. \code{Default=5}.
#' @param window_size Length of the running window (numeric, > 0, in the depth/time unit of \code{data}, several
#' sampling intervals long). A spectrum is computed for a window centred on every sample; at both ends the record is
#' extended by half a window of AR(1) noise with the mean, standard deviation and lag-1 autocorrelation of the first
#' and last half window. Must be given. \code{Default=NULL}.
#' @param run_multicore Compute the windows in parallel on several cores (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param n_cores Number of cores used when \code{run_multicore = TRUE} (integer; \code{NULL}: an automatic, safe
#'   number: physical cores minus one, at most 4, and limited by the free memory). Set
#'   \code{options(waverider.n_cores = ...)} to change the default. A value > 1 also switches on parallel computing.
#'   Results do not depend on the number of cores. \code{Default=NULL}.
#' @param genplot Plot the matrix chosen with \code{plot_res} against depth and frequency, with the mean spectrum and
#' a colour key (logical, TRUE or FALSE); needs \code{freq_min} and \code{freq_max}. \code{Default=FALSE}.
#' @param x_lab Label of the x-axis (depth/time) of the plot (character). \code{Default="depth (m)"}.
#' @param y_lab Label of the y-axis (frequency) of the plot (character). \code{Default="frequency cycle/metre"}.
#' @param plot_res Which matrix of the result is plotted (integer, 1 to 9): 1 amplitude, 2 power, 3 phase, 4 frequency
#' uncertainty, 5 AR(1) confidence level (percent), 6 AR(1) fit, 7 AR(1) 90\% power level, 8 AR(1) 95\% power level,
#' 9 AR(1) 99\% power level. \code{Default=1}.
#' @param perc_vis Lower cut-off of the colour scale as a quantile of the plotted values (numeric, between 0 and 1):
#' values below this quantile get the lowest colour, e.g. 0.5 shows only the upper half of the values in colour.
#' \code{Default=0}.
#' @param freq_max Highest frequency shown in the plot (numeric, in cycles per depth/time unit, at most the Nyquist
#' frequency); required when \code{genplot = TRUE}. \code{Default=NULL}.
#' @param freq_min Lowest frequency shown in the plot (numeric, in cycles per depth/time unit, > 0); required when
#' \code{genplot = TRUE}. \code{Default=NULL}.
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
#' @param keep_editable Keep the graphical parameters and layout after plotting, so that extra features can be added
#' (logical, TRUE or FALSE); \code{FALSE} restores the previous settings. \code{Default=FALSE}.
#' @param verbose Show a progress bar of the window computations (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param dev_new Open a new plot window (via \code{wr_open_device()}) before plotting (logical, TRUE or FALSE); a
#' separate window gives a "nice" looking plot, but in an R markdown document the plot may then not appear.
#' \code{Default=FALSE}.
#'
#' @author
#'Based on the \link[astrochron]{periodogram}
#'function of the 'astrochron' R package.
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#'
#'@examples
#'\donttest{
#'#Conduct a windowed ftt on the magnetic susceptibility record
#'#of the Sullivan core of Pas et al., (2018).
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
#'                    palette_name ="rainbow",
#'                    color_brewer= "grDevices",
#'                    keep_editable=FALSE,
#'                    verbose=FALSE,
#'                    dev_new=FALSE)
#'}
#'
#' @return
#'Returns a list which contains 10 elements
#'element 1: Amplitude matrix
#'element 2: Power matrix
#'element 3: Phase matrix
#'element 4: AR1_CL matrix
#'element 5: AR1_Fit matrix
#'element 6: AR1_90_power matrix
#'element 7: AR1_95_power matrix
#'element 8: AR1_99_power matrix
#'element 9: depth
#'element 10: y_axis
#'If genplot is \code{Default=TRUE} then a plot of one of the elements 1:8 is plotted
#'
#' @param colour_brewer British spelling of \code{color_brewer}, with the same options; if given (not \code{NULL}) it overrides
#' \code{color_brewer}. \code{Default=NULL}.
#' @export
#' @importFrom stats quantile
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats lm
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics text
#' @importFrom graphics box
#' @importFrom graphics polygon
#' @importFrom grDevices rgb
#' @importFrom graphics layout
#' @importFrom stats cor
#' @importFrom stats pchisq
#' @importFrom stats qchisq
#' @importFrom astrochron periodogram
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors

win_fft <- function(data = NULL,
    padfac = 5,
    window_size = NULL,
    run_multicore = FALSE,
    n_cores = NULL,
    genplot = FALSE,
    x_lab = c("depth (m)"),
    y_lab = c("frequency cycle/metre"),
    plot_res = 1,
    perc_vis = 0,
    freq_max = NULL,
    freq_min = NULL,
    palette_name = "viridis",
    color_brewer = "viridis",
    keep_editable = FALSE,
    verbose=FALSE,
    dev_new=FALSE,
    colour_brewer = NULL) {
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym

  n.levels = 100
  dat <- data
  dat <- na.omit(dat)
  d <- data.frame(dat)
  dt <- d[2, 1] - d[1, 1]
  Nyq <- 1 / (2 * dt)
  xmax <- Nyq

  timesteps_data <- round(((window_size / dt) / 2), 0)

  phi_x_up <-  wr_autocorrelation(dat[1:timesteps_data, 2])
  if (phi_x_up >= 1) {
    phi_x_up <- 0.99
  }
  x_up <- wr_coloured_noise(
      timesteps = timesteps_data,
      mean = mean(dat[1:timesteps_data, 2]),
      sd = sd(dat[1:timesteps_data, 2]),
      phi = phi_x_up
  )

  phi_x_down <-
      wr_autocorrelation(dat[(nrow(dat) - timesteps_data):nrow(dat), 2])

  if (phi_x_down >= 1) {
    phi_x_down <- 0.99
  }

  x_down <- wr_coloured_noise(
      timesteps = timesteps_data,
      mean = mean(dat[(nrow(dat) - timesteps_data):nrow(dat), 2]),
      sd = sd(dat[(nrow(dat) - timesteps_data):nrow(dat), 2]),
      phi = phi_x_down
  )

  depth_up <-
      sort(seq(
          from = d[1, 1] - dt,
          by = -dt,
          length.out = timesteps_data
      ))
  depth_down <-
      seq(from = (d[nrow(d), 1] + dt),
          by = dt,
          length.out = timesteps_data)

  colnames(d) <- c("a", "b")
  pad_up <- cbind(depth_up, x_up)
  colnames(pad_up) <- c("a", "b")
  pad_down <- cbind(depth_down, x_down)
  colnames(pad_down) <- c("a", "b")

  d <- rbind(pad_up, d, pad_down)




  simulations <- nrow(dat)




  i <- 1 # needed to assign 1 to ijk to avoid note
  npts <- round(((window_size / dt)), 0)

  fit <-  wr_lapply(seq_len(simulations), function(i) {
    d_subsel <- d[i:(i + (window_size / dt - 1)),]

    demean = TRUE
    if (demean == TRUE) {
      dave <- colMeans(d_subsel[2])
      d_subsel[2] <- d_subsel[2] - dave
    }

    Nyq <- 1 / (2 * dt)
    Ray <- 1 / (npts * dt)
    pad <- as.numeric(d_subsel[, 2])

    if (padfac > 1) {
      pad <- append(pad, rep(0, (npts * padfac - npts)))
    }
    if ((npts * padfac) %% 2 != 0) {
      pad <- append(pad, 0)
    }


    nf = length(pad)
    df <- 1 / (nf * dt)
    freq <- double(nf)
    ijk <- seq(1, nf, by = 1)
    ft <- fft(pad)
    nrm = 1
    if (nrm == 1) {
      ft = ft / npts
    }
    pwr <- Mod(ft) ^ 2
    amp <- sqrt(pwr)
    phase <- atan2(Im(ft), Re(ft))
    freq <- df * (ijk - 1)

    fft.out <- data.frame(cbind(freq, amp, pwr, phase))


    fft.out[,5] <-   uncertainty_freq <- freq/(2*pi*(window_size/(1/freq)))

    fft.out <- fft.out[fft.out[, 1] <= Nyq ,]
    fft.out <- fft.out[fft.out[, 1] > 0,]


    colnames(fft.out) <- c("Frequency", "Amplitude", "Power",
        "Phase","uncertainty")

    lag0 <- d[1:(npts - 1), 2]
    lag1 <- d[2:npts, 2]
    rho <- cor(lag0, lag1)
    So = mean(fft.out[, 3])
    AR = So * (1 - (rho ^ 2)) / (1 - (2 * rho * cos(pi * fft.out[,
        1] / Nyq)) + (rho ^
        2))
    dofAR = 2
    chiAR <- (fft.out[, 3] / AR) * dofAR
    chiCLAR <- pchisq(chiAR, df = dofAR)
    AR1_90 <- AR * qchisq(0.9, df = dofAR) / dofAR
    AR1_95 <- AR * qchisq(0.95, df = dofAR) / dofAR
    AR1_99 <- AR * qchisq(0.99, df = dofAR) / dofAR
    fft.out <- data.frame(
        cbind(
            fft.out[, 1],
            fft.out[,
                2],
            fft.out[, 3],
            fft.out[, 4],
            fft.out[, 5],
            chiCLAR * 100,
            AR,
            AR1_90,
            AR1_95,
            AR1_99
        )
    )
    colnames(fft.out) <- c(
        "Frequency",
        "Amplitude",
        "Power",
        "Phase",
        "uncertainty",
        "AR1_CL",
        "AR1_Fit",
        "AR1_90_power",
        "AR1_95_power",
        "AR1_99_power"
    )

    return(fft.out)
  }, n_cores = n_cores, run_multicore = run_multicore, combine = "list", progress = verbose)



  fit2 <- fit


  Amplitude_mat <-
      matrix(data = NA,
          ncol = nrow(dat),
          nrow = nrow(fit2[[1]]))
  Power_mat <-
      matrix(data = NA,
          ncol = nrow(dat),
          nrow = nrow(fit2[[1]]))

  Phase_mat <-
      matrix(data = NA,
          ncol = nrow(dat),
          nrow = nrow(fit2[[1]]))
  uncer_mat <-
      matrix(data = NA,
          ncol = nrow(dat),
          nrow = nrow(fit2[[1]]))


  AR1_CL_mat <-
      matrix(data = NA,
          ncol = nrow(dat),
          nrow = nrow(fit2[[1]]))

  AR1_Fit_mat <-
      matrix(data = NA,
          ncol = nrow(dat),
          nrow = nrow(fit2[[1]]))


  AR1_90_power_mat <-
      matrix(data = NA,
          ncol = nrow(dat),
          nrow = nrow(fit2[[1]]))

  AR1_95_power_mat <-
      matrix(data = NA,
          ncol = nrow(dat),
          nrow = nrow(fit2[[1]]))

  AR1_99_power_mat <-
      matrix(data = NA,
          ncol = nrow(dat),
          nrow = nrow(fit2[[1]]))

  for (kk in 1:length(fit2)) {
    extract <- as.data.frame(fit2[[kk]])
    Amplitude_mat[, kk] <- extract[, 2]
    Power_mat[, kk] <- extract[, 3]
    Phase_mat[, kk] <- extract[, 4]
    uncer_mat[,kk] <- extract[, 5]
    AR1_CL_mat[, kk] <- extract[, 6]
    AR1_Fit_mat[, kk] <- extract[, 7]
    AR1_90_power_mat[, kk] <- extract[, 8]
    AR1_95_power_mat[, kk] <- extract[, 9]
    AR1_99_power_mat[, kk] <- extract[, 10]
  }


  depth <- (dat[, 1])
  y_axis <- unlist(fit2[[kk]][1])

  results <- list(
      Amplitude_mat = Amplitude_mat,
      Power_mat = Power_mat,
      Phase_mat = Phase_mat,
      uncer_mat = uncer_mat,
      AR1_CL_mat = AR1_CL_mat,
      AR1_Fit_mat = AR1_Fit_mat,
      AR1_90_power_mat = AR1_90_power_mat,
      AR1_95_power_mat = AR1_95_power_mat,
      AR1_99_power_mat = AR1_99_power_mat,
      depth = dat[, 1],
      y_axis = unlist(fit2[[kk]][1])
  )



  if (genplot == TRUE) {
    y_axis <- as.numeric(unlist(results$y_axis))
    sel_cols_up <- max(which(y_axis < freq_max))
    sel_cols_down <- min(which(y_axis > freq_min))

    bottom_perc <- perc_vis

    pmax_avg_sel <- t(results[[plot_res]])
    pmax_avg_sel <- pmax_avg_sel[, sel_cols_down:sel_cols_up]

    if(dev_new==TRUE){
      wr_open_device(width = 14, height = 7)
    }

    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }

    layout.matrix <- matrix(c(3, 1, 2), nrow = 1, ncol = 3)
    layout(
        mat = layout.matrix,
        heights = c(1, 1, 1),
        # Heights of the two rows
        widths = c(8, 2, 2)
    )

    par(mar = c(4, 4, 2, 3))
    power_max_mat.levels = quantile(pmax_avg_sel,
        probs = seq(
            from = bottom_perc,
            to = 1,
            length.out = n.levels + 1
        ))


    image.plt = par()$plt

    key.cols <- .wr_key_cols(palette_name, color_brewer, n.levels)

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
  return(results)
}
