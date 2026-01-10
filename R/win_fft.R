#' @title Windowed fft based spectral analysis
#'
#' @description The \code{\link{win_fft}} function conducts a windowed spectral analysis based using the fft
#'
#'@param data Input data set  should consist of a matrix with 2 columns with
#'first column being depth and the second column being a proxy
#'@param padfac Pad record with zero, zero padding smooths out the spectra
#'@param window_size size of the running window
#'@param run_multicore Run function using multiple cores \code{Default="FALSE"}
#'@param genplot Generate plot \code{Default="FALSE"}
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
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#' @param verbose Print text \code{Default=FALSE}.
#' @param dev_new Opens a new plotting window to plot the plot, this
#' guarantees  a "nice" looking plot however when plotting in an R markdown
#'document the plot might not plot  \code{Default=FALSE}
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
#' @export
#' @importFrom Matrix rowMeans
#' @importFrom stats quantile
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom foreach foreach
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
#' @importFrom foreach %dopar%
#' @importFrom graphics layout
#' @importFrom parallel stopCluster
#' @importFrom colorednoise autocorrelation
#' @importFrom colorednoise colored_noise
#' @importFrom stats cor
#' @importFrom stats pchisq
#' @importFrom stats qchisq
#' @importFrom astrochron periodogram
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

win_fft <- function(data = NULL,
                    padfac = 5,
                    window_size = NULL,
                    run_multicore = FALSE,
                    genplot = FALSE,
                    x_lab = c("depth (m)"),
                    y_lab = c("frequency cycle/metre"),
                    plot_res = 1,
                    perc_vis = 0,
                    freq_max = NULL,
                    freq_min = NULL,
                    palette_name ="rainbow",
                    color_brewer= "grDevices",
                    keep_editable = FALSE,
                    verbose=FALSE,
                    dev_new=FALSE) {

  n.levels = 100
  dat <- data
  dat <- na.omit(dat)
  d <- data.frame(dat)
  dt <- d[2, 1] - d[1, 1]
  Nyq <- 1 / (2 * dt)
  xmax <- Nyq

  timesteps_data <- round(((window_size / dt) / 2), 0)

  phi_x_up <-  autocorrelation(dat[1:timesteps_data, 2])
  if (phi_x_up >= 1) {
    phi_x_up <- 0.99
  }
  x_up <- colorednoise::colored_noise(
    timesteps = timesteps_data,
    mean = mean(dat[1:timesteps_data, 2]),
    sd = sd(dat[1:timesteps_data, 2]),
    phi = phi_x_up
  )

  phi_x_down <-
    autocorrelation(dat[(nrow(dat) - timesteps_data):nrow(dat), 2])

  if (phi_x_down >= 1) {
    phi_x_down <- 0.99
  }

  x_down <- colorednoise::colored_noise(
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


  if (run_multicore == TRUE) {
    numCores <- detectCores()
    cl <- parallel::makeCluster(numCores - 2)
    registerDoSNOW(cl)
  } else{
    numCores <- 1
    cl <- makeCluster(numCores)
    registerDoSNOW(cl)
  }


  simulations <- nrow(dat)

  if (verbose==TRUE){
    pb <- txtProgressBar(max = simulations, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)}else{opts=NULL}



  i <- 1 # needed to assign 1 to ijk to avoid note
  npts <- round(((window_size / dt)), 0)

  fit <-  foreach (i = 1:simulations, .options.snow   = opts) %dopar% {
    d_subsel <- d[i:(i + (window_size / dt - 1)),]

    demean = TRUE
    if (demean == TRUE) {
      dave <- colMeans(d_subsel[2])
      d_subsel[2] <- d_subsel[2] - dave
    }

    Nyq <- 1 / (2 * dt)
    Ray <- 1 / (npts * dt)
    pad <- as.numeric(d_subsel[, 2])

    if (padfac > 1)
      pad <- append(pad, rep(0, (npts * padfac - npts)))
    if ((npts * padfac) %% 2 != 0)
      pad <- append(pad, 0)


    nf = length(pad)
    df <- 1 / (nf * dt)
    freq <- double(nf)
    ijk <- seq(1, nf, by = 1)
    ft <- fft(pad)
    nrm = 1
    if (nrm == 1)
      ft = ft / npts
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
  }


  stopCluster(cl)

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
    dev.new(width = 14, height = 7)}

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
