#'@title Create an age model using minimal tuning
#'
#' @description Create an age model using the minimal tuning technique.
#' This means that the distance between 2 peaks of an extracted cycle are set
#' to duration of the interpreted astronomical cycle
#'
#'@param data Input is an cycle extracted filtered in the depth domain
#' @param pts The pts parameter specifies how many points to the left/right up/down the peak detect algorithm goes in detecting
#'a peak. The peak detecting algorithm works by comparing the values left/right up/down of it, if the values are both higher or lower
#'then the value a peak. To deal with error produced by this algorithm the pts parameter can be changed which can
#'aid in peak detection. Usually increasing the pts parameter means more peak certainty, however it also means that minor peaks might not be
#'picked up by the algorithm \code{Default=5}
#'@param cycle duration in kyr of the filtered/extracted cycle
#'@param tune_opt tuning options "min", "max" and "minmax" use minima, maxima or both
#' of the cyclic signal to create the age model \code{Default="max"}
#'@param output #'The output depends on the output setting
#'If output = 0 output is a matrix of with 4 columns being; depth,proxy,sedimentation rate and time
#'If output = 1 output is a matrix of with 2 columns being; depth and sedimentation rate
#'#'If output = 2 output is a matrix of with 2 columns being; depth and time
#'@param genplot Keep option to add extra features after plotting  \code{Default=FALSE}
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#'@author
#' Part of the code is based on the \link[astrochron]{sedrate2time}
#' function of the 'astrochron' R package
#'
#'@references
#'Routines for astrochronologic testing, astronomical time scale construction, and
#'time series analysis <doi:10.1016/j.earscirev.2018.11.015>
#'
#'@return
#'The output depends on the output setting
#'If output = 0 output is a matrix of with 4 columns being (depth,proxy,sedimentation rate and time)
#'If genplot = TRUE 4 plots are generated; depth vs proxy, depth vs sedimentation rate, depth vs time and time vs proxy
#'If output = 1 output is a matrix of with 2 columns being (depth and sedimentation rate )
#'If genplot = TRUE a plot of depth vs sedimentation rate is generated
#'If output = 2 output is a matrix of with 2 columns being (depth and time)
#'If genplot = TRUE a plot of depth vs time is generated
#'
#'@examples
#'\donttest{
#'# Extract the 405kyr eccentricity cycle from the wavelet scalogram
#'# from the magnetic susceptibility record f the Sullivan core
#'# of Pas et al., (2018) and then create a age model using minimal tuning
#'# (e.g.) set the distance between peaks to 405 kyr
#'
#'mag_wt <- analyze_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = FALSE,
#' omega_nr = 10)
#'
#'
#'mag_405 <- extract_signal_stable_V2(
#'  wavelet = mag_wt,
#'  period_max = 4,
#'  period_min = 2,
#'  add_mean = FALSE,
#'  plot_residual = FALSE,
#'  keep_editable = FALSE
#')
#'
#'mag_405_min_tuning <- minimal_tuning(data = mag_405,
#'pts = 5,
#'cycle = 405,
#'tune_opt = "max",
#'output = 0,
#'genplot = FALSE,
#'keep_editable = FALSE)
#'
#'
#'}
#'
#' @export
#' @importFrom stats na.omit
#' @importFrom graphics par
#' @importFrom graphics plot.new
#' @importFrom grDevices dev.size

minimal_tuning <- function(data = NULL,
                              pts = 5,
                              cycle = 405,
                              tune_opt = "max",
                              output = 0,
                              genplot = FALSE,
                              keep_editable = FALSE){


astro_mindetect <- as.data.frame(data)
astro_mindetect$min <- 0
for (i in pts:(nrow(data) - pts)) {
  if ((data[i, 2] - data[(i + pts), 2] < 0) &
      (data[i, 2] - data[(i - (pts - 1)), 2] < 0))
  {
    astro_mindetect[i, 3] <- 1
  }
}

astro_mindetect_error_corr <- astro_mindetect
astro_mindetect_error_corr <- astro_mindetect_error_corr[astro_mindetect_error_corr$min == 1 ,]

astro_maxdetect <- as.data.frame(data)
astro_maxdetect$max <- 0
for (i in pts:(nrow(data) - pts)) {
  if ((data[i, 2] - data[(i + pts), 2] > 0) &
      (data[i, 2] - data[(i - (pts - 1)), 2]  > 0))
  {
    astro_maxdetect[i, 3] <- 1
  }
}

astro_maxdetect_error_corr <- astro_maxdetect
astro_maxdetect_error_corr <-
  astro_maxdetect_error_corr[astro_maxdetect_error_corr$max == 1 ,]

max <- astro_maxdetect_error_corr
colnames(max) <- c("A", "B", "C")
min <- astro_mindetect_error_corr
colnames(min) <- c("A", "B", "C")

min[, 3] <- -1
peaks <- rbind(max, min)
peaks <- peaks[order(peaks[, 1]), ]
i <- 1
res_rownr <- nrow(peaks)

while (i < res_rownr) {
  if ((i < res_rownr) & (peaks[i, 3] == peaks[(i + 1), 3])) {
    if ((i < res_rownr) & (peaks[i, 3]  == 1 & peaks[(i + 1), 3] == 1) &
        (peaks[i, 2] > peaks[(i + 1), 2])) {
      peaks[(i + 1), ] <- NA
      peaks <- na.omit(peaks)
      res_rownr <- res_rownr - 1
    }
    if ((i < res_rownr) &
        (peaks[i, 3]  == 1 & peaks[(i + 1), 3] == 1) &
        (peaks[i, 2] < peaks[(i + 1), 2])) {
      peaks[i, ] <- NA
      peaks <- na.omit(peaks)
      res_rownr <- res_rownr - 1
    }
    if ((i < res_rownr) &
        (peaks[i, 3] == -1 & peaks[(i + 1), 3] == -1) &
        (peaks[i, 2] < peaks[(i + 1), 2])) {
      peaks[(i + 1), ] <- NA
      peaks <- na.omit(peaks)
      res_rownr <- res_rownr - 1
    }
    if ((i < res_rownr) &
        (peaks[i, 3] == -1 & peaks[(i + 1), 3] == -1) &
        (peaks[i, 2] > peaks[(i + 1), 2])) {
      peaks[i, ] <- NA
      peaks <- na.omit(peaks)
      res_rownr <- res_rownr - 1

    }
  }
  if ((peaks[i, 3] != peaks[(i + 1), 3]) |
      is.na(peaks[i, 3] != peaks[(i + 1), 3])) {
    i <- i + 1
  }

}


if (tune_opt == "min") {
  peaks_min <- peaks[peaks[, 3] > 0,]
  dist <-
    peaks_min[2:(nrow(peaks_min)),] - peaks_min[1:(nrow(peaks_min) - 1),]
  sed_rate <- (dist[, 1] * 100) / cycle
  sed_rate <-
    cbind(sed_rate, peaks_min[1:(nrow(peaks_min) - 1), 1], peaks_min[2:(nrow(peaks_min)), 1])

}

if (tune_opt == "max") {
  peaks_min <- peaks[peaks[, 3] > 0,]
  dist <-
    peaks_min[2:(nrow(peaks_min)),] - peaks_min[1:(nrow(peaks_min) - 1),]
  sed_rate <- (dist[, 1] * 100) / cycle
  sed_rate <-
    cbind(sed_rate, peaks_min[1:(nrow(peaks_min) - 1), 1], peaks_min[2:(nrow(peaks_min)), 1])

}
if (tune_opt == "minmax") {
  peaks_min <- peaks
  dist <-
    peaks_min[2:(nrow(peaks_min)),] - peaks_min[1:(nrow(peaks_min) - 1),]
  sed_rate <- (dist[, 1] * 100) / (cycle / 2)
  sed_rate <-
    cbind(sed_rate, peaks_min[1:(nrow(peaks_min) - 1), 1], peaks_min[2:(nrow(peaks_min)), 1])

}


top <- c(sed_rate[1, 1], data[1, 1], sed_rate[1, 2])
bot <-
  c(sed_rate[nrow(sed_rate), 1], sed_rate[nrow(sed_rate), 3], data[nrow(data), 1])
sed_rate <- rbind(top, sed_rate, bot)

data[, 3] <- NA
p <- 1

for (i  in 1:nrow(data)) {
  if (data[i, 1] < sed_rate[p, 2]) {
    data[i, 3] <- sed_rate[p, 1]
  }
  if (data[i, 1] == sed_rate[p, 2] & p + 1 <= nrow(sed_rate)) {
    data[i, 3] <- (sed_rate[p, 1] + sed_rate[(p + 1), 1]) / 2
  }
  if (p > nrow(sed_rate)) {
    p <- nrow(sed_rate)
  }
  if (p == nrow(sed_rate)) {
    data[i, 3] <- sed_rate[nrow(sed_rate), 1]
  }
  if (data[i, 1] > sed_rate[p, 2]) {
    data[i, 3] <- sed_rate[p, 1]
  }
  if (data[i, 1] == sed_rate[p, 3] & p + 1 <= nrow(sed_rate)) {
    p <- p + 1
  }
}

tracked_cycle_curve <- data[, c(1, 3)]

sedrates <- data.frame(tracked_cycle_curve)
dat <- as.matrix(tracked_cycle_curve)
dat <- na.omit(dat)
dat <- dat[order(dat[, 1], na.last = NA, decreasing = F),]
npts <- length(dat[, 1])
start <- dat[1, 1]
end <- dat[length(dat[, 1]), 1]
x1 <- dat[1:(npts - 1), 1]
x2 <- dat[2:(npts), 1]
dx = x2 - x1
dt = median(dx)
xout <- seq(start, end, by = dt)
npts <- length(xout)
interp <- approx(dat[, 1], dat[, 2], xout, method = "linear",
                 n = npts)
sedrates <- as.data.frame(interp)

npts <- length(sedrates[, 1])
sedrates[1] = sedrates[1] * 100
sedrates[2] = 1 / sedrates[2]
dx = sedrates[2, 1] - sedrates[1, 1]
midptx = (sedrates[2:npts, 1] + sedrates[1:(npts - 1), 1]) / 2
slope = (sedrates[2:npts, 2] - sedrates[1:(npts - 1), 2]) / dx
yint = sedrates[2:npts, 2] - (slope * sedrates[2:npts, 1])
midpty = (slope * midptx) + yint
hsum = cumsum(midpty * dx)
hsum = append(0, hsum)
out = data.frame(cbind(sedrates[, 1] / 100, hsum))

data[, 4] <- out[, 2]
colnames(data) <-  c("depth", "proxy", "cm/kyr", "time")

if (output == 0) {
  data <- data
  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    layout.matrix <- matrix(c(1, 2, 3,4 ), nrow = 4 , ncol = 1)
    graphics::layout(
      mat = layout.matrix,
      heights = c(1),
      # Heights of the two rows
      widths = c(1)
    ) # Widths of the two columns
    par(mar = c(4, 4, 1, 1))
    plot(
      x = data[, 1],
      y = data[, 2],
      type = "l",
      main = "Data depth domain",
      xlab = "meters",
      ylab = "proxy"
    )

    plot(
      x = data[, 1],
      y = data[, 3],
      type = "l",
      xlab = "meters",
      ylab = "cm/kyr (ka)",
      main = "sedimentation rate plot"
    )

    plot(
      data[, 1],
      data[, 4],
      type = "l",
      xlab = "meters",
      ylab = "Time (ka)",
      main = "Data time domain"
    )

    plot(
      data[, 4],
      data[, 2],
      type = "l",
      xlab = "time (ka)",
      ylab = "proxy",
      main = "Data time domain"
    )

  }
}

if (output == 1) {
  data <- data[, c(1, 3)]
  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    layout.matrix <- matrix(c(1), nrow = 1 , ncol = 1)
    graphics::layout(mat = layout.matrix,
                     heights = c(1),
                     # Heights of the two rows
                     widths = c(1)) # Widths of the two columns
    par(mar = c(4, 4, 1, 1))

    plot(
      x = data[, 1],
      y = data[, 2],
      type = "l",
      xlab = "meters",
      ylab = "cm/kyr (ka)",
      main = "sedimentation rate plot"
    )

  }

}



if (output == 2) {
  data <- data[, c(1, 4)]
  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    layout.matrix <- matrix(c(1), nrow = 1 , ncol = 1)
    graphics::layout(mat = layout.matrix,
                     heights = c(1),
                     # Heights of the two rows
                     widths = c(1)) # Widths of the two columns
    par(mar = c(4, 4, 1, 1))
    plot(
      data[, 1],
      data[, 2],
      type = "l",
      xlab = "meters",
      ylab = "Time (ka)",
      main = "Data time domain"
    )

  }

}

return(data)

}



