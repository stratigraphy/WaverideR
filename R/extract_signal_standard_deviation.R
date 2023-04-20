#' @title Extract a signal  using standard deviation
#'
#' @description Extract signal from a
#' wavelet spectra in the depth domain using a the standard deviation of the omega (number of cycles)
#' as boundaries. The uncertainty is based on the Gabor uncertainty principle applied to the
#' continuous wavelet transform using a Morlet wavelet. The calculated uncertainty is the underlying
#' analytical uncertainty which is the result of applying the Gabor uncertainty principle to the
#' continuous wavelet transform using a Morlet wavelet.
#'
#'
#' @param wavelet Wavelet object created using the \code{\link{analyze_wavelet}} function.
#' @param tracked_cycle_curve Curve of the cycle tracked using the
#' \code{\link{track_period_wavelet}} function. Any input (matrix or dataframe)
#'  in which the first column is depth or time and the second column is period should work.
#' @param multi multiple of the standard deviation to be used as boundaries for the cycle extraction
#'  \code{Default=1}.
#' @param extract_cycle Period of the cycle to be extracted.
#' @param tracked_cycle_period Period of the tracked cycle.
#' @param add_mean Add mean to the extracted cycle \code{Default=TRUE}.
#' @param tune Tune data set using the \code{Default=tracked_cycle_curve} curve \code{Default=FALSE}.
#' @param genplot_uncertainty_wt Generate a wavelet spectra plot with the tracked curve and its
#' analytical uncertainty based the Gabor uncertainty principle applied
#' continuous wavelet transform using a Morlet wavelet on superimposed on top of it.
#' In the plot the red curve and blue curves are the upper and lower bounds
#'based on the \code{multi} parameter which x-times the standard deviation of uncertainty.
#'The black curve is the \code{Default=tracked_cycle_curve} curve.
#' @param genplot_extracted Generates a plot with the data set and
#' the extracted cycle on top \code{Default=TRUE} of it.
#'
#' @references
#' Gabor, Dennis. "Theory of communication. Part 1: The analysis of information."
#' Journal of the Institution of Electrical Engineers-part III: radio and
#' communication engineering 93, no. 26 (1946): 429-441.
#'
#'Russell, Brian, and Jiajun Han. "Jean Morlet and the continuous wavelet transform.
#'" CREWES Res. Rep 28 (2016): 115.
#'
#'@examples
#'\donttest{
#'#Extract the 405 kyr eccentricity cycle from the the magnetic susceptibility \cr
#'#record of the Sullivan core and use the Gabor uncertainty principle to define \cr
#'# the mathematical uncertainty of the analysis and use a factor of that standard \cr
#'#  deviation to define boundaries
#'
#'# perform the CWT
#'mag_wt <- analyze_wavelet(data = mag,
#' dj = 1/100,
#' lowerPeriod = 0.1,
#' upperPeriod = 254,
#' verbose = TRUE,
#' omega_nr = 10)
#'
#'#Track the 405 kyr eccentricity cycle in a wavelet spectra
#'
#'#mag_track <- track_period_wavelet(astro_cycle = 405,
#'#                                   wavelet=mag_wt,
#'#                                   n.levels = 100,
#'#                                   periodlab = "Period (metres)",
#'#                                   x_lab = "depth (metres)")
#'
#'#Instead of tracking, the tracked solution data set \code{\link{mag_track_solution}} is used \cr
#'mag_track <- mag_track_solution
#'
#' mag_track_complete <- completed_series(
#'   wavelet = mag_wt,
#'   tracked_curve = mag_track,
#'   period_up = 1.2,
#'   period_down = 0.8,
#'   extrapolate = TRUE,
#'   genplot = TRUE
#' )
#'
#'# smooth the tracking of the 405 kyr eccentricity cycle
#' mag_track_complete <- loess_auto(time_series = mag_track_complete,
#' genplot = TRUE, print_span = TRUE)
#'
#'# extract the 405 kyr eccentricty cycle from the wavelet spectrum and use \cr
#'# the Gabor uncertainty principle to define the mathematical uncertainty of \cr
#'# the analysis and use a multiple of the derived standard deviation to define boundaries
#'
#'mag_405_ecc <- extract_signal_standard_deviation(
#'wavelet = mag_wt,
#'tracked_cycle_curve = mag_track_complete,
#'multi = 1,
#'extract_cycle = 405,
#'tracked_cycle_period = 405,
#'add_mean = TRUE,
#'tune = FALSE,
#'genplot_uncertainty_wt = TRUE,
#'genplot_extracted = TRUE
#')
#'}
#' @return Signal extracted from the wavelet spectra.
#' Output is a matrix with the first column being depth/time
#'and the second column is the astronomical cycle extracted from the proxy record
#'
#'If \code{genplot_uncertainty_wt=TRUE} then a wavelet spectra will be plotted
#'with the uncertainty superimposed on top of it. In the plot the red curve and
#' blue curves are the upper and lower bounds
#'based on the \code{multi} parameter.The black curve is the \code{Default=tracked_cycle_curve} curve.
#'If \code{genplot_extracted=TRUE} plot with the data set and
#'the extracted cycle on top of it will be plotted.
#'
#' @export
#' @importFrom Hmisc approxExtrap
#' @importFrom stats na.omit
#' @importFrom graphics par
#' @importFrom graphics hist
#' @importFrom graphics lines
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom DescTools Closest


extract_signal_standard_deviation <- function(
wavelet=NULL,
tracked_cycle_curve= NULL,
multi = 1,
extract_cycle = NULL,
tracked_cycle_period = NULL,
add_mean=TRUE,
tune=FALSE,
genplot_uncertainty_wt=TRUE,
genplot_extracted=TRUE
){my.w <- wavelet
  my.data <- cbind(wavelet$x,wavelet$y)
  filtered_cycle <- my.data[,1]
  filtered_cycle <- as.data.frame(filtered_cycle)
  filtered_cycle$value <- NA

  completed_series <- na.omit(tracked_cycle_curve)
  completed_series[,2] <- completed_series[,2]*(extract_cycle/tracked_cycle_period)
  app <- approxExtrap(x=completed_series[,1],y=completed_series[,2],xout=my.data[,1],
                      method="linear")
  interpolated <- cbind(app$x,app$y)

  Wave = my.w$Wave
  Power = my.w$Power

  nc = my.w$nc
  nr = my.w$nr
  dt = my.w$dt
  dj = my.w$dj

  Scale = my.w$Scale
  Period = my.w$Period
  loess.span = my.w$loess.span
  rec.waves = matrix(0, nrow = nr, ncol = nc)


  for (s.ind in seq_len(nr)) {
    rec.waves[s.ind, ] = (Re(Wave[s.ind, ])/sqrt(Scale[s.ind])) *
      dj * sqrt(dt)/(pi^(-1/4))}


  interpolated <- as.data.frame(interpolated)

  ncycles <- my.w$omega_nr
  b <- (2*sqrt(2*log(2)))
  a <- ((8*log(2)/(2*pi)))
  k <- (ncycles/(8*log(2)))*2



  interpolated$f0 <- (1/(interpolated[,2]))
  interpolated$df <- (a*interpolated$f0)/k
  interpolated$sd_morlet <- interpolated$df/b


  fact_high <- 1/(interpolated$f0-(interpolated$sd_morlet*multi)
  )
  fact_low <-  1/(interpolated$f0+(interpolated$sd_morlet*multi)
  )

  fact_high[fact_high>max(my.w$Period)] <- max(my.w$Period)
  fact_low[fact_low<min(my.w$Period)] <- min(my.w$Period)


  interpolated[,3] <-fact_high
  interpolated[,4]<- fact_low

  for (i in 1:nrow(filtered_cycle)){
    row_nr_high <- Closest(Period[],interpolated[i,3],which=TRUE)
    row_nr_low <- Closest(Period[],interpolated[i,4],which=TRUE)

    row_nr_high <-row_nr_high[1]
    row_nr_low <- row_nr_low[1]

    value <- rec.waves[c(row_nr_low:row_nr_high),i]
    value  <- sum(value, na.rm = T)
    value <- as.numeric(value)
    filtered_cycle[i,2] <- value
  }

  rec_value  <- colSums(rec.waves, na.rm = T)

  filtered_cycle[,2] <- (filtered_cycle[,2]) * sd(my.data[,2])/sd(rec_value)

  if(add_mean==TRUE){
    filtered_cycle[,2] <- filtered_cycle[,2] + mean(my.data[,2])
  }


  if (genplot_uncertainty_wt==TRUE & tune==FALSE){
    plot_wavelet(
      wavelet = wavelet,
      plot.COI = TRUE,
      n.levels = 100,
      color.palette = "rainbow(n.levels, start = 0, end = 0.7)",
      useRaster = TRUE,
      periodlab = "Period (metres)",
      x_lab = "depth (metres)")


    combined_sedrate <-cbind(interpolated[,1],1/interpolated[,2],interpolated$sd_morlet*multi)

    xcords <- c(combined_sedrate[,1],sort(combined_sedrate[,1],decreasing = TRUE))
    xcords
    data_sort1 <- combined_sedrate[order(combined_sedrate[,1],decreasing = TRUE), ]
    ycords <- c(1/(combined_sedrate[,2]+combined_sedrate[,3]),1/(data_sort1[,2]-data_sort1[,3]))

    polygon(x = xcords,
            y = log2(ycords),
            col=rgb(0.5, 0.5, 0.5,0.5),
            border = "black")

    lines(interpolated[,1],log2(interpolated[,2]),lwd=2)
    lines(interpolated[,1],log2(interpolated[,3]),col="red",lwd=2)
    lines(interpolated[,1],log2(interpolated[,4]),col="blue",lwd=2)

  }


  if (genplot_uncertainty_wt==TRUE & tune==TRUE){
    data_set <- cbind(wavelet$x,wavelet$y)
    data_set_time <- curve2tune(
      data = data_set,
      tracked_cycle_curve = tracked_cycle_curve,
      tracked_cycle_period = tracked_cycle_period,
      genplot = TRUE
    )

    dat <- as.matrix(data_set_time)
    dat <- na.omit(dat)
    dat <- dat[order(dat[, 1], na.last = NA, decreasing = F),]
    npts <- length(dat[, 1])
    start <- dat[1, 1]
    end <- dat[length(dat[, 1]), 1]
    x1 <- dat[1:(npts - 1), 1]
    x2 <- dat[2:(npts), 1]
    dx = x2 - x1
    dt = median(dx)
    steps_size <- dt
    xout <- seq(start, end, by = dt)
    npts <- length(xout)
    interp <- approx(dat[, 1], dat[, 2], xout, method = "linear",
                     n = npts)
    data_set_time <- as.data.frame(interp)

    data_set_time_wt <-
      analyze_wavelet(
        data = data_set_time,
        dj = 1/200,
        lowerPeriod = steps_size,
        upperPeriod =  data_set_time[nrow(data_set_time),1]-data_set_time[1,1],
        verbose = FALSE,
        omega_nr = wavelet$omega_nr)

  plot_wavelet(
      wavelet = data_set_time_wt,
      plot.COI = TRUE,
      n.levels = 100,
      color.palette = "rainbow(n.levels, start = 0, end = 0.7)",
      useRaster = TRUE,
      periodlab = "Period (kyr)",
      x_lab = "depth (metres)")


interpolated_time <- curve2tune(
    data = interpolated[,c(1,2)],
    tracked_cycle_curve = tracked_cycle_curve,
    tracked_cycle_period = tracked_cycle_period,
    genplot = FALSE)


combined_sedrate <-cbind(interpolated_time[,1],1/(interpolated[,2]/(tracked_cycle_curve[,2]*tracked_cycle_period))
                         ,interpolated$sd_morlet*multi*(tracked_cycle_curve[,2]*tracked_cycle_period))
    xcords <- c(combined_sedrate[,1],sort(combined_sedrate[,1],decreasing = TRUE))
    data_sort1 <- combined_sedrate[order(combined_sedrate[,1],decreasing = TRUE), ]
    ycords <- c(1/(combined_sedrate[,2]+combined_sedrate[,3]),1/(data_sort1[,2]-data_sort1[,3]))
    polygon(x = xcords,
            y = log2(1/ycords),
            col=rgb(0.5, 0.5, 0.5,0.5),
            border = "black")

    lines(interpolated_time[,1],log2(interpolated[,2]),lwd=2)
    lines(interpolated_time[,1],log2((combined_sedrate[,2]+combined_sedrate[,3])),col="red",lwd=2)
    lines(interpolated_time[,1],log2((combined_sedrate[,2]-combined_sedrate[,3])),col="blue",lwd=2)

  }



  if(tune==TRUE){
    filtered_cycle <- curve2tune(
      data = filtered_cycle,
      tracked_cycle_curve = tracked_cycle_curve,
      tracked_cycle_period = tracked_cycle_period,
      genplot = FALSE)
    }


  if(genplot_extracted==TRUE){
    dev.new(width=7,height=4,noRStudioGD = TRUE)
    if(tune==TRUE){
      my.data_time <- curve2tune(
        data = my.data,
        tracked_cycle_curve = tracked_cycle_curve,
        tracked_cycle_period = tracked_cycle_period,
        genplot = FALSE)
      plot(my.data_time,type="l")
    }else{
      plot(my.data[,1], my.data[,2],type="l")
      }

    if(add_mean==TRUE){
      lines(filtered_cycle[,1],filtered_cycle[,2],col="blue")
    }else{
      lines(filtered_cycle[,1],filtered_cycle[,2]+ mean(my.data[,2]),col="blue")
      }


  }

  filtered_cycle
}






