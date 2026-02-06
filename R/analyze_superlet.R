#' @title Conduct the superlet transform on a time series or signal
#'
#' @description
#' Compute the superlet transform using a complex Morlet wavelet following
#' the approach of Moca et al. (2021). The superlet transform increases
#' time frequency resolution by combining multiple wavelet orders at each
#' frequency. Both the frequency sampling and the scaling of superlet order
#' with frequency can be explicitly controlled.
#'
#' @param data Input data as a matrix or data frame. The first column must
#' represent depth or time, and the second column the signal or proxy record.
#'
#' @param upperPeriod Maximum period to be analysed. Controls the lowest
#' analysed frequency.
#'
#' @param lowerPeriod Minimum period to be analysed. Controls the highest
#' analysed frequency.
#'
#' @param Nf Number of frequencies used to construct the superlet spectrum.
#'
#' @param freq.scale Character string defining how frequencies are sampled.
#' One of \code{"log2"}, \code{"linear"}, or \code{"sqrt"}.
#' Default is \code{"log2"}, following the original superlet formulation.
#'
#' @param c1 Base number of cycles of the Morlet wavelet. Acts as the
#' fundamental wavelet width.
#'
#' @param o numeric vector of length two defining the minimum and
#' maximum superlet order.
#'
#' @param mult Logical. If TRUE, the number of cycles increases
#' multiplicatively with superlet order. If FALSE, cycles increase
#' additively.
#'
#' @param scaling Character string defining how the number of wavelet
#' cycles varies with scale. One of
#' \code{"log2"}, \code{"linear"}, \code{"sqrt"}, \code{"quadratic"},
#' or \code{"power"}.
#'
#' @param alpha Numeric. Exponent used when \code{scaling = "power"}.
#' Values greater than one emphasize high frequency sharpening, whereas
#' values smaller than one emphasize low frequency sharpening.
#'
#' @param verbose Logical. If TRUE, print progress and interpolation
#' information.
#'
#' @return
#'A list with class \code{"analyze.superlet"} containing
#'Power: time frequency power spectrum
#'dt: sampling interval after interpolation
#'dj: number of frequencies
#'Power.avg: average spectral power
#'Period: physical periods corresponding to frequencies
#'nc: number of columns in the input data
#'nr: number of frequency levels
#'axis.1: x axis values (time or depth)
#'axis.2: y axis values (log2 scaled periods)
#'c1: base number of wavelet cycles
#'o: numeric vector of length two defining the minimum and
#' maximum superlet order. If NULL, all frequencies are analysed using
#' order one.
#' freq.scale: scaling of the frequencies
#' order.mode: order.mode ("none", "log2", "linear", "sqrt", "power")
#'x: interpolated x values
#'y: interpolated signal values
#'
#' @author
#' Adapted from MATLAB implementations by V. V. Moca et al. (2021)
#'
#' @references
#' Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., and Mureșan, R. C. (2021).
#' Time frequency super resolution with superlets.
#' Nature Communications, 12, 337.
#' \url{https://doi.org/10.1038/s41467-020-20539-9}
#'
#' @examples
#' \donttest{
#' ## Example 1. Using the Total Solar Irradiance data set of Steinhilber et al. (2012)
#' TSI_sl <-
#'   analyze_superlet(
#'     data = TSI,
#'     upperPeriod = 8192,
#'     lowerPeriod = 16,
#'     Nf = 100,
#'     c1 = 2,
#'     o = c(1, 5),
#'     mult = TRUE,
#'     order_scaling = "log2",
#'     order_alpha = 1,
#'     verbose = FALSE
#'   )
#'
#'
#' ## Example 2. Using the magnetic susceptibility data set of Pas et al. (2018)
#' mag_sl <-
#'   analyze_superlet(
#'     data = mag,
#'     upperPeriod = 254,
#'     lowerPeriod = 0.1,
#'     Nf = 100,
#'     c1 = 2,
#'     o = c(1, 5),
#'     mult = TRUE,
#'     order_scaling = "log2",
#'     order_alpha = 1,
#'     verbose = FALSE
#'   )
#'
#'
#' ## Example 3. Using the greyscale data set of Zeeden et al. (2013)
#' grey_sl <-
#'   analyze_superlet(
#'     data = grey,
#'     upperPeriod = 256,
#'     lowerPeriod = 0.02,
#'     Nf = 100,
#'     c1 = 2,
#'     o = c(2,5),
#'     mult = TRUE,
#'     order_scaling = "log2",
#'     order_alpha = 1,
#'     verbose = FALSE)
#' }
#'
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom stats fft

analyze_superlet <- function(data,
                             upperPeriod,
                             lowerPeriod,
                             Nf,
                             c1,
                             o = NULL,
                             mult = TRUE,
                             order_scaling = "log2",
                             order_alpha = 1,
                             verbose = FALSE) {

  order_scaling <- match.arg(order_scaling)


  ## 1. Basic setup and data preparation


  Fs <- 1
  lower_freq <- 1 / upperPeriod
  upper_freq <- 1 / lowerPeriod

  dat <- as.matrix(data)
  dat <- na.omit(dat)
  dat <- dat[order(dat[, 1], na.last = NA, decreasing = FALSE), ]

  npts <- length(dat[, 1])
  dx <- dat[2:npts, 1] - dat[1:(npts - 1), 1]
  dt <- median(dx)

  xout <- seq(from = dat[1, 1], to = dat[npts, 1], by = dt)
  interp <- approx(dat[, 1], dat[, 2], xout, method = "linear")
  dat <- as.data.frame(interp)

  if (verbose) {
    cat("Dataset interpolated to spacing:", round(dt, 10), "\n")
  }

  x_axis <- dat[, 1]
  input <- dat[, 2]
  Npoints <- length(input)


  ## 2. Frequency grid


  Fi <- c(lower_freq * dt, upper_freq * dt)
  if (Fi[1] > Fi[2]) Fi <- rev(Fi)

  log2_Freqs <- seq(log2(Fi[1]), log2(Fi[2]), length.out = Nf)
  Freqs <- 2^log2_Freqs

  ## Normalised frequency coordinate (0–1)
  u <- (log2(Freqs) - min(log2(Freqs))) /
    (max(log2(Freqs)) - min(log2(Freqs)))


  ## 3. Superlet order definition with explicit scaling


  if (is.null(o) || order_scaling == "none") {

    order_frac <- rep(ifelse(is.null(o), 1, o[1]), Nf)

  } else {

    if (order_scaling == "linear") {
      u_eff <- (Freqs - min(Freqs)) / (max(Freqs) - min(Freqs))
    } else {
      u_eff <- u  # log2-based default
    }

    if (order_scaling == "sqrt") {
      u_eff <- sqrt(u_eff)
    }

    if (order_scaling == "power") {
      u_eff <- u_eff^order_alpha
    }

    order_frac <- o[1] + u_eff * (o[2] - o[1])
  }

  order_int <- ceiling(order_frac)


  ## 4. Complex Morlet wavelet


  cxmorlet <- function(Fc, Nc, Fs) {
    sd <- (Nc / 2) * abs(1 / Fc) / 2.5
    wl <- 2 * floor((6 * sd * Fs) / 2) + 1
    off <- floor(wl / 2)
    t <- (seq_len(wl) - 1 - off) / Fs
    g <- exp(-(t^2) / (2 * sd^2)) / (sd * sqrt(2 * pi))
    w <- g * exp(2i * pi * Fc * t)
    w / sum(g)
  }


  ## 5. Precompute wavelets


  wavelets <- vector("list", Nf)
  padding <- 0
  max_wl <- 0

  for (i_freq in seq_len(Nf)) {

    wavelets[[i_freq]] <- vector("list", order_int[i_freq])

    for (i_ord in seq_len(order_int[i_freq])) {

      n_cyc <- if (mult) i_ord * c1 else i_ord + c1
      w <- cxmorlet(-Freqs[i_freq], n_cyc, Fs)

      wavelets[[i_freq]][[i_ord]] <- list(wlen = length(w), w = w)

      padding <- max(padding, floor(length(w) / 2))
      max_wl <- max(max_wl, length(w))
    }
  }


  ## 6. FFT preparation


  nfft <- 2^ceiling(log2(Npoints + 2 * padding + max_wl))
  inv_nfft <- 1 / nfft

  for (i_freq in seq_len(Nf)) {
    for (i_ord in seq_len(order_int[i_freq])) {

      w <- wavelets[[i_freq]][[i_ord]]$w
      hpad <- complex(nfft)
      hpad[seq_along(w)] <- rev(Conj(w))

      wavelets[[i_freq]][[i_ord]]$fft <- fft(hpad)

      center <- (length(w) + 1) / 2
      wavelets[[i_freq]][[i_ord]]$idx <-
        (center + padding):(center + padding + Npoints - 1)

      wavelets[[i_freq]][[i_ord]]$w <- NULL
    }
  }


  ## 7. Superlet transform


  wave_mat <- matrix(0, nrow = Nf, ncol = Npoints)

  xpad <- complex(nfft)
  xpad[(padding + 1):(padding + Npoints)] <- input
  X <- fft(xpad)

  for (i_freq in seq_len(Nf)) {

    n_wavelets <- floor(order_frac[i_freq])
    frac <- order_frac[i_freq] - n_wavelets
    geometric_mean <- rep(1 + 0i, Npoints)

    for (i_ord in seq_len(order_int[i_freq])) {

      wf <- wavelets[[i_freq]][[i_ord]]
      res <- fft(X * wf$fft, inverse = TRUE)[wf$idx] * inv_nfft

      if (i_ord <= n_wavelets) {
        geometric_mean <- geometric_mean * res
      } else if (frac > 0) {
        geometric_mean <- geometric_mean * (res^frac)
      }
    }

    wave_mat[i_freq, ] <- geometric_mean^(1 / order_frac[i_freq])
  }


  ## 8. Output


  wave_mat <- wave_mat[nrow(wave_mat):1, ]

  Periods_phys <- 1 / (Freqs / dt)
  Periods_phys <- rev(Periods_phys)

  output <- list(
    Wave = wave_mat,
    Power = t(Mod(wave_mat)^2),
    dt = dt,
    dj = 1 / Nf,
    Period = Periods_phys,
    Scale = Periods_phys,
    nc = Npoints,
    nr = Nf,
    axis.1 = x_axis,
    axis.2 = log2(Periods_phys),
    c1 = c1,
    o = o,
    order_scaling = order_scaling,
    order_alpha = order_alpha,
    x = x_axis,
    y = input
  )

  class(output) <- "analyze.superlet"
  return(output)
}




# odl function
# analyze_superlet <- function(data,
#                              upperPeriod,
#                              lowerPeriod,
#                              Nf,
#                              c1,
#                              o = NULL,
#                              mult = TRUE,
#                              verbose = FALSE) {
#   ## 1. Basic setup and frequency limits
#
#
#   # Sampling frequency is normalised to 1 because the time or depth
#   # spacing is handled explicitly via dt
#   Fs <- 1
#
#   # Convert period limits to physical frequencies
#   # upperPeriod corresponds to the lowest frequency
#   # lowerPeriod corresponds to the highest frequency
#   lower_freq <- 1 / upperPeriod
#   upper_freq <- 1 / lowerPeriod
#
#
#
#   ## 2. Data preparation and interpolation
#
#
#   # Convert input to matrix and remove missing values
#   dat <- as.matrix(data)
#   dat <- na.omit(dat)
#
#   # Ensure data are ordered in ascending time or depth
#   dat <- dat[order(dat[, 1], na.last = NA, decreasing = FALSE), ]
#
#   # Estimate sampling interval using the median spacing
#   npts <- length(dat[, 1])
#   dx <- dat[2:npts, 1] - dat[1:(npts - 1), 1]
#   dt <- median(dx)
#
#   # Interpolate data onto an evenly spaced grid
#   xout <- seq(from = dat[1, 1], to = dat[npts, 1], by = dt)
#   interp <- approx(dat[, 1], dat[, 2], xout, method = "linear")
#   dat <- as.data.frame(interp)
#
#   if (verbose) {
#     cat("Dataset interpolated to spacing:", round(dt, 10), "\n")
#   }
#
#   # Store axis and signal vectors
#   x_axis <- dat[, 1]
#   x <- dat[, 2]
#
#   # Input signal used for convolution
#   input <- dat[, 2]
#
#
#
#   ## 3. Define frequency grid
#
#
#   # Frequencies scaled by sampling interval
#   Fi <- c(lower_freq * dt, upper_freq * dt)
#
#   # Ensure correct ordering of frequency bounds
#   if (Fi[1] > Fi[2]) Fi <- rev(Fi)
#
#   # Logarithmically spaced frequencies following superlet theory
#   log2_Fi <- log2(Fi)
#   log2_Freqs <- seq(log2_Fi[1], log2_Fi[2], length.out = Nf)
#   Freqs <- 2^log2_Freqs
#
#
#
#   ## 4. Define superlet orders
#
#
#   # If an order range is supplied, smoothly vary order with frequency
#   # Otherwise, default to standard wavelet behaviour (order = 1)
#   if (!is.null(o)) {
#     order_frac <- seq(o[1], o[2], length.out = Nf)
#     order_int <- ceiling(order_frac)
#   } else {
#     order_frac <- rep(1, Nf)
#     order_int <- rep(1, Nf)
#   }
#
#
#
#   ## 5. Prepare signal buffers
#
#
#   # Convert vector input to a matrix for generality
#   if (is.vector(input)) {
#     input <- matrix(input, nrow = 1)
#   }
#
#   Nbuffers <- nrow(input)  # number of independent signals
#   Npoints <- ncol(input)   # number of samples per signal
#
#
#
#   ## 6. Define complex Morlet wavelet
#
#
#   # Generates a complex Morlet wavelet for a given frequency and number
#   # of oscillatory cycles
#   cxmorlet <- function(Fc, Nc, Fs) {
#     sd <- (Nc / 2) * abs(1 / Fc) / 2.5
#     wl <- 2 * floor((6 * sd * Fs) / 2) + 1
#     off <- floor(wl / 2)
#     t <- (seq_len(wl) - 1 - off) / Fs
#     g <- exp(-(t^2) / (2 * sd^2)) / (sd * sqrt(2 * pi))
#     w <- g * exp(2i * pi * Fc * t)
#     w / sum(g)
#   }
#
#
#
#   ## 7. Precompute wavelets and FFT padding
#
#
#   wavelets <- vector("list", Nf)
#   padding <- 0
#   max_wl <- 0
#
#   # Generate all wavelets in advance for efficiency
#   for (i_freq in seq_len(Nf)) {
#     wavelets[[i_freq]] <- vector("list", order_int[i_freq])
#
#     for (i_ord in seq_len(order_int[i_freq])) {
#
#       # Number of cycles increases with order
#       n_cyc <- if (mult) i_ord * c1 else i_ord + c1
#
#       w <- cxmorlet(-Freqs[i_freq], n_cyc, Fs)
#
#       wavelets[[i_freq]][[i_ord]] <- list(wlen = length(w), w = w)
#
#       padding <- max(padding, floor(length(w) / 2))
#       max_wl <- max(max_wl, length(w))
#     }
#   }
#
#
#
#   ## 8. FFT preparation
#
#
#   # FFT length chosen as next power of two for speed
#   nfft <- 2^ceiling(log2(Npoints + 2 * padding + max_wl))
#   inv_nfft <- 1 / nfft
#
#
#
#   ## 9. Precompute FFTs of wavelets
#
#
#   for (i_freq in seq_len(Nf)) {
#     for (i_ord in seq_len(order_int[i_freq])) {
#
#       w <- wavelets[[i_freq]][[i_ord]]$w
#       hpad <- complex(nfft)
#       hpad[seq_along(w)] <- rev(Conj(w))
#
#       wavelets[[i_freq]][[i_ord]]$fft <- fft(hpad)
#
#       center <- (length(w) + 1) / 2
#       wavelets[[i_freq]][[i_ord]]$idx <-
#         (center + padding):(center + padding + Npoints - 1)
#
#       # Free memory
#       wavelets[[i_freq]][[i_ord]]$w <- NULL
#     }
#   }
#
#
#
#   ## 10. Main superlet convolution loop
#
#
#   wtresult <- matrix(0, nrow = Nf, ncol = Npoints)
#   bufbegin <- padding + 1
#   bufend <- padding + Npoints
#
#   for (i_buf in seq_len(Nbuffers)) {
#
#     # FFT of padded signal
#     xpad <- complex(nfft)
#     xpad[bufbegin:bufend] <- input[i_buf, ]
#     X <- fft(xpad)
#
#     for (i_freq in seq_len(Nf)) {
#
#       temp <- 1
#       n_wavelets <- floor(order_frac[i_freq])
#       frac <- order_frac[i_freq] - n_wavelets
#
#       for (i_ord in seq_len(order_int[i_freq])) {
#
#         wf <- wavelets[[i_freq]][[i_ord]]
#         res <- fft(X * wf$fft, inverse = TRUE)
#         res_slice <- res[wf$idx]
#
#         mag2 <- (Re(res_slice)^2 + Im(res_slice)^2) * (inv_nfft^2)
#
#         if (i_ord <= n_wavelets) {
#           temp <- temp * (2 * mag2)
#         } else if (frac > 0) {
#           temp <- temp * (2 * mag2)^frac
#         }
#       }
#
#       wtresult[i_freq, ] <- wtresult[i_freq, ] +
#         temp^(1 / order_frac[i_freq])
#     }
#   }
#
#   wtresult <- wtresult[nrow(wtresult):1, ]
#
#   ## 11. Assemble output object
#   Periods_phys <- 1/(Freqs/dt)
#   output <- list(Power = t(wtresult/Nbuffers), dt = dt, dj = Nf,
#                  Power.avg = rev(rowMeans(wtresult/Nbuffers)), Period = Periods_phys,
#                  nc = length(x), nr = Nf, axis.1 = x_axis, axis.2 = rev(log2(Periods_phys)),
#                  c1 = c1,o=o, x = dat[, 1], y = dat[, 2])
#   class(output) <- "analyze.superlet"
#   return(invisible(output))
# }
