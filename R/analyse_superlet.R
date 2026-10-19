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
#'
#' @param c1 Base number of cycles of the Morlet wavelet. Acts as the
#' fundamental wavelet width.
#'
#' @param o numeric vector of length two defining the minimum and
#' maximum superlet order.
#'
#' @param mult Logical. If TRUE, the number of cycles increases
#' multiplicatively with superlet order (\eqn{c_i = i c_1}). If FALSE, cycles increase
#' additively (\eqn{c_i = c_1 + i - 1}, Moca et al., 2021).
#'
#' @param order_scaling Character string defining how the superlet order
#' increases from \code{o[1]} (lowest frequency) to \code{o[2]} (highest frequency). One of
#' \code{"none"} (default, same as \code{"log2"}: linear in log frequency), \code{"log2"}, \code{"linear"} (linear in
#' frequency, as in Moca et al., 2021), \code{"sqrt"}, \code{"quadratic"}
#' or \code{"power"}.
#'
#' @param order_alpha Numeric. Exponent used when \code{scaling = "power"}.
#' Values greater than one emphasize high frequency sharpening, whereas
#' values smaller than one emphasize low frequency sharpening.
#'
#' @param verbose Logical. If TRUE, print progress and interpolation
#' information.
#'
#' @param detrend Removal of a trend before the transform: \code{"linear"} (default),
#' \code{"mean"} or \code{"none"}. A mean or trend that is not removed leaks into all
#' periods and gives large edge effects, which is a problem for stratigraphic records
#' (for example grey values or counts) that do not have a zero mean.
#'
#' @param run_multicore Logical. Compute the frequencies on several cores \code{Default=FALSE}.
#' @param n_cores Number of cores used when \code{run_multicore = TRUE} (\code{NULL}: a safe automatic number;
#' see \code{\link{analyse_wavelet}}). Results do not depend on the number of cores.
#' @param standardize Logical. If TRUE (default) the record is divided by its standard deviation, as in
#' \code{\link{analyse_wavelet}}, so that the power does not depend on the units of the proxy.
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
#'order: superlet order of every period
#'coi: distance from the start and end of the record within which edge effects are important (cone of influence), for every period
#'coi.1, coi.2: the cone of influence as polygons (x and log2 period), drawn by the plot functions
#'x: interpolated x values
#'y: interpolated signal values
#'
#' @author
#' Adapted from MATLAB implementations by V. V. Moca et al. (2021)
#'
#' @references
#' Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., & Mureșan, R. C. (2021).
#' Time-frequency super-resolution with superlets.
#' Nature Communications, 12(1), 337.
#' \doi{10.1038/s41467-020-20539-9}
#'
#' @examples
#' \donttest{
#' ## Example 1. Using the Total Solar Irradiance data set of Steinhilber et al. (2012)
#' TSI_sl <-
#'   analyse_superlet(
#'     data = TSI,
#'     upperPeriod = 8192,
#'     lowerPeriod = 16,
#'     Nf = 20,
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
#'   analyse_superlet(
#'     data = mag,
#'     upperPeriod = 254,
#'     lowerPeriod = 0.1,
#'     Nf = 20,
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
#'   analyse_superlet(
#'     data = grey,
#'     upperPeriod = 256,
#'     lowerPeriod = 0.02,
#'     Nf = 20,
#'     c1 = 2,
#'     o = c(2,5),
#'     mult = TRUE,
#'     order_scaling = "log2",
#'     order_alpha = 1,
#'     verbose = FALSE)
#' }
#'
#' @param wavelet Mother wavelet, given as a name (PyWavelets style names with the order appended, such as
#' \code{"db4"}, \code{"sym8"}, \code{"cgau3"}, \code{"gaus2"} or \code{"cmor1.5-1"}, also work) or as a function
#' \code{psi_hat(omega)} that returns the Fourier transform of a user defined mother wavelet.
#' Analytic (complex) wavelets, which give phase and amplitude and can be synchrosqueezed: \code{"morlet"}
#' (default; aliases \code{"gabor"}, \code{"amor"}), \code{"morlet0"} (exactly zero mean Morlet, admissible for a
#' small omega_nr; alias \code{"morlex"}), \code{"paul"} (alias \code{"cauchy"}), \code{"morse"} (generalised Morse
#' wavelet of Olhede and Walden, 2002; aliases \code{"gmw"}, \code{"airy"}), \code{"bump"}, \code{"cmor"} (complex
#' Morlet with bandwidth and centre frequency as in Matlab/PyWavelets), \code{"cgau"} (analytic Gaussian derivative;
#' \code{"cmhat"} is the complex Mexican hat), \code{"shannon"} (complex Shannon, box spectrum), \code{"fbsp"}
#' (frequency B-spline), \code{"meyer"}, \code{"loggabor"} (alias \code{"lognormal"}), \code{"gammatone"} (asymmetric in
#' time: sharp onset, slow decay) and \code{"hann"} (raised cosine band).
#' Real wavelets, which give a real transform (no phase or envelope): \code{"dog"} (derivative of Gaussian;
#' \code{"mexh"}, \code{"ricker"}, \code{"marr"} and \code{"sombrero"} are the Mexican hat, order 2), \code{"morl"} (real
#' Morlet), \code{"hhat"} (Hermitian hat), \code{"poisson"}, and the orthogonal wavelets \code{"haar"}, \code{"db"}
#' (Daubechies), \code{"sym"} (Symlets) and \code{"battle"} (Battle-Lemarie spline wavelets), which are computed in the
#' frequency domain from their filters. The period axis of the orthogonal wavelets is approximate because their
#' spectra are broad and not symmetric. A narrow wavelet (Paul, low-order DOG, Shannon, Haar) resolves position
#' well, a wide one (Morse with a large beta, Morlet with a high omega_nr, bump) resolves period well.
#' @param wavelet_order Parameter(s) of the wavelet (see \code{\link{analyse_wavelet}}); when it is NULL the wavelet with the same relative bandwidth as a Morlet wavelet with the superlet cycles is used at every frequency (cycles are converted as omega_0 = 2 pi cycles / 5). Parameters: \code{paul}: order m (default 4); \code{dog}: derivative m
#' (default 2); \code{morse}: \code{c(gamma, beta)} or beta alone with gamma = 3 (default \code{c(3, 20)});
#' \code{bump}: \code{c(mu, sigma)} (default \code{c(5, 0.6)}); \code{cmor}: \code{c(bandwidth, centre frequency)}
#' (default \code{c(1.5, 1)}); \code{cgau}: derivative m (default 4); \code{shannon}: \code{c(bandwidth, centre
#' frequency)} (default \code{c(2, 1.5)}); \code{fbsp}: \code{c(m, bandwidth, centre frequency)} (default
#' \code{c(2, 2, 1.5)}); \code{loggabor}: \code{c(centre angular frequency, log-width)} (default \code{c(6, 0.25)});
#' \code{meyer} and \code{hhat} have no parameter; \code{morlet0}: omega_0 (default omega_nr); \code{morl}: omega_0 (default 5);
#' \code{poisson}: order n (default 2); \code{gammatone}: \code{c(n, centre angular frequency, bandwidth)} (default
#' \code{c(4, 6, 1.5)}); \code{hann}: \code{c(centre angular frequency, half width)} (default \code{c(6, 3)});
#' \code{db}: number of vanishing moments 1 to 20 (default 4, \code{"haar"} = \code{"db1"}); \code{sym}: 2 to 20
#' (default 8); \code{battle}: spline order 1, 3, 5 or 7 (default 3). For the Morlet wavelet use \code{omega_nr}.
#'
#' @param standardise British spelling of \code{standardize}; if given it overrides \code{standardize}.
#' @export
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats fft
#' @importFrom stats fft

analyse_superlet <- function(data,
                             upperPeriod,
                             lowerPeriod,
                             Nf,
                             c1,
                             o = NULL,
                             mult = TRUE,
                             order_scaling = c("none", "log2", "linear", "sqrt", "quadratic", "power"),
                             order_alpha = 1,
                             verbose = FALSE,
                             detrend = c("linear", "mean", "none"),
                             standardize = TRUE,
                             run_multicore = FALSE,
                             n_cores = NULL,
  standardise = NULL,
                             wavelet = "morlet",
                             wavelet_order = NULL) {
  if (!is.null(standardise)) standardize <- standardise  # British spelling synonym

  order_scaling <- match.arg(order_scaling)
  detrend <- match.arg(detrend)

  ## 1. Data preparation: sort, remove missing values, interpolate to the median spacing
  dat <- as.matrix(data)
  dat <- na.omit(dat)
  dat <- dat[order(dat[, 1], na.last = NA, decreasing = FALSE), , drop = FALSE]
  npts <- nrow(dat)
  dt <- median(diff(dat[, 1]))
  xout <- seq(from = dat[1, 1], to = dat[npts, 1], by = dt)
  dat <- as.data.frame(approx(dat[, 1], dat[, 2], xout, method = "linear"))
  if (verbose) cat("Dataset interpolated to spacing:", round(dt, 10), "\n")
  if (lowerPeriod < 2 * dt) {
    warning("lowerPeriod is shorter than twice the sampling interval (", signif(2 * dt, 4),
            "): periods below the Nyquist period cannot be resolved")
  }

  x_axis <- dat[, 1]
  input <- dat[, 2]
  Npoints <- length(input)
  # a mean or trend leaks into all periods and produces large edge effects: remove it
  if (detrend == "mean") input <- input - mean(input)
  if (detrend == "linear") input <- input - stats::fitted(stats::lm(input ~ seq_len(Npoints)))
  if (standardize && sd(input) > 0) input <- input / sd(input)

  ## 2. Frequency grid (cycles per sample), log spaced
  Fs <- 1
  Fi <- sort(c(1 / upperPeriod, 1 / lowerPeriod) * dt)
  Freqs <- 2^seq(log2(Fi[1]), log2(Fi[2]), length.out = Nf)

  ## 3. Superlet order per frequency (Moca et al. 2021: the order grows with frequency from o[1] to o[2])
  if (is.null(o)) o <- c(1, 1)
  if (length(o) == 1) o <- c(o, o)
  u_lin <- (Freqs - min(Freqs)) / (max(Freqs) - min(Freqs))
  u_log <- (log2(Freqs) - min(log2(Freqs))) / (max(log2(Freqs)) - min(log2(Freqs)))
  if (Nf == 1) u_lin <- u_log <- 0
  u_eff <- switch(order_scaling,
                  none = u_log,
                  log2 = u_log,
                  linear = u_lin,
                  sqrt = sqrt(u_lin),
                  quadratic = u_lin^2,
                  power = u_lin^order_alpha)
  order_frac <- o[1] + u_eff * (o[2] - o[1])
  order_int <- ceiling(order_frac)

  ## 4. Complex Morlet wavelet with Nc cycles (Moca et al. 2021)
  wl_of <- function(Fc, Nc) 2 * floor((6 * (Nc / 2) * abs(1 / Fc) / 2.5 * Fs) / 2) + 1
  cxmorlet <- function(Fc, Nc) {
    sd <- (Nc / 2) * abs(1 / Fc) / 2.5
    wl <- 2 * floor((6 * sd * Fs) / 2) + 1
    off <- floor(wl / 2)
    t <- (seq_len(wl) - 1 - off) / Fs
    g <- exp(-(t^2) / (2 * sd^2)) / (sd * sqrt(2 * pi))
    g * exp(2i * pi * Fc * t) / sum(g)
  }
  cycles_of <- function(i_ord) if (mult) i_ord * c1 else c1 + i_ord - 1
  use_morlet <- identical(.wr_canonical_name(wavelet), "morlet") && is.null(wavelet_order)
  # Morlet wavelet of Nc cycles: sd = Nc / 2 / 2.5 periods, i.e. omega_0 = 2 * pi * Nc / 5
  omega_of <- function(Nc) 2 * pi * Nc / 5
  wl_cache <- list()
  wdef_of <- function(Nc) {
    key <- as.character(round(omega_of(Nc), 2))
    if (is.null(wl_cache[[key]])) wl_cache[[key]] <<- .wr_wavelet_equiv(wavelet, omega_of(Nc), wavelet_order, ngrid = 4001)
    wl_cache[[key]]
  }

  ## 5. Padding and FFT length depend on the longest wavelet (only the lengths are needed here, the wavelets
  ## themselves are made per frequency so that they do not have to be kept in memory)
  padding <- 0
  max_wl <- 0
  for (i_freq in seq_len(Nf)) {
    for (i_ord in seq_len(order_int[i_freq])) {
      wl <- wl_of(Freqs[i_freq], cycles_of(i_ord))
      if (!use_morlet) wl <- wl * 2   # the frequency domain wavelets are not truncated: leave generous room
      padding <- max(padding, floor(wl / 2))
      max_wl <- max(max_wl, wl)
    }
  }
  nfft <- 2^ceiling(log2(Npoints + 2 * padding + max_wl))
  inv_nfft <- 1 / nfft

  # cone of influence: e-folding distance of the Gaussian envelope of the (fractional order) widest wavelet; this
  # varies smoothly with the period because the cycles follow the fractional order
  cycles_eff <- if (mult) order_frac * c1 else c1 + order_frac - 1
  coi <- sqrt(2) * ((cycles_eff / 2) * (1 / Freqs) / 2.5) * dt
  if (!use_morlet) {
    # edge zone relative to the Morlet wavelet of the same bandwidth: time width = decay time * scale
    # = decay time * period / Fourier factor
    nc_mid <- mean(range(cycles_eff))
    wd <- wdef_of(nc_mid)
    wm <- .wr_wavelet("morlet", omega_nr = omega_of(nc_mid), ngrid = 4001)
    coi <- coi * (.wr_decay_time(wd) / wd$fourier_factor) / (.wr_decay_time(wm) / wm$fourier_factor)
  }

  ## 6. Superlet transform: geometric mean of the wavelet amplitudes, circular mean of the phases
  xpad <- complex(nfft)
  xpad[(padding + 1):(padding + Npoints)] <- input
  X <- fft(xpad)
  tiny <- .Machine$double.xmin
  omega_k <- 2 * pi * c(0:(nfft / 2), -((nfft / 2 - 1):1)) / nfft     # angular frequency per sample
  one_frequency <- function(i_freq) {
    n_wavelets <- floor(order_frac[i_freq])
    frac <- order_frac[i_freq] - n_wavelets
    log_amp <- numeric(Npoints)            # sum of log amplitudes (a plain product can underflow)
    phasor <- complex(Npoints)
    for (i_ord in seq_len(order_int[i_freq])) {
      wt <- if (i_ord <= n_wavelets) 1 else frac        # the last, fractional, order has a fractional weight
      if (wt > 0) {
        if (use_morlet) {
          w <- cxmorlet(Freqs[i_freq], cycles_of(i_ord))
          hpad <- complex(nfft)
          hpad[seq_along(w)] <- rev(Conj(w))
          center <- (length(w) + 1) / 2
          idx <- (center + padding):(center + padding + Npoints - 1)
          res <- fft(X * fft(hpad), inverse = TRUE)[idx] * inv_nfft
        } else {
          # frequency domain wavelet, amplitude calibrated like the Morlet version (a unit sinusoid at the
          # peak frequency has an amplitude of 1/2 in the complex coefficient)
          wd <- wdef_of(cycles_of(i_ord))
          scale_i <- (1 / Freqs[i_freq]) / wd$fourier_factor
          res <- fft(X * wd$shape(scale_i, omega_k), inverse = TRUE)[(padding + 1):(padding + Npoints)] * inv_nfft
        }
        amp <- pmax(Mod(res), tiny)
        log_amp <- log_amp + wt * log(amp)
        phasor <- phasor + wt * res / amp
      }
    }
    pm <- Mod(phasor)
    unit <- ifelse(pm > 0, phasor / pmax(pm, tiny), 1 + 0i)
    exp(log_amp / order_frac[i_freq]) * unit
  }
  rows <- wr_lapply(seq_len(Nf), one_frequency, n_cores = n_cores, run_multicore = run_multicore,
                    combine = "list", progress = verbose)
  wave_mat <- matrix(unlist(rows), nrow = Nf, ncol = Npoints, byrow = TRUE)

  ## 8. Output (lowest frequency = longest period first, as before)
  wave_mat <- wave_mat[nrow(wave_mat):1, , drop = FALSE]
  Periods_phys <- rev(1 / (Freqs / dt))
  coi <- rev(coi)
  # the adaptive superlet uses many cycles at short periods and few at long ones, so the widest wavelet (and the
  # edge zone) can be shorter at a long period than at a shorter one, which gave an hourglass-shaped cone; the cone
  # is made monotonic in the period (conservative: never narrower than at any shorter period)
  coi <- cummax(coi)

  # cone of influence as one polygon: the edge zones at the start and end of the record
  lp <- log2(Periods_phys)
  xl <- x_axis[1] - 0.5 * dt
  xr <- x_axis[Npoints] + 0.5 * dt
  # one polygon: the left and right edge zones meet in the middle of the record where they would overlap
  cc <- pmin(coi, (xr - xl) / 2)
  coi.1 <- c(xl, xl + cc, xr - rev(cc), xr, xr, xl)
  coi.2 <- c(lp[1], lp, rev(lp), lp[1], lp[Nf], lp[Nf])

  output <- list(
    Wave = wave_mat,
    Power = t(Mod(wave_mat)^2),
    dt = dt,
    dj = 1 / Nf,
    Period = Periods_phys,
    Power.avg = rowMeans(Mod(wave_mat)^2),
    Scale = Periods_phys,
    nc = Npoints,
    nr = Nf,
    axis.1 = x_axis,
    axis.2 = log2(Periods_phys),
    c1 = c1,
    o = o,
    order_scaling = order_scaling,
    order_alpha = order_alpha,
    mult = mult,
    detrend = detrend,
    standardize = standardize,
    wavelet = .wr_canonical_name(wavelet),
    wavelet_order = wavelet_order,
    order = rev(order_frac),
    coi = coi,
    coi.1 = coi.1,
    coi.2 = coi.2,
    x = x_axis,
    y = input
  )
  class(output) <- "analyze.superlet"
  return(output)
}
